"""The DesignRecord schema and its I/O (plan 3).

"Every stage reads and writes this and nothing else."

Invariants enforced here rather than trusted (plan 3):
  * `scores` is APPEND-ONLY and namespaced by stage ("tmol.ddE_ts_sub").
    `add_scores` refuses to overwrite another stage's key.
  * No large arrays in the record -- write them to the design dir and store
    the path. `add_scores` rejects non-scalar values.
  * Resumability (plan 0.3): `has_stage` lets a stage skip designs that already
    carry its output, so a killed job can be re-run as-is.
  * Failure is recorded, never swallowed (plan 0.2): `fail` sets status and
    appends the exception text, and the batch continues.
"""

from __future__ import annotations

import hashlib
import json
import os
import tempfile
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Any, Iterator, Mapping

SCALAR = (int, float, bool, str, type(None))


def design_id(cell: str, backbone_seed: int, seq_seed: int) -> str:
    """Stable 12-hex-char id from (cell, backbone_seed, seq_seed) (plan 3).

    Stable across runs and machines, so a design keeps its identity when a
    stage is re-run or a run is resumed.
    """
    key = f"{cell}|{backbone_seed}|{seq_seed}"
    return hashlib.sha256(key.encode()).hexdigest()[:12]


@dataclass
class DesignRecord:
    """One design's complete state. Serialised as one JSON object per line."""

    design_id: str
    run_id: str
    cell: dict = field(default_factory=dict)
    provenance: dict = field(default_factory=dict)
    paths: dict = field(default_factory=dict)
    sequence: str | None = None
    length: int | None = None
    catalytic_residues: dict = field(default_factory=dict)  # 1-indexed
    scores: dict = field(default_factory=dict)
    gates: dict = field(default_factory=dict)
    status: str = "ok"
    errors: list = field(default_factory=list)

    # --- scores: append-only, namespaced, scalar ------------------------
    def add_scores(self, stage: str, values: Mapping[str, Any]) -> None:
        """Add `stage`-namespaced scores. Refuses overwrites and non-scalars."""
        for key, value in values.items():
            full = key if key.startswith(f"{stage}.") else f"{stage}.{key}"
            if not isinstance(value, SCALAR):
                raise TypeError(
                    f"score {full!r} is {type(value).__name__}; scores must be "
                    "scalar -- write arrays to the design dir and store the path"
                )
            if full in self.scores and self.scores[full] != value:
                raise ValueError(
                    f"score {full!r} already set to {self.scores[full]!r}; "
                    "scores are append-only, never overwrite another stage's key"
                )
            self.scores[full] = value

    def set_gate(self, stage: str, passed: bool, reasons: list[str] | None = None) -> None:
        self.gates[stage] = "pass" if passed else "fail"
        self.gates[f"{stage}_reasons"] = list(reasons or [])

    def fail(self, stage: str, exc: BaseException | str) -> None:
        """Record a stage failure and carry on with the rest of the batch."""
        self.status = "failed"
        text = (
            f"{stage}: {type(exc).__name__}: {exc}"
            if isinstance(exc, BaseException)
            else f"{stage}: {exc}"
        )
        self.errors.append(text)

    def has_stage(self, stage: str) -> bool:
        """True if `stage` already produced scores or a gate -- skip it."""
        return stage in self.gates or any(
            k.startswith(f"{stage}.") for k in self.scores
        )

    # --- serialisation --------------------------------------------------
    def to_dict(self) -> dict:
        return asdict(self)

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> "DesignRecord":
        known = {f for f in cls.__dataclass_fields__}
        unknown = set(data) - known
        if unknown:
            # Loud: an unknown key means a schema change we have not handled.
            raise ValueError(f"DesignRecord has unknown keys: {sorted(unknown)}")
        return cls(**{k: v for k, v in data.items() if k in known})

    def write(self, design_dir: str | Path) -> Path:
        """Write record.json into the design dir, atomically."""
        d = Path(design_dir)
        d.mkdir(parents=True, exist_ok=True)
        return _atomic_write_json(d / "record.json", self.to_dict())

    @classmethod
    def read(cls, path: str | Path) -> "DesignRecord":
        with open(path) as fh:
            return cls.from_dict(json.load(fh))


def _atomic_write_json(path: Path, obj: Any) -> Path:
    """Write JSON via a temp file + rename, so a kill cannot truncate it.

    Jobs here get preempted (plan 0.3); a half-written record.json would be an
    unresumable run.
    """
    path = Path(path)
    fd, tmp = tempfile.mkstemp(dir=str(path.parent), suffix=".tmp")
    try:
        with os.fdopen(fd, "w") as fh:
            json.dump(obj, fh, indent=2, sort_keys=True)
            fh.flush()
            os.fsync(fh.fileno())
        os.replace(tmp, path)
    except BaseException:
        Path(tmp).unlink(missing_ok=True)
        raise
    return path


# --- manifest ----------------------------------------------------------


def append_manifest(manifest_path: str | Path, record: DesignRecord) -> None:
    """Append one record as a JSONL line.

    Append-mode writes of a single short line are atomic enough for our use,
    and keep the manifest readable while a run is in flight.
    """
    p = Path(manifest_path)
    p.parent.mkdir(parents=True, exist_ok=True)
    with open(p, "a") as fh:
        fh.write(json.dumps(record.to_dict(), sort_keys=True) + "\n")
        fh.flush()
        os.fsync(fh.fileno())


def read_manifest(manifest_path: str | Path) -> Iterator[DesignRecord]:
    """Yield every record in a manifest.

    A corrupt line raises with its line number rather than being skipped --
    silently dropping designs would corrupt every pass-rate denominator.
    """
    p = Path(manifest_path)
    if not p.exists():
        return
    with open(p) as fh:
        for n, line in enumerate(fh, start=1):
            line = line.strip()
            if not line:
                continue
            try:
                yield DesignRecord.from_dict(json.loads(line))
            except Exception as exc:
                raise ValueError(f"{p}:{n}: corrupt manifest line: {exc}") from exc


def latest_by_design_id(manifest_path: str | Path) -> dict[str, DesignRecord]:
    """Collapse an append-only manifest to the newest record per design.

    Later lines win, so a resumed or re-run stage supersedes its earlier entry.
    """
    out: dict[str, DesignRecord] = {}
    for rec in read_manifest(manifest_path):
        out[rec.design_id] = rec
    return out
