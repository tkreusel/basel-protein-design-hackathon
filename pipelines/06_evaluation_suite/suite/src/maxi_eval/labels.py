"""Activity labels for the evaluated proteins.

Sources (later wins): `entries:`/`extra_entries:` in config/entries.yaml, then config/labels.csv.
labels.csv columns:  entry, kcat_KM  [, kcat, KM, class, note]
  entry   protein id WITHOUT the `_boltz` suffix (e.g. 7K4Q, HG3_17, design_042); matched case-insensitively
  kcat_KM catalytic efficiency in M^-1 s^-1 (blank / NA if unknown or inactive-below-detection; then give `class`)
  class   optional free label, e.g. active / inactive, used for AUROC when kcat_KM is missing
"""
from __future__ import annotations
import csv
from pathlib import Path
import yaml


def load_labels(config_yaml: str = "config/entries.yaml", labels_csv: str = "config/labels.csv") -> dict[str, dict]:
    out: dict[str, dict] = {}
    cfg = yaml.safe_load(Path(config_yaml).read_text())
    for src in (cfg.get("entries") or {}, cfg.get("extra_entries") or {}):
        for k, v in src.items():
            if v.get("kcat_KM") is not None:
                out[k.lower()] = {"kcat_KM": float(v["kcat_KM"]), "class": v.get("class"), "source": v.get("activity_source")}
    p = Path(labels_csv)
    if p.exists():
        for row in csv.DictReader(p.open()):
            if not row["entry"] or row["entry"].lstrip().startswith("#"):
                continue
            k = row["entry"].strip().lower()
            kc = (row.get("kcat_KM") or "").strip()
            out[k] = {"kcat_KM": float(kc) if kc and kc.upper() != "NA" else None,
                      "class": (row.get("class") or "").strip() or None, "source": row.get("note") or "labels.csv"}
    return out


def label_for(entry: str, labels: dict) -> dict | None:
    key = entry.lower()
    for suffix in ("_boltz",):
        if key.endswith(suffix):
            key = key[: -len(suffix)]
    return labels.get(key)
