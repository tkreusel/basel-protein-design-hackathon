"""Fetch an MSA from the public ColabFold server for one or more sequences.

The organisers confirmed this is allowed: the "use only the provided computing
resources" rule in GUIDE.md is about GPUs, not about an MSA service that could
equally run locally -- and the evaluation itself will generate MSAs this way,
so matching the protocol matters.

Expectation to test, not assume: the organisers' own verified example reports
`"kemp_design": {"depth": 1, "nonquery_rows": 0}` for a de novo design against
the full 183 GB UniRef90. If de novo Kemp designs genuinely have no natural
homologues, a ColabFold MSA will also come back query-only and will change
nothing downstream. This script exists to check that rather than guess.

Politeness and privacy:
  * one request at a time, with backoff; this is a free shared service
  * a generic User-Agent -- no personal or project identifying details
  * sequences are sent as form fields, never interpolated into a URL

Usage:
  python -I colabfold_msa.py <out_dir> <name>:<sequence> [<name>:<seq> ...]
  python -I colabfold_msa.py <out_dir> --fasta <file.fasta> [--limit N]
"""

from __future__ import annotations

import subprocess
import sys
import tarfile
import time
from pathlib import Path

API = "https://api.colabfold.com"
# Generic client string. Deliberately carries no user, project or host detail.
UA = "kemp-design-pipeline/1.0"
MODE = "env"          # UniRef + environmental, ColabFold's default
POLL_SECONDS = 15
MAX_WAIT_SECONDS = 1800


def _curl(args: list[str]) -> tuple[int, str]:
    p = subprocess.run(["curl", "-sS", "-A", UA, *args],
                       capture_output=True, text=True)
    return p.returncode, (p.stdout or p.stderr)


def submit(name: str, seq: str) -> str:
    """Submit one sequence, return the ticket id."""
    fasta = f">{name}\n{seq}\n"
    # The API wants application/x-www-form-urlencoded, NOT multipart. A
    # multipart POST returns HTTP 400 "http: no such file" (the Go server
    # looking for a field it cannot find). --data-urlencode also means the
    # sequence is encoded rather than interpolated raw.
    rc, out = _curl(["-X", "POST",
                     "--data-urlencode", f"q={fasta}",
                     "--data-urlencode", f"mode={MODE}",
                     f"{API}/ticket/msa"])
    if rc != 0:
        raise RuntimeError(f"submit failed (curl rc={rc}): {out[:300]}")
    import json
    try:
        data = json.loads(out)
    except Exception as exc:
        raise RuntimeError(f"submit: unparseable response: {out[:300]}") from exc
    if "id" not in data:
        raise RuntimeError(f"submit: no ticket id in response: {data}")
    return data["id"]


def wait(ticket: str) -> str:
    """Poll until the ticket completes. Returns the final status."""
    import json
    waited = 0
    last = ""
    while waited < MAX_WAIT_SECONDS:
        rc, out = _curl([f"{API}/ticket/{ticket}"])
        if rc == 0:
            try:
                st = json.loads(out).get("status", "")
            except Exception:
                st = ""
            if st and st != last:
                print(f"    status: {st}")
                last = st
            if st == "COMPLETE":
                return st
            if st in ("ERROR", "UNKNOWN", "RATELIMIT"):
                raise RuntimeError(f"ticket {ticket} ended as {st}")
        time.sleep(POLL_SECONDS)
        waited += POLL_SECONDS
    raise TimeoutError(f"ticket {ticket} still pending after {MAX_WAIT_SECONDS}s")


def download(ticket: str, dest: Path) -> Path:
    dest.mkdir(parents=True, exist_ok=True)
    tgz = dest / f"{ticket}.tar.gz"
    rc, out = _curl(["-o", str(tgz), f"{API}/result/download/{ticket}"])
    if rc != 0 or not tgz.exists() or tgz.stat().st_size == 0:
        raise RuntimeError(f"download failed for {ticket}: {out[:300]}")
    # Downloaded archive = untrusted data: extract into its own empty dir.
    extract = dest / ticket
    extract.mkdir(exist_ok=True)
    with tarfile.open(tgz) as tf:
        for m in tf.getmembers():
            if m.isdir() or m.isfile():
                p = Path(m.name)
                if p.is_absolute() or ".." in p.parts:
                    raise RuntimeError(f"unsafe path in archive: {m.name}")
        tf.extractall(extract)
    return extract


def summarise(a3m: Path) -> dict:
    """Depth and non-query row count, matching the organisers' reporting."""
    rows, seq = [], None
    for line in a3m.read_text().splitlines():
        if line.startswith("#"):
            continue
        if line.startswith(">"):
            if seq is not None:
                rows.append(seq)
            seq = ""
        elif seq is not None:
            # Keep only alignment characters. A3M files can carry blank or
            # whitespace-only lines, which otherwise inflate query_length.
            seq += "".join(c for c in line.strip() if c.isalpha() or c in "-.")
    if seq is not None:
        rows.append(seq)
    rows = [r for r in rows if r]
    non_query = sum(1 for r in rows[1:] if any(c.isupper() for c in r))
    return {"depth": len(rows), "nonquery_rows": non_query,
            "query_length": len(rows[0].replace("-", "")) if rows else 0}


def read_fasta(p: Path) -> list[tuple[str, str]]:
    out, name, buf = [], None, []
    for line in p.read_text().splitlines():
        if line.startswith(">"):
            if name:
                out.append((name, "".join(buf)))
            name, buf = line[1:].split()[0], []
        elif line.strip():
            buf.append(line.strip())
    if name:
        out.append((name, "".join(buf)))
    return out


def main(argv: list[str]) -> None:
    out_dir = Path(argv[0])
    jobs: list[tuple[str, str]] = []
    if len(argv) > 1 and argv[1] == "--fasta":
        jobs = read_fasta(Path(argv[2]))
        if "--limit" in argv:
            jobs = jobs[: int(argv[argv.index("--limit") + 1])]
    else:
        for a in argv[1:]:
            n, _, s = a.partition(":")
            jobs.append((n, s))

    out_dir.mkdir(parents=True, exist_ok=True)
    import json
    stats = {}
    for name, seq in jobs:
        print(f"\n[{name}] {len(seq)} aa -> ColabFold ({MODE})")
        try:
            t = submit(name, seq)
            print(f"    ticket: {t}")
            wait(t)
            ex = download(t, out_dir)
            a3ms = sorted(ex.rglob("*.a3m"))
            if not a3ms:
                raise RuntimeError(f"no .a3m in result for {name}")
            target = out_dir / f"{name}.a3m"
            target.write_text(a3ms[0].read_text())
            s = summarise(target)
            stats[name] = s
            print(f"    depth={s['depth']}  nonquery_rows={s['nonquery_rows']}"
                  f"  query_length={s['query_length']}")
            print(f"    saved: {target}")
        except Exception as exc:  # noqa: BLE001 - report and continue
            print(f"    FAILED: {type(exc).__name__}: {exc}")
            stats[name] = {"error": f"{type(exc).__name__}: {exc}"}
        time.sleep(2)

    (out_dir / "msa_stats.json").write_text(json.dumps(stats, indent=2))
    print(f"\nwrote {out_dir / 'msa_stats.json'}")
    print("\nReference, organisers' verified example:")
    print("  designed sequence   depth=1   nonquery_rows=0   (no homologues)")
    print("  natural control     depth=301 nonquery_rows=300")


if __name__ == "__main__":
    main(sys.argv[1:])
