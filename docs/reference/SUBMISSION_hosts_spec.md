# Design submission

Submit one FASTA file to S3. Works from the hackathon clusters or from your own
machine; it requires only your team AWS sign-in.

Bucket: `perlin-pdh-2026-paris` (eu-west-3, Paris)
Key:    `submissions/<YOUR_ACCOUNT_ID>/designs.fasta`

---

## What to submit

One file, named `designs.fasta`.

| Rule | Value |
|---|---|
| Designs per team | at most **20** |
| Sequence length | **50-260** residues |
| Alphabet | the 20 standard amino acids, uppercase |
| Design name | letters, digits, `_`, `-`; unique; max 64 characters |
| Catalytic annotation | **`base=` and `oxy=` on every header** |
| Duplicates | identical sequences are scored once and waste a slot |

```
>my_design_01 base=19 oxy=24
MKAYLVTGGAGFIGSHLVDRLLEQGHEVIVLDNLSTGKRENVPAGATLVEGDLGDAALLD
>my_design_02 base=19 oxy=41
MTAYLITGGAGFIGSHLCERLLADGHDVVVLDNLSTGSRANLDPSASLVVGDIRDGALVR
```

Sequence only. Do not submit the ligand; every design is cofolded against
**6NT**, the same ligand for every team.

## Declaring your catalytic residues

Design is unindexed, so scoring cannot identify your catalytic residues from
the sequence alone. Declare them in the header:

| Field | Meaning | Allowed residues |
|---|---|---|
| `base=` | the catalytic base that deprotonates C3 | **ASP** or **GLU** (`D`, `E`) |
| `oxy=` | the oxyanion-hole donor to the leaving oxygen | **GLN, ASN, SER, THR, TYR, LYS, ARG** (`Q N S T Y K R`) |

- Both are **1-based** indices into the sequence you submitted — not HG3.17
  numbering, not 0-based. The first residue is 1.
- The validator checks that each index exists and that the residue there is of
  the allowed type. A `base=` pointing at an alanine is rejected, not scored:
  an off-by-one measures the distance to an unrelated residue and returns a
  plausible but wrong number.
- The two must be different residues.
- `oxyanion=` is accepted as a synonym for `oxy=`, and spaces around `=` are fine.
- Nothing else may appear in the header. Free text after the name is rejected.
- **HIS cannot be declared as the donor.** It donates from NE2 or ND1 depending on
  tautomer and the header has no way to say which, so a His-base dyad cannot be
  expressed in this format. Design around it.

## Eligibility

Your design must be novel. Two rules, both as announced at the kickoff:

- **Backbones with a TM-score > 0.7 to published structures do not qualify.**
- **Sequences with more than 70% identity to published enzymes do not qualify.**

This is a novelty requirement, not a quality one. HG3.17 itself does not qualify.
Starting from a known Kemp eliminase and mutating it will not produce an eligible
entry — the design has to be far enough away on both sequence and fold.

## What we are looking for

Each design is cofolded with the ligand in **Boltz-2** and ranked against the
other submissions by the same procedure. Four criteria:

1. **Catalytic contact distances.** Your declared base close to the atom it
   deprotonates, your declared donor close to the atom that develops negative
   charge. This carries the most weight of anything here.
2. **Approach geometry.** Distance alone is not enough: each contact is also
   judged on the angle it approaches the ring at and on whether it lies in the
   ring plane. A donor coming in nearly end-on scores poorly even at a perfect
   distance. The angle and the in-plane terms do not average — the worse of the
   two decides.
3. **Interface confidence (iPDE).** The mean predicted distance error over the
   protein-ligand interface of the cofolded complex. Lower is better; minimise
   it.
4. **Active-site packing.** The catalytic base enclosed in a pocket, not
   solvent-exposed.
5. **Pose consistency.** The ligand placed in the same position across repeated
   predictions. Variance between predictions counts against the design.

Geometry — items 1 and 2 together — outweighs everything else combined. Strength
in one term does not compensate for failure in another: the terms multiply.

**There is an optimum, not a direction.** Each contact distance is scored against
a target taken from crystal structures of working Kemp eliminases and from
textbook hydrogen-bond geometry, which agree. Too short is a clash and is
penalised harder than the same error too long.

Designs are folded **with MSAs**, so a design on a natural scaffold is not
handicapped relative to a de novo one.

## Get the two tools

Both are in the shared bucket, readable by every team account. Download them
into the same directory; `submit-designs.sh` expects the validator beside
itself:

```
cd /mnt/work/projects/<your project>        # on the cluster: /mnt is read-only
aws s3 cp s3://perlin-pdh-2026-paris/release/submission/validate-submission.py . --region eu-west-3
aws s3 cp s3://perlin-pdh-2026-paris/release/submission/submit-designs.sh      . --region eu-west-3
chmod +x submit-designs.sh
```

Off the cluster, any writable directory works. If you upload by hand, you need
only `validate-submission.py`.

## Check it before you upload

```
python3 validate-submission.py designs.fasta
```

Python 3 only: no pip install, no AWS, no network. It lists every problem and
exits non-zero if the file would be rejected. **Run it early** — a file that
fails validation at the deadline cannot be scored.

## Upload

**From anywhere:**

```
./submit-designs.sh designs.fasta
```

Validates, resolves your account id, and refuses to upload a file that would be
rejected.

**By hand with the CLI:**

```
aws s3 cp designs.fasta \
  s3://perlin-pdh-2026-paris/submissions/$(aws sts get-caller-identity --query Account --output text)/designs.fasta \
  --region eu-west-3
```

**From the browser:** sign in to your team AWS account, open your folder using
the direct link in your team packet (`CONNECT.txt`), then use *Upload*. Opening the bucket
from the S3 home page shows *Access denied*; this is expected, as you can see
only your own folder.

## Rules of the road

- You may re-upload as often as you like. **The file present at the deadline is
  the one scored**; each upload replaces the previous one.
- You cannot see, read or overwrite another team's submission, and they cannot
  see yours.
- After the deadline, S3 refuses uploads. There is no grace period.
