#!/usr/bin/env python3
"""Validate a hackathon design submission before it is uploaded.

Pure Python 3 standard library: no pip install, no AWS, no network.
Runs anywhere, including on a team's own laptop or cluster.

    python3 validate-submission.py designs.fasta

Exit code 0 means the file is accepted by the scoring pipeline.
Any other exit code means it would be rejected; every problem is listed.
"""

import argparse
import re
import sys

# --- Submission rules. Organisers: these four numbers are the whole policy. ---
# MAX_DESIGNS is settled at 20 and is the dial that makes judging fit its
# window: 12 teams x 20 designs is the cofold budget the scoring hardware was
# sized against. Raising it does not make scoring a little slower, it overruns
# the window, so change it only together with the scoring plan.
MAX_DESIGNS = 20
MIN_LENGTH = 50
MAX_LENGTH = 260          # competition limit; the design tooling enforces it too
MAX_ID_LENGTH = 64

AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")
ID_PATTERN = re.compile(r"^[A-Za-z0-9_-]+$")
ATTR_PATTERN = re.compile(r"^([A-Za-z_]+)=(.*)$")

# Catalytic annotation. The design is unindexed, so scoring cannot find these
# residues by alignment; the team has to declare them.
REQUIRED_ATTRS = ("base", "oxy")
# Catalytic base: the carboxylate that deprotonates C3 of the benzisoxazole.
BASE_RESIDUES = {"D": "ASP", "E": "GLU"}
# Oxyanion-hole donor: H-bond donor to the developing phenolate oxygen.
# HIS is excluded on purpose: it donates from NE2 or ND1 depending on tautomer and a
# two-field header cannot say which, so a His-base dyad cannot be expressed here.
OXY_RESIDUES = {
    "Q": "GLN", "N": "ASN", "S": "SER", "T": "THR", "Y": "TYR",
    "K": "LYS", "R": "ARG",
}
ATTR_SPEC = {
    "base": ("catalytic base", BASE_RESIDUES),
    "oxy": ("oxyanion-hole donor", OXY_RESIDUES),
}
# The scorer's header regex accepts 'oxyanion' as a synonym for 'oxy'; mirror that here
# so a file this script passes is a file the scorer parses, and vice versa.
ATTR_ALIASES = {"oxyanion": "oxy"}

# Ambiguity codes the scoring pipeline cannot cofold; named so the error can
# explain itself rather than just refusing the character.
AMBIGUOUS = {
    "X": "unknown residue",
    "B": "asparagine/aspartate ambiguity",
    "Z": "glutamine/glutamate ambiguity",
    "J": "leucine/isoleucine ambiguity",
    "U": "selenocysteine",
    "O": "pyrrolysine",
}


def parse_fasta(text):
    """Yield (line_number, header, sequence). Tolerates CRLF and blank lines."""
    records = []
    header = None
    start_line = 0
    chunks = []
    for number, raw in enumerate(text.splitlines(), start=1):
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                records.append((start_line, header, "".join(chunks)))
            header = line[1:].strip()
            start_line = number
            chunks = []
        else:
            if header is None:
                raise ValueError(
                    "line %d: sequence data appears before the first '>' header" % number
                )
            chunks.append(line)
    if header is not None:
        records.append((start_line, header, "".join(chunks)))
    return records


def split_header(header):
    """Split 'name base=12 oxy=50' into (name, {attr: raw_value}, problems)."""
    problems = []
    # The scorer's regex tolerates spaces around '='; collapse them so 'base = 12'
    # parses here too rather than being reported as stray free text.
    header = re.sub(r"\s*=\s*", "=", header)
    fields = header.split()
    name = fields[0] if fields else ""
    attrs = {}
    for field in fields[1:]:
        match = ATTR_PATTERN.match(field)
        if not match:
            problems.append(
                "'%s' is not a 'key=value' annotation. The header is the design name "
                "followed by 'base=' and 'oxy=' only, e.g. "
                "'>my_design_01 base=127 oxy=50'." % field
            )
            continue
        key, value = match.group(1).lower(), match.group(2)
        key = ATTR_ALIASES.get(key, key)
        if key not in ATTR_SPEC:
            problems.append(
                "unknown annotation '%s='. Expected 'base=' and 'oxy=', exactly those "
                "spellings." % match.group(1)
            )
        elif key in attrs:
            problems.append("'%s=' is given more than once." % key)
        else:
            attrs[key] = value
    return name, attrs, problems


def check_catalytic(attrs, sequence, where, problems):
    """Validate base= and oxy= against the submitted sequence."""
    parsed = {}
    for key in REQUIRED_ATTRS:
        label, allowed = ATTR_SPEC[key]
        if key not in attrs:
            problems.append(
                "%s: the header is missing '%s=', the 1-based index of the %s. Every "
                "design must declare it, e.g. '>my_design_01 base=127 oxy=50'."
                % (where, key, label)
            )
            continue

        raw = attrs[key]
        try:
            index = int(raw)
        except ValueError:
            problems.append(
                "%s: '%s=%s' is not a whole number. Give the 1-based position of the "
                "%s in this sequence." % (where, key, raw, label)
            )
            continue

        parsed[key] = index

        if not sequence:
            continue

        if not 1 <= index <= len(sequence):
            problems.append(
                "%s: '%s=%d' is outside this sequence, which is %d residues long. "
                "Indices are 1-based and count into the sequence you submitted, not "
                "HG3.17 numbering." % (where, key, index, len(sequence))
            )
            continue

        residue = sequence[index - 1]
        if residue not in allowed:
            problems.append(
                "%s: '%s=%d' points at '%s', but the %s must be one of %s. Check for "
                "an off-by-one: a wrong index silently measures the distance to an "
                "unrelated residue."
                % (
                    where,
                    key,
                    index,
                    residue,
                    label,
                    ", ".join(sorted(allowed.values())),
                )
            )

    # The scorer requires the two positions to be different residues.
    if parsed.get("base") is not None and parsed.get("base") == parsed.get("oxy"):
        problems.append(
            "%s: 'base' and 'oxy' both point at residue %d. They must be two different "
            "residues: one donates the proton, the other stabilises the leaving oxygen."
            % (where, parsed["base"])
        )


def validate(text):
    """Return a list of human-readable problems. Empty list means valid."""
    problems = []

    try:
        records = parse_fasta(text)
    except ValueError as exc:
        return ["%s" % exc, "The file must be FASTA: a '>' header line, then the sequence."]

    if not records:
        return ["The file contains no FASTA records (no line starting with '>')."]

    if len(records) > MAX_DESIGNS:
        problems.append(
            "%d designs found, but the limit is %d. Submit your best %d."
            % (len(records), MAX_DESIGNS, MAX_DESIGNS)
        )

    seen_ids = {}
    seen_sequences = {}

    for line_number, header, sequence in records:
        design_id, attrs, header_problems = split_header(header)
        where = "line %d (>%s)" % (line_number, design_id or "<empty>")

        for item in header_problems:
            problems.append("%s: %s" % (where, item))

        if not design_id:
            problems.append(
                "%s: the header is empty; use '>my_design_01 base=127 oxy=50'." % where
            )
        elif len(design_id) > MAX_ID_LENGTH:
            problems.append(
                "%s: the design name is %d characters; the limit is %d."
                % (where, len(design_id), MAX_ID_LENGTH)
            )
        elif not ID_PATTERN.match(design_id):
            problems.append(
                "%s: the design name may only contain letters, digits, '_' and '-'. "
                "Only the 'base=' and 'oxy=' annotations may follow it." % where
            )

        if design_id in seen_ids:
            problems.append(
                "%s: duplicate design name, already used on line %d. Names must be unique."
                % (where, seen_ids[design_id])
            )
        else:
            seen_ids[design_id] = line_number

        if not sequence:
            problems.append("%s: no sequence follows this header." % where)
            check_catalytic(attrs, "", where, problems)
            continue

        upper = sequence.upper()
        if upper != sequence:
            problems.append("%s: sequence must be uppercase." % where)

        bad = sorted(set(upper) - AMINO_ACIDS)
        if bad:
            detail = ", ".join(
                "'%s' (%s)" % (ch, AMBIGUOUS[ch]) if ch in AMBIGUOUS else "'%s'" % ch
                for ch in bad
            )
            problems.append(
                "%s: sequence contains characters that are not one of the 20 standard "
                "amino acids: %s. The scoring pipeline cannot cofold these." % (where, detail)
            )

        if not (MIN_LENGTH <= len(upper) <= MAX_LENGTH):
            problems.append(
                "%s: sequence is %d residues; it must be between %d and %d."
                % (where, len(upper), MIN_LENGTH, MAX_LENGTH)
            )

        check_catalytic(attrs, upper, where, problems)

        if upper in seen_sequences:
            problems.append(
                "%s: this sequence is identical to '%s'. Duplicate designs are scored once "
                "and waste one of your %d slots."
                % (where, seen_sequences[upper], MAX_DESIGNS)
            )
        else:
            seen_sequences[upper] = design_id

    return problems


def main():
    parser = argparse.ArgumentParser(
        description="Validate a design submission FASTA before uploading it."
    )
    parser.add_argument("fasta", help="path to your designs.fasta")
    parser.add_argument(
        "--quiet", action="store_true", help="print nothing; use the exit code only"
    )
    args = parser.parse_args()

    try:
        with open(args.fasta, "r", encoding="utf-8", errors="replace") as handle:
            text = handle.read()
    except OSError as exc:
        if not args.quiet:
            print("FAILED: cannot read %s: %s" % (args.fasta, exc))
        return 2

    problems = validate(text)
    records = []
    try:
        records = parse_fasta(text)
    except ValueError:
        pass

    if problems:
        if not args.quiet:
            print("REJECTED: %s" % args.fasta)
            print("%d problem(s) must be fixed before the deadline:\n" % len(problems))
            for item in problems:
                print("  - %s" % item)
            print("\nFix these and run this script again.")
        return 1

    if not args.quiet:
        lengths = [len(sequence) for _, _, sequence in records]
        print("ACCEPTED: %s" % args.fasta)
        print(
            "  %d design(s), %d-%d residues (limit %d designs, %d-%d residues)"
            % (
                len(records),
                min(lengths),
                max(lengths),
                MAX_DESIGNS,
                MIN_LENGTH,
                MAX_LENGTH,
            )
        )
        print("  base= and oxy= declared and checked on every design.")
        print("  This file will be accepted by the scoring pipeline.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
