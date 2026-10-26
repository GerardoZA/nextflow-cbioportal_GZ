#!/usr/bin/env python3
"""Generate a template cBioPortal data_clinical_sample.txt for a genomic-only study.

Only values the pipeline actually knows are written:
  PATIENT_ID     subject_id from the genomic samplesheet
  SAMPLE_ID      sample_id from the genomic samplesheet
  ONCOTREE_CODE  --cancer-type (the cancer_type param), upper-cased

Every other column is written empty. Nothing is inferred: the empty columns
must be filled in by hand before the study is imported into cBioPortal.

Format: https://docs.cbioportal.org/file-formats/#clinical-data
"""

import argparse
import csv
import re
import sys

# (attribute name, display name, description, datatype, priority)
COLUMNS = [
    ("PATIENT_ID", "Patient Identifier", "Identifier to uniquely specify a patient.", "STRING", "1"),
    ("SAMPLE_ID", "Sample Identifier", "A unique sample identifier.", "STRING", "1"),
    ("ONCOTREE_CODE", "Oncotree Code", "Oncotree Code", "STRING", "1"),
    ("CANCER_TYPE", "Cancer Type", "Cancer Type", "STRING", "1"),
    ("CANCER_TYPE_DETAILED", "Cancer Type Detailed", "Cancer Type Detailed", "STRING", "1"),
    ("SAMPLE_TYPE", "Sample Type", "The type of sample (i.e., normal, primary, met, recurrence).", "STRING", "1"),
]

# cBioPortal: IDs may only contain numbers, letters, points, underscores and hyphens
INVALID_ID_CHARACTERS = re.compile(r"[^A-Za-z0-9._-]")


def read_linking(path):
    """Return (subject_id, sample_id) pairs from the linking TSV, sorted by sample_id."""
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        missing = {"subject_id", "sample_id"} - set(reader.fieldnames or [])
        if missing:
            sys.exit(f"ERROR: {path} is missing column(s): {', '.join(sorted(missing))}")
        pairs = [((row["subject_id"] or "").strip(), (row["sample_id"] or "").strip()) for row in reader]

    if not pairs:
        sys.exit(f"ERROR: no samples found in {path}")

    seen = set()
    for subject, sample in pairs:
        if not subject or not sample:
            sys.exit(f"ERROR: blank subject_id or sample_id in {path}: '{subject}' / '{sample}'")
        if sample in seen:
            sys.exit(f"ERROR: sample_id defined twice in {path}: {sample}")
        seen.add(sample)
        for value in (subject, sample):
            if INVALID_ID_CHARACTERS.search(value):
                print(
                    f"WARNING: '{value}' contains characters cBioPortal does not allow in IDs "
                    "(only letters, numbers, points, underscores and hyphens)",
                    file=sys.stderr,
                )

    return sorted(pairs, key=lambda pair: pair[1])


def main():
    parser = argparse.ArgumentParser(description="Create a template cBioPortal data_clinical_sample.txt.")
    parser.add_argument("--linking", required=True, help="TSV with subject_id and sample_id columns")
    parser.add_argument("--cancer-type", default="", help="OncoTree code from the cancer_type param (may be empty)")
    parser.add_argument("-o", "--output", default="data_clinical_sample.txt", help="Output path (default: data_clinical_sample.txt)")
    args = parser.parse_args()

    pairs = read_linking(args.linking)
    oncotree_code = args.cancer_type.strip().upper()

    with open(args.output, "w", newline="") as fh:
        # Four '#' rows describing the attributes, then the attribute names
        fh.write("#" + "\t".join(col[1] for col in COLUMNS) + "\n")
        fh.write("#" + "\t".join(col[2] for col in COLUMNS) + "\n")
        fh.write("#" + "\t".join(col[3] for col in COLUMNS) + "\n")
        fh.write("#" + "\t".join(col[4] for col in COLUMNS) + "\n")
        fh.write("\t".join(col[0] for col in COLUMNS) + "\n")
        for subject, sample in pairs:
            values = {"PATIENT_ID": subject, "SAMPLE_ID": sample, "ONCOTREE_CODE": oncotree_code}
            fh.write("\t".join(values.get(col[0], "") for col in COLUMNS) + "\n")

    empty = [col[0] for col in COLUMNS if col[0] not in ("PATIENT_ID", "SAMPLE_ID")
             and not (col[0] == "ONCOTREE_CODE" and oncotree_code)]
    print(f"Wrote {len(pairs)} samples to {args.output}", file=sys.stderr)
    print(f"TEMPLATE: columns left empty, to be filled in by hand: {', '.join(empty)}", file=sys.stderr)


if __name__ == "__main__":
    main()
