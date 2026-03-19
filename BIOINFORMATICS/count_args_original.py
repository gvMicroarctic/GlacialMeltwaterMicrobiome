#!/usr/bin/env python3

import sys
import csv
from collections import Counter

def main():
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <input_file.tsv> <output_prefix>")
        sys.exit(1)

    input_file = sys.argv[1]
    prefix = sys.argv[2]

    columns = [
        "ARG",
        "AMR_category",
        "AMR_sub_class",
        "Resistance_mechanism"
    ]

    counters = {col: Counter() for col in columns}

    with open(input_file, newline='', encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")

        # ensure required columns exist
        for col in columns:
            if col not in reader.fieldnames:
                print(f"Error: column '{col}' not found in input file.")
                sys.exit(1)

        for row in reader:
            for col in columns:
                val = row[col].strip()
                counters[col].update([val])

    # write output files
    for col in columns:
        out_file = f"{prefix}_{col}_counts.tsv"
        with open(out_file, "w", encoding="utf-8") as out:
            out.write(f"{col}\tCount\n")
            for value, count in counters[col].most_common():
                out.write(f"{value}\t{count}\n")


if __name__ == "__main__":
    main()
