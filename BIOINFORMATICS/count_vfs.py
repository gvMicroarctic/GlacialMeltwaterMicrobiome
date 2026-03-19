#!/usr/bin/env python3

import sys
import csv
from collections import Counter

def main():
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} <input_file.tsv>")
        sys.exit(1)

    input_file = sys.argv[1]

    pathogenic_count = 0
    confidence_counts = Counter()

    with open(input_file, newline='', encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")

        # Check columns exist
        for col in ["Virulence_prediction", "Virulence_confidence_level"]:
            if col not in reader.fieldnames:
                print(f"Error: column '{col}' not found in input file.")
                sys.exit(1)

        # Count values
        for row in reader:
            if row["Virulence_prediction"] == "pathogenic":
                pathogenic_count += 1

            level = row["Virulence_confidence_level"].strip()
            confidence_counts[level] += 1

    # ---- print results to STDOUT ----
    print(f"Virulence_prediction == pathogenic\t{pathogenic_count}")
    print("\nVirulence_confidence_level counts:")
    print("Level\tCount")

    for level, count in confidence_counts.most_common():
        print(f"{level}\t{count}")


if __name__ == "__main__":
    main()
