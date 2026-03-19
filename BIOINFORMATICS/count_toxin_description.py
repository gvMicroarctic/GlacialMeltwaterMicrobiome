#!/usr/bin/env python3

import sys
import csv
from collections import Counter

def main():
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <input_file.tsv> <output_file.tsv>")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    counts = Counter()

    with open(input_file, newline='', encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")

        if "Description" not in reader.fieldnames:
            print("Error: No 'Description' column found in input file.")
            sys.exit(1)

        for row in reader:
            desc = row["Description"].strip()
            counts[desc] += 1

    # Write results to output file
    with open(output_file, "w", newline='', encoding="utf-8") as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow(["Description", "Count"])
        for desc, count in counts.most_common():
            writer.writerow([desc, count])


if __name__ == "__main__":
    main()
