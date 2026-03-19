#!/usr/bin/env python3

import glob

# Path pattern for all MAG prediction files
pattern = "/storage/varliero/icevirome_prj/MAG/pathofact/PathoFact_MAGs_*/PathoFact_report/PathoFact_bin_*_predictions.tsv"

# Output file
outfile = "/storage/varliero/icevirome_prj/MAG/pathofact/PathoFact_predictions.tsv"

files = sorted(glob.glob(pattern))

if not files:
    print("No prediction files found.")
    exit(1)

print(f"Found {len(files)} files. Merging into:")
print(outfile)

with open(outfile, "w") as out:
    header_written = False

    for f in files:
        with open(f, "r") as infile:
            for i, line in enumerate(infile):
                # Write header only once
                if i == 0:
                    if not header_written:
                        out.write(line)
                        header_written = True
                else:
                    out.write(line)

print("Done.")
