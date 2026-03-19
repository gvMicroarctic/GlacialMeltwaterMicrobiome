#!/usr/bin/env python3

import glob

# Path pattern for all MAG prediction files
pattern = "/storage/varliero/icevirome_prj/MAG/pathofact/PathoFact_MAGs_*/PathoFact_intermediate/Prodigal/bin_*.gff"

# Output file
outfile = "/storage/varliero/icevirome_prj/MAG/pathofact/PathoFact.gff"

files = sorted(glob.glob(pattern))

if not files:
    print("No prediction files found.")
    exit(1)

print(f"Found {len(files)} files. Merging into:")
print(outfile)

with open(outfile, "w") as out:

    for f in files:
        with open(f, "r") as infile:
            for i, line in enumerate(infile):
                out.write(line)

print("Done.")
