#!/usr/bin/env python3

import sys

def read_fasta(path):
    seqs = {}
    header = None
    seq_lines = []

    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header:
                    seqs[header] = "".join(seq_lines)
                header = line[1:].split()[0]  # up to first whitespace
                seq_lines = []
            else:
                seq_lines.append(line)

        if header:
            seqs[header] = "".join(seq_lines)

    return seqs


def write_fasta(seqs, out_file):
    with open(out_file, "w") as out:
        for h, s in seqs.items():
            out.write(f">{h}\n")
            for i in range(0, len(s), 60):
                out.write(s[i:i+60] + "\n")


def main():
    if len(sys.argv) != 4:
        sys.exit(
            "Usage:\n"
            "  python compare_fasta_headers.py file1.fasta file2.fasta output_prefix"
        )

    f1, f2, prefix = sys.argv[1:]

    seqs1 = read_fasta(f1)
    seqs2 = read_fasta(f2)

    headers1 = set(seqs1)
    headers2 = set(seqs2)

    only_1 = {h: seqs1[h] for h in headers1 - headers2}
    only_2 = {h: seqs2[h] for h in headers2 - headers1}
    both   = {h: seqs1[h] for h in headers1 & headers2}

    write_fasta(only_1, f"{prefix}_one.fasta")
    write_fasta(only_2, f"{prefix}_two.fasta")
    write_fasta(both,   f"{prefix}_both.fasta")

    print(f"Only file1: {len(only_1)} sequences → {prefix}_one.fasta")
    print(f"Only file2: {len(only_2)} sequences → {prefix}_two.fasta")
    print(f"Both files: {len(both)} sequences → {prefix}_both.fasta")


if __name__ == "__main__":
    main()
