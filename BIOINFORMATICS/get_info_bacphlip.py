#!/usr/bin/env python3

import argparse

def main():
    parser = argparse.ArgumentParser(
        description="Parse Bacphlip output and classify viral lifestyle"
    )
    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Bacphlip output file (.bacphlip)"
    )
    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Output file with lifestyle classification"
    )
    parser.add_argument(
        "-t", "--threshold",
        type=float,
        default=0.80,
        help="Probability threshold for classification (default: 0.80)"
    )

    args = parser.parse_args()

    count_virulent = 0
    count_temperate = 0
    count_none = 0

    with open(args.output, "w") as out_file, open(args.input, "r") as in_file:
        in_file.readline()  # skip header
        for line in in_file:
            line = line.strip()
            info = line.split("\t")

            if float(info[1]) >= args.threshold:  # virulent
                out_file.write(f"{info[0]}\tvirulent\n")
                count_virulent += 1
            elif float(info[2]) >= args.threshold:  # temperate
                out_file.write(f"{info[0]}\ttemperate\n")
                count_temperate += 1
            else:
                out_file.write(f"{info[0]}\tnone\n")
                count_none += 1

    print(f"Virulent viruses are {count_virulent}.")
    print(f"Temperate viruses are {count_temperate}.")
    print(f"Unclassified viruses are {count_none}.")


if __name__ == "__main__":
    main()
