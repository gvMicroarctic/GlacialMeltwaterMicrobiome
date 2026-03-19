#!/usr/bin/env python3

import pandas as pd
import argparse
import re

parser = argparse.ArgumentParser(description="Filter iPHoP host predictions")
parser.add_argument("-i", "--input", required=True, help="Input CSV file")
parser.add_argument("-o", "--output", required=True, help="Output CSV file")

args = parser.parse_args()

df = pd.read_csv(args.input)

# Convert column to numeric
df["Confidence score"] = pd.to_numeric(df["Confidence score"], errors="coerce")

# Filter rows
df = df[df["Confidence score"] >= 90]

def filter_methods(method_string):
    if pd.isna(method_string):
        return method_string

    pairs = re.findall(r'([A-Za-z\-]+);([0-9.]+)', method_string)

    kept = []
    for method, score in pairs:
        if float(score) >= 90:
            kept.append(f"{method};{score}")

    return " ".join(kept) if kept else None

df["List of methods"] = df["List of methods"].apply(filter_methods)

df.to_csv(args.output, index=False)
