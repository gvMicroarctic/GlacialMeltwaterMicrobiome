#!/usr/bin/env python3

import pandas as pd
import re
import argparse

parser = argparse.ArgumentParser(description="Filter iPHoP host predictions and add method code")
parser.add_argument("-i", "--input", required=True, help="Input CSV file")
parser.add_argument("-o", "--output", required=True, help="Output CSV file")
args = parser.parse_args()

# Read CSV
df = pd.read_csv(args.input, skipinitialspace=True)

# Remove repeated headers inside file
df = df[df["Virus"] != "Virus"]

# Convert Confidence score to numeric
df["Confidence score"] = pd.to_numeric(df["Confidence score"], errors="coerce")

# Keep only rows with confidence >= 90
df = df[df["Confidence score"] >= 90]

# Filter methods column (remove methods with score < 90)
def filter_methods(method_string):
    if pd.isna(method_string):
        return None
    pairs = re.findall(r'([A-Za-z\-]+);([0-9.]+)', method_string)
    kept = [f"{m};{s}" for m, s in pairs if float(s) >= 90]
    return " ".join(kept) if kept else None

df["List of methods"] = df["List of methods"].apply(filter_methods)

# Assign Method_code
def assign_method_code(method_string):
    if pd.isna(method_string):
        return None
    
    methods = [m.split(";")[0] for m in method_string.split()]
    unique_methods = set(methods)
    
    if unique_methods == {"CRISPR"}:
        return 2
    elif unique_methods == {"blast"}:
        return 3
    elif unique_methods == {"RaFAH"}:
        return 4
    elif "iPHoP-RF" in unique_methods or len(unique_methods) > 1:
        return 1
    else:
        return 1  # fallback for any other cases

df["Method_code"] = df["List of methods"].apply(assign_method_code)

# Write output
df.to_csv(args.output, index=False, header=True)
