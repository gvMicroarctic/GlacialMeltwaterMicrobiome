#!/usr/bin/env python3
"""
plot_transposons_combined_cat_color.py

Plot all sequences from a GFF (with optional FASTA) in one figure.
- Taller gene boxes
- Larger labels
- Colors based on GFF attribute (cat=)
- Tracks ordered by MAG ID
- Titles show MAG, contig, coordinates, and taxonomy

Requires:
pip install biopython dna-features-viewer matplotlib bcbio-gff
"""

import os
import sys
import re
from BCBio import GFF
from Bio import SeqIO
from dna_features_viewer import GraphicFeature, GraphicRecord
import matplotlib.pyplot as plt

plt.rcParams["svg.fonttype"] = "none"

# -------------------------
# Default file paths
# -------------------------
GFF_FILE = "filtered_genes_megaclusters_fasta_eggnog_description_cat.gff"
FASTA_FILE = "transposon.fasta"
META_FILE = "specifics_corrected_drep_new.txt"

# -------------------------
# Load bin → MAG + taxonomy
# -------------------------
def load_bin_metadata(meta_file):
    meta = {}
    with open(meta_file) as fh:
        header = fh.readline().rstrip().split("\t")
        idx = {name: i for i, name in enumerate(header)}

        def get(parts, key):
            i = idx.get(key)
            if i is None or i >= len(parts) or parts[i] == "":
                return "NA"
            return parts[i]

        for line in fh:
            if not line.strip():
                continue

            parts = line.rstrip().split("\t")
            bin_id = parts[idx["File"]]

            meta[bin_id] = {
                "mag": get(parts, "Mag"),
                "phylum": get(parts, "phylum"),
                "class": get(parts, "class"),
                "family": get(parts, "family"),
            }

    return meta


# -------------------------
# Load GFF records
# -------------------------
def load_records(gff_path, fasta_path=None):
    if fasta_path and os.path.exists(fasta_path):
        seq_dict = SeqIO.to_dict(SeqIO.parse(fasta_path, "fasta"))
        with open(gff_path) as gff_handle:
            records = list(GFF.parse(gff_handle, base_dict=seq_dict))
    else:
        with open(gff_path) as gff_handle:
            records = list(GFF.parse(gff_handle))

    if not records:
        raise ValueError("No records parsed from GFF.")
    return records

# -------------------------
# Feature helpers
# -------------------------
def guess_label(feature):
    q = feature.qualifiers
    for k in ("gene", "Name", "product", "ID"):
        if k in q and q[k]:
            if q[k][0] == "Unclassified":
                return None
            return q[k][0]
    return feature.type

def pick_color(feature):
    cat = feature.qualifiers.get("cat", ["unknown"])[0]

    if "transposase" in cat:
        return "red"
        
    # Information Storage & Processing
    elif "K" in cat:
        return "yellow"
    elif "L" in cat:
        return "orange"
        
    # Cellular Processes & Signaling
    elif "M" in cat:
        return "cyan"
    elif "T" in cat:
        return "blue"
        
    # Metabolism
    elif "P" in cat:
        return "green"
    elif "C" in cat:
        return "lime"
    elif "H" in cat:
        return "yellowgreen"
    elif "S" in cat or "-" in cat or "Unclassified" in cat:
        return "white"
    else:
        return "#cccccc"

# -------------------------
# Convert features for plotting
# -------------------------
def features_from_record(record, pad=50, box_height=1.2):
    features = []
    starts = [int(f.location.start) for f in record.features]
    ends = [int(f.location.end) for f in record.features]

    if not starts or not ends:
        return [], 0

    min_c, max_c = min(starts), max(ends)

    for f in record.features:
        try:
            start = int(f.location.start) - min_c + pad
            end = int(f.location.end) - min_c + pad
        except Exception:
            continue

        features.append(
            GraphicFeature(
                start=start,
                end=end,
                strand=f.location.strand if f.location.strand in (1, -1, 0) else 0,
                color=pick_color(f),
                label=guess_label(f),
                height=box_height
            )
        )

    return features, max_c - min_c + 2 * pad

# -------------------------
# Record start–end (1-based)
# -------------------------
def record_start_end(record):
    starts = [int(f.location.start) for f in record.features]
    ends = [int(f.location.end) for f in record.features]

    if not starts or not ends:
        return None, None

    return min(starts) + 1, max(ends)

# -------------------------
# MAG sorting key
# -------------------------
def mag_sort_key(record, meta):
    info = record.id.split("_")[:-2]
    bin_id = info[0] + "_" + info[1]

    if bin_id not in meta:
        return (float("inf"), bin_id)

    mag = meta[bin_id]["mag"]
    m = re.search(r"(\d+)", mag)
    mag_num = int(m.group(1)) if m else float("inf")

    return (mag_num, mag)

# -------------------------
# Main
# -------------------------
def main():
    gff_path = sys.argv[1] if len(sys.argv) > 1 else GFF_FILE
    fasta_path = sys.argv[2] if len(sys.argv) > 2 else (
        FASTA_FILE if os.path.exists(FASTA_FILE) else None
    )

    records = load_records(gff_path, fasta_path)
    meta = load_bin_metadata(META_FILE)

    records.sort(key=lambda r: mag_sort_key(r, meta))

    fig, axes = plt.subplots(
        len(records), 1,
        figsize=(20, 4 * len(records)),
        squeeze=False
    )

    for i, rec in enumerate(records):
        features, seq_len = features_from_record(rec)
        if seq_len == 0:
            continue

        gr = GraphicRecord(sequence_length=seq_len, features=features)
        ax = axes[i, 0]
        gr.plot(ax=ax, figure_width=20)

        for t in ax.texts:
            t.set_fontsize(14)
            t.set_weight("bold")

        info = rec.id.split("_")[:-2]
        bin_id = info[0] + "_" + info[1]
        contig_id = info[2] + "_" + info[3]

        md = meta.get(bin_id, {})
        mag = md.get("mag", bin_id)
        phylum = md.get("phylum", "NA")
        class_ = md.get("class", "NA")
        family = md.get("family", "NA")

        start, end = record_start_end(rec)

        if start is not None:
            title = (
                f"{mag} - {contig_id} - "
                f"{start:,}–{end:,} bp - "
                f"{phylum}, {class_}, {family}"
            )
        else:
            title = f"{mag} - {contig_id} - {phylum}, {class_}, {family}"

        ax.set_title(title, fontsize=15, weight="bold")

    plt.tight_layout()
    plt.savefig("all_transposons_combined_cat_color.png", dpi=300, bbox_inches="tight")
    plt.savefig("all_transposons_combined_cat_color.svg", bbox_inches="tight")

    print("Saved figure with requested title format.")

if __name__ == "__main__":
    main()
