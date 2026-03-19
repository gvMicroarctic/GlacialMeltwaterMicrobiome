#!/usr/bin/env python3

def parse_attributes(attr_string):
    """Parse GFF attribute column into a dictionary"""
    attrs = {}
    for item in attr_string.strip().split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            attrs[key] = value
    return attrs


def parse_dram_gff(gff_path):
    """
    Parse DRAM genes.gff
    Returns dict:
    {
        protein_id: {
            contig,
            start,
            end,
            new_protein_name
        }
    }
    """
    proteins = {}

    with open(gff_path) as f:
        for line in f:
            if line.startswith("#"):
                continue

            cols = line.rstrip().split("\t")
            contig_full, _, feature, start, end, *_ , attributes = cols

            if feature != "CDS":
                continue

            # contig formatting: bin_112_NODE_1_bin_112 -> NODE_1_bin_112
            contig = "_".join(contig_full.split("_")[2:])

            attr_dict = parse_attributes(attributes)
            protein_id = attr_dict["ID"]

            new_name = f"{contig}_{start}_{end}"

            proteins[protein_id] = {
                "contig": contig,
                "start": int(start),
                "end": int(end),
                "new_protein_name": new_name
            }

    return proteins

def parse_pathofact_gff(gff_path):
    """
    Parse PathoFact.gff
    Returns dict:
    {
        protein_id: {
            contig,
            start,
            end,
            new_protein_name
        }
    }
    """
    proteins = {}

    with open(gff_path) as f:
        for line in f:
            if line.startswith("#"):
                continue

            cols = line.rstrip().split("\t")
            contig, _, feature, start, end, *_ , attributes = cols

            if feature != "CDS":
                continue

            attr_dict = parse_attributes(attributes)
            protein_id = contig + "_" + attr_dict["ID"]

            new_name = f"{contig}_{start}_{end}"

            proteins[protein_id] = {
                "contig": contig,
                "start": int(start),
                "end": int(end),
                "new_protein_name": new_name
            }

    return proteins

# ----------- Run parsing -----------

dram_gff = "/storage/varliero/icevirome_prj/MAG/dram_drep/genes.gff"
pathofact_gff = "/storage/varliero/icevirome_prj/MAG/pathofact/PathoFact.gff"

dram_dict = parse_dram_gff(dram_gff)
pathofact_dict = parse_pathofact_gff(pathofact_gff)

# ----------- Load AMRFinder results with annotation -----------

def load_amrfinder_amr_ids(amrfinder_path, dram_dict):
    """
    Read AMRFinder output and return dict:
    {
        protein_id: {
            'contig', 'start', 'end', 'new_protein_name',
            'Class', 'Subclass'
        }
    }
    Only where Type == 'AMR'
    """
    matched = {}

    with open(amrfinder_path) as f:
        header = f.readline().rstrip().split("\t")

        type_idx = header.index("Type")
        class_idx = header.index("Class")
        subclass_idx = header.index("Subclass")

        for line in f:
            if not line.strip():
                continue

            fields = line.rstrip().split("\t")
            if fields[type_idx] != "AMR":
                continue

            protein_id = fields[0]

            if protein_id in dram_dict:
                # Add AMRFinder annotation into DRAM dict
                info = dram_dict[protein_id].copy()
                info["Class"] = fields[class_idx]
                info["Subclass"] = fields[subclass_idx]
                matched[protein_id] = info

    return matched

amrfinder_file = "/storage/varliero/icevirome_prj/MAG/amrfinder.txt"
dram_dict_arg = load_amrfinder_amr_ids(amrfinder_file, dram_dict)

print(f"Total DRAM proteins: {len(dram_dict)}")
print(f"AMR proteins matched in DRAM: {len(dram_dict_arg)}")

# ----------- Load PathoFact ARG predictions -----------

def load_pathofact_arg_predictions(pred_file, pathofact_dict):
    """
    Load PathoFact ARG predictions.
    Only keep proteins with ARG != "-"
    Returns dictionary of matched proteins with original info from pathofact_dict
    Adds ARG annotations:
        ARG, AMR_category, AMR_sub_class, Resistance_mechanism
    """
    pathofact_arg_dict = {}

    with open(pred_file) as f:
        header = f.readline().rstrip().split("\t")
        orf_idx = header.index("ORF")
        contig_idx = header.index("Contig")
        arg_idx = header.index("ARG")
        amr_cat_idx = header.index("AMR_category")
        amr_sub_idx = header.index("AMR_sub_class")
        res_mech_idx = header.index("Resistance_mechanism")

        for line in f:
            if not line.strip():
                continue

            fields = line.rstrip().split("\t")
            if fields[arg_idx] == "-":
                continue

            orf = fields[orf_idx]
            contig = fields[contig_idx]
            protein_id = f"{contig}_{orf}"

            if protein_id in pathofact_dict:
                info = pathofact_dict[protein_id].copy()
                # Add ARG annotations
                info["ARG"] = fields[arg_idx]
                info["AMR_category"] = fields[amr_cat_idx]
                info["AMR_sub_class"] = fields[amr_sub_idx]
                info["Resistance_mechanism"] = fields[res_mech_idx]
                pathofact_arg_dict[protein_id] = info

    return pathofact_arg_dict


pathofact_pred_file = "/storage/varliero/icevirome_prj/MAG/pathofact/PathoFact_predictions.tsv"
pathofact_dict_arg = load_pathofact_arg_predictions(pathofact_pred_file, pathofact_dict)

print(f"Total PathoFact proteins: {len(pathofact_dict)}")
print(f"PathoFact ARG proteins matched: {len(pathofact_dict_arg)}")

# ----------- Find common proteins based on new_protein_name -----------

dram_names = set(info["new_protein_name"] for info in dram_dict_arg.values())
pathofact_names = set(info["new_protein_name"] for info in pathofact_dict_arg.values())

common_names = dram_names.intersection(pathofact_names)

print(f"Number of common proteins: {len(common_names)}")
print("Common proteins:", list(common_names))


# ----------- Merge DRAM and PathoFact ARGs based on new_protein_name -----------

output_file = "/storage/varliero/icevirome_prj/MAG/pathofact/merged_amrfinder_pathofact_ARGs.tsv"

# Define all possible fields for output
output_fields = [
    "new_protein_name",
    "AMRFinder_protein_id",
    "AMRFinder_contig",
    "AMRFinder_start",
    "AMRFinder_end",
    "AMRFinder_Class",
    "AMRFinder_Subclass",
    "PathoFact_protein_id",
    "PathoFact_contig",
    "PathoFact_start",
    "PathoFact_end",
    "PathoFact_ARG",
    "PathoFact_AMR_category",
    "PathoFact_AMR_sub_class",
    "PathoFact_Resistance_mechanism"
]

merged_rows = []

# Create mapping: new_protein_name -> DRAM info / PathoFact info
dram_name_to_id = {info["new_protein_name"]: pid for pid, info in dram_dict_arg.items()}
patho_name_to_id = {info["new_protein_name"]: pid for pid, info in pathofact_dict_arg.items()}

all_names = set(dram_name_to_id.keys()).union(patho_name_to_id.keys())

for name in all_names:
    # DRAM info
    if name in dram_name_to_id:
        dram_id = dram_name_to_id[name]
        dram_info = dram_dict_arg[dram_id]
        dram_row = [
            name,
            dram_id,
            dram_info.get("contig", "NA"),
            dram_info.get("start", "NA"),
            dram_info.get("end", "NA"),
            dram_info.get("Class", "NA"),
            dram_info.get("Subclass", "NA")
        ]
    else:
        dram_row = [name] + ["NA"]*6

    # PathoFact info
    if name in patho_name_to_id:
        patho_id = patho_name_to_id[name]
        patho_info = pathofact_dict_arg[patho_id]
        patho_row = [
            patho_id,
            patho_info.get("contig", "NA"),
            patho_info.get("start", "NA"),
            patho_info.get("end", "NA"),
            patho_info.get("ARG", "NA"),
            patho_info.get("AMR_category", "NA"),
            patho_info.get("AMR_sub_class", "NA"),
            patho_info.get("Resistance_mechanism", "NA")
        ]
    else:
        patho_row = ["NA"]*8

    merged_rows.append(dram_row + patho_row)

# Write to TSV
with open(output_file, "w") as out:
    out.write("\t".join(output_fields) + "\n")
    for row in merged_rows:
        out.write("\t".join(str(x) for x in row) + "\n")

print(f"Merged file written to: {output_file}")
print(f"Total merged proteins: {len(merged_rows)}")