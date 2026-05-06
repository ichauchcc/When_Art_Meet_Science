import pandas as pd

input_file = "clinvar_SCN5A_raw.tsv"
output_file = "clinvar_SCN5A.tsv"

df = pd.read_csv(
    input_file,
    sep="\t",
    header=None,
    names=[
        "chr",
        "pos",
        "ref",
        "alt",
        "gene_raw",
        "consequence_raw",
        "classification",
        "conflicting_classification_detail"
    ]
)

# Clean gene: SCN5A:6331 -> SCN5A
df["gene"] = df["gene_raw"].astype(str).str.split(":").str[0]

# Clean consequence: SO:0001583|missense_variant -> missense_variant
df["consequence"] = (
    df["consequence_raw"]
    .astype(str)
    .str.split("|")
    .str[-1]
)

# Keep only SCN5A
df = df[df["gene"] == "SCN5A"]

# Replace missing conflict detail "." with empty string
df["conflicting_classification_detail"] = (
    df["conflicting_classification_detail"]
    .replace(".", "")
)

df_final = df[
    [
        "chr",
        "pos",
        "ref",
        "alt",
        "gene",
        "consequence",
        "classification",
        "conflicting_classification_detail"
    ]
]

df_final.to_csv(output_file, sep="\t", index=False)

print(f"Saved {len(df_final)} SCN5A ClinVar variants to {output_file}")
