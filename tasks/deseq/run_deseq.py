import pandas as pd
import numpy as np
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
import os
from goatools.base import download_go_basic_obo
from goatools.obo_parser import GODag
from goatools.goea.go_enrichment_ns import GOEnrichmentStudyNS
from goatools.associations import dnld_assc

# Set working directory
os.chdir("/app")

# Read count data from multiple files
count_files = [
    "data/processing/7_mapping/SRR1278968_ReadsPerGene.out.tab",
    "data/processing/7_mapping/SRR1278969_ReadsPerGene.out.tab",
    "data/processing/7_mapping/SRR1278970_ReadsPerGene.out.tab",
    "data/processing/7_mapping/SRR1278971_ReadsPerGene.out.tab",
    "data/processing/7_mapping/SRR1278972_ReadsPerGene.out.tab",
    "data/processing/7_mapping/SRR1278973_ReadsPerGene.out.tab",
]

sample_names = [
    "SRR1278968",
    "SRR1278969",
    "SRR1278970",
    "SRR1278971",
    "SRR1278972",
    "SRR1278973",
]

# Read and combine count data
count_data = []
gene_ids = None

for i, file_path in enumerate(count_files):
    df = pd.read_csv(file_path, sep="\t", header=0)
    if gene_ids is None:
        gene_ids = df.iloc[:, 0]  # First column contains gene IDs
    count_data.append(df.iloc[:, 3])  # Fourth column contains counts

# Create count matrix
count_matrix = pd.concat(count_data, axis=1)
count_matrix.columns = sample_names
count_matrix.index = gene_ids

# Remove first 3 rows (summary statistics from STAR)
count_matrix = count_matrix.iloc[3:]

# Create sample metadata
sample_info = pd.DataFrame(
    {
        "sample": sample_names,
        "condition": [
            "Plankton",
            "Plankton",
            "Plankton",
            "Biofilm",
            "Biofilm",
            "Biofilm",
        ],
    }
)
sample_info = sample_info.set_index("sample")

# Transpose count matrix to have samples as rows and genes as columns
count_matrix_t = count_matrix.T  # Transpose: samples as rows, genes as columns

print("Count matrix shape (samples x genes):", count_matrix_t.shape)
print("Sample info shape:", sample_info.shape)
print("Count matrix index (samples):", count_matrix_t.index.tolist())
print("Sample info index (samples):", sample_info.index.tolist())

# Run DESeq2 analysis
dds = DeseqDataSet(
    counts=count_matrix_t,
    metadata=sample_info,
    design_factors="condition",
    refit_cooks=True,
)

dds.deseq2()

# Get results - comparing Biofilm vs Plankton
stat_res = DeseqStats(dds, contrast=("condition", "Biofilm", "Plankton"))
stat_res.summary()
results_df = stat_res.results_df

# Filter for significantly upregulated genes in biofilm
# log2FoldChange > 2 and padj < 0.01
upreg_mask = (
    (results_df["log2FoldChange"] > 2)
    & (results_df["padj"] < 0.01)
    & (~results_df["padj"].isna())
)
upreg_in_biofilm = results_df[upreg_mask].index.tolist()

# Save upregulated genes as CSV with statistics
upreg_results = results_df[upreg_mask][["log2FoldChange", "pvalue", "padj"]].copy()
upreg_results.index.name = "gene_id"
upreg_results.to_csv("results/up_regulated_genes.csv", index=True)

print(f"Found {len(upreg_in_biofilm)} upregulated genes in biofilm condition")

print("Analysis complete!")