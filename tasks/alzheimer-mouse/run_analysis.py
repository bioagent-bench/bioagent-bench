import os
import pandas as pd
import numpy as np
import mygene
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests
import gseapy as gp
import requests
import networkx as nx
from io import StringIO
import logging


def setup_logging():
    """Configure logging settings"""
    os.makedirs("outputs", exist_ok=True)
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler('outputs/analysis.log'),
            logging.StreamHandler()
        ]
    )

def annotate_genes(df):
    logging.info(f"Starting gene annotation for {len(df)} genes")
    mg = mygene.MyGeneInfo()
    ensembl_ids = df.index.str.split(".").str[0].tolist()
    gene_info = mg.querymany(ensembl_ids, scopes='ensembl.gene', fields='symbol', species='mouse')
    gene_df = pd.DataFrame(gene_info).dropna(subset=['symbol']).drop_duplicates(subset='query')
    df.index = ensembl_ids
    df['Gene_Name'] = df.index.map(gene_df.set_index('query')['symbol'])
    cols = ['Gene_Name'] + [col for col in df.columns if col != 'Gene_Name']
    logging.info(f"Completed gene annotation. Found symbols for {len(gene_df)} genes")
    return df[cols]

def filter_normalize(data, min_cpm=0.7, min_samples=2):
    logging.info(f"Starting filtering and normalization. Initial genes: {len(data)}")
    gene_names = data.iloc[:, 0]
    raw_counts = data.iloc[:, 1:]
    cpm = raw_counts.apply(lambda x: (x / x.sum()) * 1e6, axis=0)
    mask = (cpm > min_cpm).sum(axis=1) >= min_samples
    filtered_counts = raw_counts[mask]
    filtered_gene_names = gene_names[mask]
    logging.info(f"After filtering: {len(filtered_counts)} genes remaining")
    
    geometric_means = filtered_counts.apply(lambda row: np.exp(np.log(row[row > 0]).mean()), axis=1)
    size_factors = filtered_counts.div(geometric_means, axis=0).median(axis=0)
    normalized_counts = filtered_counts.div(size_factors, axis=1)
    logging.info("Normalization completed")
    return pd.concat([filtered_gene_names, normalized_counts], axis=1)

def perform_dea(df, treated_columns, control_columns):
    logging.info(f"Starting differential expression analysis with {len(df)} genes")
    logging.info(f"Treatment samples: {len(treated_columns)}, Control samples: {len(control_columns)}")
    
    results = []
    for gene in df.index:
        treated = pd.to_numeric(df.loc[gene, treated_columns], errors='coerce').dropna()
        control = pd.to_numeric(df.loc[gene, control_columns], errors='coerce').dropna()
        if treated.empty or control.empty:
            continue
        log2fc = np.log2((treated.mean() + 1) / (control.mean() + 1))
        t_stat, p_val = ttest_ind(treated, control)
        results.append({
            "gene": gene,
            "Gene_Name": df.loc[gene, "Gene_Name"],
            "log2fc": log2fc,
            "t_stat": t_stat,
            "p_val": p_val
        })
    
    results_df = pd.DataFrame(results)
    results_df['p_val'] = pd.to_numeric(results_df['p_val'], errors='coerce')
    results_df = results_df.dropna(subset=['p_val'])
    results_df['p_adj'] = multipletests(results_df['p_val'], method='fdr_bh')[1]
    results_df['abs_log2fc'] = results_df['log2fc'].abs()
    
    final_results = results_df[(results_df['p_adj'] < 0.075) & (results_df['abs_log2fc'] > 0.75)]
    logging.info(f"DEA completed. Found {len(final_results)} significant genes")
    return final_results

def fetch_string_ppi(gene_list, species='10090'):
    logging.info(f"Fetching STRING PPI data for {len(gene_list)} genes")
    base_url = "https://string-db.org/api/tsv/network"
    genes = "\n".join(gene_list)
    params = {'identifiers': genes, 'species': species, 'limit': 1}
    try:
        response = requests.post(base_url, data=params)
        if response.status_code == 200:
            logging.info("Successfully retrieved PPI data from STRING")
            return response.text
        else:
            logging.error(f"Failed to fetch PPI data. Status code: {response.status_code}")
            return None
    except Exception as e:
        logging.error(f"Error fetching PPI data: {str(e)}")
        return None

def main():
    setup_logging()
    logging.info("Starting Alzheimer's mouse model analysis")
    
    os.makedirs("outputs", exist_ok=True)
    
    # 1. Load and annotate datasets
    logging.info("Loading 5xFAD dataset")
    countlist_5xFAD = pd.read_csv("./data/GSE168137_countList.txt", sep="\t", index_col=0)
    countlist_5xFAD = countlist_5xFAD.loc[:, countlist_5xFAD.columns.str.contains('cortex_8mo')]
    countlist_5xFAD.columns = [
        "5xFAD_cortex_8mon_Female_295", "5xFAD_cortex_8mon_Female_312", "5xFAD_cortex_8mon_Female_339", "5xFAD_cortex_8mon_Female_341", "5xFAD_cortex_8mon_Female_342",
        "5xFAD_cortex_8mon_Male_299", "5xFAD_cortex_8mon_Male_300", "5xFAD_cortex_8mon_Male_307", "5xFAD_cortex_8mon_Male_387", "5xFAD_cortex_8mon_Male_390",
        "BL6_cortex_8mon_Female_322", "BL6_cortex_8mon_Female_338", "BL6_cortex_8mon_Female_340", "BL6_cortex_8mon_Female_348", "BL6_cortex_8mon_Female_351",
        "BL6_cortex_8mon_Male_389", "BL6_cortex_8mon_Male_396", "BL6_cortex_8mon_Male_399", "BL6_cortex_8mon_Male_410", "BL6_cortex_8mon_Male_412"
    ]
    countlist_5xFAD = annotate_genes(countlist_5xFAD)

    logging.info("Loading 3xTG_AD dataset")
    countlist_3xTG_AD = pd.read_csv("./data/GSE161904_Raw_gene_counts_cortex.txt", sep="\t", index_col=0)
    countlist_3xTG_AD = countlist_3xTG_AD.loc[:, countlist_3xTG_AD.columns.str.contains('G3')]
    countlist_3xTG_AD.columns = ['3xTG_AD_Cortex_R1', '3xTG_AD_Cortex_R3', '3xTG_AD_Cortex_R4', 'WT_Cortex_R10', 'WT_Cortex_R17', 'WT_Cortex_R9']
    countlist_3xTG_AD = annotate_genes(countlist_3xTG_AD)

    logging.info("Loading PS3O1S dataset")
    DEA_PS3O1S = pd.read_csv("./data/DEA_PS3O1S.csv")

    # 2. Normalize
    logging.info("Starting data normalization")
    norm_5xFAD = filter_normalize(countlist_5xFAD, min_cpm=0.7)
    norm_3xTG_AD = filter_normalize(countlist_3xTG_AD, min_cpm=0.7, min_samples=1)

    # 3. DEA
    logging.info("Performing differential expression analysis")
    deg_5xFAD = perform_dea(norm_5xFAD,
                            treated_columns=[col for col in norm_5xFAD.columns if "5xFAD" in col],
                            control_columns=[col for col in norm_5xFAD.columns if "BL6" in col])
    deg_3xTG_AD = perform_dea(norm_3xTG_AD,
                              treated_columns=['3xTG_AD_Cortex_R1','3xTG_AD_Cortex_R3','3xTG_AD_Cortex_R4'],
                              control_columns=['WT_Cortex_R10','WT_Cortex_R17','WT_Cortex_R9'])
    DEA_PS3O1S['abs_log2fc'] = DEA_PS3O1S['log2fc'].abs()
    deg_PS3O1S = DEA_PS3O1S[(DEA_PS3O1S['abs_log2fc'] > 0.75) & (DEA_PS3O1S['pval'] < 0.075)]

    deg_5xFAD.to_csv("outputs/DEG_5xFAD.csv", index=False)
    deg_3xTG_AD.to_csv("outputs/DEG_3xTG_AD.csv", index=False)
    deg_PS3O1S.to_csv("outputs/DEG_PS3O1S.csv", index=False)

    # 4. GO Enrichment
    logging.info("Starting GO enrichment analysis")
    for label, gene_list in zip(
        ["5xFAD", "3xTG_AD", "PS3O1S"],
        [deg_5xFAD['Gene_Name'], deg_3xTG_AD['Gene_Name'], deg_PS3O1S['gene_name']]
    ):
        logging.info(f"Running GO enrichment for {label}")
        try:
            enr = gp.enrichr(gene_list=gene_list.dropna().tolist(), gene_sets=["GO_Biological_Process_2018"], organism='mouse')
            enr.results.to_csv(f"outputs/GO_BP_{label}.csv", index=False)
        except Exception as e:
            logging.error(f"Error in GO enrichment for {label}: {str(e)}")

    # 5. Pathway Enrichment
    for label, gene_list in zip(
        ["5xFAD", "3xTG_AD", "PS3O1S"],
        [deg_5xFAD['Gene_Name'], deg_3xTG_AD['Gene_Name'], deg_PS3O1S['gene_name']]
    ):
        enr = gp.enrichr(gene_list=gene_list.dropna().tolist(), gene_sets=["KEGG_2016"], organism='mouse')
        enr.results.to_csv(f"outputs/KEGG_{label}.csv", index=False)

    # 6. PPI Analysis
    gene_list_all = deg_5xFAD['Gene_Name'].dropna().tolist() + \
                    deg_3xTG_AD['Gene_Name'].dropna().tolist() + \
                    deg_PS3O1S['gene_name'].dropna().tolist()

    ppi_data = fetch_string_ppi(gene_list_all)
    ppi_df = pd.read_csv(StringIO(ppi_data), sep="\t")
    ppi_df_filtered = ppi_df[ppi_df['score'] > 0.7]
    ppi_df_filtered.to_csv("outputs/STRING_PPI_filtered.csv", index=False)

    G = nx.Graph()
    for _, row in ppi_df_filtered.iterrows():
        G.add_edge(row[0], row[1], weight=row['score'])

    degree_df = pd.DataFrame(sorted(nx.degree_centrality(G).items(), key=lambda x: x[1], reverse=True)[:10],
                             columns=["Protein", "Degree_Centrality"])
    degree_df.to_csv("outputs/Top10_Degree_Centrality.csv", index=False)

    betweenness_df = pd.DataFrame(sorted(nx.betweenness_centrality(G).items(), key=lambda x: x[1], reverse=True)[:10],
                                  columns=["Protein", "Betweenness_Centrality"])
    betweenness_df.to_csv("outputs/Top10_Betweenness_Centrality.csv", index=False)

    logging.info("Analysis complete! All outputs saved to the outputs/ directory")

if __name__ == "__main__":
    main()
