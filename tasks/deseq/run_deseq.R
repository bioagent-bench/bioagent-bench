options(repos = c(CRAN = "https://cloud.r-project.org"))

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "3.18")

# Install required packages
BiocManager::install(c(
    "DESeq2",
    "BiocParallel",
    "clusterProfiler",
    "GenomeInfoDbData",
    "GO.db"
))

if (!require("ggpubr"))
    install.packages("ggpubr")

library(DESeq2)
library(ggpubr)
library(clusterProfiler)
library(GO.db)


setwd(".")
data <- NULL
for (f in c("data/processing/7_mapping/SRR1278968_ReadsPerGene.out.tab",
            "data/processing/7_mapping/SRR1278969_ReadsPerGene.out.tab",
            "data/processing/7_mapping/SRR1278970_ReadsPerGene.out.tab",
            "data/processing/7_mapping/SRR1278971_ReadsPerGene.out.tab", 
            "data/processing/7_mapping/SRR1278972_ReadsPerGene.out.tab", 
            "data/processing/7_mapping/SRR1278973_ReadsPerGene.out.tab")) {
  file <- read.table(f, header = TRUE, sep = "\t")
  data <- cbind(data, file[,4])
}

rownames(data) <- file[,1]
colnames(data) <- c("SRR1278968", "SRR1278969", "SRR1278970",
                    "SRR1278971", "SRR1278972", "SRR1278973")
data <- data[-c(1,2,3),]    

groups <- factor(c(rep("Plankton", 3), rep("Biofilm", 3)), levels=c("Plankton","Biofilm"))
colData <- DataFrame(condition=groups)
pre_dds <- DESeqDataSetFromMatrix(data, colData, design = ~condition)
dds <- DESeq(pre_dds)
res <- results(dds, contrast = c("condition", "Biofilm","Plankton"))

ggmaplot(res, main = "Differential expression in Biofilm vs Planktonic form",
    fdr = 0.01, fc = 4, size = 3,
    palette = c("red", "blue", "darkgray"),
    legend = "top", top = 0)
ggsave("data/processing/8_deseq/maplot.png")


norm_data<-vst(dds)
plotPCA(norm_data,ntop=5000)
ggsave("data/processing/8_deseq/pca.png")

upreg_in_biofilm<-row.names(res[!is.na(res$padj) & res$log2FoldChange > 2 & res$padj<0.01,])
write.table(upreg_in_biofilm, file="data/processing/8_deseq/up_regulated_genes.txt", sep = "\t", quote = FALSE, col.names = F, row.names = F )

all_go_tems<-as.data.frame(GOTERM)
all_go_tems<-all_go_tems[!(all_go_tems$go_id %in% c("GO:0003674","GO:0008150","GO:0005575")),]
BP_terms<-all_go_tems$go_id[all_go_tems$Ontology=="BP"]

all_lines <- readLines("data/processing/9_gorich/gene_association.cgd")
data_lines <- all_lines[!grepl("^!", all_lines)]
filtered_lines <- data_lines[grepl("578454", data_lines)]
cpar_go <- read.table(text = filtered_lines, sep = "\t", stringsAsFactors = FALSE)
cpar_go$V2 <- gsub("\\|.*$", "", cpar_go$V11)
BP_universe <- cpar_go[cpar_go$V5 %in% BP_terms, c("V5", "V2", "V9")]
colnames(BP_universe) <- c("V1", "V2", "V3")

goterms <- Term(GOTERM)
a <- as.data.frame(goterms)
go_names <- data.frame(
  row.names(a),
  a$goterms,
  stringsAsFactors = FALSE
)
colnames(go_names) <- c("row.names(a)", "goterms")

upreg_in_biofilm <- read.table("data/processing/8_deseq/up_regulated_genes.txt",header = FALSE)
upreg_in_biofilm <- as.character(upreg_in_biofilm$V1)


ego<-enricher(upreg_in_biofilm, universe = as.character(BP_universe$V2), TERM2GENE = BP_universe,TERM2NAME = go_names)

dotplot(ego, showCategory = 15)
ggsave("data/processing/9_gorich/up_reg_go_terms.png")