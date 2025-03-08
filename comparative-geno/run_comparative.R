library(DECIPHER)

setwd(".")
seqs_dir <- 'data/seqs'

DBPATH <- tempfile()

# Get all .fna files directly from the seqs directory
fna_files <- list.files(seqs_dir, pattern = "\\.fna$", full.names = TRUE)
fna_files <- fna_files[1:min(2, length(fna_files))]

GeneCalls <- vector('list', length = length(fna_files))

for (i in seq_along(fna_files)){
    
    fna_file <- fna_files[i]
    
    Seqs2DB(seqs = fna_file,
            type = "FASTA", 
            dbFile = DBPATH,
            identifier = as.character(i),
            verbose = TRUE)
    
    dnaGenome <- readDNAStringSet(fna_file)
    geneLocs <- FindGenes(dnaGenome)
    
    data("NonCodingRNA_Bacteria")
    ncRNA <- NonCodingRNA_Bacteria
    geneticRegions <- FindNonCoding(ncRNA, dnaGenome)

    GeneCalls[[i]] <- geneLocs
}
names(GeneCalls) <- seq_along(GeneCalls)

Syn <- FindSynteny(dbFile = DBPATH,
                   verbose = TRUE)

Overlaps <- NucleotideOverlap(SyntenyObject = Syn,
                               GeneCalls = GeneCalls,
                               Verbose = TRUE)

Pairs <- PairSummaries(SyntenyLinks = Overlaps,
                     GeneCalls = GeneCalls,
                     DBPATH = DBPATH,
                     PIDs = FALSE, 
                     Score = FALSE,
                     Verbose = TRUE)
