library(DECIPHER)
library(SynExtend)

setwd('.')
COGExampleDir <- 'data/micrococcus_wgff'

DBPATH <- tempfile()

genomedirs <- list.files(COGExampleDir, full.names = TRUE)
genomedirs <- genomedirs[grep('json', genomedirs, fixed=T, invert=T)]

GeneCalls <- vector('list', length=length(genomedirs))

for (i in seq_along(genomedirs)){
  subfiles <- list.files(genomedirs[i], full.names = TRUE)
  
  fna_file <- subfiles[which(grepl('.*fna$', subfiles))]
  gff_file <- subfiles[which(grepl('.*gff$', subfiles))]
  
  Seqs2DB(seqs = fna_file,
          type = "FASTA",
          dbFile = DBPATH,
          identifier = as.character(i),
          verbose = TRUE)
  
  GeneCalls[[i]] <- gffToDataFrame(GFF = gff_file,
                                    Verbose = TRUE)
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
                     PIDs = FALSE, # Set to TRUE for better accuracy (slower) 
                     Score = FALSE, # Set to TRUE for better accuracy (slower)
                     Verbose = TRUE)

# These methods only work if we set PIDs and Score to TRUE
# Unfortunately we don't have time in this workshop to use these
# Feel free to try them out on your own with a larger dataset!

# P02 <- BlockExpansion(Pairs = P01,
#                       DBPATH = DBPATH,
#                       Verbose = TRUE,
#                       NewPairsOnly = FALSE)
# P03 <- BlockReconciliation(Pairs = P02,
#                            PIDThreshold = 0.75,
#                            SCOREThreshold = 200,
#                            Verbose = TRUE)
# Pairs <- P03[P03$PID > 0.4, ]

head(Pairs)

## Finding COGs
COGSets <- DisjointSet(Pairs = Pairs,
                        Verbose = TRUE)

COGSets[1:3]

# Extract sequences for COGs with at least 5 orthologs
Sequences <- ExtractBy(x = Pairs,
                       y = DBPATH,
                       z = COGSets[lengths(COGSets) >= 5],
                       Verbose = TRUE)

# These come back in different orders, so let's match them up
allnames <- lapply(Sequences, names)
COGMapping <- sapply(COGSets, function(x) {
                       which(sapply(allnames, function(y) setequal(x, y)))
                     })
COGMapping <- COGMapping[sapply(COGMapping, function(x) length(x) > 0)]

MatchedCOGSets <- COGSets[names(COGMapping)]
MatchedSequences <- Sequences[unlist(COGMapping)]
names(MatchedSequences) <- names(COGMapping)


MatchedCOGSets[1:3]
MatchedSequences[1:3]