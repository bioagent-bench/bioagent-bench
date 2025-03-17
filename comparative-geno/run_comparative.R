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

P01 <- PairSummaries(SyntenyLinks = Overlaps,
                     GeneCalls = GeneCalls,
                     DBPATH = DBPATH,
                     PIDs = TRUE, 
                     Score = TRUE,
                     Verbose = TRUE,
                     processors=NULL)

P02 <- BlockExpansion(Pairs = P01,
                      DBPATH = DBPATH,
                      Verbose = TRUE,
                      NewPairsOnly = FALSE)
P03 <- BlockReconciliation(Pairs = P02,
                           PIDThreshold = 0.75,
                           SCOREThreshold = 200,
                           Verbose = TRUE)
Pairs <- P03[P03$PID > 0.4, ]

## Finding COGs
COGSets <- DisjointSet(Pairs = Pairs,
                        Verbose = TRUE)

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

# Build phylogenetic trees for each COG
# Create a list to store all the trees
COGTrees <- vector("list", length=length(MatchedSequences))
names(COGTrees) <- names(MatchedSequences)

# Loop through each COG and build a tree
for (i in seq_along(MatchedSequences)) {
  # Get the current COG sequences
  currentCOG <- MatchedSequences[[i]]
  
  # Since these are coding regions, align using AlignTranslation
  alignedCOG <- AlignTranslation(currentCOG)
  
  # Build a Maximum Likelihood tree with ancestral state reconstruction
  # Using a reasonable time limit to balance accuracy and runtime
  COGTrees[[i]] <- TreeLine(alignedCOG, 
                           method = "ML",
                           reconstruct = TRUE,
                           maxTime = 0.05,
                           processors=NULL)
  
  # Print progress
  cat("Completed tree", i, "of", length(MatchedSequences), "\n")
}

# Save the trees and sequences for later use
save(MatchedSequences, COGTrees, file = "COG_trees_and_sequences.RData")

# Now I need to loop through each COG and translate it and annotate it functionally
# Create a list to store all the annotations
CogsAnnot <- vector("list", length=length(MatchedSequences))
names(CogsAnnot) <- names(MatchedSequences)

# Translate the sequences
geneSeqs <- lapply(MatchedSequences, translate)

# Load the training set for taxonomic classification
load('data/Actinobacteria.RData')  # This loads 'trainingSet'

# Loop through each COG and annotate it
for (i in seq_along(geneSeqs)) {
  # Get the current COG protein sequences
  currentProtein <- geneSeqs[[i]]
  CogsAnnot[[i]] <- IdTaxa(currentProtein, trainingSet, processors=NULL)
  
  # Print progress
  cat("Completed annotation", i, "of", length(geneSeqs), "\n")
}

# Save the annotations for later use
save(MatchedSequences, COGTrees, CogsAnnot, file = "COG_trees_sequences_annotations.RData")

# Subsetting COGs based on specific criteria
# Get assembly identifiers for each COG
truncCOGs <- lapply(MatchedSequences, function(x) sort(as.integer(gsub('^([0-9]+)_.*', '\\1', names(x)))))

# Find COGs without paralogs (each genome appears at most once)
noParas <- sapply(truncCOGs, function(x) length(x) == length(unique(x)))

# Get genes in 4 or more organisms
inFourOrMore <- sapply(truncCOGs, function(x) length(unique(x)) >= 4)

# Make sure COGs are coding elements
codingCOGs <- sapply(CogsAnnot, function(x) is(x, 'Taxa'))

# At least one high confidence annotation
highConf <- sapply(CogsAnnot, function(x) 
                   if(is(x, 'Taxa')) 
                     max(sapply(x, function(y) 
                                y$confidence[length(y$confidence)])) > 50
                   else FALSE
                   )

# Apply all filters
FilteredCOGs <- noParas & inFourOrMore & codingCOGs & highConf

# Subset our data
FilteredSequences <- MatchedSequences[FilteredCOGs]
FilteredAnnots <- CogsAnnot[FilteredCOGs]
FilteredTrees <- COGTrees[FilteredCOGs]

# Save the filtered data
save(FilteredSequences, FilteredTrees, FilteredAnnots, file = "Filtered_COG_data.RData")

# Print summary of filtering
cat("Original COGs:", length(MatchedSequences), "\n")
cat("COGs after filtering:", length(FilteredSequences), "\n")
cat("Breakdown of filters:\n")
cat("  No paralogs:", sum(noParas), "\n")
cat("  In 4+ organisms:", sum(inFourOrMore), "\n")
cat("  Coding elements:", sum(codingCOGs), "\n")
cat("  High confidence annotation:", sum(highConf), "\n")

consAnnots <- vector('character', length=length(FilteredAnnots))
for (i in seq_along(FilteredAnnots)) {
  taxaentry <- FilteredAnnots[[i]]
  
  # If no annotation, it's a noncoding gene
  if (!is(taxaentry, 'Taxa'))
    consAnnots[i] <- 'NONCODING'
  # Otherwise it's a coding gene
  else {
    # Grab all the annotations aside from "Unclassified"
    annots <- sapply(taxaentry, function(y) y$taxon[length(y$taxon)])
    annots <- annots[annots != 'unclassified_Root']
    
    # If we only have "Unclassified", just mark it as uncharacterized
    if (length(annots) == 0)
      consAnnots[i] <- 'Uncharacterized'
    
    # Otherwise take the most common annotation
    else
      consAnnots[i] <- names(sort(table(annots), decreasing=TRUE))[1]
  }
}

# Create EvoWeaver object
pw <- EvoWeaver(FilteredTrees)

# Make predictions of functional associations
preds <- predict(pw)
print(preds)

# Find clusters of coevolving COGs using igraph
library(igraph)
set.seed(123) # For reproducibility

adjMatrix <- as.matrix(preds)
g <- graph_from_adjacency_matrix(adjMatrix, weighted=TRUE,
                               mode='undirected', diag=FALSE)

clusters <- cluster_louvain(g)

# Getting the clusters & identifying COGs by consensus annotation
clusterLabels <- vector('list', length(clusters))
for (i in seq_along(clusterLabels)) {
  cluster <- communities(clusters)[[i]]
  labs <- consAnnots[as.integer(cluster)]
  clusterLabels[[i]] <- labs[order(sapply(labs, function(x) 
                                         ifelse(grepl(" ", x), 
                                                strsplit(x, ' ')[[1]][3], 
                                                x)))]
}

# Save the clustering results
save(preds, clusters, clusterLabels, consAnnots, file = "EvoWeaver_results.RData")

# Print information about the clusters
cat("Number of clusters identified:", length(clusters), "\n")
for (i in seq_along(clusterLabels)) {
  cat("\nCluster", i, "contains", length(clusterLabels[[i]]), "COGs\n")
  if (length(clusterLabels[[i]]) <= 10) {
    cat("Annotations:\n")
    for (j in seq_along(clusterLabels[[i]])) {
      cat("  ", clusterLabels[[i]][j], "\n")
    }
  } else {
    cat("First 10 annotations:\n")
    for (j in 1:10) {
      cat("  ", clusterLabels[[i]][j], "\n")
    }
    cat("  ... and", length(clusterLabels[[i]]) - 10, "more\n")
  }
}

cat("Analysis complete. Results saved to EvoWeaver_results.RData\n")