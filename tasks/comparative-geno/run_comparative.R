library(DECIPHER)
library(SynExtend)

# Enable debugging output
cat("=== COMPARATIVE GENOMICS ANALYSIS DEBUG MODE ===\n")
cat("Script started at:", as.character(Sys.time()), "\n\n")

# Set working directory and check data
setwd('.')
COGExampleDir <- 'data'

cat("Step 1: Checking data directory...\n")
if (!dir.exists(COGExampleDir)) {
  stop("ERROR: Data directory '", COGExampleDir, "' does not exist!")
}

cat("Data directory found:", COGExampleDir, "\n")
cat("Creating temporary database...\n")
DBPATH <- tempfile()
cat("Database path:", DBPATH, "\n")

cat("\nStep 2: Finding genome directories...\n")
genomedirs <- list.files(COGExampleDir, full.names = TRUE)
genomedirs <- genomedirs[grep('json', genomedirs, fixed=T, invert=T)]
cat("Found", length(genomedirs), "genome directories:\n")
for (i in seq_along(genomedirs)) {
  cat("  ", i, ":", basename(genomedirs[i]), "\n")
}

if (length(genomedirs) == 0) {
  stop("ERROR: No genome directories found!")
}

cat("\nStep 3: Processing genome files and building database...\n")
GeneCalls <- vector('list', length=length(genomedirs))

for (i in seq_along(genomedirs)){
  cat("Processing genome", i, "of", length(genomedirs), ":", basename(genomedirs[i]), "\n")
  
  tryCatch({
    subfiles <- list.files(genomedirs[i], full.names = TRUE)
    cat("  Found", length(subfiles), "files in directory\n")
    
    fna_file <- subfiles[which(grepl('.*fna$', subfiles))]
    gff_file <- subfiles[which(grepl('.*gff$', subfiles))]
    
    if (length(fna_file) == 0) {
      stop("No .fna file found in ", genomedirs[i])
    }
    if (length(gff_file) == 0) {
      stop("No .gff file found in ", genomedirs[i])
    }
    
    cat("  FNA file:", basename(fna_file), "\n")
    cat("  GFF file:", basename(gff_file), "\n")
    
    cat("  Adding sequences to database...\n")
    Seqs2DB(seqs = fna_file,
            type = "FASTA",
            dbFile = DBPATH,
            identifier = as.character(i),
            verbose = TRUE)
    
    cat("  Parsing GFF file...\n")
    GeneCalls[[i]] <- gffToDataFrame(GFF = gff_file,
                                      Verbose = TRUE)
    cat("  Found", nrow(GeneCalls[[i]]), "gene calls\n")
    
  }, error = function(e) {
    cat("ERROR processing genome", i, ":", e$message, "\n")
    stop("Failed at genome processing step")
  })
}
names(GeneCalls) <- seq_along(GeneCalls)
cat("Successfully processed all genomes!\n")

cat("\nStep 4: Finding synteny...\n")
tryCatch({
  Syn <- FindSynteny(dbFile = DBPATH,
                     verbose = TRUE)
  cat("Synteny analysis completed successfully\n")
}, error = function(e) {
  cat("ERROR in FindSynteny:", e$message, "\n")
  stop("Failed at synteny step")
})

cat("\nStep 5: Finding nucleotide overlaps...\n")
tryCatch({
  Overlaps <- NucleotideOverlap(SyntenyObject = Syn,
                                 GeneCalls = GeneCalls,
                                 Verbose = TRUE)
  cat("Nucleotide overlap analysis completed successfully\n")
}, error = function(e) {
  cat("ERROR in NucleotideOverlap:", e$message, "\n")
  stop("Failed at nucleotide overlap step")
})

cat("\nStep 6: Generating pair summaries...\n")
tryCatch({
  P01 <- PairSummaries(SyntenyLinks = Overlaps,
                       GeneCalls = GeneCalls,
                       DBPATH = DBPATH,
                       PIDs = TRUE, 
                       Score = TRUE,
                       Verbose = TRUE,
                       processors=NULL)
  cat("Pair summaries completed successfully, found", nrow(P01), "pairs\n")
}, error = function(e) {
  cat("ERROR in PairSummaries:", e$message, "\n")
  stop("Failed at pair summaries step")
})

cat("\nStep 7: Block expansion...\n")
tryCatch({
  P02 <- BlockExpansion(Pairs = P01,
                        DBPATH = DBPATH,
                        Verbose = TRUE,
                        NewPairsOnly = FALSE)
  cat("Block expansion completed successfully, now have", nrow(P02), "pairs\n")
}, error = function(e) {
  cat("ERROR in BlockExpansion:", e$message, "\n")
  stop("Failed at block expansion step")
})

cat("\nStep 8: Block reconciliation...\n")
tryCatch({
  P03 <- BlockReconciliation(Pairs = P02,
                             PIDThreshold = 0.75,
                             SCOREThreshold = 200,
                             Verbose = TRUE)
  cat("Block reconciliation completed successfully\n")
  
  Pairs <- P03[P03$PID > 0.4, ]
  cat("After PID filtering (>0.4):", nrow(Pairs), "pairs remain\n")
}, error = function(e) {
  cat("ERROR in BlockReconciliation:", e$message, "\n")
  stop("Failed at block reconciliation step")
})

cat("\nStep 9: Finding COGs (Clusters of Orthologous Groups)...\n")
tryCatch({
  COGSets <- DisjointSet(Pairs = Pairs,
                          Verbose = TRUE)
  cat("Found", length(COGSets), "COG sets\n")
  
  # Print size distribution
  sizes <- lengths(COGSets)
  cat("COG size distribution:\n")
  cat("  Min:", min(sizes), "Max:", max(sizes), "Mean:", round(mean(sizes), 2), "\n")
  cat("  COGs with ≥4 members:", sum(sizes >= 4), "\n")
}, error = function(e) {
  cat("ERROR in DisjointSet:", e$message, "\n")
  stop("Failed at COG finding step")
})

cat("\nStep 10: Extracting sequences for large COGs...\n")
tryCatch({
  # Extract sequences for COGs with at least 5 orthologs
  largeCOGs <- COGSets[lengths(COGSets) >= 4]
  cat("Extracting sequences for", length(largeCOGs), "COGs with ≥4 members\n")
  
  Sequences <- ExtractBy(x = Pairs,
                         y = DBPATH,
                         z = largeCOGs,
                         Verbose = TRUE)
  cat("Successfully extracted sequences for", length(Sequences), "COGs\n")
}, error = function(e) {
  cat("ERROR in ExtractBy:", e$message, "\n")
  stop("Failed at sequence extraction step")
})

cat("\nStep 11: Matching COG sets with sequences...\n")
tryCatch({
  # These come back in different orders, so let's match them up
  allnames <- lapply(Sequences, names)
  COGMapping <- sapply(COGSets, function(x) {
                         which(sapply(allnames, function(y) setequal(x, y)))
                       })
  COGMapping <- COGMapping[sapply(COGMapping, function(x) length(x) > 0)]
  
  MatchedCOGSets <- COGSets[names(COGMapping)]
  MatchedSequences <- Sequences[unlist(COGMapping)]
  names(MatchedSequences) <- names(COGMapping)
  
  cat("Successfully matched", length(MatchedSequences), "COG sets with sequences\n")
}, error = function(e) {
  cat("ERROR in COG matching:", e$message, "\n")
  stop("Failed at COG matching step")
})

cat("\nStep 12: Building phylogenetic trees...\n")
tryCatch({
  # Build phylogenetic trees for each COG
  # Create a list to store all the trees
  COGTrees <- vector("list", length=length(MatchedSequences))
  names(COGTrees) <- names(MatchedSequences)
  
  cat("Building trees for", length(MatchedSequences), "COGs\n")
  
  # Loop through each COG and build a tree
  for (i in seq_along(MatchedSequences)) {
    cat("Building tree", i, "of", length(MatchedSequences), "\n")
    
    tryCatch({
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
      
    }, error = function(e) {
      cat("WARNING: Failed to build tree for COG", i, ":", e$message, "\n")
      COGTrees[[i]] <- NULL
    })
  }
  
  # Remove failed trees
  successful_trees <- !sapply(COGTrees, is.null)
  cat("Successfully built", sum(successful_trees), "trees out of", length(COGTrees), "attempted\n")
  
}, error = function(e) {
  cat("ERROR in tree building:", e$message, "\n")
  stop("Failed at tree building step")
})

cat("\nStep 13: Saving trees and sequences...\n")
tryCatch({
  save(MatchedSequences, COGTrees, file = "COG_trees_and_sequences.RData")
  cat("Successfully saved trees and sequences\n")
}, error = function(e) {
  cat("ERROR saving trees and sequences:", e$message, "\n")
  stop("Failed at saving step")
})

cat("\nStep 14: Functional annotation...\n")
tryCatch({
  # Create a list to store all the annotations
  CogsAnnot <- vector("list", length=length(MatchedSequences))
  names(CogsAnnot) <- names(MatchedSequences)
  
  # Translate the sequences
  cat("Translating sequences...\n")
  geneSeqs <- lapply(MatchedSequences, translate)
  
  # Load the training set for taxonomic classification
  cat("Loading training set...\n")
  if (!file.exists('data/Actinobacteria.RData')) {
    stop("Training set file 'data/Actinobacteria.RData' not found!")
  }
  load('data/Actinobacteria.RData')  # This loads 'trainingSet'
  cat("Training set loaded successfully\n")
  
  # Loop through each COG and annotate it
  cat("Annotating", length(geneSeqs), "COGs...\n")
  for (i in seq_along(geneSeqs)) {
    cat("Annotating COG", i, "of", length(geneSeqs), "\n")
    
    tryCatch({
      # Get the current COG protein sequences
      currentProtein <- geneSeqs[[i]]
      CogsAnnot[[i]] <- IdTaxa(currentProtein, trainingSet, processors=NULL)
      
    }, error = function(e) {
      cat("WARNING: Failed to annotate COG", i, ":", e$message, "\n")
      CogsAnnot[[i]] <- NULL
    })
  }
  
  successful_annotations <- !sapply(CogsAnnot, is.null)
  cat("Successfully annotated", sum(successful_annotations), "COGs out of", length(CogsAnnot), "attempted\n")
  
}, error = function(e) {
  cat("ERROR in functional annotation:", e$message, "\n")
  stop("Failed at functional annotation step")
})

cat("\nStep 15: Saving annotations...\n")
tryCatch({
  save(MatchedSequences, COGTrees, CogsAnnot, file = "COG_trees_sequences_annotations.RData")
  cat("Successfully saved annotations\n")
}, error = function(e) {
  cat("ERROR saving annotations:", e$message, "\n")
  stop("Failed at saving annotations step")
})

cat("\nStep 16: Filtering COGs...\n")
tryCatch({
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
  
  # Print summary of filtering
  cat("Original COGs:", length(MatchedSequences), "\n")
  cat("COGs after filtering:", length(FilteredSequences), "\n")
  cat("Breakdown of filters:\n")
  cat("  No paralogs:", sum(noParas), "\n")
  cat("  In 4+ organisms:", sum(inFourOrMore), "\n")
  cat("  Coding elements:", sum(codingCOGs), "\n")
  cat("  High confidence annotation:", sum(highConf), "\n")
  
  if (length(FilteredSequences) == 0) {
    stop("ERROR: No COGs passed all filtering criteria!")
  }
  
}, error = function(e) {
  cat("ERROR in COG filtering:", e$message, "\n")
  stop("Failed at COG filtering step")
})

cat("\nStep 17: Saving filtered data...\n")
tryCatch({
  save(FilteredSequences, FilteredTrees, FilteredAnnots, file = "Filtered_COG_data.RData")
  cat("Successfully saved filtered data\n")
}, error = function(e) {
  cat("ERROR saving filtered data:", e$message, "\n")
  stop("Failed at saving filtered data step")
})

cat("\nStep 18: Creating consensus annotations...\n")
tryCatch({
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
  cat("Created consensus annotations for", length(consAnnots), "COGs\n")
}, error = function(e) {
  cat("ERROR creating consensus annotations:", e$message, "\n")
  stop("Failed at consensus annotation step")
})

cat("\nStep 19: Creating EvoWeaver object...\n")
tryCatch({
  pw <- EvoWeaver(FilteredTrees)
  cat("EvoWeaver object created successfully\n")
}, error = function(e) {
  cat("ERROR creating EvoWeaver object:", e$message, "\n")
  stop("Failed at EvoWeaver creation step")
})

cat("\nStep 20: Making predictions...\n")
tryCatch({
  preds <- predict(pw)
  cat("Predictions completed successfully\n")
  print(preds)
}, error = function(e) {
  cat("ERROR making predictions:", e$message, "\n")
  stop("Failed at prediction step")
})

cat("\nStep 21: Clustering analysis...\n")
tryCatch({
  # Find clusters of coevolving COGs using igraph
  library(igraph)
  set.seed(123) # For reproducibility
  
  adjMatrix <- as.matrix(preds)
  g <- graph_from_adjacency_matrix(adjMatrix, weighted=TRUE,
                                 mode='undirected', diag=FALSE)
  
  clusters <- cluster_louvain(g)
  cat("Clustering completed, found", length(clusters), "clusters\n")
  
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
  
  cat("Cluster labels created successfully\n")
}, error = function(e) {
  cat("ERROR in clustering analysis:", e$message, "\n")
  stop("Failed at clustering step")
})

cat("\nStep 22: Saving final results...\n")
tryCatch({
  save(preds, clusters, clusterLabels, consAnnots, file = "EvoWeaver_results.RData")
  cat("Final results saved successfully\n")
}, error = function(e) {
  cat("ERROR saving final results:", e$message, "\n")
  stop("Failed at saving final results step")
})

cat("\nStep 23: Generating final report...\n")
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

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("Script completed successfully at:", as.character(Sys.time()), "\n")
cat("Results saved to EvoWeaver_results.RData\n")