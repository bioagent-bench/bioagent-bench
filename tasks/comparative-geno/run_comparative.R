library(DECIPHER)
library(SynExtend)

# Enable debugging output
cat("=== COMPARATIVE GENOMICS ANALYSIS DEBUG MODE ===\n")
cat("Script started at:", as.character(Sys.time()), "\n\n")

# For automated debugging runs, uncomment one of these lines:
# resume_from_checkpoint <- TRUE  # Skip to Step 9.5 (filtering and tree building)
# resume_from_prediction <- TRUE  # Skip to Step 14 (consensus annotations and predictions)

# Check for intermediate checkpoint to resume from
intermediate_file <- "./outputs/intermediate_step9.RData"
prediction_checkpoint <- "./outputs/COG_trees_sequences_annotations.RData"
resume_from_checkpoint <- FALSE
resume_from_prediction <- FALSE

if (file.exists(prediction_checkpoint)) {
  cat("Found prediction checkpoint file:", prediction_checkpoint, "\n")
  if (!resume_from_prediction) {
    cat("Do you want to resume from Step 15 (EvoWeaver predictions)? [Y/n]: ")
    # For automated runs, you can comment out the readline and set resume_from_prediction = TRUE
    user_input <- readline()
    if (tolower(trimws(user_input)) %in% c("", "y", "yes")) {
      resume_from_prediction <- TRUE
    }
  }
  
  if (resume_from_prediction) {
    cat("Resuming from prediction checkpoint...\n")
    load(prediction_checkpoint)
    cat("Loaded prediction checkpoint successfully.\n")
    cat("Variables loaded: FilteredSequences, FilteredAnnots, COGTrees\n")
    cat("Skipping to Step 14: Creating consensus annotations...\n\n")
  } else {
    cat("Checking for earlier checkpoint...\n")
  }
}

if (!resume_from_prediction && file.exists(intermediate_file)) {
  cat("Found intermediate checkpoint file:", intermediate_file, "\n")
  cat("Do you want to resume from Step 9.5 (filtering and tree building)? [Y/n]: ")
  # For automated runs, you can comment out the readline and set resume_from_checkpoint = TRUE
  user_input <- readline()
  if (tolower(trimws(user_input)) %in% c("", "y", "yes")) {
    resume_from_checkpoint <- TRUE
    cat("Resuming from checkpoint...\n")
    load(intermediate_file)
    cat("Loaded intermediate results successfully.\n")
    cat("Variables loaded: GeneCalls, Syn, Overlaps, Pairs, COGSets,\n")
    cat("                  MatchedCOGSets, MatchedSequences, Sequences, genomedirs\n")
    cat("Skipping to Step 9.5: Early COG filtering...\n\n")
  } else {
    cat("Starting fresh analysis from the beginning...\n")
  }
} else if (!resume_from_prediction) {
  cat("No intermediate checkpoint found. Starting fresh analysis...\n")
}

# Set working directory and check data (skip if resuming)
if (!resume_from_checkpoint && !resume_from_prediction) {
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
genomedirs <- genomedirs[grep('RData', genomedirs, fixed=T, invert=T)]
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
  Pairs <- PairSummaries(SyntenyLinks = Overlaps,
                       GeneCalls = GeneCalls,
                       DBPATH = DBPATH,
                       PIDs = TRUE, 
                       Score = TRUE,
                       Verbose = TRUE,
                       processors=NULL)

  P02 <- BlockExpansion(Pairs = Pairs,
                      DBPATH = DBPATH,
                      Verbose = TRUE,
                      NewPairsOnly = FALSE)
  P03 <- BlockReconciliation(Pairs = P02,
                           PIDThreshold = 0.75,
                           SCOREThreshold = 200,
                           Verbose = TRUE)
  Pairs <- P03[P03$PID > 0.4, ]
  
  cat("Pair summaries completed successfully, found", nrow(Pairs), "pairs\n")
}, error = function(e) {
  cat("ERROR in PairSummaries:", e$message, "\n")
  stop("Failed at pair summaries step")
})

cat("\nStep 7: Finding COGs (Clusters of Orthologous Groups)...\n")
tryCatch({
  COGSets <- DisjointSet(Pairs = Pairs,
                          Verbose = TRUE)
  cat("Found", length(COGSets), "COG sets\n")
  
  # Print size distribution
  sizes <- lengths(COGSets)
  cat("COG size distribution:\n")
  cat("  Min:", min(sizes), "Max:", max(sizes), "Mean:", round(mean(sizes), 2), "\n")
  cat("  COGs with ≥5 members:", sum(sizes >= 5), "\n")
}, error = function(e) {
  cat("ERROR in DisjointSet:", e$message, "\n")
  stop("Failed at COG finding step")
})

cat("\nStep 8: Extracting sequences for large COGs...\n")
tryCatch({
  # Extract sequences for COGs with at least 5 orthologs
  largeCOGs <- COGSets[lengths(COGSets) >= 5]
  cat("Extracting sequences for", length(largeCOGs), "COGs with ≥5 members\n")
  
  Sequences <- ExtractBy(x = Pairs,
                         y = DBPATH,
                         z = largeCOGs,
                         Verbose = TRUE)
  cat("Successfully extracted sequences for", length(Sequences), "COGs\n")
}, error = function(e) {
  cat("ERROR in ExtractBy:", e$message, "\n")
  stop("Failed at sequence extraction step")
})

cat("\nStep 9: Matching COG sets with sequences...\n")
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

cat("\nStep 9.1: Saving intermediate results (before filtering)...\n")
tryCatch({
  # Create outputs directory if it doesn't exist
  if (!dir.exists("./outputs")) {
    dir.create("./outputs", recursive = TRUE)
  }
  
  # Save all intermediate results up to this point (before filtering)
  save(GeneCalls, Syn, Overlaps, Pairs, COGSets, MatchedCOGSets, 
       MatchedSequences, Sequences, genomedirs,
       file = "./outputs/intermediate_step9.RData")
  
  cat("Successfully saved intermediate results to ./outputs/intermediate_step9.RData\n")
  cat("You can resume from this point by loading this file if needed.\n")
}, error = function(e) {
  cat("ERROR saving intermediate results:", e$message, "\n")
  stop("Failed at saving intermediate results step")
})

} # End of if (!resume_from_checkpoint) block

# After Step 9 (before tree building):
if (!resume_from_prediction) {
cat("\nStep 9.5: Early COG filtering (sequence-based)...\n")
tryCatch({
  # Get assembly identifiers for each COG
  truncCOGs <- lapply(MatchedSequences, function(x) sort(as.integer(gsub('^([0-9]+)_.*', '\\1', names(x)))))

  
  noParas <- sapply(truncCOGs, function(x) length(x) == length(unique(x)))
  inFiveOrMore <- sapply(truncCOGs, function(x) length(unique(x)) <= 4)
  
  # Early filtering (combine all desired filters)
  EarlyFiltered <- inFiveOrMore
  
  # Subset data before expensive tree building
  PreFilteredSequences <- MatchedSequences[EarlyFiltered]
  
  cat("COGs before early filtering:", length(MatchedSequences), "\n")
  cat("COGs after early filtering:", length(PreFilteredSequences), "\n")
  cat("Breakdown of early filters:\n")
  cat("  In 4 or less organisms:", sum(inFiveOrMore), "\n")
  cat("  Combined filters passed:", sum(EarlyFiltered), "\n")
}, error = function(e) {
  cat("ERROR in early filtering:", e$message, "\n")
  stop("Failed at early filtering step")
})

# Move functional annotation before tree building for better efficiency

cat("\nStep 10: Functional annotation...\n")
tryCatch({
  # Create a list to store all the annotations
  CogsAnnot <- vector("list", length=length(PreFilteredSequences))
  names(CogsAnnot) <- names(PreFilteredSequences)
  
  # Translate the sequences
  cat("Translating sequences...\n")
  geneSeqs <- lapply(PreFilteredSequences, translate)
  
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

cat("\nStep 11: Filtering COGs (annotation-based)...\n")
tryCatch({
  # Apply annotation-based filters
  
  # Make sure COGs are coding elements
  codingCOGs <- sapply(CogsAnnot, function(x) is(x, 'Taxa'))
  
  # At least one high confidence annotation
  highConf <- sapply(CogsAnnot, function(x) 
                     if(is(x, 'Taxa')) 
                       max(sapply(x, function(y) 
                                  y$confidence[length(y$confidence)])) > 50
                     else FALSE
                     )

  # Check for valid functional annotations (no NAs or poor annotations)
  validAnnots <- sapply(CogsAnnot, function(x) {
    if(!is(x, 'Taxa')) return(FALSE)  # Must be a Taxa object
    
    # Get all annotations
    annots <- sapply(x, function(y) y$taxon[length(y$taxon)])
    annots <- annots[annots != 'unclassified_Root']
    
    # Must have at least one good annotation
    return(length(annots) > 0)
  })
  
  # Apply annotation-based filters
  FilteredCOGs <- codingCOGs & highConf & validAnnots
  
  # Subset our data for final filtered set
  FilteredSequences <- PreFilteredSequences[FilteredCOGs]
  FilteredAnnots <- CogsAnnot[FilteredCOGs]
  
  # Print summary of filtering
  cat("COGs after early filtering (Step 9.5):", length(PreFilteredSequences), "\n")
  cat("COGs after annotation filtering:", length(FilteredSequences), "\n")
  cat("Breakdown of annotation filters:\n")
  cat("  Coding elements:", sum(codingCOGs), "\n")
  cat("  High confidence annotation:", sum(highConf), "\n")
  cat("  Valid annotations:", sum(validAnnots), "\n")
  
  if (length(FilteredSequences) == 0) {
    stop("ERROR: No COGs passed all filtering criteria!")
  }
  
}, error = function(e) {
  cat("ERROR in COG filtering:", e$message, "\n")
  stop("Failed at COG filtering step")
})

cat("\nStep 12: Creating consensus annotations...\n")
tryCatch({
  consAnnots <- vector('character', length=length(FilteredAnnots))
  for (i in seq_along(FilteredAnnots)) {
    taxaentry <- FilteredAnnots[[i]]
    
    # If no annotation, it's a noncoding gene
    if (!is(taxaentry, 'Taxa')) {
      consAnnots[i] <- 'NONCODING'
    }
    # Otherwise it's a coding gene
    else {
      # Extract all annotations with better error handling
      all_annots <- character(0)
      
      tryCatch({
        # Get annotations from each sequence in the taxa entry
        for (j in seq_along(taxaentry)) {
          if (length(taxaentry[[j]]$taxon) > 0) {
            # Get the most specific annotation (last element)
            annotation <- taxaentry[[j]]$taxon[length(taxaentry[[j]]$taxon)]
            if (!is.null(annotation) && !is.na(annotation) && annotation != "") {
              all_annots <- c(all_annots, annotation)
            }
          }
        }
      }, error = function(e) {
        cat("Warning: Error extracting annotations for COG", i, ":", e$message, "\n")
      })
      
      # Filter annotations - keep unclassified only if it has meaningful text
      filtered_annots <- character(0)
      for (annot in all_annots) {
        # Keep annotation if:
        # 1. It's not unclassified_Root (too generic)
        # 2. It's meaningful text (more than just "unclassified")
        # 3. It's not empty or NA
        
        if (!is.na(annot) && annot != "" && annot != "unclassified_Root") {
          # Allow unclassified annotations if they have additional meaningful text
          if (grepl("unclassified", annot, ignore.case = TRUE)) {
            # Keep if it has more than just "unclassified" (e.g., "unclassified Bacteria")
            if (nchar(annot) > 12 && !grepl("^unclassified$", annot, ignore.case = TRUE)) {
              filtered_annots <- c(filtered_annots, annot)
            }
          } else {
            # Keep all non-unclassified annotations
            filtered_annots <- c(filtered_annots, annot)
          }
        }
      }
      
      # Assign consensus annotation
      if (length(filtered_annots) == 0) {
        consAnnots[i] <- 'Uncharacterized'
      } else {
        # Take the most common meaningful annotation
        annot_table <- table(filtered_annots)
        consAnnots[i] <- names(sort(annot_table, decreasing=TRUE))[1]
      }
      
      # Additional debugging for first few COGs
      if (i <= 3) {
        cat("COG", i, "debug:\n")
        cat("  Raw annotations found:", length(all_annots), "\n")
        cat("  Filtered annotations:", length(filtered_annots), "\n")
        cat("  Final annotation:", consAnnots[i], "\n")
        if (length(filtered_annots) > 0) {
          cat("  Top annotations:", paste(head(names(sort(table(filtered_annots), decreasing=TRUE)), 3), collapse=", "), "\n")
        }
      }
    }
  }
  
  # Summary of annotation results
  cat("Created consensus annotations for", length(consAnnots), "COGs\n")
  annotation_summary <- table(consAnnots)
  cat("Annotation summary:\n")
  for (annot_type in names(sort(annotation_summary, decreasing=TRUE))) {
    cat("  ", annot_type, ":", annotation_summary[annot_type], "\n")
  }
  
}, error = function(e) {
  cat("ERROR creating consensus annotations:", e$message, "\n")
  stop("Failed at consensus annotation step")
})

cat("\nStep 13: Building phylogenetic trees (final filtered set)...\n")
tryCatch({
  # Build phylogenetic trees for each COG (using final filtered set)
  # Create a list to store all the trees
  COGTrees <- vector("list", length=length(FilteredSequences))
  names(COGTrees) <- names(FilteredSequences)
  
  cat("Building trees for", length(FilteredSequences), "COGs (final filtered set)\n")
  
  # Loop through each COG and build a tree
  for (i in seq_along(FilteredSequences)) {
    cat("Building tree", i, "of", length(FilteredSequences), "\n")
    
    tryCatch({
      # Get the current COG sequences
      currentCOG <- FilteredSequences[[i]]
      
      # Since these are coding regions, align using AlignTranslation
      alignedCOG <- AlignTranslation(currentCOG)
      
      # Build a tree with ancestral state reconstruction
      # Using a reasonable time limit to balance accuracy and runtime
      distMatrix <- DistanceMatrix(alignedCOG)
      COGTrees[[i]] <- TreeLine(myDistMatrix=distMatrix, 
                               myXStringSet=alignedCOG,
                               method = "UPGMA",
                               reconstruct = TRUE,
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

cat("\nStep 14: Saving trees, sequences and annotations...\n")
tryCatch({
  save(FilteredSequences, FilteredAnnots, COGTrees, consAnnots, file = "./outputs/COG_trees_sequences_annotations.RData")
  cat("Successfully saved trees, sequences and annotations to output/\n")
}, error = function(e) {
  cat("ERROR saving trees, sequences and annotations:", e$message, "\n")
  stop("Failed at saving step")
})

} # End of if (!resume_from_prediction) block

cat("\nStep 15: Creating EvoWeaver object...\n")
tryCatch({
  pw <- EvoWeaver(COGTrees)
  cat("EvoWeaver object created successfully\n")
}, error = function(e) {
  cat("ERROR creating EvoWeaver object:", e$message, "\n")
  stop("Failed at EvoWeaver creation step")
})

cat("\nStep 16: Making predictions...\n")
tryCatch({
  preds <- predict(pw, Method="Jaccard")
  cat("Predictions completed successfully\n")
  print(preds)
}, error = function(e) {
  cat("ERROR making predictions:", e$message, "\n")
  stop("Failed at prediction step")
})

cat("\nStep 17: Clustering analysis...\n")
tryCatch({
  # Find clusters of coevolving COGs using igraph
  library(igraph)
  
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

cat("\nStep 18: Saving final results...\n")
tryCatch({
  # Create results directory if it doesn't exist
  if (!dir.exists("./results")) {
    dir.create("./results", recursive = TRUE)
  }
  
  save(preds, clusters, clusterLabels, consAnnots, file = "./results/EvoWeaver_results.RData")
  cat("Final results saved successfully\n")
}, error = function(e) {
  cat("ERROR saving final results:", e$message, "\n")
  stop("Failed at saving final results step")
})

cat("\nStep 19: Creating COG-to-cluster CSV mapping...\n")
tryCatch({
  # Create a data frame with COG information and cluster assignments
  cog_cluster_df <- data.frame(
    COG_ID = character(0),
    Cluster_Number = integer(0),
    Consensus_Annotation = character(0),
    stringsAsFactors = FALSE
  )
  
  # Get cluster membership for each COG
  cluster_membership <- membership(clusters)
  
  # Create the mapping
  for (i in seq_along(cluster_membership)) {
    cog_id <- names(FilteredSequences)[i]
    cluster_num <- cluster_membership[i]
    annotation <- consAnnots[i]
    
    cog_cluster_df <- rbind(cog_cluster_df, data.frame(
      COG_ID = cog_id,
      Cluster_Number = cluster_num,
      Consensus_Annotation = annotation,
      stringsAsFactors = FALSE
    ))
  }
  
  # Sort by cluster number then by COG ID
  cog_cluster_df <- cog_cluster_df[order(cog_cluster_df$Cluster_Number, cog_cluster_df$COG_ID), ]
  
  # Write to CSV
  write.csv(cog_cluster_df, file = "./results/COG_cluster_mapping.csv", row.names = FALSE)
  cat("Successfully created COG-to-cluster CSV mapping with", nrow(cog_cluster_df), "COGs\n")
  cat("CSV saved to /results/COG_cluster_mapping.csv\n")
}, error = function(e) {
  cat("ERROR creating COG-to-cluster CSV:", e$message, "\n")
  stop("Failed at CSV creation step")
})

cat("\nStep 20: Generating final report...\n")
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
cat("Results saved to /results/EvoWeaver_results.RData\n")
cat("COG-to-cluster mapping saved to /results/COG_cluster_mapping.csv\n")