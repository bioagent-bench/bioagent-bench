#!/bin/bash

./scan10/nextflow ./scan10/src/scan10.nf \
    -profile docker \
    -resume \
    --fastq "/home/dev/benchmark/bio-agent-benchmark/scan/data/cells-dataset/filtered" \
    --quantif cellranger \
    --version 106 \
    --species mouse \
    --chemistry V2 \
    --filtergtf \
    --min_feature_RNA 500 \
    --min_cells 3 \
    --max_percent_mito 5 \