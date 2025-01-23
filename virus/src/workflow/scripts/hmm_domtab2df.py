import sys
import numpy as np
import pandas as pd
from Bio import SearchIO
import argparse
import multiprocessing
import pickle

# Create an argument parser object
parser = argparse.ArgumentParser(description='Example script')

# Add the -in and -out options
parser.add_argument('-in', dest='input_file', help='Input file', required=True)
parser.add_argument('-out', dest='output_file', help='Output file', required=True)
# Add the --threads option
parser.add_argument('-t', dest='num_threads', type=int, help='Number of threads', required=True)

# Parse the command line arguments
args = parser.parse_args()

# Get the input file name and output file name from the arguments
input_file = args.input_file
output_file = args.output_file
# Get the number of threads from the arguments
num_threads = args.num_threads

d = []
with open(input_file, 'r') as fin:
    for record in SearchIO.parse(fin, "hmmsearch3-domtab"):
        for hit in record.hits:
            for hsp in hit.hsps:
                d.append(
                    {
                        'contig_id': hit.id,
                        'hmm': record.id ,
                        # 'evalue': hit.bitscore,
                        'hit_start': hsp.hit_start,
                        'hit_end': hsp.hit_end,
                        # 'hit_length': (hsp.hit_end - hsp.hit_start)*3
                    })

hmm_df = pd.DataFrame(d)
hmm_df['contig_id'] = hmm_df['contig_id'].map(lambda x: str(x)[:-2])

contig_hmm_dict = {}
for index, row in hmm_df.iterrows():
    if row['contig_id'] not in contig_hmm_dict:
             contig_hmm_dict[row['contig_id']] = [[row['hit_start'], row['hit_end']]]
    elif row['contig_id'] in contig_hmm_dict:
             contig_hmm_dict[row['contig_id']].append([row['hit_start'], row['hit_end']])

all_contig_hmm_dict = list([contig_hmm_dict])
contig_hmm_dict

merge_dict = {}

for index in contig_hmm_dict:
    contig = index
    arr = contig_hmm_dict[index]
    # print(index,": ", arr)
    arr.sort(key = lambda x: x[0])

    # array to hold the merged intervals
    m = []
    s = -10000
    max = -100000
    for i in range(len(arr)):
        a = arr[i]
        if a[0] > max:
            if i != 0:
                m.append([s,max])
            max = a[1]
            s = a[0]
        else:
            if a[1] >= max:
                max = a[1]

    #'max' value gives the last point of
    # that particular interval
    # 's' gives the starting point of that interval
    # 'm' array contains the list of all merged intervals

    if max != -100000 and [s, max] not in m:
        # interval = max - s
        m.append([s, max])
    # print("\n", index, ":", length, end = " ")
    mergeinterval = 0
    for i in range(len(m)):
        # print(m[i], end = " ")
        interval = m[i][1]-m[i][0]
        # print(interval, end = " ")
        mergeinterval += interval
    # print("mergeinterval", mergeinterval, ratio, end = " ")
    # times 3 for converting the amino acid seq lengths to nucleotide seq lengths
    merge_dict[contig] = mergeinterval*3
    
# hmm_merge_hit_length_df = pd.DataFrame.from_dict(merge_dict, orient = 'index', columns = ['merged_hit_length'])
# hmm_merge_hit_length_df.to_csv(output_file, sep='\t')

with open(output_file, 'wb') as f:
    pickle.dump(merge_dict, f)

