import sys
import os
import csv

if len(sys.argv) < 3:
    #No input file specified
    print("ERROR! Usage: python3 gene_freq.py <input_file> <number_of_samples> [frequency_threshold]")
elif len(sys.argv) == 3: 
    #Set input file
    infile_path = sys.argv[1]

    #Set number of samples
    n = int(sys.argv[2])
    
    #Default to 0.8 if unspecified
    freq_thresh = 0.8
else:
    #Set input file
    infile_path = sys.argv[1]

    #Set number of samples
    n = int(sys.argv[2])

    try:
        #Set to second argument
        freq_thresh = float(sys.argv[3])
        #Default to 0.8 if invalid
        if freq_thresh > 1:
            freq_thresh = 0.8
    except:
        #Default to 0.8 if invalid
        freq_thresh = 0.8

inpath = infile_path.split('/')[0]
infile = infile_path.split('/')[1]

r = csv.reader(open(infile_path, 'r'))
gene_counts = {rows[0]: rows[1] for rows in r}
for k in gene_counts:
    gene_counts[k] = int(gene_counts[k])

freq_genes = {key: gene_counts[key] for key in gene_counts if gene_counts[key] / n >= freq_thresh}

outfile = inpath + "/gene_thresh_" + str(freq_thresh) + ".csv"

w = csv.writer(open(outfile, 'w'))
for k, v in freq_genes.items():
    w.writerow([k, v])
