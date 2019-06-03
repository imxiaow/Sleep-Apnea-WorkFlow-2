import sys
import os
from random import sample
import time

#Check for correct usage
try:
    len(sys.argv) == 3
except:
    print("ERROR! Usage: python3 gene_list_subset.py <input_filename> <proportionate_size_of_subset>")
    sys.exit()

#Save input filename and subset proportion
infile = sys.argv[1]
prop = sys.argv[2]

#Initialize full list
genes = []

#Open input file in read mode
with open(infile, 'r') as f:
    #Iterate through input file lines
    for line in f:
        #Remove newlines
        curr = line[:-1]

        #Append line to gene list
        genes.append(curr)

#Close input file
f.close()

#Separate infile into filename and extension
inname = infile[:-4]
inext = infile[-4:]

#Compute number of items, m, in subset
n = len(genes)
m = round(n * float(prop))

#Create subset list
subset = sample(genes, m)

#Create output directory if it doesn't exist
outdir = "samples/"
if not os.path.exists(outdir):
    os.mkdir(outdir)

#Build output file path and name
#   (Timestamp to ensure unique name, regardless of proportion overlap)
outfile = outdir + inname + "_" + prop + "subset_" + str(time.time()) + inext

#Open output file
with open(outfile, 'w') as f:
    #Iterate through genes in subset and write to output file
    for gene in subset:
        f.write("%s\n" % gene)

#Close output file
f.close()
