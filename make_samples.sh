#!/bin/bash

# Check for required input file argument
if [[ $# < 2 ]]; then
    echo "ERROR! Usage: ./make_samples <input_filename> <replicates> [list_of_subset_proportions]"
    exit 1
fi

#Get gene list file
gene_list=$1

# Set number of replicates at each subset size
reps=$2

# Create samples
for i in $(seq 1 $reps); do
	# If subset proportion(s) are specified
	if [[ $# > 2 ]]; then
		# Create subset for each proportion
		for arg in ${@:3}; do
			python3 gene_list_subset.py $gene_list $arg
		done
	# If no subset proportion is specified
	else
		# Default to 1.0
		python3 gene_list_subset.py $gene_list 1.0
	fi
done
