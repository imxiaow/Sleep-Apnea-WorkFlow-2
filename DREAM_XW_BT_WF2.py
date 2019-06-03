import sys
import os
import csv
import shutil
from BioLink.biolink_client import BioLinkWrapper
import pandas as pd
from os import makedirs
from html3.html3 import XHTML
from Modules.Mod1A_functional_sim import FunctionalSimilarity
from Modules.Mod1B1_phenotype_similarity import PhenotypeSimilarity
from Modules.Mod1E_interactions import GeneInteractions

#======================================================================================================
pyptha = sys.executable.split('/')
pyptha[-2]= 'lib'
pypth='/'.join(pyptha) + '*/site-packages'

# Hack to get around problematic updating of distutils installed PyYAML and a slightly older pandas requiring a compatible numpy
shutil.rmtree(pypth + '/PyYAML*', ignore_errors=True)
shutil.rmtree(pypth + '/numpy*', ignore_errors=True)

sys.path.append("../mvp-module-library")


#======================================================================================================
#Set frequency threshold
if len(sys.argv) == 3:
    #Default to 0.8 if unspecified
    freq_thresh = 0.8
    # Default
    threshold_Mod1A = 0.7
    # Default
    threshold_Mod1B = 0.7
else:
    try:
        #Set to first argument
        freq_thresh = float(sys.argv[1])
        #Default to 0.8 if invalid
        if freq_thresh > 1:
            freq_thresh = 0.8

        threshold_Mod1A = float(sys.argv[2])
        if threshold_Mod1A > 1:
            threshold_Mod1A = 0.7

        threshold_Mod1B = float(sys.argv[3])
        if threshold_Mod1B > 1:
            threshold_Mod1B = 0.7

    except:
        #Default to 0.8 if invalid
        freq_thresh = 0.8
        # Default
        threshold_Mod1A = 0.7
        # Default
        threshold_Mod1B = 0.7


#===== DEFINED FUNCTIONS ===============================================================================
def output_file(tag, title, ext):
    basepath = "./Tidbit/" + tag
    filename = title.replace(" ", "_")
    filepath = basepath + "/" + filename + "." + ext
    makedirs(basepath, exist_ok=True)
    output = open(filepath, "w+")
    output.info = {'tag': tag, 'title': title}
    return output


def dump_html(output, body):
    title = output.info['title'] + " for " + output.info['tag']

    doc = XHTML()

    doc.head.title(title)
    doc.body.h1(title)
    doc.body.p(body.to_html())

    output.write(str(doc))


def load_genes(model, data, threshold):
    # Module specification
    inputParameters = {
        'input': data,
        'parameters': {
            'taxon': 'human',
            'threshold': threshold,
        },
    }

    # Load the computation parameters
    model.load_input_object(inputParameters)
    model.load_gene_set()


def similarity(model, data, threshold, input_disease_symbol, module, title):
    # Initialize
    load_genes(model, data, threshold)
    model.load_associations()

    # Perform the comparison
    results = model.compute_similarity()

    # Process the results
    results_table = pd.DataFrame(results)
    results_table = results_table[
        ~results_table['hit_id'].isin(alternative_hit_ID_check)].sort_values('score', ascending=False)
    results_table['module'] = module

    # save the gene list to a file under the "Tidbit" subdirectory
    output = output_file(input_disease_symbol, title, "html")
    dump_html(output, results_table)
    output.close()

    return results_table


def interactions(model, data, threshold, input_disease_symbol, module, title):
    # Initialize
    load_genes(model, data, threshold)

    # Perform the comparison
    results = model.get_interactions()

    # Process the results
    results_table = pd.DataFrame(results)
    results_table = results_table[
        ~results_table['hit_id'].isin(alternative_hit_ID_check)].sort_values('score', ascending=False)
    results_table['module'] = module

    # save the gene list to a file under the "Tidbit" subdirectory
    output = output_file(input_disease_symbol, title, "html")
    dump_html(output, results_table)
    output.close()

    return results_table


#======================================================================================================
#Set output directory name
input_disease_symbol = "Sleep_Apnea_Output"

#Initialize lists to hold sample symbol/id pairs
samples = []
#Initialize lists to hold sample ids
sample_ids = []

#Set sample directory
sampdir = os.fsencode("samples/")

#======================================================================================================
for sampfile in os.listdir(sampdir):
    #Initialize lists to hold current sample information
    currsamp = []
    currsamp_ids = []

    #Build filename
    fname = "samples/" + os.fsdecode(sampfile)

    #Open file for reading
    with open(fname, 'r') as f:
        #Iterate through genes in sample
        for line in f:
            #Split into symbol and id
            sym_id = line[:-1].split(',')

            #Add gene to current sample lists
            currgene = {'hit_id': sym_id[1], 'hit_symbol': sym_id[0]}
            currsamp.append(currgene)
            currsamp_ids.append(sym_id[1])

    #Close file
    f.close()

    #Append sample info to sample lists
    samples.append(currsamp)
    sample_ids.append(currsamp_ids)


#=== Repeat Mod1A and Mod1B for each sample file =========================================================
outgenes = []
#Initialize dictionary of counts for each gene
gene_counts = {}

#Iterate through sample lists
for i in range(len(samples)):
    print("Sample " + str(i + 1) + "/" + str(len(samples)) + "...")

    #Set input set
    input_curie_set = samples[i]
    #Set ID check
    alternative_hit_ID_check = sample_ids[i]
    #===========================
    # Functional Similarity using Jaccard index threshold
    func_sim_human = FunctionalSimilarity()
    Mod1A_results = similarity( func_sim_human, input_curie_set, threshold_Mod1A, input_disease_symbol, 'Mod1A', "Functionally Similar Genes" )
    #Get list of IDs from Mod1E_results
    Mod1A_ids = Mod1A_results['hit_id'].values.tolist()
    print("\n" + "Mod1A Functional Similarity results: (threshold of " + str(threshold_Mod1A) + ")\n" + str(Mod1A_results) + "\n")

    #===========================
    # Phenotypic simulatiry using OwlSim calculation threshold
    pheno_sim_human = PhenotypeSimilarity()
    Mod1B_results = similarity( pheno_sim_human, input_curie_set, threshold_Mod1B, input_disease_symbol, 'Mod1B', "Phenotypically Similar Genes" )
    #Get list of IDs from Mod1E_results
    Mod1B_ids = Mod1B_results['hit_id'].values.tolist()
    print("\n" + "Mod1B Phenotype Similarity results: (threshold of " + str(threshold_Mod1B) + ")\n" + str(Mod1B_results) +"\n")

    #Add each unique ID to output genes list
    outgenes.append(list(set(Mod1A_ids)))
    outgenes.append(list(set(Mod1B_ids)))

    #===========================
    # Checking if there's intersections between Mod1A_results AND Mod1B_results.
    # If there's a gene result from both module, get its hit_id and input_ID etc.
    intersect_Mod1A_Mod1B = set(Mod1A_ids) & set(Mod1B_ids)
    print("\n" + "Intersect between Mod1A and Mod1B: \n" +str(intersect_Mod1A_Mod1B) + "\n")

'''
    if (intersect_Mod1A_Mod1B != None):
        if(len(intersect_Mod1A_Mod1B) ==1):
            intersect_input_id = Mod1A_results['input_id'].values.tolist()[Mod1A_results['hit_id'].values.tolist().index(intersect_Mod1A_Mod1B)]
            intersect_input_symbol = Mod1A_results['input_symbol'].values.tolist()[Mod1A_results['hit_id'].values.tolist().index(intersect_Mod1A_Mod1B)]
            print("\n" + "Input ID for intesect hit: \n" + str(intersect_input_id) "\n" + "Input Symbol for intesect hit: \n" + str(intersect_input_symbol)+ "\n")

        if(len(intersect_Mod1A_Mod1B) !=1):
            intersect_input_id_lst = []
            intersect_input_symbol_lst=[]
            for item in intersect_Mod1A_Mod1B:
                intersect_input_id = Mod1A_results['input_id'].values.tolist()[Mod1A_results['hit_id'].values.tolist().index(intersect_Mod1A_Mod1B)]
                intersect_input_id_lst.append(intersect_input_id)
                intersect_input_symbol = Mod1A_results['input_symbol'].values.tolist()[Mod1A_results['hit_id'].values.tolist().index(intersect_Mod1A_Mod1B)]
                intersect_input_symbol_lst.append(intersect_input_symbol)
            print("\n" + "Input ID for intesect hit: \n" + str(intersect_input_id_lst) + "\n" + "Input Symbol for intesect hit: \n" + str(intersect_input_id_lst) + "\n")
'''

    #===========================
#    # Gene-gene interactions
#    gene_gene_interaction = GeneInteractions()
#    Mod1E_results = interactions( gene_gene_interaction, input_curie_set, 0, input_disease_symbol, 'Mod1E', "Genes-Gene interactions" )
#    #Get list of IDs from Mod1E_results
#    Mod1E_ids = Mod1E_results['hit_id'].values.tolist()
    #print("Mod1E Genes-Gene interactions results: (threshold N.A.)\n"+str(Mod1E_results) +"\n")




#======================================================================================================
#Iterate through output gene samples
for samp in outgenes:
    #Iterate through genes in current sample
    for gene in samp:
        #Add gene to gene_counts dictionary or increment count if already exists
        if gene in gene_counts:
            gene_counts[gene] += 1
        else:
            gene_counts[gene] = 1

#Get only genes that appear in sample outputs with a frequency at or above the threshold
freq_genes = {key: gene_counts[key] for key in gene_counts if gene_counts[key] / len(samples) >= freq_thresh}

#Create output directory if it doesn't already exist
if not os.path.exists(input_disease_symbol):
    os.mkdir(input_disease_symbol)

#Build output filename
outfile_all = input_disease_symbol + "/full_gene_counts.csv"
outfile_freq = input_disease_symbol + "/gene_thresh_" + str(freq_thresh) + ".csv"

#Write frequent gentes to output CSV
f = csv.writer(open(outfile_all, 'w'))
for k, v in gene_counts.items():
    f.writerow([k, v])

g = csv.writer(open(outfile_freq, 'w'))
for k, v in freq_genes.items():
    g.writerow([k, v])



#======================================================================================================
# END
