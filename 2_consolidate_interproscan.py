# Need to run in a python env
import sys
from collections import defaultdict
import re

# Input and output filenames
input_filename = "pangenome_interpro.tsv"  # Replace with your actual filename
output_filename = "consolidated_pangenome_interpro.tsv"

# Dictionary to store gene names and their GO terms
gene_go_dict = defaultdict(set)

# Open the input file for reading
with open(input_filename, 'r') as infile:
    for line in infile:
        columns = line.strip().split('\t')  # Split the line by tab
        gene_name = columns[0]  # First column is the gene name
        
        # Extract and clean GO terms, removing content inside parentheses
        for col in columns:
            if 'GO:' in col:
                go_terms = re.findall(r'GO:\d+', col)  # Extract only GO terms without extra text
                gene_go_dict[gene_name].update(go_terms)

# Write the output to a file
with open(output_filename, 'w') as outfile:
    outfile.write("Gene\tGO_Terms\n")  # Write the header
    for gene, go_terms in gene_go_dict.items():
        go_terms_combined = "|".join(sorted(go_terms))  # Join GO terms with a pipe (|) and sort them
        outfile.write(f"{gene}\t{go_terms_combined}\n")

