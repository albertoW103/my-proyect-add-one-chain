#!/bin/python3

import sys
import os

# Uncomment the line below if you need to use custom functions from a specific directory.
from my_functions import *

#######################################################################
# Inputs:
# This script takes several command-line arguments:
# 1. input_filename: The name of the input file containing the protein structure.
# 2. output_filename: The base name for the output files.
# 3. terminal_position: The terminal position of the protein, either 'Nter' (N-terminus) or 'Cter' (C-terminus).
# 4. nconfs: A list of integers representing the number of configurations to generate.
# 5. bridge: The sequence of the bridge to be used in the construction.
# 6. peptide (optional): The peptide sequence, if any. If not provided, it's set to None.

input_filename = sys.argv[1]
output_filename = sys.argv[2]
terminal_position = sys.argv[3]
nconfs = [int(x) for x in sys.argv[4].split()]
bridge = sys.argv[5]
if len(sys.argv) > 6:
    peptide = sys.argv[6]
else:
    peptide = None

#######################################################################
# Generate the complete sequence:
# If a peptide sequence is provided, it is appended to the bridge sequence.
if peptide:
    sequence = bridge + peptide
else:
    sequence = bridge

# Load filtered polymers:
# The function `get_blocks_polymer` reads in polymers from the 'filtered_polymers.xyz' file.
filtered_polymers = get_blocks_polymer('filtered_polymers.xyz')

# Determine the terminal position in the protein structure:
# The 'Nter' option corresponds to the N-terminus (index 0),
# and the 'Cter' option corresponds to the C-terminus (index -1).
if terminal_position == 'Nter':
    position = 0
elif terminal_position == 'Cter':
    position = -1
else:
    print('An option must be provided for terminal_position')
    exit(1)

# Obtain the protein structure oriented based on the selected terminal position.
protein = get_protein(input_filename, position)

# Merge the protein with all filtered polymers:
# This merges the protein with the filtered polymers, positioning them based on the terminal position and sequence.
protein_polymer_all = get_merge_all(protein, filtered_polymers, terminal_position, sequence)
print(f'protein and polymer merged all = {len(protein_polymer_all)}')

# Save the merged protein-polymer structure to an XYZ file.
get_block2xyz(protein_polymer_all, 'protein_polymer_all.xyz')

# Generate random configurations:
# For each value in `nconfs`, generate a random configuration by merging the protein with the polymers.
for nconf in nconfs:
    protein_merge_random = get_merge_random(protein_polymer_all, nconf)
    print(f'protein and polymer merged random = {len(protein_merge_random)}')

    # Save the random configuration results.
    save_results(terminal_position, protein_merge_random, output_filename, bridge, peptide, nconf)

exit()

