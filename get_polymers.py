#!/bin/python3

import sys
import os
import subprocess
import sys

# Uncomment the line below if you need to use custom functions from a specific directory.
from my_functions import *

#######################################################################
# Inputs:
# This script takes several command-line arguments:
# 1. input_filename: The name of the input file containing the protein structure.
# 2. output_filename: The base name for the output files.
# 3. terminal_position: The terminal position of the protein, either 'Nter' (N-terminus) or 'Cter' (C-terminus).
# 4. bridge: The sequence of the bridge to be used in the construction.
# 5. peptide (optional): The peptide sequence, if any. If not provided, it's set to None.

input_filename    = sys.argv[1]
output_filename   = sys.argv[2]
terminal_position = sys.argv[3]
sequence          = sys.argv[4]

# Cutoff distance for filtering polymers based on proximity to the protein structure.
cutoff_distance = 0.38













#######################################################################
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

# obtain the protein structure oriented based on the selected terminal position.
protein = get_protein_relative(input_filename, position)

# generate the polymers based on the sequence.
polymers = get_blocks_from_polymers_relative(sequence)
print(f'npolymers = {len(polymers)}')

# save the generated polymers to an XYZ file.
get_block2xyz(polymers, 'polymers.xyz')

# filter the polymers based on their distance to the protein using the specified cutoff distance.
filtered_polymers = get_filtered_polymers(protein, polymers, cutoff_distance)
print(f'filtered npolymers = {len(filtered_polymers)}')

# Save the filtered polymers to a separate XYZ file.
get_block2xyz(filtered_polymers, 'filtered_polymers.xyz')

exit()

