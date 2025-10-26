import os
import numpy as np
import copy
import random

# This mapping is for the peptide sequence to insert
aa_mapping = {
    "A": "Ala",
    "R": "Arg",
    "N": "Asn",
    "D": "Asp",
    "C": "Cys",
    "Q": "Gln",
    "E": "Glu",
    "G": "Gly",
    "H": "His",
    "I": "Ile",
    "L": "Leu",
    "K": "Lys",
    "M": "Met",
    "F": "Phe",
    "P": "Pro",
    "S": "Ser",
    "T": "Thr",
    "W": "Trp",
    "Y": "Tyr",
    "V": "Val"
}

# Inverted amino acid mapping
inverted_aa_mapping = {v: k for k, v in aa_mapping.items()}

#############################################################################
# functions:
#############################################################################

def get_protein_relative(input_filename, position):
    """
    Reads a protein structure from a file (XYZ),
    repositions all atoms relative to a specified bead (Nter or Ceter),
    and returns the relative coordinates.

    Parameters:
    - input_filename (str): Path to the input file containing protein coordinates.
    - position (int): Index of the bead to use as the reference point for repositioning.

    Returns:
    - protein (list): List of repositioned coordinates for the protein.
    """
        
    ##############################
    # read file:
    ##############################
    skiprows = 2                            # skip first to rows, typical for xyz from molecular theory
    with open(input_filename, "r") as file:
        lines = file.readlines()[skiprows:] # real all files
    
        # get coords from file:
        block = []     # empty list to store configurations
        for line in lines:
            line_xyz = line.strip().split()  # line: restype, x, y, z
            
            # only modified the coods:
            line_xyz[1] = float(line_xyz[1]) # make float X
            line_xyz[2] = float(line_xyz[2]) # make float Y
            line_xyz[3] = float(line_xyz[3]) # make float Z
            block.append(line_xyz)           # append the line
    
    protein = copy.deepcopy(block)              # copy and define name of block
    terminal = copy.deepcopy(protein[position]) # protein[position] mush be 0 or -1
    
    ##############################
    # relative the coords:
    ##############################
    for line in protein:
        # only modified the coods:
        line[1] = round(line[1] - terminal[1], 4)  # x
        line[2] = round(line[2] - terminal[2], 4)  # y
        line[3] = round(line[3] - terminal[3], 4)  # z
    
    # return protein:
    return protein




def get_blocks_from_polymer_xyz(input_filename):
    """
    Reads polymer structures from a file (XYZ) and returns a list of blocks, where each block contains
    coordinates of atoms belonging to a single polymer segment.

    Parameters:
    - input_filename (str): Path to the input file containing polymer coordinates.

    Returns:
    - blocks (list): List of blocks, each containing coordinates of a polymer segment.
    """
    
    # read file:
    with open(input_filename, 'r') as file:

        # elimina espacios vacios y saltos en cada linea (\n), elimina la segunda linea de polymers:
        lines = []
        for line in file:                   # recorre línea por línea el archivo
            if line.strip():                # strip() elimina espacios y saltos si 'line' no está vacía
                lines.append(line.strip())  # guarda la línea ya 'limpia'
        
        # create empty list:
        blocks = []                         # list where all polymer conformations (blocks) will be stored
        current_block = []                  # temporary collector for a single polymer conformation
        
        # main loop that read all lines:
        for line in lines:
        
            # definimos una variable:
            line_value = line.strip().split()[0]
            
            ##################################################
            # # primera linea con de "numero de segmentos"
            ##################################################
            #if line_value != 'C' and not current_block:
            if line_value != 'C' and len(current_block) == 0:
                #current_block = []
                continue  # salta a la siguiente linea
            
            ##################################################
            # n linea de "numero de segmentos"
            ##################################################
            #elif line_value != 'C' and current_block:
            elif line_value != 'C' and len(current_block) > 0:
                blocks.append(current_block)           
                current_block = []
            
            ##################################################
            # append coords:
            ##################################################                                            
            else:
                # append coords with line_value = 'C'
                # linea de "coordenada xyz"
                # C   -2.8238004418837050       0.41636354969003631        2.5085438199303338
                line_xyz = line.strip().split()
                line_xyz[1] = float(line_xyz[1])
                line_xyz[2] = float(line_xyz[2])
                line_xyz[3] = float(line_xyz[3]) 
                current_block.append(line_xyz)
        
        # esta linea guarda el ultimo bloque:
        if current_block:
            blocks.append(current_block)
    
    return blocks




def get_blocks_from_protein_xyz(input_filename, comment):
    """
    Reads protein structures from a file and returns a list of blocks, where each block contains
    coordinates of atoms belonging to a single protein segment.

    Parameters:
    - input_filename (str): Path to the input file containing protein coordinates.
    - comment (str): The comment line used to separate blocks in the file.

    Returns:
    - blocks (list): List of blocks, each containing coordinates of a protein segment.
    """
    
    with open(input_filename, 'r') as file:  # read a file
        # definimos una variable:
        #lines = [line for line in file if not line.split()[0][0].isdigit()]
        lines = []
        for line in file:
            stripped = line.strip()
            if not stripped:
                continue # ignora lineas vacias
            
            parts = stripped.split()
            # si el primer carácter del primer token no es un dígito
            if not parts[0][0].isdigit():
                lines.append(stripped)
        
        ##############################################
        # get 
        ##############################################
        # get empty lists:
        blocks        = []
        current_block = []
        
        # main loop that read all lines:
        for line in lines:
            
            ##################################################
            # crea un nuevo bloque
            ##################################################
            #if comment in line and not current_block:
            if comment in line and len(current_block) == 0:
                #current_block = []
                continue
            
            ##################################################
            # crea un nuevo bloque
            # cuando encuentra una linea que no es 'C' 
            ##################################################
            #elif comment in line and current_block:
            elif commebt in line and len(current_block) > 0:
                blocks.append(current_block)                   
                current_block = []
            
            ##################################################
            # guarda las coordenadas
            ##################################################                          
            else:
                line_xyz = line.strip().split()
                line_xyz[1] = float(line_xyz[1])
                line_xyz[2] = float(line_xyz[2])
                line_xyz[3] = float(line_xyz[3]) 
                current_block.append(line_xyz)
        
        ##########################
        # esta linea guarda el ultimo bloque
        ##########################
        if current_block:
            blocks.append(current_block)
    
    return blocks  




def a_to_nm(block):
    """
    Converts atomic coordinates from angstroms to nanometers for a given block of coordinates.

    Parameters:
    - block (list): List of coordinates to be converted.

    Returns:
    - block (list): List of coordinates converted to nanometers.
    """
    for line in block:
        line[1] = round(line[1]/10,4)
        line[2] = round(line[2]/10,4)
        line[3] = round(line[3]/10,4)
    
    return block




def get_blocks_from_polymers_relative(sequence):
    """
    Generates polymer chains based on a given sequence and returns them as a list of blocks.

    Parameters:
    - sequence (str): The sequence of monomers to generate the polymer.

    Returns:
    - polymers (list): List of blocks, each containing coordinates of a polymer segment.
    """
    
    # get the lenght of the chain to generate:
    chain_length = len(sequence) + 1   # +1 in order to concatenate the te chain
    
    # get random:
    os.system('cd random/ && bash compi.sh')
    os.system(f'cd random/ && echo {chain_length} | ./polymer.x')
    
    # get blocks of polymers:
    polymers = get_blocks_from_polymer_xyz('random/polymer.xyz')
    
    # eliminate the firt residue from lopymer:
    # remeber that the firt elemnt is start in zero ('O')
    for polymer in polymers:
        a_to_nm(polymer)   # convert from amstron to nm
        polymer.pop(0)     # eilinate the first line

    return polymers




def get_filtered_polymers(protein, polymers, cutoff_distance):
    """
    Filters polymer segments by removing those that overlap with the protein based on a cutoff distance.

    Parameters:
    - protein (list): List of coordinates for the protein.
    - polymers (list): List of polymer segments to filter.
    - cutoff_distance (float): Minimum allowed distance between polymer and protein coordinates.

    Returns:
    - filtered_polymers (list): List of polymer segments that do not overlap with the protein.
    """
    
    # create a empty list to store the filtered polymers:
    filtered_polymers = []
    
    # loop for all polymers:
    for polymer in polymers:
        
        ###############################################################
        # for one polymer get all distances betwen all residues on protein and polymer 
        ###############################################################
        
        # get empty list:
        dxyz_list = []
        
        for coord_protein in protein:                  # over each residue of the protein
            coord1 = np.array(coord_protein[1:])       # select only xyz
            
            for coord_polymer in polymer:              # over each residue of the polymer
                coord2 = np.array(coord_polymer[1:])   # select only xyz
                dxyz = np.linalg.norm(coord1 - coord2) # get distane
                dxyz_list.append(dxyz)                 # append
        
        ######################################
        # condition to store the polymer configuration:
        #if all(element >= cutoff_distance for element in dxyz_list):
        #    filtered_polymers.append(polymer)
        # si algun elemento del bloque no satisface la condition no se agrega
        append_block = True
        for distance in dxyz_list:
            if distance < cutoff_distance:
                append_block = False
                break
        
        if append_block: 
            filtered_polymers.append(polymer)
        
    # return:
    return filtered_polymers



def add_sequence(molecule, sequence, start_sequence, end_sequence):
    """
    Adds a peptide sequence to a given molecule.

    Parameters:
    - molecule (list): List of coordinates representing the molecule.
    - sequence (str): Peptide sequence to add to the molecule.
    - start_sequence (int): Starting index in the molecule to begin adding the sequence.
    - end_sequence (int): Ending index in the molecule to stop adding the sequence.

    Returns:
    - molecule (list): Molecule with the peptide sequence added.
    """
    
    for i in range(start_sequence, end_sequence):
        index = i - start_sequence
        molecule[i][0] = aa_mapping[sequence[index]]

    return molecule



def add_terminals(molecule):
    """
    Adds N-terminal and C-terminal labels to the first and last beads of a molecule.

    Parameters:
    - molecule (list): List of coordinates representing the molecule.

    Returns:
    - molecule (list): Molecule with N-terminal and C-terminal labels added.
    """
    molecule[0][0] = inverted_aa_mapping[molecule[0][0]] + '_Nt'
    molecule[-1][0] = inverted_aa_mapping[molecule[-1][0]] + '_Ct'

    return molecule



def get_merge_all(protein, filtered_polymers, terminal_position, sequence):
    """
    Merges protein and polymer segments into a single structure, adding peptide sequences
    and terminal labels as necessary.

    Parameters:
    - protein (list): List of coordinates representing the protein.
    - filtered_polymers (list): List of filtered polymer segments.
    - terminal_position (str): Position of the terminal to use ('Nter' or 'Cter').
    - sequence (str): Peptide sequence to add to the merged structure.

    Returns:
    - molecule_list (list): List of merged structures.
    """
    
    # list to store the configurations:
    molecule_list = []     
    
    # conditional for Nter:
    if terminal_position == 'Nter':
        for polymer in filtered_polymers:
            polymer = copy.deepcopy(polymer)[::-1]
            polymer_protein = copy.deepcopy(polymer) + copy.deepcopy(protein)
            polymer_protein = add_sequence(polymer_protein, sequence[::-1], 0, len(sequence))
            polymer_protein = add_terminals(polymer_protein)
            molecule_list.append(polymer_protein)
    
    # conditional for Cter: 
    elif terminal_position == 'Cter':
        for polymer in filtered_polymers:
            protein_polymer = copy.deepcopy(protein) + copy.deepcopy(polymer)
            protein_polymer = add_sequence(protein_polymer, sequence, len(protein), len(protein) + len(sequence))
            protein_polymer = add_terminals(protein_polymer)
            molecule_list.append(protein_polymer)
                                             
    else:
        print('It must be provided an option to include a peptide')

    return molecule_list



def get_merge_random(protein_merge_all, nconf):
    """
    Randomly selects a specified number of merged protein-polymer structures from a list.

    Parameters:
    - protein_merge_all (list): List of merged protein-polymer structures.
    - nconf (int): Number of structures to select.

    Returns:
    - protein_merge_random (list): Randomly selected structures.
    """
    seed_value = 0
    random.seed(seed_value)
    protein_merge_random = random.sample(protein_merge_all, nconf)
         
    return protein_merge_random


       
def save_results(terminal_position, molecules, output_filename, bridge, peptide, nconf):
    """
    Saves the coordinates and sequences of molecules to .xyz and .seq files.

    Parameters:
    - terminal_position (str): Position of the terminal to use ('Nter' or 'Cter').
    - molecules (list): List of molecules to save.
    - output_filename (str): Base name for the output files.
    - bridge (str): Bridge sequence (if any) used in the molecule.
    - peptide (str): Peptide sequence (if any) used in the molecule.
    - nconf (int): Configuration number for the output files.
    """
    if peptide and bridge:
        output_filename_xyz = f'{output_filename}_{terminal_position}_bridge-{bridge}_peptide-{peptide}_conf-{nconf}.xyz'
        output_filename_seq = f'{output_filename}_{terminal_position}_bridge-{bridge}_peptide-{peptide}_conf-{nconf}.seq'
    elif peptide:
        output_filename_xyz = f'{output_filename}_{terminal_position}_peptide-{peptide}_conf-{nconf}.xyz'
        output_filename_seq = f'{output_filename}_{terminal_position}_peptide-{peptide}_conf-{nconf}.seq'
    elif bridge:
        output_filename_xyz = f'{output_filename}_{terminal_position}_bridge-{bridge}_conf-{nconf}.xyz'
        output_filename_seq = f'{output_filename}_{terminal_position}_bridge-{bridge}_conf-{nconf}.seq'
    else:
        print('An option must be provided')

    with open(output_filename_xyz, 'w') as file:
        n_beads = len(molecules[0])
        for molecule in molecules:
            file.write(f'{n_beads}\nresultado_salida\n')
            for line in molecule:
                file.write('\t'.join(map(str, line)) + '\n')

    with open(output_filename_seq, 'w') as file:
        file.write('resultado_salida\n')
        index = 0
        for element in molecule:
            amino_acid = element[0]
            index += 1
            file.write(f'{index}\t{amino_acid}\n')



def get_block2xyz(blocks, output_filename):
    """
    Saves a list of blocks to an .xyz file.

    Parameters:
    - blocks (list): List of blocks to save.
    - output_filename (str): Name of the output .xyz file.
    """
    
    with open(output_filename, 'w') as file:
        for block in blocks:
            nseg = len(block)
            file.write(f'{nseg}\nresultado_salida\n') # add first and second line
            for line in block:
                file.write(f'{line[0]}\t{line[1]}\t{line[2]}\t{line[3]}\n')



###########################################################################
# function definitions:
def splitxyz(input_filename, comment, output_basename):
    """
    Splits a file containing multiple .xyz structures into individual files.

    Parameters:
    - input_filename (str): Path to the input file containing multiple .xyz structures.
    - comment (str): Comment line that separates the structures.
    - output_basename (str): Base name for the output files.

    Returns:
    - None
    """
    blocks = []
    current_block = []
    
    with open(input_filename, 'r') as file:
        lines = [line for line in file if len(line.split()) != 1]
        
        for line in lines:
            pattern2_to_test = line.strip()
            if pattern2_to_test == comment and not current_block:
                current_block = []
            elif pattern2_to_test == comment and current_block:
                blocks.append(current_block)                   
                current_block = []
            else:
                line_read = line.strip().split()
                current_block.append(line_read[:4])

        if current_block:
            blocks.append(current_block)

    for count, block in enumerate(blocks):
        output = f'{output_basename}-{count:02d}.xyz'
        with open(output, 'w') as file:
            file.write(f'{len(block)}\n{comment}\n')
            for line in block:
                file.write(f'{line[0]}\t{line[1]}\t{line[2]}\t{line[3]}\n')



def get_thetas(input_filename, comment, n1, n2):
    """
    Calculates the theta angle between two points for each block in a protein file.

    Parameters:
    - input_filename (str): Path to the input file containing protein coordinates.
    - comment (str): Comment line that separates the structures.
    - n1 (int): Index of the first point.
    - n2 (int): Index of the second point.

    Returns:
    - nblocks (int): Number of blocks processed.
    - thetas (list): List of calculated theta angles.
    """
    blocks = get_blocks_from_protein_xyz(input_filename, comment)
    thetas = []
    
    for block in blocks:
        delx = block[n2][1] - block[n1][1]
        dely = block[n2][2] - block[n1][2]
        delz = block[n2][3] - block[n1][3]
        delr = np.sqrt(delx**2 + dely**2 + delz**2)
        theta = np.arccos(delz / delr)
        thetas.append(theta)
                
    return len(blocks), thetas 



def get_costhetas(input_filename, comment, n1, n2):
    """
    Calculates the cosine of the theta angle between two points for each block in a protein file.

    Parameters:
    - input_filename (str): Path to the input file containing protein coordinates.
    - comment (str): Comment line that separates the structures.
    - n1 (int): Index of the first point.
    - n2 (int): Index of the second point.

    Returns:
    - nblocks (int): Number of blocks processed.
    - costhetas (list): List of calculated cosine theta values.
    """
    blocks = get_blocks_from_protein_xyz(input_filename, comment)
    costhetas = []
    
    for block in blocks:
        delx = block[n2][1] - block[n1][1]
        dely = block[n2][2] - block[n1][2]
        delz = block[n2][3] - block[n1][3]
        delr = np.sqrt(delx**2 + dely**2 + delz**2)
        costheta = np.float64(delz) / np.float64(delr)
        costhetas.append(costheta)
                
    return len(blocks), costhetas 



def get_phis(input_filename, comment, n1, n2):
    """
    Calculates the phi angle between two points for each block in a protein file.

    Parameters:
    - input_filename (str): Path to the input file containing protein coordinates.
    - comment (str): Comment line that separates the structures.
    - n1 (int): Index of the first point.
    - n2 (int): Index of the second point.

    Returns:
    - nblocks (int): Number of blocks processed.
    - phis (list): List of calculated phi angles.
    """
    blocks = get_blocks_from_protein_xyz(input_filename, comment)
    phis = []
    
    for block in blocks:
        delx = block[n2][1] - block[n1][1]
        dely = block[n2][2] - block[n1][2]
        delz = block[n2][3] - block[n1][3]
        phi = np.arctan2(dely, delx)
        phis.append(phi)
                
    return len(blocks), phis 



def calculate_distances_polymer(input_filename, n1, n2):
    """
    Calculates the distance between two points for each block in a polymer file.

    Parameters:
    - input_filename (str): Path to the input file containing polymer coordinates.
    - n1 (int): Index of the first point.
    - n2 (int): Index of the second point.

    Returns:
    - nblocks (int): Number of blocks processed.
    - distances (list): List of calculated distances.
    """
    blocks = get_blocks_from_polymer_xyz(input_filename)
    distances = []
    
    for block in blocks:
        position1 = block[n1][1:]
        position2 = block[n2][1:]
        distance = np.linalg.norm(np.array(position2) - np.array(position1))
        distances.append(distance)

    return len(blocks), distances



def calculate_distances_protein(input_filename, comment, n1, n2):
    """
    Calculates the distance between two points for each block in a protein file.
    
    Parameters:
    - input_filename (str): Path to the input file containing protein coordinates.
    - comment (str): Comment line that separates the structures.
    - n1 (int): Index of the first point.
    - n2 (int): Index of the second point.

    Returns:
    - nblocks (int): Number of blocks processed.
    - distances (list): List of calculated distances.
    """
    blocks = get_blocks_from_protein_xyz(input_filename, comment)
    distances = []
    
    for block in blocks:
        position1 = block[n1][1:]
        position2 = block[n2][1:]
        distance = np.linalg.norm(np.array(position2) - np.array(position1))
        distances.append(distance)

    return len(blocks), distances









