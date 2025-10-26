#!/bin/bash

# inputs:
input_filename='1J05.xyz'         # the name of the structure
output_filename='1J05'            # sequence
terminal_positions='Nter Cter'    # it can be Nter or Cter
bridges='GGGSLPET
         GGGSGGGSLPET
         GGGSGGGSGGGSGGGSLPET
         GGGSGGGSGGGSGGGSGGGSGGGSLPET
         GGGSGGGSGGGSGGGSGGGSGGGSGGGSGGGSLPET
         GGGSGGGSGGGSGGGSGGGSGGGSGGGSGGGSGGGSGGGSLPET
         GGGSGGGSGGGSGGGSGGGSGGGSGGGSGGGSGGGSGGGSGGGSGGGSLPET
         GGGSGGGSGGGSGGGSGGGSGGGSGGGSGGGSGGGSGGGSGGGSGGGSGGGSGGGSLPET
         GGGSGGGSGGGSGGGSGGGSGGGSGGGSGGGSGGGSGGGSGGGSGGGSGGGSGGGSGGGSGGGSLPET
         GGGSGGGSGGGSGGGSGGGSGGGSGGGSGGGSGGGSGGGSGGGSGGGSGGGSGGGSGGGSGGGSGGGSGGGSLPET
         GGGSGGGSGGGSGGGSGGGSGGGSGGGSGGGSGGGSGGGSGGGSGGGSGGGSGGGSGGGSGGGSGGGSGGGSGGGSGGGSLPET'
peptide='DSARGFKKPGKR'            # histag sequence
nconfs='50 100 150 200 250 300 350 400 450 500 550 600 650 700 750 800 850 900 950 1000 1500 2000 2500 3000 4000 5000' # number of configurations
#DSAQGFQQPGQQ
#DSARGFKKPGKR

# 1, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20
###################################################################################
# Check if the input file exists
if [ ! -f "$input_filename" ]; then
    echo "Error: Input file '$input_filename' not found!"
    exit 1
fi

# loop for Nter and Cter:
for terminal_position in $terminal_positions; do

    # loop for different bridge:
    for bridge in $bridges; do

        # create directory for each construction:
        if [ -n "$peptide" ]; then
            mkdir ${output_filename}_position-${terminal_position}_bridge-${bridge}_peptide-${peptide}/
            cd ${output_filename}_position-${terminal_position}_bridge-${bridge}_peptide-${peptide}/
        else
            mkdir ${output_filename}_position-${terminal_position}_bridge-${bridge}/
            cd ${output_filename}_position-${terminal_position}_bridge-${bridge}/
        fi

        # copy fortran chain generator script:
        cp -r ../random/ .

        # run scripts:
        # for polymer
        python3 ../get_polymers.py ../$input_filename $output_filename $terminal_position $bridge $peptide
        
        # for polymer-protein
        python3 ../get_protein-polymer.py ../$input_filename $output_filename $terminal_position "$nconfs" $bridge $peptide
        
        # go back:
        cd ../
        
    done
done
###################################################################################

exit

