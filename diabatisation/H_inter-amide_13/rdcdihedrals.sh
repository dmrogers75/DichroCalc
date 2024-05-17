#!/bin/bash

# Input $1 is pdb file

# Path to DichroCalc scripts in ~/bin dir to execute dihedrals script
/home/$USER/bin/DichroCalc/scripts/dihedrals $1 > $1.dihedrals.txt

# Extract phi and psi angles only
sed '1,5d' $1.dihedrals.txt | awk '{print $3, $4}' | sed '$d' > $1.dihedralsPhiPsi.txt
