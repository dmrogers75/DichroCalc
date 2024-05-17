#!/bin/bash

# ~/bin/dcPhiPsiHdc13.py has input arguments:
# (phi, psi) dihedral angles file; 
# DichroCalc Hamiltonian file;
# Euclidean cut-off distance, default is to use 0.0 to ignore;
# Inter-amide couplings for diamide 2d (-135 deg, 135 deg), default is 1 for 2ds4.

 pdbs=('1ppn')

for i in "${pdbs[@]}"; do
 echo "$i"
 python3 ~/bin/dcPhiPsiHdc13.py $i.pdb.dihedralsPhiPsi.txt $i.ham 0.0 1
 echo "$i done"
done

exit 0

