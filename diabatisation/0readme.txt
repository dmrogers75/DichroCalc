Description of 'diabatisation' directory contents.


coupling.f90 - Fortran source code to perform diabatisation.

matrix.cpp - Modified DichroCalc C++ code to read in a Hamiltonian matrix from file 'new.mat'. Replace file in DichroCalc source code at path "DichroCalc/dichrocalc/src/matrix.cpp". 

rplcHam.py - Python script to replace nearest neighbour inter-amide interactions in Hamiltonian.


coupling_TDM_input.txt - Example input file (for diamide 2a) for executable of coupling.f90.

Execute coupling.exe using command: 
./coupling.exe coupling_TDM_input.txt


av_20A_p180p180.ham - Example DichroCalc Hamiltonian (for model peptide 2a).
diabat_2a.txt - Inter-amide couplings from diabatisation for diamide 2a.
av_20A_p180p180.ham.Rep.txt - Modified DichroCalc Hamiltonian for model peptide 2a (output from rplcHam.py).

Execute rplcHam.py using command:
python rplcHam.py av_20A_p180p180.ham diabat_2a.txt
