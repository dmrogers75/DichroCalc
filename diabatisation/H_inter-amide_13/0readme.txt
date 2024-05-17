Files to modify a DichroCalc Hamiltonian using inter-amide couplings from diabatisation of monomer and dimer transition properties for 13 diamide geometries (2a to 2m).


For example for papain (PDB entry: 1ppn).


Script rdcdihedrals.sh executes DichroCalc's dihedrals script to generate phi and psi dihedral angles for an inputed PDB file.

bash rdcdihedrals.sh 1ppn.pdb

Generates 1ppn.pdb.dihedrals.txt and 1ppn.pdb.dihedralsPhiPsi.txt.


Script rdcphipsiHdc.sh executes dcPhiPsiHdc13.py to modify a DichroCalc Hamiltonian (e.g. 1ppn.ham) using the (phi, psi) angles (e.g. 1ppn.pdb.dihedralsPhiPsi.txt) and the set of 13 inter-amide couplings.

bash rdcphipsiHdc.sh > rdcphipsiHdc.1ppn.out

Generates 1ppn.ham.Hdc13_0.0_2ds4.txt (the modified Hamiltonian).
