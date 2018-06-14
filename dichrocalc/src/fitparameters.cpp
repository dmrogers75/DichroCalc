// #################################################################################################
//
//  Program:      fitparameters.cpp
//
//  Function:     Part of DichroCalc:
//                Fits the respective parameters (monopoles and dipole moments) to the position of
//                each chromophore
//
//  Author:       Benjamin M. Bulheller
//
//  Affiliation:  Research group Jonathan D. Hirst
//                University of Nottingham, UK
//
//                Research group Shaul Mukamel
//                University of California at Irvine, USA
//
//  Version:      $Revision: 4620 $, $Date: 2009-07-06 09:09:38 +0100 (Mon, 06 Jul 2009) $
//
//  Date:         May 2009
//
// #################################################################################################


#include "../include/dichrocalc.h"


// ================================================================================


int Dichro::FitParameters ( void )
{
	int  AtomNumParSet, AtomNumGroup, CurAtomIndex, AtomIndex, GroupTransNum;
	int  k, l, Atom, Group, Type, Trans, Mono, Coord, FirstTransition;
	bool ChargeTransfer, There;
	unsigned int Index;
	vector<int>  AtomIndices; // to keep track of which atoms have already been added to the vector
	
	int    NumberOfGroups = DC_Input.Chromophores.Type.size();
	double AverageFitDistance = 0;
	double LargestFitDistance = 0;
	double FitDistance = 0;
	int    AtomNumber = 0;
	string ParSetName, Origin;
	
	if (DC_Verbose) Dichro::NewTask ( "Fitting Parameters" );
	
	if (DC_Debug > 2) {
		Dichro::NewFileTask (DC_DbgFile, "Fitting Parameters to Chromophore Atoms");
		fprintf (DC_DbgFile, "Number of chromophores: %d\n", NumberOfGroups);
	}
	
	// only needed if DC_PrintXyzFiles is set true in the constructor
	FILE* CoordinatesParSet;
	FILE* CoordinatesGroup;
	
	// DEBUG OUTPUT
	if (DC_PrintXyzFiles) {
		string CoordsParSet   = DC_InFileBaseName + ".par.xyz";
		string CoordsGroup    = DC_InFileBaseName + ".pdb.xyz";
		CoordinatesParSet = fopen (CoordsParSet.c_str(), "w");
		CoordinatesGroup  = fopen (CoordsGroup.c_str(), "w");
	}
	
	// initialize some household variables in DC_System
	DC_System.NumberOfAtoms = 0;
	DC_System.NumberOfGroups = NumberOfGroups;
	DC_System.NumberOfTransitions = 0;
	DC_System.MatrixDimension = 0; // this is the same as the number of transitions
	
	// run over each chromophore in the system
	for (Group = 0; Group < NumberOfGroups; Group++) {
		// the matrices are (row, colum) order like so:
		// atom1x, atom1y, atom1z
		// atom2x, atom2y, atom2z
		// ...
		
		// the current chromophore type, e.g. 1 for a peptide group
		Type = DC_Input.Chromophores.Type.at(Group);
		
		ParSetName = DC_ParSets.at(Type).Name;
		
		// initialize the vector for this particular group
		SystemGroup CurGroup;
		
		
		// --------------------------------------------------------------------------------
		// Prepare the parameter set atoms
		// --------------------------------------------------------------------------------
		
		// the number of atoms in the parameter set group and the chromophore
		AtomNumParSet = DC_ParSets.at(Type).Atoms.size();
		CurGroup.NumberOfAtoms = AtomNumParSet;
		CurGroup.ParameterSet  = ParSetName;

		// this variable counts all parameter set atoms, i.e. not the number of atoms in the PDB!
		DC_System.NumberOfAtoms += AtomNumParSet;
		
		// the position vector for moving the parameter set atoms to the origin
		vector<double> PosVecParSet = DC_ParSets.at(Type).Reference;
		
		// create a vector with the coordinates of the parameter set atoms
		vector< vector<double> > CoordParSet;
		vector< vector<double> > CoordParSetOrigin (AtomNumParSet, vector<double>(3, 0.0));
		
		for (Atom = 0; Atom < AtomNumParSet; Atom++) {
			CoordParSet.push_back (DC_ParSets.at(Type).Atoms.at(Atom).Coord);
			
			// Translate each coordinate about the respective reference coordinate.
			// This moves the parameter set atoms (as well as later on the chromophore atoms)
			// to the origin, this is needed for the fitting.
			for (Coord = 0; Coord < 3; Coord++)
				CoordParSetOrigin.at(Atom).at(Coord) = CoordParSet.at(Atom).at(Coord)
				                                         - PosVecParSet.at(Coord);
		}
		
		// DEBUG OUTPUT
		// printf ("\n\n%s\n", ParSetName.c_str());
		
		if (DC_Debug > 4) {
			Dichro::NewFileTask (DC_FitFile, "Chromophore " + tostring(Group));
			
			fprintf (DC_FitFile, "   Original position of the %d atoms in the %s parameter set:",
			                     AtomNumParSet, DC_ParSets.at(Type).Name.c_str() );
			for ( Atom = 0; Atom < AtomNumParSet; Atom++ ) {
				fprintf (DC_FitFile, "\n      Atom %2d:  ", Atom);
				
				for ( Coord = 0; Coord < 3; Coord++ )
					fprintf (DC_FitFile, "    %12.6f", CoordParSet.at(Atom).at(Coord));
			}
			
			fprintf (DC_FitFile, "\n\n   Position vector: %12.6f    %12.6f    %12.6f\n",
			         PosVecParSet.at(0), PosVecParSet.at(1), PosVecParSet.at(2));
			fprintf (DC_FitFile, "\n");
			
			fprintf (DC_FitFile, "   Translated to the origin:");
			
			for ( Atom = 0; Atom < AtomNumParSet; Atom++ ) {
				fprintf (DC_FitFile, "\n      Atom %2d:  ", Atom);
				
				for ( Coord = 0; Coord < 3; Coord++ )
					fprintf (DC_FitFile, "    %12.6f", CoordParSetOrigin.at(Atom).at(Coord));
			}
			
			fprintf (DC_FitFile, "\n\n");
		}
		
		
		// --------------------------------------------------------------------------------
		// Prepare the chromophore atoms
		// --------------------------------------------------------------------------------
		
		AtomNumGroup  = DC_Input.Chromophores.Atoms.at(Group).size();
		
		vector<int> GroupAtomIndices = DC_Input.Chromophores.Atoms.at(Group);
		
		// The atom indices are the array indeces of the group's atoms in the DC_Coordinate
		// array. In DC_System they are only needed to check for overlapping groups quickly
		// (rather than using the atom coordinates).
		CurGroup.AtomIndices = GroupAtomIndices;
//		PrintVector (&GroupAtomIndices);
//		PrintVector (&CurGroup.AtomIndices);
		
		// the position vector for moving the group atoms to the origin
		vector<double> PosVecGroup (3, 0.0);  // initialize it with three zeros
		
		// create a vector with the coordinates of the parameter set atoms
		vector< vector<double> > CoordGroup;
		vector< vector<double> > CoordGroupOrigin (AtomNumGroup, vector<double>(3, 0.0));
		
		// add up the xyz coordinates
		for (Atom = 0; Atom < AtomNumGroup; Atom++) {
			// the atom index of the group, e.g. the "3" in "3  4  6"
			AtomIndex = DC_Input.Chromophores.Atoms.at(Group).at(Atom);
			CoordGroup.push_back (DC_Input.Coordinates.Groups.at(AtomIndex));
			
			for (Coord = 0; Coord < 3; Coord++)
				PosVecGroup.at(Coord) += DC_Input.Coordinates.Groups.at(AtomIndex).at(Coord);
		}
		
		// divide every coordinate by the number of atoms (similar to the center of mass)
		for (Coord = 0; Coord < 3; Coord++)
			PosVecGroup.at(Coord) /= AtomNumGroup;
		
		if (DC_Debug > 4) {
			fprintf (DC_FitFile, "\n   Original position of the %d chromophore atoms:", AtomNumGroup);
			for ( Atom = 0; Atom < AtomNumGroup; Atom++ ) {
				fprintf (DC_FitFile, "\n      Atom %2d:  ", Atom);
				
				for ( Coord = 0; Coord < 3; Coord++ )
					fprintf (DC_FitFile, "    %12.6f", CoordGroup.at(Atom).at(Coord));
			}
		}
		
		for (Atom = 0; Atom < AtomNumGroup; Atom++) {
			// Translate each coordinate about the respective reference coordinate.
			// This moves the parameter set atoms (as well as later on the chromophore atoms)
			// to the origin, this is needed for the fitting.
			for (Coord = 0; Coord < 3; Coord++)
				CoordGroupOrigin.at(Atom).at(Coord) = CoordGroup.at(Atom).at(Coord)
				                                      - PosVecGroup.at(Coord);
		}
		
		if (DC_Debug > 4) {
			fprintf (DC_FitFile, "\n\n   Position vector: %12.6f    %12.6f    %12.6f\n\n",
			         PosVecGroup.at(0), PosVecGroup.at(1), PosVecGroup.at(2));
			fprintf (DC_FitFile, "   Translated to the origin:");
			for ( Atom = 0; Atom < AtomNumGroup; Atom++ ) {
				fprintf (DC_FitFile, "\n     Atom %2d:   ", Atom);
				
				for ( Coord = 0; Coord < 3; Coord++ )
					fprintf (DC_FitFile, "    %12.6f", CoordGroupOrigin.at(Atom).at(Coord));
			}
			fprintf (DC_FitFile, "\n\n");
		}
		
		// --------------------------------------------------------------------------------
		
		if (AtomNumParSet != AtomNumGroup) {
			cerr << "\nERROR: The number of atoms in chromophore " << Group + 1 << " ("
			     << AtomNumGroup << ") does not match\n"
			     << "       the number of atoms in the chromophore type "
			     << DC_Input.Parameters.Name.at(Group) << " (" << AtomNumParSet << ").\n\n";
			
			DC_Error = "Wrong number of atoms in assigned parameter set.";
			DC_ErrorCode = 140;
			return 140;
		}
		
		// the matrix with the coordinates of the parameter set atoms
		Matrix ParSetMatrix (AtomNumParSet, 3);   ParSetMatrix = 0.0;  // initialize to zero
		// the matrix with the coordinates of the chromophore read from the input file
		Matrix GroupMatrix  (AtomNumGroup, 3);    GroupMatrix = 0.0;  // initialize to zero
		
		for (Atom = 0; Atom < AtomNumParSet; Atom++) { // the rows are the atoms
			for (Coord = 0; Coord < 3; Coord++) {          // the cols are the xyz coordinates
				// M(r,c) starts counting at 1, M.element(r,c) starts counting at 0
				ParSetMatrix.element(Atom, Coord) = CoordParSetOrigin.at(Atom).at(Coord);
				GroupMatrix.element(Atom, Coord)  = CoordGroupOrigin.at(Atom).at(Coord);
			}
		}
		
		Matrix RotMatrixUnitary (3, 3);
		Matrix RotMatrixNonUnitary (3, 3);
		RotationMatrix (Group, ParSetMatrix, GroupMatrix, &RotMatrixNonUnitary, &RotMatrixUnitary);
		
		ParSetName = DC_ParSets.at(Type).Name;                // parset name as string
		ChargeTransfer = DC_ParSets.at(Type).ChargeTransfer;  // CT? - true/false
		CurGroup.ChargeTransfer = ChargeTransfer;
		
		GroupTransNum = DC_Input.Parameters.Trans.at(Type);         // number of transitions
		CurGroup.NumberOfTransitions = GroupTransNum;
		// count all transitions of the whole system
		DC_System.NumberOfTransitions += GroupTransNum;
		
		if (DC_Debug > 4) {
			Dichro::OutputFileSeparator (DC_FitFile);
			Dichro::OutputFileHeadline (DC_FitFile, "Performing fitting of the parameters:");
			fprintf (DC_FitFile, "   Name of parameter set:   %s\n", ParSetName.c_str());
			fprintf (DC_FitFile, "   Number of transitions:   %d\n", GroupTransNum);
			fprintf (DC_FitFile, "\n   Rotation Matrix (non-unitary):\n");
			FilePrintMatrix (DC_FitFile, &RotMatrixNonUnitary);
			fprintf (DC_FitFile, "\n   Rotation Matrix (unitary):\n");
			FilePrintMatrix (DC_FitFile, &RotMatrixUnitary);
			fprintf (DC_FitFile, "\n");
		}
		
		int State;             // excitation from ground state to pi* is state 0
		int StateTransNum;     // the number of transitions to be read from the current state
		
		ParSet CurParSet;                // a shortcut to the current parameter set
		ParSetTrans CurTransParSet;      // a shortcut to the current transition
		
		CurParSet = DC_ParSets.at(Type);
		
		AtomNumber += CurParSet.NumberOfAtoms; // count the atoms of all groups
		
		// Rotation of the atoms. This is not needed for the actual calculation but
		// to calculate the fitting accuracy/fitting errors/fitting deviations (pick one).
		
		// // rotate the group position vector
		for (k = 0; k < 3; k++){
			for (l = 0; l < 3; l++) {
				PosVecGroup.at(k) = PosVecGroup.at(k) -
				                  ( RotMatrixUnitary.element(l, k) * CurParSet.Reference.at(l) );
			}
		}
		
		// initialize a vector element for this group of atoms
		vector< vector<double> > GroupAtoms;
		CurGroup.Atoms = GroupAtoms;
		
		vector<SystemTransition> GroupTrans;
		CurGroup.Trans = GroupTrans;
		
		for (Atom = 0; Atom < CurParSet.NumberOfAtoms; Atom++) {
			// create a fresh vector for each atom
			vector<double> AtomCoords (3, 0.0);
			// rotate the atom around the origin
			Rotate (&CurParSet.Atoms.at(Atom).Coord, &AtomCoords, &RotMatrixUnitary);
			
			for (Coord = 0; Coord < 3; Coord++)
				// move it to the position of the chromophore
				AtomCoords.at(Coord) += PosVecGroup.at(Coord);
			
			FitDistance = PointDistance (&AtomCoords, &CoordGroup.at(Atom));
			AverageFitDistance = AverageFitDistance + FitDistance;
			
			if (FitDistance > LargestFitDistance) LargestFitDistance = FitDistance;
			
			// add the atom to the vector containing all atoms of the group
			CurGroup.Atoms.push_back (AtomCoords);
			
			CurAtomIndex = DC_Input.Chromophores.Atoms.at(Group).at(Atom);
			
			There = false;
			for (Index = 0; Index < AtomIndices.size(); Index++)
				if (AtomIndices.at(Index) == CurAtomIndex) There = true;
			
			if (not There) {
				DC_System.Atoms.push_back (AtomCoords);
				AtomIndices.push_back (CurAtomIndex);
			}
		}
		
		if (DC_Debug > 4) {
			fprintf (DC_FitFile, "   Group atoms to be matched:\n");
			for (Atom = 0; Atom < CurParSet.NumberOfAtoms; Atom++) {
				FilePrintCoord (DC_FitFile,
				                &DC_Input.Coordinates.Groups.at( GroupAtomIndices.at(Atom) ));
				
				if (DC_PrintXyzFiles)
					FilePrintCoord (CoordinatesGroup,
					                &DC_Input.Coordinates.Groups.at( GroupAtomIndices.at(Atom) ));
			}
			
			fprintf (DC_FitFile, "\n   Position vector:   \n");
			FilePrintCoord (DC_FitFile, &PosVecGroup);
			
			fprintf (DC_FitFile, "\n   Parameter set atoms before:\n");
			for (Atom = 0; Atom < CurParSet.NumberOfAtoms; Atom++)
				FilePrintCoord (DC_FitFile, &CurParSet.Atoms.at(Atom).Coord);
			
			fprintf (DC_FitFile, "\n   Parameter set atoms after:\n");
			for (Atom = 0; Atom < CurParSet.NumberOfAtoms; Atom++) {
				FilePrintCoord (DC_FitFile, &CurGroup.Atoms.at(Atom));
				
				if (DC_PrintXyzFiles)
					FilePrintCoord (CoordinatesParSet, &CurGroup.Atoms.at(Atom));
			}
			
			fprintf (DC_FitFile, "\n");
		}
		
		/*
		Pattern of how many transitions of each state are read in and what they represent.
		Example for NMA4FIT2, peptide bond with four transitions
		    4 trans in S0   3 in S1   2 in S2   1 in S3
		        n -> pi*    -\ -\ -\
		       pi -> pi*    -/  |  |    -\ -\
		      pib -> pi*       -/  |    -/  |     -\
		       n' -> pi*          -/       -/     -/
		
		Charge-Transfer Set:
		The four monomer transitions at the beginning are not necessarily in that order!
		
		 8 trans in S0    7 trans in S1      6 in S2     5 in S3    4 in S4    3 S5   2 S6  1 S7
		     n1 -> pi*    -\-\-\-\-\-\-\
		    pi1 -> pi*    -/ | | | | | |  -\-\-\-\-\-\
		     n2 -> pi*      -/ | | | | |  -/ | | | | |  -\-\-\-\-\
		    pi2 -> pi*        -/ | | | |    -/ | | | |  -/ | | | |  -\-\-\-\
		                         | | | |       | | | |     | | | |   | | | |
		    CT1                 -/ | | |      -/ | | |    -/ | | |  -/ | | |  -\-\-\
		    CT2                   -/ | |        -/ | |      -/ | |    -/ | |  -/ | |  -\-\
		    CT3                     -/ |          -/ |        -/ |      -/ |    -/ |  -/ |  -\
		    CT4                       -/            -/          -/        -/      -/    -/  -/
		*/
		
		bool Permanent = false;
		int CTBonus = 0;
		if (ChargeTransfer) GroupTransNum = GroupTransNum + 4;
		
		// GroupTransNum are required from state 0 (e.g. n->pi*, pi->pi*)
		// For N transitions, N+1 ("<=") states need to be processed (one more for the perm. moments)
		for (State = 0; State <= GroupTransNum; State++) {
			if (State < GroupTransNum + CTBonus)  { // if it is a transition
				// from each state one transition less is needed than from the state before
				StateTransNum = GroupTransNum - State;
				
				// if a specific backbone transition was requested and this is the first parameter
				// set type (i.e. the backbone parameters)
				if (DC_Input.Configuration.BBTrans > -1 and Type == 0) {
					// run from the requested transition and read only one transition (+1)
					FirstTransition = DC_Input.Configuration.BBTrans;
					StateTransNum = FirstTransition + 1;
				}
				// if a specific charge-transfer transition was requested and this is a CT chromophore
				else if (DC_Input.Configuration.CTTrans > -1 and ChargeTransfer) {
					// run from the requested transition and read only one transition (+1)
					FirstTransition = DC_Input.Configuration.CTTrans;
					StateTransNum = FirstTransition + 1;
				}
				else {
					FirstTransition = 0;
				}
			} // of if (State < GroupTransNum)  { // if it is a transition
			else { // if it is the permanent moments
				Permanent = true;
				// jump to the last state, which are the permanent moments
				State = CurParSet.States.size() - 1;
				
				// if a specific backbone transition was requested and this is the first parameter set
				// type (i.e. the backbone parameters)
				if (DC_Input.Configuration.BBTrans > -1 and Type == 0) {
					// read the first permanent moment (only one)
					FirstTransition = 0;
					StateTransNum   = 1;
				}
				// if a specific charge-transfer transition was requested and this is a CT chromophore
				else if (DC_Input.Configuration.CTTrans > -1 and ChargeTransfer) {
					// read the first permanent moment (only one)
					FirstTransition = 0;
					StateTransNum   = 1;
				}
				else {
					FirstTransition = 0;
					// for each transition the permanent moment is needed
					StateTransNum = GroupTransNum;
				}
			}
			
			// DEBUG OUTPUT
			// if (State == 0) printf ("\n");
			// printf ("State %d - reading %d transitions\n", State, StateTransNum);
			
			for (Trans = FirstTransition; Trans < StateTransNum; Trans++) {
				// create a string for the debug file to describe where the transition originated from
				// this string is created before the 'masking' of the transition number for CT groups
				Origin = ParSetName + " - State "      + tostring (State) +
				                      " - Transition " + tostring (Trans);
				
				// a shortcut to the original data from the parset
				CurTransParSet = CurParSet.States.at(State).at(Trans);
				// instantiate a new vector for the chromophore
				SystemTransition CurTransGroup;
				
				// DEBUG OUTPUT
				// if (not Permanent)
				// 	printf ("State %d,   Transition %d   Energy %8.3f  Wavelength %8.3f\n",
				// 	         State, Trans, CurTransParSet.Energy, CurTransParSet.Wavelength);
				// else
				// 	printf ("Permanent Moments,   Transition %d\n", Trans);
				// printf ("     %s\n", Origin.c_str());
				
				CurTransGroup.Origin     = Origin;
				CurTransGroup.Permanent  = Permanent;
				CurTransGroup.Energy     = CurTransParSet.Energy;
				CurTransGroup.Wavelength = CurTransParSet.Wavelength;
				
				// Note that there is no .at(Trans).
				// This takes the state separation out of the equation, i.e. State.Trans for 3 trans
				// 1.1, 1.2, 1.3 - 2.1, 2.2 - 3.1 - 4.1, 4.2, 4.3 (4. are the permanent moments)
				// is changed into sequentially numbered transitions:
				// 1    2    3     4    5     6     7    8    9   (7, 8, 9 are the perm. moments)
				// and these are accessed later by [Trans * (Trans+1)]/2 + 1 = 7
				
				// rotation of the transition dipole moments
				vector<double> EDM;
				Rotate (&CurTransParSet.EDM, &EDM, &RotMatrixUnitary);
				CurTransGroup.EDM = EDM;
				
				if (DC_Debug > 4) {
					fprintf (DC_FitFile, "   Elec. dipole moment before:  ");
					FilePrintCoord (DC_FitFile, &CurTransParSet.EDM, true);
					fprintf (DC_FitFile, "   Elec. dipole moment after:   ");
					FilePrintCoord (DC_FitFile, &EDM, true);
					fprintf (DC_FitFile, "\n");
				}
				
				if (not Permanent) {
					vector<double> MDM;
					Rotate (&CurTransParSet.MDM, &MDM, &RotMatrixUnitary);
					CurTransGroup.MDM = MDM;
					
					if (DC_Debug > 4) {
						fprintf (DC_FitFile, "   Mag. dipole moment before:   ");
						FilePrintCoord (DC_FitFile, &CurTransParSet.MDM, true);
						fprintf (DC_FitFile, "   Mag. dipole moment after:    ");
						FilePrintCoord (DC_FitFile, &MDM, true);
						fprintf (DC_FitFile, "\n");
					}
				}
				
				// all monopoles for this transition
				vector<ParSetMonopole> Monopoles;
				
				// rotation and translation of the monopoles
				for (Mono = 0; Mono < CurParSet.States.at(State).at(Trans).NumberOfMonopoles; Mono++) {
					ParSetMonopole Monopole;
					
					Rotate (&CurTransParSet.Monopoles.at(Mono).Coord, &Monopole.Coord, &RotMatrixUnitary);
					
					for (Coord = 0; Coord < 3; Coord++)
						Monopole.Coord.at(Coord) += PosVecGroup.at(Coord);
					
					Monopole.Charge = CurTransParSet.Monopoles.at(Mono).Charge;
					Monopoles.push_back (Monopole);
				}
				
				// the position (or reference) vector is later needed for the CD calculation
				CurGroup.Reference = PosVecGroup;
				
				CurTransGroup.Monopoles = Monopoles;
				CurTransGroup.NumberOfMonopoles =
						CurParSet.States.at(State).at(Trans).NumberOfMonopoles;
				
				if (not Permanent)
					CurGroup.Trans.push_back (CurTransGroup);
				else
					CurGroup.Perm.push_back (CurTransGroup);
			} // of for (Trans = 0; Trans < GroupTransNum; Trans++)
			
			// DEBUG OUTPUT
			// Dichro::OutputSystemGroupClass (&CurGroup);
		} // of for (State = 0; State < GroupTransNum; State++)
		
		if (DC_Debug > 2) {
			Dichro::OutputFileHeadline (DC_DbgFile, "   Chromophore " + tostring(Group));
			
			for (Atom = 0; Atom < CurGroup.NumberOfAtoms; Atom++)
				fprintf (DC_DbgFile, "      Atom %2d:       %12.4f %12.4f %12.4f\n",
				         Atom, CurGroup.Atoms.at(Atom).at(0), CurGroup.Atoms.at(Atom).at(1),
				               CurGroup.Atoms.at(Atom).at(2));
				
				fprintf (DC_DbgFile, "\n      Reference point: %10.4f   %10.4f   %10.4f\n",
				    CurGroup.Reference.at(0), CurGroup.Reference.at(1), CurGroup.Reference.at(2));
				fprintf (DC_DbgFile,   "      Fit distance:     %9.4f\n", FitDistance);
				
				fprintf (DC_DbgFile, "\n      Rotation Matrix (non-unitary):\n");
				FilePrintMatrix (DC_DbgFile, &RotMatrixNonUnitary);
				
				fprintf (DC_DbgFile, "\n      Rotation Matrix (unitary):\n");
				FilePrintMatrix (DC_DbgFile, &RotMatrixUnitary);
				fprintf (DC_DbgFile, "\n");
				
				if (Group < NumberOfGroups-1) // don't print separator after the last group
					Dichro::OutputFileSeparator (DC_DbgFile, 3);
		}
		
		DC_System.Groups.push_back (CurGroup);
	} // of for (Group = 0; Group < NumberOfGroups; Group++)
	
	// not really necessary, just for a better readability at some points
	DC_System.MatrixDimension = DC_System.NumberOfTransitions;
	
	// complete DC_System class is printed after the diagonalization
	// Dichro::OutputSystemClass ();
	
	AverageFitDistance = AverageFitDistance / AtomNumber;
	
	if (DC_Verbose) {
		printf ("   Average fit distance:  %8.3f\n",   AverageFitDistance);
		printf ("   Largest fit distance:  %8.3f\n", LargestFitDistance);
	}
	
	if (DC_PrintXyzFiles) {
		fclose (CoordinatesGroup);
		fclose (CoordinatesParSet);
	}
	
	return 0;
} // of Dichro::FitParameters


// ================================================================================


int Dichro::RotationMatrix ( int Group, Matrix ParSetMatrix, Matrix GroupMatrix,
                             Matrix* RotMatrixNonUnitary, Matrix* RotMatrixUnitary  )
// calculated the rotation matrix to turn the second set of atoms into the first one
{
	int Atom;
	int Cols = ParSetMatrix.ncols();
	int Rows = ParSetMatrix.nrows();
	
	//	printf ("RotationMatrix for Group %d\n", Group);
	
	// adds an additional atom, if a planar system is found
	Dichro::CheckPlanar ( &ParSetMatrix, &GroupMatrix );
	
	// ================================================================================
	
	if ( (GroupMatrix.ncols() != ParSetMatrix.ncols()) or
	      (GroupMatrix.nrows() != ParSetMatrix.nrows() ) ) {
		printf ("\nERROR: Matrix dimensions do not match:\n\n");
		printf (  "       Chromophore:   %2d x %-2d\n", GroupMatrix.nrows(), GroupMatrix.ncols());
		printf (  "       Parameter set: %2d x %-2d\n", ParSetMatrix.nrows(), ParSetMatrix.ncols());
		printf ("\n\n");
		DC_Error = "Error during parameter fitting";
		DC_ErrorCode = 143;
		return 143;
	}
	
	// cout << "Group matrix\n";
	// cout << setw(15) << GroupMatrix << "\n\n";
	// cout << "Parameter set matrix\n";
	// cout << setw(15) << ParSetMatrix << "\n\n";
	
	Matrix U (Rows, Cols);
	Matrix V (Cols, Cols);
	DiagonalMatrix W (Cols);
	DiagonalMatrix Winv (Cols);
	
	SVD (ParSetMatrix, W, U, V);
	 
	for (Atom = 0; Atom < 3; ++Atom) {
		if (W.element(Atom) == 0) {
			cerr << "\nERROR: The W matrix contains one or more 0 values, cannot perform the\n"
			     << "       the determination of the rotation matrix.\n\n";
			DC_Error = "SVD error during group fitting.";
			DC_ErrorCode = 146;
			return 146;
		}
		
		Winv.element(Atom) = 1/W.element(Atom);
	}
	
	Matrix PseudoInverse (Cols, Rows);
	PseudoInverse = (V * Winv) * U.t();
	
// cout << "Pseudo inverse = ( V.t * Winv )  *  U.t\n";
// cout << setw(15) << PseudoInverse << "\n\n";
	
	// the original matrix is (m x n) = (atoms x 3)
	// the pseudo inverse is (m x n), its transpose (n x m) = (3 x atoms)
	Matrix RNonUnitary (3, 3);
	
	RNonUnitary = GroupMatrix.t() * PseudoInverse.t();
	
// cout << "Non-unitary RotMatrix = matrix of atoms * pseudo inverse\n";
// cout << setw(15) << RNonUnitary << "\n\n";
	
	// RNonUnitary is (atoms x atoms)
	Matrix U_rot (Rows, Rows);    U_rot = 0.0;  // initialize to zero
	Matrix V_rot (Rows, Rows);    V_rot = 0.0;  // initialize to zero
	DiagonalMatrix W_rot (Rows);  W_rot = 0.0;  // initialize to zero
	
	SVD (RNonUnitary, W_rot, U_rot, V_rot);
	
	Matrix RUnitary (3, 3);
	RUnitary = U_rot * V_rot.t();
	
// cout << "R unitary:\n";
// cout << setw(15) << RUnitary << "\n\n";
	
	*RotMatrixNonUnitary = RNonUnitary.t();
	*RotMatrixUnitary = RUnitary.t();
	
	return 0;
} // of Dichro::RotationMatrix


// ================================================================================


void Dichro::CheckPlanar ( Matrix* ParSetMatrix, Matrix* GroupMatrix )
{
	int Atom;
	int AtomNumParSet = ParSetMatrix->nrows();
	int AtomNumGroup  = GroupMatrix->nrows();
	bool Planar = false;
	
	/*
	From http://mathworld.wolfram.com/SingularValueDecomposition.html:
	Note that there are several conflicting notational conventions in use in the literature.
	Press et al. (1992) define    U as (m x n), W as (n x n), and V as (n x n).
	Mathematica [& Wikip.] define U as (m x m), W as (m x n), and V as (n x n).
	*/
	Matrix U (AtomNumParSet, 3); U = 0.0;  // initialize to zero, m x n, unitary
	Matrix V (3, 3);             V = 0.0;  // initialize to zero, n x n, unitary
	DiagonalMatrix W (3);        W = 0.0;  // initialize to zero, n x n, diagonal
	
	// the product U x W x V
	Matrix UWV (AtomNumParSet, 3);   UWV = 0.0;  // initialize to zero
	
	// the pseudoinverse (U x 1/W x V)
	Matrix PseudoInverse (3, AtomNumParSet);   PseudoInverse = 0;  // initialize to zero
	
	if (DC_Debug > 4) {
		fprintf (DC_FitFile, "\n   Matrix of chromophore atom coordinates:\n");
		FilePrintMatrix (DC_FitFile, GroupMatrix);
		
		fprintf (DC_FitFile, "\n   Matrix of the parameter set coordinates:\n");
		FilePrintMatrix (DC_FitFile, ParSetMatrix);
	}
	
	// perform a singular value decomposition
	SVD (*ParSetMatrix, W, U, V);
	
	if (DC_Debug > 4) {
		fprintf (DC_FitFile, "\n   Matrix U (parameter set):\n");
		FilePrintMatrix (DC_FitFile, &U);
		
		fprintf (DC_FitFile, "\n   Check, if it is unitary:\n");
		Matrix U_uni (AtomNumParSet, AtomNumParSet);
		U_uni = U * U.t();
		FilePrintMatrix (DC_FitFile, &U_uni);
		
		fprintf (DC_FitFile, "\n   Matrix V (parameter set):\n");
		FilePrintMatrix (DC_FitFile, &V);
		
		fprintf (DC_FitFile, "\n     Check, if it is unitary:\n");
		Matrix V_uni (3, 3);
		V_uni = V * V.t();
		FilePrintMatrix (DC_FitFile, &V_uni);
		
		fprintf (DC_FitFile, "\n   Matrix W (parameter set):\n");
		FilePrintMatrix (DC_FitFile, &W);
	}
	
	// IMPORTANT: a planar system would run into a division by zero
	// W has always 3 columns as the initial matrix has 3 (x, y, z)
	for (Atom = 0; Atom < 3; Atom++)
		if (W.element(Atom) == 0) Planar = true;
	
	if (Planar) {
		if (DC_Debug > 4) {
			fprintf (DC_FitFile, "\n   Planar system found, adding virtual atom\n");
			fprintf (DC_FitFile,   "   ----------------------------------------\n\n");
		}
		
		// ================================================================================
		// Add an additional point to the parameter set atoms
		// ================================================================================
		
		// vectors needed for planar systems
		Matrix Diff12 (1, 3);
		Matrix Diff32 (1, 3);
		Matrix Cross  (1, 3);
		
		// calculate the difference between the 1st and 2nd point
		Diff12 = ParSetMatrix->row(1) - ParSetMatrix->row(2);
		// calulate the difference between the 2nd and 3rd point
		Diff32 = ParSetMatrix->row(3) - ParSetMatrix->row(2);
		// calculate the cross product between the difference vectors
		Cross  = crossproduct (Diff12, Diff32);
		Cross  = Cross + ParSetMatrix->row(2);
		
		if (DC_Debug > 4) {
			fprintf (DC_FitFile, "   ParSet diff 1-2:     ");
			FilePrintMatrix (DC_FitFile, &Diff12);
			fprintf (DC_FitFile, "   ParSet diff 2-3:     ");
			FilePrintMatrix (DC_FitFile, &Diff32);
			fprintf (DC_FitFile, "   ParSet crossproduct: ");
			FilePrintMatrix (DC_FitFile, &Cross);
		}
		
		// add another row to the matrix
		++AtomNumParSet;  // the virtual atom
		// (V is still 3x3, 3 being the columns of the original matrix, i.e. x,y,z)
		ParSetMatrix->resize_keep (AtomNumParSet, 3);  // an additional row
		Matrix U (AtomNumParSet, 3); U = 0.0;
		PseudoInverse.resize_keep (3, AtomNumParSet); // an additional column
		
		// add the cross product as 'virtual' additional point to the array
		ParSetMatrix->row(AtomNumParSet) = Cross;
		
		
		// ================================================================================
		// Add an additional point to the chromophore atoms
		// ================================================================================
		
		// calculate the difference between the 1st and 2nd point
		Diff12 = GroupMatrix->row(1) - GroupMatrix->row(2);
		// calulate the difference between the 2nd and 3rd point
		Diff32 = GroupMatrix->row(3) - GroupMatrix->row(2);
		// calculate the cross product between the difference vectors
		Cross  = crossproduct (Diff12, Diff32);
		Cross  = Cross + GroupMatrix->row(2);
		
		if (DC_Debug > 4) {
			fprintf (DC_FitFile, "\n\n");
			fprintf (DC_FitFile, "   Group diff 1-2:      ");
			FilePrintMatrix (DC_FitFile, &Diff12);
			fprintf (DC_FitFile, "   Group diff 2-3:      ");
			FilePrintMatrix (DC_FitFile, &Diff32);
			fprintf (DC_FitFile, "   Group crossproduct:  ");
			FilePrintMatrix (DC_FitFile, &Cross);
		}
		
		// add another row to the matrix
		GroupMatrix->resize_keep (AtomNumGroup+1, 3);
		
		// add the cross product as 'virtual' additional point to the array
		GroupMatrix->row(AtomNumParSet) = Cross;
		
		++AtomNumGroup;  // the virtual atom
		
		if (DC_Debug > 4) {
			fprintf (DC_FitFile, "\n   Extended matrix of chromophore atom coordinates:\n");
			FilePrintMatrix (DC_FitFile, GroupMatrix);
		}
	} // of if (Planar)
} // of Dichro::CheckPlanar


// ================================================================================


void Dichro::Rotate ( vector<double> *In, vector<double> *Out, Matrix *RotMatrix )
// rotate vector using a given rotation matrix
{
	int k, l;
	Out->clear(); // clear all elements
	
	// THIS WAS VERIFIED TO PRODUCE IDENTICAL OUTPUT AS THE FORTRAN VERSION
	for (k = 0; k < 3; k++) {
		Out->push_back (0.0);  // set all three elements to 0
		
		for (l = 0; l < 3; l++) {
			Out->at(k) = Out->at(k) + (RotMatrix->element(l, k) * In->at(l) );
		}
	}
	
	return;
} // of Dichro::Rotate


// ================================================================================

