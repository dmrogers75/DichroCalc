// #################################################################################################
//
//  Program:      matrix.cpp
//
//  Function:     Part of DichroCalc:
//                Routines to create the Hamiltonian matrix
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


int  Dichro::HamiltonianMatrix ( void )
{
	// The diagonal are the excitation energies, each transition is represented by
	// one element on the diagonal. That is, the total number of transitions in the
	// system is equal to the matrix dimension (as it is a square matrix).
	int MatrixDimension     = DC_System.MatrixDimension;
	int NumberOfGroups      = DC_System.NumberOfGroups;
	int NumberOfTransitions = DC_System.NumberOfTransitions;
	
	int Group, Trans, iGroup, jGroup, iTrans, jTrans, row, col;
	vector<int> GroupSeq, TransSeq;
	vector<string> ParSetSeq;
	SystemGroup* CurGroup;
	double Interaction;
	
	if (DC_Verbose) Dichro::NewTask ( "Setting up Hamiltonian Matrix" );
	
	if (DC_Debug > 2)
		Dichro::NewFileTask (DC_DbgFile, "Setting up Hamiltonian Matrix");
	
	SymmetricMatrix Hamiltonian (MatrixDimension);
	
	Hamiltonian = 0.0;
	// PrintMatrix (&Hamiltonian);
	
	// The matrix is filled by looping over each group and each transition on it and
	// calculate the interaction with each transition on other groups. This is slightly
	// complicated by the fact that the number of transitions is not constant in all
	// groups. The following nested loops create two vectors with the sequence of groups
	// and their transitions along the diagonal
	for (Group = 0; Group < NumberOfGroups; Group++) {
		CurGroup = &DC_System.Groups.at(Group);
		
		// initialize a vector element for the current group (needed in CD_Calculation)
		ResultsGroup NewGroup;
		DC_Results.Groups.push_back (NewGroup);
		
		for (Trans = 0; Trans < CurGroup->NumberOfTransitions; Trans++) {
			GroupSeq.push_back (Group);
			TransSeq.push_back (Trans);
			ParSetSeq.push_back (DC_System.Groups.at(Group).ParameterSet);
			
			DC_Results.Trans.EDM.push_back (CurGroup->Trans.at(Trans).EDM);
			DC_Results.Trans.MDM.push_back (CurGroup->Trans.at(Trans).MDM);
			DC_Results.Trans.MDMconv.push_back (CurGroup->Trans.at(Trans).MDM);
			DC_Results.Trans.Energy.push_back (CurGroup->Trans.at(Trans).Energy);
			DC_Results.Trans.Reference.push_back (CurGroup->Reference);
			
			DC_Results.Trans.GroupSequence.push_back (Group);
			DC_Results.Trans.TransSequence.push_back (Trans);
			DC_Results.Trans.ParSetSequence.push_back (DC_System.Groups.at(Group).ParameterSet);
		}
	}
	
	// DEBUG OUTPUT
	// printf ("\n\nDC_Results - Initial EDM\n\n");
	// for (Trans = 0; Trans < MatrixDimension; Trans++)
	// 	printf ("   %2d     %12.6f   %12.6f   %12.6f\n", Trans, DC_Results.EDM.at(Trans).at(0),
	//         DC_Results.EDM.at(Trans).at(1), DC_Results.EDM.at(Trans).at(2));
	//
	//printf ("\n\nDC_Results - Initial MDM\n\n");
	//for (Trans = 0; Trans < MatrixDimension; Trans++)
	//	printf ("   %2d     %12.6f   %12.6f   %12.6f\n", Trans, DC_Results.MDM.at(Trans).at(0),
	//	        DC_Results.MDM.at(Trans).at(1),DC_Results.MDM.at(Trans).at(2));
	
	if (DC_Debug > 2) {
		int i = 0;
		
		fprintf (DC_DbgFile, "\n   Sequence of groups and transitions along the diagonal:\n\n");
		fprintf (DC_DbgFile, "   row , col   =   group - trans\n");
		
		for (i = 0; i < NumberOfTransitions; i++) {
			fprintf (DC_DbgFile, "   %3d , %-3d   =     %3d - %-3d\n",
			                       i, i, GroupSeq.at(i), TransSeq.at(i));
		}
		fprintf (DC_DbgFile, "\n\n");
	}
	
	if (DC_Debug > 2) {
		fprintf (DC_DbgFile,
		         "\n   Interaction terms for the lower triangle of the Hamiltonian:\n\n");
		fprintf (DC_DbgFile, "   row , col   =   group - trans / group - trans\n");
		
		for (row = 0; row < NumberOfTransitions; row++) {
			for (col = 0; col <= row; col++) {
				iGroup = GroupSeq.at(row);   iTrans = TransSeq.at(row);
				jGroup = GroupSeq.at(col);   jTrans = TransSeq.at(col);
				
				fprintf (DC_DbgFile, "   %3d , %-3d   =     %3d - %-3d   /   %3d - %-3d\n",
				                        row,  col,      iGroup,  iTrans, jGroup, jTrans);
			}
		}
		fprintf (DC_DbgFile, "\n\n");
	}
	
	if (DC_Debug > 2)
		fprintf (DC_DbgFile,
			"                  row , col   =   group - trans / group - trans        Interaction\n");
	
	// the two for loops run over the lower triangle of the matrix (including the diagonal)
	for (row = 0; row < NumberOfTransitions; row++) {
		for (col = 0; col <= row; col++) {
			iGroup = GroupSeq.at(row);
			iTrans = TransSeq.at(row);
			
			jGroup = GroupSeq.at(col);
			jTrans = TransSeq.at(col);
			
			Interaction = 0.0;
			
			// if it is a diagonal element (this is equal to row == col)
			if (iGroup == jGroup && iTrans == jTrans) {
				if (DC_System.Groups.at(iGroup).ChargeTransfer) iTrans = iTrans + 4;
				if (DC_System.Groups.at(jGroup).ChargeTransfer) jTrans = jTrans + 4;
				
				Hamiltonian.element (row, col) =
					DC_System.Groups.at(iGroup).Trans.at(iTrans).Energy;
			}
			
			// if it is an interaction between transitions on the same group
			else if ( GroupsOverlap (iGroup, jGroup) ) {
				Interaction = Dichro::SameGroupInteraction (iGroup, iTrans, jGroup, jTrans );
				Hamiltonian.element (row, col) = Interaction;
//				printf ("+ Interaction = %12.8f\n", Interaction);
				
				if (DC_Debug > 2)
					fprintf (DC_DbgFile,
					"   Overlap:       %3d , %-3d   =     %3d - %-3d   /   %3d - %-3d   =   %14.6f\n",
					                   row,  col,      iGroup,  iTrans , jGroup, jTrans , Interaction);
			}
			
			// if it is an interaction between transitions on different groups
			else if ( not GroupsOverlap (iGroup, jGroup) ) {
				if (DC_System.Groups.at(iGroup).ChargeTransfer) iTrans = iTrans + 4;
				if (DC_System.Groups.at(jGroup).ChargeTransfer) jTrans = jTrans + 4;
				
				Interaction = Dichro::DifferentGroupInteraction
				                      (iGroup, iTrans, jGroup, jTrans, false );
				Hamiltonian.element (row, col) = Interaction;
//				printf ("+ Interaction = %12.8f\n", Interaction);
				
				if (DC_Debug > 2)
					fprintf (DC_DbgFile,
					"   Non-overlap:   %3d , %-3d   =     %3d - %-3d   /   %3d - %-3d   =   %12.4f\n",
					                    row,  col,      iGroup,  iTrans , jGroup, jTrans , Interaction);
			}
			
			else {
				cerr << "\nERROR\n";
			}
		}
	}
	
	// printf ("\nHamiltonian matrix in J:\n\n");
	// PrintMatrix (&Hamiltonian);
	
	for (row = 0; row < NumberOfTransitions; row++) {
		for (col = 0; col < row; col++) {
			// convert the off-diagonal elements from Joule to cm-1
			Hamiltonian.element (row, col) *= 5036.0;
		}
	}

	//---Hainam's hack. Reading Hamiltonian from file
	fstream input;
	input.open("new.mat", ios::in);
	cout << "Open file new.mat for reading Hamiltonian " << endl;
	if (input.is_open()) {
          for (row = 0; row < NumberOfTransitions; row++) {
            for (col = 0; col < NumberOfTransitions; col++) {
	      input >> Hamiltonian.element (row, col);
	      //cout << row+1 << " " << col+1 << " " << Hamiltonian.element (row, col) << endl;
	    }
	  }
	}
	
	if (DC_Verbose) printf ("   Diagonalizing\n");
	
	Matrix Eigenvectors (MatrixDimension, MatrixDimension);
	SymmetricMatrix WorkSpace (MatrixDimension);
	DiagonalMatrix Eigenvalues (MatrixDimension);
	
	// diagonalize the Hamiltonian using the Jacobi mechanism (extremly reliable but slower)
	Jacobi (Hamiltonian, Eigenvalues, WorkSpace, Eigenvectors);
	// diagonalize using the householder mechanism (faster)
	// eigenvalues (Hamiltonian, Eigenvalues, Eigenvectors);
	
	row = 0;
	col = 0;
	
	int CurRow, CurCol;
	int StartCol = 0;
	
	// copy the submatrix of this group out of the complete Hamiltonian
	for (Group = 0; Group < NumberOfGroups; Group++) {
		int TransNumber = DC_System.Groups.at(Group).NumberOfTransitions;
		SymmetricMatrix Submatrix (TransNumber);
		
		// cut out the matrix part belonging to the current group
		for (CurRow = 0; CurRow < TransNumber; CurRow++) {
			col = StartCol; // reset the column to the start column of the submatrix
			
			for (CurCol = 0; CurCol <= CurRow; CurCol++) {
				Submatrix.element(CurRow, CurCol) = Hamiltonian.element(row, col);
				col++;
			}
			row++;
		}
		StartCol += TransNumber;
		DC_Results.Groups.at(Group).Submatrix = Submatrix;
	}
	
	// the actual results of the diagonalization
	DC_Results.Hamiltonian         = Hamiltonian;
	DC_Results.Eigenvalues         = Eigenvalues;
	DC_Results.Eigenvectors        = Eigenvectors;
	
	// copy some information on the system to DC_Results
	DC_Results.MatrixDimension     = DC_System.MatrixDimension;
	DC_Results.NumberOfAtoms       = DC_System.NumberOfAtoms;
	DC_Results.NumberOfGroups      = DC_System.NumberOfGroups;
	DC_Results.NumberOfTransitions = DC_System.NumberOfTransitions;

	if (DC_PrintMat) {
		if (DC_Verbose) printf ("      Output written to %s\n", DC_MatFilename.c_str());
		
		// the Hamiltonian matrix
		Dichro::OutputFileHeadline (DC_MatFile, "$DC_Results.Hamiltonian: Hamiltonian Matrix", false);
		FilePrintMatrix (DC_MatFile, &DC_Results.Hamiltonian);
		
		// the eigenvectors
		Dichro::OutputFileHeadline (DC_MatFile, "$DC_Results.Eigenvectors: Eigenvectors");
		FilePrintMatrix (DC_MatFile, &DC_Results.Eigenvectors);
		
		// the eigenvalues
		Dichro::OutputFileHeadline (DC_MatFile, "$DC_Results.Eigenvalues: Eigenvalues");
		FilePrintMatrix (DC_MatFile, &DC_Results.Eigenvalues, false );
	}
	
	return 0;
} // of Dichro::HamiltonianMatrix


// ================================================================================


bool Dichro::GroupsOverlap ( int iGroup, int jGroup )
// checks if two groups overlap (i.e. any of their atoms if it's about CT groups for example)
{
	if (iGroup == jGroup) return true;
	
	unsigned int iAtom, jAtom;
	
	// if it is a charge-transfer chromophore, the group numbers will be different and
	// the single atoms have to be checked in order to be sure
	for (iAtom = 0; iAtom < DC_System.Groups.at(iGroup).AtomIndices.size(); iAtom++)
		for (jAtom = 0; jAtom < DC_System.Groups.at(jGroup).AtomIndices.size(); jAtom++)
			// if the atom index in group i is equal to the atom index in group j, return true
			if (DC_System.Groups.at(iGroup).AtomIndices.at(iAtom) ==
				 DC_System.Groups.at(jGroup).AtomIndices.at(jAtom))
					return true;
	
	return false;
} // of Dichro::GroupsOverlap


// ================================================================================


double Dichro::DifferentGroupInteraction (int iGroup, int iTrans, int jGroup, int jTrans, bool Perm)
// calculates the interaction of transitions on different groups
{
	vector<ParSetMonopole>* iMonopoles = &DC_System.Groups.at(iGroup).Trans.at(iTrans).Monopoles;
	vector<ParSetMonopole>* jMonopoles;
	
	if (not Perm)
		jMonopoles = &DC_System.Groups.at(jGroup).Trans.at(jTrans).Monopoles;
	else
		jMonopoles = &DC_System.Groups.at(jGroup).Perm.at(jTrans).Monopoles;
	
	double DistanceThreshold = 0.01;
	
	unsigned int iMono, jMono;
	double Distance, TempInt, Interaction;
	
	// DEBUG OUTPUT
	// printf ("iGroup %2d  iTrans %2d   (Mono %14.8f) - jGroup %2d  jTrans %2d   (Mono %14.8f)\n",
	//          iGroup+1,     iTrans+1, iMonopoles->at(0).Charge, jGroup+1,     jTrans+1, jMonopoles->at(0).Charge);
	
	Interaction = 0.0;
	
	for (iMono = 0; iMono < iMonopoles->size(); iMono++) {
		TempInt = 0.0;
		
		for (jMono = 0; jMono < jMonopoles->size(); jMono++) {
			Distance = PointDistance ( &iMonopoles->at(iMono).Coord, &jMonopoles->at(jMono).Coord );
			// printf ("jMono %2d  Charge %8.3f  Distance %8.3f\n",
			//         jMono, jMonopoles->at(jMono).Charge, Distance);
			
			if (Distance < DistanceThreshold) {
				printf ("WARNING: Monopole distance below %f Angstrom for the calculation of\n",
				        DistanceThreshold);
				printf ("         the interaction on different groups. Skipped.\n");
				printf ("         Group 1 = %4d, Trans 1 = %4d  -  Group 2 = %4d, Trans 2 = %4d\n",
				        iGroup, iTrans, jGroup, jTrans);
				printf ("         Distance %8.3f Angstrom\n", Distance);
			}
			else
				TempInt = TempInt + ( jMonopoles->at(jMono).Charge / Distance );
		}
		
		Interaction = Interaction + ( iMonopoles->at(iMono).Charge * TempInt );
		// printf ("Interaction = %12.8f    TempInt = %12.8f\n", Interaction, TempInt);
		// printf ("iGroup %2d  iTrans %2d   jGroup %2d  jTrans %2d  Int %14.3f\n",
		//          iGroup,     iTrans,      jGroup,     jTrans,     Interaction);
		// printf ("iMono.Charge  %8.3f  *  TempInt %14.3f =  Interaction %14.3f\n\n",
		//         iMonopoles->at(iMono).Charge, TempInt, iMonopoles->at(iMono).Charge * TempInt);
	}
	
	return Interaction;
} // of Dichro::DifferentGroupInteraction


// ================================================================================


double Dichro::SameGroupInteraction ( int iGroup, int iTrans, int jGroup, int jTrans )
// calculates the interaction of transitions on the same group
{
	SystemGroup* iCurGroup = &DC_System.Groups.at(iGroup);
	SystemGroup* jCurGroup = &DC_System.Groups.at(jGroup);
	double Interaction;
	int Group, Trans;
	
	// catch overlapping CT-groups (different group numbers but sharing a peptide bond)
	if (iGroup != jGroup and
	           (iCurGroup->ChargeTransfer and jCurGroup->ChargeTransfer) ) return 0.0;
	
	Interaction = 0.0;

	if (iCurGroup->ChargeTransfer) {
		Group = iGroup;
		Trans = iCurGroup->NumberOfTransitions + 4;
	}
	else if (jCurGroup->ChargeTransfer) {
		Group = jGroup;
		Trans = iCurGroup->NumberOfTransitions + 4; // no mistake, it's iGroup
	}
	else { // non-CT case, iGroup == jGroup
		Group = iGroup;
		Trans = iCurGroup->NumberOfTransitions;
	}

	if (iCurGroup->ChargeTransfer) iTrans = iTrans + 4;
	if (jCurGroup->ChargeTransfer) jTrans = jTrans + 4;
	
	
	// the "+1" and --Transitions was necessary to deal with the fact that counting starts
	// at 0 (this wasn't needed in the Fortran version)
	int MinTrans    = min (iTrans+1, jTrans+1);
	int MaxTrans    = max (iTrans+1, jTrans+1);
	int Transitions = MinTrans * Trans - (MinTrans * (MinTrans+1)) / 2 + MaxTrans;
	--Transitions;
	
	// printf ("Min %d    Max %d   Transitions %d\n", MinTrans, MaxTrans, Transitions);
	// printf ("Same group:  "); //iGroup - iTrans / jGroup - jTrans");
	// printf ("  %3d - %-3d  /  %3d - %-3d  =>  %3d Transitions\n",
	//        iGroup+1, iTrans+1, jGroup+1, jTrans+1, Transitions);
	
	// Now run over all other groups (not the same one of the actual transition, despite
	// the function name) and interact the GS density with higher states (e.g. 2->3)
	for (jGroup = 0; jGroup < DC_System.NumberOfGroups; jGroup++) {
		
		if (DC_System.Groups.at(jGroup).ChargeTransfer) {
			// CT groups are ignored here since the monomers will be considered later to
			// represent the ground state
			// if (DC_Debug > 3)
				// fprintf (DC_DbgFile, "   Interaction (same group, CT):     %4d - %4d = %8.3f\n",
				//                        iGroup, jGroup, 0.0);
		}
		else {
			if (not GroupsOverlap (iGroup, jGroup)) {
				Interaction = Interaction +
				              Dichro::DifferentGroupInteraction (Group, Transitions, jGroup, 0, true);

				// DEBUG OUTPUT
				// printf ("& iGroup %2d  iTrans %2d  (Mono %8.3f)  -  jGroup %2d  jTrans %2d (Mono %8.3f)  =>  %12.8f (MinTrans = %d, MaxTrans = %d)\n",
				//         Group+1, Transitions+1, DC_System.Groups.at(Group).Trans.at(Transitions).Monopoles.at(0).Charge,
				// 		  jGroup+1, 0+1,           DC_System.Groups.at(jGroup).Perm.at(0).Monopoles.at(0).Charge,
				// 		  Interaction, MinTrans, MaxTrans);
			}
		}
	}
	
	//printf ("iGroup %2d  iTrans %2d   jGroup %2d  Int %14.3f\n",
	//         iGroup,     iTrans,      jGroup,     Interaction);
	
	// printf ("* Interaction = %12.8f\n", Interaction);
	return Interaction;
} // of Dichro::SameGroupInteraction


// ================================================================================


