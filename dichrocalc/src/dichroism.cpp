// #################################################################################################
//
//  Program:      dichroism.cpp
//
//  Function:     Part of DichroCalc:
//                Routines to calculate circular and linear dichroism
//
//  Author:       Benjamin M. Bulheller
//
//  Affiliation:  Research group Jonathan D. Hirst
//                University of Nottingham, UK
//
//                Research group Shaul Mukamel
//                University of California at Irvine, USA
//
//  Version:      $Revision: 4584 $, $Date: 2009-06-15 22:13:23 -0700 (Mon, 15 Jun 2009) $
//
//  Date:         May 2009
//
// #################################################################################################


#include "../include/dichrocalc.h"


// ================================================================================


int  Dichro::CD_Calculation ( void )
{
	const double MagicNumber = 3.3879E-6;
	
	int Trans, iGroup, iTrans, jGroup, jTrans, Coord, iCount, jCount;
	double RotationalStrength, DipoleStrength, Wavelength;
	double RotationalSum = 0;
	
	SystemGroup* iCurGroup;
	SystemGroup* jCurGroup;
	SystemTransition* jCurTrans;
	
	DiagonalMatrix* Eigenvalues  = &DC_Results.Eigenvalues;
	Matrix*         Eigenvectors = &DC_Results.Eigenvectors;
	
	int NumberOfTransitions = DC_System.NumberOfTransitions;
	
	// define and initialize the required vectors
	vector<double> EDM (3, 0.0);
	vector<double> MDM (3, 0.0);
	
	if (DC_Verbose) {
		Dichro::NewTask ( "Dichroism Calculation" );
		printf ("   Calculating CD\n");
	}
	
	for (Trans = 0; Trans < NumberOfTransitions; Trans++) {
		DC_Results.Trans.MDMconv.at(Trans).at(0) += (
		   MagicNumber * DC_Results.Trans.Energy.at(Trans) *
		      ( DC_Results.Trans.Reference.at(Trans).at(1) * DC_Results.Trans.EDM.at(Trans).at(2) -
		          DC_Results.Trans.Reference.at(Trans).at(2) * DC_Results.Trans.EDM.at(Trans).at(1) )
		);
		
		DC_Results.Trans.MDMconv.at(Trans).at(1) += (
		   MagicNumber * DC_Results.Trans.Energy.at(Trans) *
		      ( DC_Results.Trans.Reference.at(Trans).at(2) * DC_Results.Trans.EDM.at(Trans).at(0) -
		          DC_Results.Trans.Reference.at(Trans).at(0) * DC_Results.Trans.EDM.at(Trans).at(2) )
		);
		
		DC_Results.Trans.MDMconv.at(Trans).at(2) += (
		   MagicNumber * DC_Results.Trans.Energy.at(Trans) *
		      ( DC_Results.Trans.Reference.at(Trans).at(0) * DC_Results.Trans.EDM.at(Trans).at(1) -
		          DC_Results.Trans.Reference.at(Trans).at(1) * DC_Results.Trans.EDM.at(Trans).at(0) )
		);
	}
	
	
	if (DC_Debug > 2) {
		Dichro::NewFileTask (DC_DbgFile, "Circular Dichroism Calculation");
		
		fprintf (DC_DbgFile, "\n\n   $DC_Results - Magnetic Transition Dipole Moment Conversion\n\n");
		fprintf (DC_DbgFile, "                             initial moments");
		fprintf (DC_DbgFile, "                            converted moments\n");
		fprintf (DC_DbgFile, "    trans            x              y              z");
		fprintf (DC_DbgFile, "              x              y              z\n");
		
		for (Trans = 0; Trans < NumberOfTransitions; Trans++)
			fprintf (DC_DbgFile,
			       "      %2d     %12.6f   %12.6f   %12.6f   %12.6f   %12.6f   %12.6f\n", Trans,
			       DC_Results.Trans.MDM.at(Trans).at(0),
			       DC_Results.Trans.MDM.at(Trans).at(1),
			       DC_Results.Trans.MDM.at(Trans).at(2),
			       DC_Results.Trans.MDMconv.at(Trans).at(0),
			       DC_Results.Trans.MDMconv.at(Trans).at(1),
			       DC_Results.Trans.MDMconv.at(Trans).at(2));
		
		fprintf (DC_DbgFile, "\n");
	}
	
	iCount = 0;
	
	if (DC_Verbose and DC_PrintCdl)
		printf ("      Output written to %s\n", DC_CdlFilename.c_str());
	
	for (iGroup = 0; iGroup < DC_System.NumberOfGroups; iGroup++) {
		iCurGroup = &DC_System.Groups.at(iGroup);
		
		for (iTrans = 0; iTrans < iCurGroup->NumberOfTransitions; iTrans++) {
			
			RotationalStrength = 0.0;
			DipoleStrength     = 0.0;
			
			for (Coord = 0; Coord < 3; Coord++) {
				EDM.at(Coord) = 0.0;
				MDM.at(Coord) = 0.0;
			}
			
			jCount = 0;
			
			for (jGroup = 0; jGroup < DC_System.NumberOfGroups; jGroup++) {
				jCurGroup = &DC_System.Groups.at(jGroup);
				
				for (jTrans = 0; jTrans < jCurGroup->NumberOfTransitions; jTrans++) {
					jCurTrans = &jCurGroup->Trans.at(jTrans);
					// DEBUG OUTPUT
					// printf ("iGroup %2d  iTrans %2d  -  jGroup %2d  jTrans %2d\n",
					//                 iGroup,     iTrans,        jGroup,     jTrans);
					
					// DEBUG OUTPUT
					// printf ("Input-EDM: %8.3f %8.3f %8.3f\n",
					//        jCurTrans->EDM.at(0),  jCurTrans->EDM.at(1),  jCurTrans->EDM.at(2));
					// printf ("Input-MDM: %8.3f %8.3f %8.3f\n",
					//        DC_Results.Trans.MDMconv.at(jCount).at(0),
					//        DC_Results.Trans.MDMconv.at(jCount).at(1),
					//        DC_Results.Trans.MDMconv.at(jCount).at(2));
					
					for (Coord = 0; Coord < 3; Coord++) {
						MDM.at(Coord) +=
						   ( Eigenvectors->element(jCount,iCount) *
						            DC_Results.Trans.MDMconv.at(jCount).at(Coord) );
						
						EDM.at(Coord) +=
						   ( Eigenvectors->element(jCount,iCount) * jCurTrans->EDM.at(Coord) *
						        jCurTrans->Energy / Eigenvalues->element(iCount) );
					}
					
					++jCount;
				} // of for (jTrans = 0; jTrans < jCurGroup->NumberOfTransitions; jTrans++)
			} // of for (jGroup = 0; jGroup < DC_System.NumberOfGroups; jGroup++)
			
			// DEBUG OUTPUT
			// printf ("Final EDM:  %8.3f %8.3f %8.3f\n", EDM.at(0), EDM.at(1), EDM.at(2));
			// printf ("Final MDM:  %8.3f %8.3f %8.3f\n", MDM.at(0), MDM.at(1), MDM.at(2));
			
			for (Coord = 0; Coord < 3; Coord++) {
				// Rosenfeld equation:
				// calculate the dot product of the elec. and mag. dipole moment:
				RotationalStrength += (EDM.at(Coord) * MDM.at(Coord));
				
				// calculate the 'square' dot product of the elec. dipole moment
				DipoleStrength += ( EDM.at(Coord) * EDM.at(Coord) );
			}
			
			RotationalSum += RotationalStrength;
			
			// convert the wavelength to nanometer
			Wavelength = 1E7 / Eigenvalues->element(iCount);
			
			if (DC_PrintCdl)
				fprintf (DC_CdlFile, "%14.8f %14.8f\n", Wavelength, RotationalStrength);
			
			// Add all results of the calculation for this transition to DC_Results. The first
			// three assignments actually update the initial input data (i.e. it replaces the
			// monomer data with the result for the interacting system).
			DC_Results.Trans.EDM.at(iCount) = EDM;
			DC_Results.Trans.MDM.at(iCount) = MDM;
			DC_Results.Trans.Energy.at(iCount) = Eigenvalues->element(iCount);
			DC_Results.Trans.RotationalStrength.push_back (RotationalStrength);
			DC_Results.Trans.DipoleStrength.push_back (DipoleStrength);
			DC_Results.Trans.Wavelength.push_back (Wavelength);
			DC_Results.Trans.Reference.push_back (iCurGroup->Reference);
			 
			// add the results for access via the group count (the vectors have been initialized
			// in Dichro::HamiltonianMatrix during setting up the group/transition sequence)
			DC_Results.Groups.at(iGroup).EDM.push_back (EDM);
			DC_Results.Groups.at(iGroup).MDM.push_back (MDM);
			DC_Results.Groups.at(iGroup).Energy.push_back (Eigenvalues->element(iCount));
			DC_Results.Groups.at(iGroup).Wavelength.push_back (Wavelength);
			DC_Results.Groups.at(iGroup).RotationalStrength.push_back (RotationalStrength);
			DC_Results.Groups.at(iGroup).DipoleStrength.push_back (DipoleStrength);
			DC_Results.Groups.at(iGroup).Reference = iCurGroup->Reference;
			DC_Results.Groups.at(iGroup).ParameterSet = iCurGroup->ParameterSet;
			DC_Results.Groups.at(iGroup).NumberOfTransitions = iCurGroup->NumberOfTransitions;
			DC_Results.Groups.at(iGroup).ChargeTransfer = iCurGroup->ChargeTransfer;
			
			++iCount;
		} // of for (iTrans = 0; iTrans < CurGroup->NumberOfTransitions; Trans++)
	}	// of for (iGroup = 0; iGroup < DC_System.NumberOfGroups; iGroup++)
	
	return 0;
} // of Dichro::CD_Calculation


// ================================================================================


int  Dichro::LD_Calculation ( void )
{
	int iGroup, iTrans, jGroup, jTrans;
	double Wavelength, Energy, Coefficient, TotalPolarization, TotalOscillatorStrength;
	int i, j, k, iCount, jCount, Coord;
	SystemGroup* iCurGroup;
	SystemGroup* jCurGroup;
	
	vector< vector< vector<double> > > PolTensor; // the polarization tensor
	
	if (DC_Verbose) printf ("   Calculating LD\n");
	if (DC_Debug > 2) Dichro::NewFileTask (DC_DbgFile, "Linear Dichroism Calculation");
	
	// determine the maximum number of transitions on a group as this is needed for he
	// dimension of the polarization tensor
	int MaxNumberOfTransitions = 0;
	
	for (iGroup = 0; iGroup < DC_System.NumberOfGroups; iGroup++)
		if (MaxNumberOfTransitions < DC_System.Groups.at(iGroup).NumberOfTransitions)
			MaxNumberOfTransitions = DC_System.Groups.at(iGroup).NumberOfTransitions;
	
	//PolTensor = &DC_Results.PolTensor;
	
	// Used in publications PCCP 2002, 4, 4051-4057 and PCCP 2007, 9, 2020-2035
	
	// initialize the polarization tensor
	// CAUTION: k was changed to the maximum number of transitions on one group, before it was
	// limited to only 2 transitions in general
	for (k = 0; k < MaxNumberOfTransitions+1; k++) {
		// create a vector of 3 vectors, all initialized with MaxNumberOfTransitions+2 of zeroes
		vector< vector<double> > Dim2 (3, vector<double>(MaxNumberOfTransitions+2, 0.0));
		// push it into the polarization tensor to create the third dimension
		PolTensor.push_back (Dim2);
	}
	
	vector< vector<double> > OscillatorStrength
	           (DC_System.NumberOfGroups * MaxNumberOfTransitions, // this is a 'safe' value
	                   vector<double>(MaxNumberOfTransitions, 0.0));
	
	if (DC_PrintPol) {
		if (DC_Verbose) printf ("      Output written to %s\n", DC_PolFilename.c_str());
		
		fprintf (DC_PolFile, "\nTotal Transition Dipole Moments\n");
		fprintf (DC_PolFile,   "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n");
		fprintf (DC_PolFile, "   Trans Wavelength        x            y            z      ");
		
		for (iTrans = 0; iTrans < MaxNumberOfTransitions; iTrans++)
			fprintf (DC_PolFile, "  trans %d    ", iTrans);
		
		fprintf (DC_PolFile, "    Total     Norm\n");
	}
	
	iCount = 0;
	
	for (iGroup = 0; iGroup < DC_System.NumberOfGroups; iGroup++) {
		iCurGroup = &DC_System.Groups.at(iGroup);
		
		for (iTrans = 0; iTrans < iCurGroup->NumberOfTransitions; iTrans++) {
			Wavelength = DC_Results.Trans.Wavelength.at(iCount);
			
			// Planck      h = 6.62608 E-34 Js = m^2 kg s^-1
			// lightspeed  c = 2.997924 E+8 m s^-1
			// Joule       J = 2.29271276 E+17 Hartree
			
			// Energy in a.u.
			// Energy  = 45.563421 / Wavelength
			// Energy  = ->Energy.at(Trans);
			Energy  = 2.29371276 * 6.62608 * 2.997924 /  Wavelength;
			// printf ("Energy: %12.6f\n", Energy);
			
			// First get the total transition dipole moment for the Trans-th transition.
			// For the moment these are kept separated into the individual contributions
			// 
			// 
			
			// initialize the polarization vector for each transition of this group
			vector< vector<double> > Dxyz (iCurGroup->MaxNumberOfTransitions*3, vector<double>(3, 0.0)); // DMR 20180601 NumberOfTransitions to MaxNumberOfTransitions
			
			jCount = 0;
			
			for (jGroup = 0; jGroup < DC_System.NumberOfGroups; jGroup++) {
				jCurGroup = &DC_System.Groups.at(jGroup);
				
				for (jTrans = 0; jTrans < jCurGroup->NumberOfTransitions; jTrans++) {
					// convert from Debye to a.u.
					// 1 a.u. = 8.4784E-30 Cm = 2.5417 D
					Coefficient = DC_Results.Eigenvectors.element(jCount, iCount) / 2.5417477;
					
					// DEBUG OUTPUT
					// printf ("Coefficient:  %12.6f\n", Coefficient);
					// printf ("EDM:          %12.6f %12.6f %12.6f\n",
					//         jCurGroup->Trans.at(jTrans).EDM.at(0),
					//         jCurGroup->Trans.at(jTrans).EDM.at(1),
					//         jCurGroup->Trans.at(jTrans).EDM.at(2));
					
					for (Coord = 0; Coord < 3; Coord++)
						Dxyz.at(jTrans).at(Coord) +=
						                (Coefficient * jCurGroup->Trans.at(jTrans).EDM.at(Coord));
					
					jCount++;
				} // of for (jTrans = 0; jTrans < jCurGroup->NumberOfTransitions; jTrans++)
			} // of for (jGroup = 0; jGroup < DC_System.Groups.size(); jGroup++)
			
			// initialize a fresh vector for every new transition
			vector<double> PolVec (3, 0.0);
			
			// PolVec contains the polarization vector of the iCount-th transition
			for (jTrans = 0; jTrans < iCurGroup->NumberOfTransitions; jTrans++) {
				if (DC_Debug > 2) {
					fprintf (DC_DbgFile, "   jTrans: %2d  %12.6f%12.6f%12.6f\n", jTrans,
					       Dxyz.at(jTrans).at(0), Dxyz.at(jTrans).at(1), Dxyz.at(jTrans).at(2));
				}
				
				for (Coord = 0; Coord < 3; Coord++)
					PolVec.at(Coord) += Dxyz.at(jTrans).at(Coord);
			}
			
			if (DC_Debug > 2) {
				fprintf (DC_DbgFile, "   iCount: %2d  %12.6f%12.6f%12.6f\n\n", iCount,
				       PolVec.at(0), PolVec.at(1), PolVec.at(2));
			}
			
			// 
			// TotalPolarization = total
			// Dtot.at(0) = nPi*
			// Dtot.at(1) = PiPi*
			
			vector<double> Dtot (jCurGroup->NumberOfTransitions, 0.0);
			
			for (jTrans = 0; jTrans < jCurGroup->NumberOfTransitions; jTrans++)
				Dtot.at(jTrans) = VectorNorm (&Dxyz.at(jTrans));
				
			TotalPolarization = VectorNorm (&PolVec);
			
			// The total TDM's for the iCount'th transition are in Dtot, these are also broken
			// down into the components from each type of transition. Now we are in a position
			// to add the contribution from this state to the polarization tensor.
			for (j = 0; j < 3; j++)
				for (i = 0; i < 3; i++)
					for (k = 1; k < MaxNumberOfTransitions+1; k++)
						PolTensor.at(k).at(i).at(j) =
						       PolTensor.at(k).at(i).at(j)
						     + ( 2 * Dxyz.at(k).at(i) * Dxyz.at(k).at(j) / Energy );
			
			for (j = 0; j < 3; j++)
				for (i = 0; i < 3; i++)
						PolTensor.at(0).at(i).at(j) =
						       PolTensor.at(0).at(i).at(j)
						     + ( 2 * PolVec.at(i) * PolVec.at(j) / Energy );
			
			for (jTrans = 0; jTrans < jCurGroup->NumberOfTransitions; jTrans++)
				OscillatorStrength.at(iCount).at(jTrans) = 2 * Dtot.at(jTrans) / (3 * Energy);
			
			TotalOscillatorStrength = 2 * TotalPolarization / (3 * Energy);
			
			// --------------------------------------------------------------------------------
			
			if (DC_PrintPol) {
				// write the data for the current transition to the .pol file
				fprintf (DC_PolFile, "  %5d %10.3f %12.6f %12.6f %12.6f",
						iCount, Wavelength, PolVec.at(0), PolVec.at(1), PolVec.at(2));
				
				for (jTrans = 0; jTrans < jCurGroup->NumberOfTransitions; jTrans++)
					fprintf (DC_PolFile, " %12.6f", OscillatorStrength.at(iCount).at(jTrans));
				
				fprintf (DC_PolFile, " %12.6f\n", TotalOscillatorStrength);
				
				if ((iCount+1) % DC_System.NumberOfGroups == 0)
					Dichro::PrintPolarizationTensor (&PolTensor, MaxNumberOfTransitions);
			}
			
			// --------------------------------------------------------------------------------
			
			if (DC_PrintVec) {
				// write the .vec file (just the information required for absobance spectra)
				fprintf (DC_VecFile, "%8.3f %12.6f %12.6f %12.6f\n",
				        Wavelength, PolVec.at(0), PolVec.at(1), PolVec.at(2));
			}
			
			DC_Results.Trans.PolarizationVector.push_back (PolVec);
			DC_Results.Trans.OscillatorStrength.push_back (TotalOscillatorStrength);
			
			DC_Results.Groups.at(iGroup).PolarizationVector.push_back (PolVec);
			DC_Results.Groups.at(iGroup).OscillatorStrength.push_back (TotalOscillatorStrength);
			
			iCount++;
		} // of for (iTrans = 0; iTrans < iCurGroup->NumberOfTransitions; iTrans++) {
	} // of for (iGroup = 0; iGroup < DC_System.Groups.size(); iGroup++)
	
	if (DC_Verbose and DC_PrintVec) printf ("      Output written to %s\n", DC_VecFilename.c_str());
	return 0;
} // of Dichro::LD_Calculation


// ================================================================================


void Dichro::PrintPolarizationTensor ( vector< vector< vector<double> > >* PolTensor, int n )
{
	int i, j;
	double Average;
	
	if (not DC_PrintPol) return;
	
	for (i = 0; i < n-1; i++) {
		fprintf (DC_PolFile, "\n");
		fprintf (DC_PolFile, "   Polarization Tensor %d:\n", n);
	
		for (j = 0; j < 3; j++)
			fprintf (DC_PolFile, "      %12.3f %12.3f %12.3f\n",
			PolTensor->at(i).at(j).at(0),PolTensor->at(i).at(j).at(1),PolTensor->at(i).at(j).at(2));
	}
	
	// --------------------------------------------------------------------------------
	
	fprintf (DC_PolFile, "\n");
	fprintf (DC_PolFile, "   Total Polarization Tensor:\n");
	
	for (j = 0; j < 3; j++)
		fprintf (DC_PolFile, "      %12.3f %12.3f %12.3f\n",
		PolTensor->at(0).at(j).at(0), PolTensor->at(0).at(j).at(1), PolTensor->at(0).at(j).at(2));
	
	// --------------------------------------------------------------------------------
	
	fprintf (DC_PolFile, "\n");
	fprintf (DC_PolFile, "   Averaged Polarizations:\n");
	
	for (i = 0; i < n; i++) {
		Average = (   PolTensor->at(i).at(0).at(0)
		            + PolTensor->at(i).at(1).at(1)
		            + PolTensor->at(i).at(2).at(2) ) / 3;
		fprintf (DC_PolFile, "    %d   %12.3f\n", i, Average);
	}
	
	// --------------------------------------------------------------------------------
	
	fprintf (DC_PolFile, "\n");
	Average = (   PolTensor->at(0).at(0).at(0)
	            + PolTensor->at(0).at(1).at(1)
	            + PolTensor->at(0).at(2).at(2) ) / 3;
	fprintf (DC_PolFile, "    Total average:   %12.3f\n\n", Average);
	
	return;
} // of Dichro::PrintPolarizationTensor


// ================================================================================


