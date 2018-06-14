// #################################################################################################
//
//  Header:       dichrocalc.h
//
//  Author:       Benjamin M. Bulheller
//
//  Version:      $Revision: 4602 $, $Date: 2009-06-27 05:16:06 -0700 (Sat, 27 Jun 2009) $
//
// #################################################################################################

#include <stdio.h>         // standard input output header file
#include <iostream>        // provides file streams
#include <string.h>        // better strings (well, strings at all!)
#include <vector>          // better arrays
#include <iomanip>         // for printf-like format switches for cout
#include <fstream>         // for using file streams to read/write files
#include <sstream>         // for using string streams to handle strings
#include <regex.h>         // for the use of regular expressions
#include <dirent.h>        // to read directory listings
#include <errno.h>         // error handling

using namespace std;       // to abbreviate e.g. std::cout

// required includes for the use of the NewMat libraries
#define WANT_MATH          // include.h will get math functions
#include <newmat.h>        // the main newmat header file
#include <newmatio.h>      // matrix input output routines
#include <newmatap.h>      // matrix input output routines


// ================================================================================
// Class definition and declarations of object functions
// ================================================================================

class Dichro {
	public:
		
		// --------------------------------------------------------------------------
		// general global variables of the object
		// --------------------------------------------------------------------------
		
		string DC_InFile;            // .inp file with the protein and chromophore information
		string DC_Params;            // directory name with parameter files
		bool   DC_Verbose;           // give output for each step of the calculations true/false
		bool   DC_PrintVec;          // print the .vec file (for absorbance/LD)
		bool   DC_PrintPol;          // print the .pol file (transition polarizations)
		bool   DC_PrintMat;          // print the .mat file (matrix, eigenvectors, eigenvalues)
		bool   DC_PrintCdl;          // print the .mat file (matrix, eigenvectors, eigenvalues)
		bool   DC_PrintXyzFiles;     // create two files with the atom coordinates
		int    DC_Debug;             // set output level, 1--5, the higher the more output
		
		FILE *DC_CdlFile;            // file for the CD line spectrum
		FILE *DC_DbgFile;            // file for debugging information
		FILE *DC_FitFile;            // file for the debugging output of the fitting
		FILE *DC_PolFile;            // file for the polarization calculation
		FILE *DC_VecFile;            // file for the output of the polarization vectors
		FILE *DC_MatFile;            // file for the matrix, eigenvectors, and eigenvalues

		string DC_DbgFilename;       // file for debugging information
		string DC_FitFilename;       // file for the debugging output of the fitting
		string DC_PolFilename;       // file for the polarization calculation
		string DC_VecFilename;       // file for the output of the polarization vectors
		string DC_MatFilename;       // file for the matrix, eigenvectors, and eigenvalues
		string DC_CdlFilename;       // file for the CD line spectrum
		
		string DC_Error;             // message of a fatal error
		int    DC_ErrorCode;         // the error code returned from the died-off function
		vector<string> Warnings;     // non-fatal problems
		string DC_InFileBaseName;    // the base name of the input file
		string DC_ParamsDefault;     // the default directory of parameter files
		
		// general limits for the calculations
		// the maximum number of transitions than can be on a group
		static const int DC_MaxGroupTransitions = 20;
		
		// --------------------------------------------------------------------------
		// classes for the data read from the input file (.inp)
		// --------------------------------------------------------------------------
		
		class InputConfiguration  {  // the $CONFIGURATION block
			public:
				int BBTrans;           // -1 = all backbone transitions, 0-3 = a specific BB transition
				int CTTrans;           // -1 = all CT transitions, 0-3 = a specific CT transition
				int Factor;            // scale factor for intensities (number of residues)
				int MinWL;             // minimum wavelength
				int MaxWL;             // maximum wavelength
		};
		
		class InputParameters {      // the $PARAMETERS block
			public:
				vector<string> Name;   // the names of the parameter sets
				vector<int>    Trans;  // the number of transition on each parameter set
		};
		
		class InputChromophores {    // the $CHROMOPHORES block
			public:
				vector<int> Type;              // the index in DC_Parameters
				vector< vector<int> > Atoms;   // the atom indices in DC_Input.Coordinates
		};
		
		class InputCoordinates {        // the $COORDINATES block
			public:
				vector< vector<double> > Groups;   // the xyz coordinates
				vector<string>  Labels;   // the PDB atom labels
				vector<int>     Atoms;    // the array indices, identical to .Chromophores vector
		};
		
		class Input { // combines all information read from the .inp file
			public:
				InputConfiguration Configuration;
				InputParameters    Parameters;
				InputChromophores  Chromophores;
				InputCoordinates   Coordinates;
		} DC_Input;
		
		
		// --------------------------------------------------------------------------
		// variables for the parameter sets
		// --------------------------------------------------------------------------
		
		class ParSetAtom {      // a single atom in a parameter set
			public:
				string Label;           // the PDB label of this atom, just used for debug output
				double Weighting;       // a weighting factor
				vector<double> Coord;   // x, y, z coordinates, usually in Angstrom
		};
		
		class ParSetMonopole {  // a single monopole in a transition
			public:
				vector <double> Coord;     // the x, y, z coordinates in Angstrom
				double Charge;             // the charge in 10^-19 esu
		};
		
		class ParSetTrans {     // a single transition in a ParSet class
			public:
				bool  Permanent;          // whether it is a permanent moment true/false
				int   NumberOfMonopoles;  // the number of monopoles of the transition
				double Energy;            // the excitation energy in cm^-1
				double Wavelength;        // the energy converted to nm
				double ScaleFactor;       // a scale factor for the dipole moments
				vector<double> EDM;       // coordinates of the electric trans. dip. mom.
				vector<double> MDM;       // coordinates of the magnetic trans. dip. mom.
				vector<ParSetMonopole> Monopoles;  // coordinates and charges of the monopoles
		};
		
		class ParSet {          // the class instantiated for each parameter set read from a .par
			public:
				string  Name;                // the name of the parameter set (also filename)
				int     NumberOfAtoms;       // the number of atoms
				int     NumberOfTransitions; // the number of transitions
				bool    ChargeTransfer;      // CT groups have four additional transitions
				vector<ParSetAtom> Atoms;    // atom coordinates and scale factors
				vector<double> Reference;    // the reference coordinate for this chromophore
				vector< vector<ParSetTrans> > States;
		};
		
		vector<ParSet> DC_ParSets;   // holds all information read from the .par files
		
		
		// --------------------------------------------------------------------------
		// class to describe the initial system
		// --------------------------------------------------------------------------
		
		class SystemTransition {  // all information describing one particular transition of a group
			public:
				bool   Permanent;         // whether it is a permanent moment true/false
				double  Energy;           // the excitation energy in cm^-1
				double  Wavelength;       // the energy converted to nm
				string Origin;            // the parameter set, state and transition number it came from
			// trans   coords
				vector<double> EDM;       // coordinates of the electric trans. dip. mom.
				vector<double> MDM;       // coordinates of the magnetic trans. dip. mom.
			// trans   monopole
				int NumberOfMonopoles;            // the number of monopoles of the transition
				vector<ParSetMonopole> Monopoles; // coordinates and charges of the monopoles
		};
		
		class SystemGroup {
			public:
				int NumberOfAtoms;              // number of atoms in this group
				int NumberOfTransitions;        // the number of transitions on this group
				bool ChargeTransfer;            // charge-transfer group true/false
				string ParameterSet;            // the name of the used parameter set
				vector<double> Reference;       // the reference vector for this group
				vector<SystemTransition> Trans;  // all transitions on this group
				vector<SystemTransition> Perm;   // all permanent moments of this group
				vector< vector<double> > Atoms; // all atoms of this group
				vector<int> AtomIndices;        // the indices in the DC_Input.Coordinates array
		};
		
		class System {       // all information necessary to describe the system
			public:
				int NumberOfAtoms;        // number of atoms (not PDB but of all parameter sets)
				int NumberOfGroups;       // the total number of groups in the system
				int NumberOfTransitions;  // the total number of transitions on all groups
				int MatrixDimension;      // technically this is the same as the NumberOfTransitions
				
			// atom    coords
				vector< vector<double> > Atoms;
			// group   trans
				vector< SystemGroup > Groups;
		} DC_System;
		
		
		// --------------------------------------------------------------------------
		// class for the results (all derived from DC_System)
		// --------------------------------------------------------------------------
		
		class ResultsGroup {   // all results for access by the group
			public:
				int NumberOfTransitions;           // the number of transitions on this group
				bool ChargeTransfer;               // charge-transfer group true/false
				string ParameterSet;               // the name of the used parameter set
				vector<double> Reference;          // the reference vector for this group
				vector<double> Energy;             // the excitation energy in cm^-1
				vector<double> Wavelength;         // the energy converted to nm
				vector<double> DipoleStrength;     //
				vector<double> RotationalStrength; //
				SymmetricMatrix Submatrix;
				
			// trans   coords
				vector< vector<double> > EDM;      // EDM coordinates for each trans. of the group
				vector< vector<double> > MDM;      // MDM coordinates for each trans. of the group
				vector< vector<double> > PolarizationVector;  // polarization vectors
				vector<double> OscillatorStrength;  //
		};
		
		class ResultsTrans {   // all results for access by the transition
			public:
				vector<double> Energy;              // the excitation energy in cm^-1
				vector<double> Wavelength;          // the energy converted to nm
				vector<double> DipoleStrength;      //
				vector<double> RotationalStrength;  //
				vector<double> OscillatorStrength;  //
				
				vector<int> GroupSequence;      // the sequence of the groups along the matrix diagonal
				vector<int> TransSequence;      // the sequence of the transitions within the groups
				vector<string> ParSetSequence;  // the name of the used parameter set
				
			// trans   coords
				vector< vector<double> > EDM;       // coordinates of the electric trans. dip. mom.
				vector< vector<double> > MDM;       // coordinates of the magnetic trans. dip. mom.
				vector< vector<double> > MDMconv;   // magnetic trans. dip. mom. before mixing
				vector< vector<double> > Reference; // the reference vector for each group
				
				vector< vector<double> > PolarizationVector; // polarization vectors
		};
		
		class Results {     // all data calculated from DC_System
			public:
				int NumberOfAtoms;        // number of atoms (not PDB but of all parameter sets)
				int NumberOfGroups;       // the total number of groups in the system
				int NumberOfTransitions;  // the total number of transitions on all groups
				int MatrixDimension;      // technically this is the same as the NumberOfTransitions
				
				SymmetricMatrix Hamiltonian;
				DiagonalMatrix  Eigenvalues;
				Matrix          Eigenvectors;
				
				vector<ResultsGroup> Groups;   // all information sorted for access by the group
				ResultsTrans Trans;            // all information sorted for access by the trans.
				
				// results of a polarization calculation
				vector< vector< vector<double> > > PolTensor; // the polarization tensor
		} DC_Results;
		
		
		// --------------------------------------------------------------------------
		// declarations of the internal functions
		// --------------------------------------------------------------------------
		
		// iolibrary.cpp
		void NewTask ( string Message );
		void NewFileTask ( FILE* File, string Message );
		void OutputHeadline ( string Headline );
		void OutputFileHeadline ( FILE* File, string Headline, bool Lines = true );
		void OutputFileSeparator ( FILE* File, int Indent = 0 );
		void OutputInputConfigurationClass ( void );
		void OutputInputParametersClass ( void );
		void OutputInputChromophoresClass ( void );
		void OutputInputCoordinatesClass ( void );
		void OutputInputClass ( void );
		void OutputParSetMonopoleClass ( ParSetMonopole* Monopole );
		void OutputParSetTransClass ( ParSetTrans* Transition, int Trans = -1,
		                              int State = -1, string Name = "#");
		void OutputParSetClass ( ParSet* CurParSet );
		void OutputSystemTransitionClass ( SystemTransition* CurTrans, int Group, int Trans );
		void OutputSystemGroupClass ( SystemGroup* CurGroup, int Group = -1 );
		void OutputSystemClass ( void );
		void OutputResultsTransClass ( void );
		void OutputResultsGroupClass ( void );
		void OutputResultsClass ( void );
		void OutputSystemData ( void );
		void DebugOutput ( vector< vector<double> >* Vector );
		
		// readinput.cpp
		int  ReadInput ( void );
		int  ReadInputSection ( ifstream *File, string CurrentBlock );
		int  CheckInputData ( void );
		int  ReadParameters ( void );
		int  ReadTransition ( vector<string> FileLines, vector< vector<string> > FileFields,
		                      string ParSetName, unsigned int *FilePos, ParSetTrans *CurTrans,
		                      int Trans, bool Permanent );
		int  ColumnError (string ParSet, string Line, int Columns);
		
		// fitparameters.cpp
		int  FitParameters ( void );
		int  RotationMatrix ( int Chrom, Matrix ParSetMatrix, Matrix GroupMatrix,
		                      Matrix *RotMatrixNonUnitary, Matrix* RotMatrixUnitary );
		void CheckPlanar ( Matrix* ParSetMatrix, Matrix* GroupMatrix );
		void Rotate ( vector<double> *In, vector<double> *Out, Matrix *RotMatrix );
		
		// matrix.cpp
		int    HamiltonianMatrix ( void );
		double SameGroupInteraction ( int iGroup, int iTrans, int jGroup, int jTrans );
		double DifferentGroupInteraction (int iGroup, int iTrans, int jGroup, int jTrans, bool Perm);
		bool   GroupsOverlap ( int iGroup, int jGroup );
		
		// dichroism.cpp
		int  CD_Calculation ( void );
		int  LD_Calculation ( void );
		void PrintPolarizationTensor ( vector< vector< vector<double> > >* PolTensor, int n );
		
	public:
		Dichro  ( string InFile, string Params, bool Verbose, int Debug,
		          bool PrintVec = false, bool PrintPol = false, bool PrintMat = false );
		~Dichro ( void );
};


// ================================================================================


#include "iolibrary.h"

