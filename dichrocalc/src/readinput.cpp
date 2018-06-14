// #################################################################################################
//
//  Program:      readinput.cpp
//
//  Function:     Part of DichroCalc:
//                Reads and parses the input data and all required parameter sets
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


// the class constructor
Dichro::Dichro ( string InFile, string Params, bool Verbose, int Debug,
                 bool PrintVec, bool PrintPol, bool PrintMat )
{
	// ===================================================================
	// Global configuration parameters
	// ===================================================================
	
	DC_PrintCdl      = true;               // whether to print the .cdl file
	DC_PrintXyzFiles = true;               // files with atom coordinates
	DC_ParamsDefault = "/bin/params";      // $HOME is added  
	
	// ===================================================================
	
	DC_InFile   = InFile;
	DC_Params   = Params;
	DC_Verbose  = Verbose;
	DC_Debug    = Debug;
	DC_PrintVec = PrintVec;
	DC_PrintPol = PrintPol;
	DC_PrintMat = PrintMat;
	
	DC_Input.Configuration.BBTrans = -1;
	DC_Input.Configuration.CTTrans = -1;
	
	DC_Error     = "";
	DC_ErrorCode = 0;
	
	DC_InFileBaseName = DC_InFile;
	
	if (DC_Params.size() == 0) {
		char* HOME = getenv ("HOME");      // read the environment variable
		DC_Params  = HOME;                 // convert from c-string to string
		DC_Params  += DC_ParamsDefault;    // set to ~/bin/params
	}
	
	if (DC_InFileBaseName.find (".inp") != string::npos) { // if the extension .inp is found
		// delete everything      from the position of .inp to the total length of the string
		DC_InFileBaseName.erase (DC_InFileBaseName.find (".inp"), DC_InFileBaseName.size());
	}
	
	if (DC_Debug > 0) {
		DC_DbgFilename = DC_InFileBaseName + ".dbg";
		DC_DbgFile = fopen (DC_DbgFilename.c_str(), "w");
	}
	
	if (DC_Debug > 4) { // write out the fitting data of FitParameters
		DC_FitFilename = DC_InFileBaseName + ".fit";
		DC_FitFile = fopen (DC_FitFilename.c_str(), "w");
	}
	
	if (DC_PrintVec) { // the polarization vectors for the absorbance/LD
		DC_VecFilename = DC_InFileBaseName + ".vec";
		DC_VecFile = fopen (DC_VecFilename.c_str(), "w");
	}
	
	if (DC_PrintPol) { // the polarizations broken into single transition contributions
		DC_PolFilename = DC_InFileBaseName + ".pol";
		DC_PolFile = fopen (DC_PolFilename.c_str(), "w");
	}
	
	if (DC_PrintMat) { // to print out the matrix, eigenvectors, and eigenvalues
		DC_MatFilename = DC_InFileBaseName + ".mat";
		DC_MatFile = fopen (DC_MatFilename.c_str(), "w");
	}
	
	if (DC_PrintCdl) { // to print out the matrix, eigenvectors, and eigenvalues
		DC_CdlFilename = DC_InFileBaseName + ".cdl";
		DC_CdlFile = fopen (DC_CdlFilename.c_str(), "w");
	}
	
	if (DC_Error == "") { ReadInput ();          }
	if (DC_Error == "") { CheckInputData ();     }
	if (DC_Error == "") { ReadParameters ();     }
	if (DC_Error == "") { FitParameters ();      }
	if (DC_Error == "") { HamiltonianMatrix ();  }
	if (DC_Error == "") { CD_Calculation ();     }
	if (DC_Error == "") { LD_Calculation ();     }
	
	if (DC_Debug > 1) Dichro::OutputSystemClass  ();
	if (DC_Debug > 0) Dichro::OutputResultsClass ();
	
	if (DC_Debug > 0) fclose (DC_DbgFile);
	if (DC_Debug > 4) fclose (DC_FitFile);
	
	if (DC_PrintCdl) fclose (DC_CdlFile);
	if (DC_PrintPol) fclose (DC_PolFile);
	if (DC_PrintVec) fclose (DC_VecFile);
	if (DC_PrintMat) fclose (DC_MatFile);
} // of Dichro::Dichro

// the class destructor
Dichro::~Dichro ( void )
{
} // of Dichro::~Dichro


// ================================================================================


int Dichro::ReadInput ( void )
// reads and parses the input file for a dichroism calculation
{
	ifstream InFile;
	vector<string> Fields;
	string Label, Line;
	int ErrorCode = 0;
	
	if (DC_Verbose) Dichro::NewTask ( "Reading Input File" );
	
	if (DC_Debug > 3) {
		Dichro::NewFileTask (DC_DbgFile, "Reading Input File");
		fprintf (DC_DbgFile, "   File: %s\n\n", DC_InFile.c_str());
	}
	
	InFile.open (DC_InFile.c_str(), ios::in);
	
	if ( not InFile ) {
		// try adding the .inp extension
		DC_InFile = DC_InFile + ".inp";
		InFile.open (DC_InFile.c_str(), ios::in);
		
		if ( not InFile ) {
			printf ("\nERROR: Could not open file %s.\n\n", DC_InFile.c_str());
			DC_Error = "Unable to open input file";
			DC_ErrorCode = 100;
			return 100;
		}
	}
	
	while ( not InFile.eof() ) {             // until the end of the file is reached
		Line = NextLine (&InFile);            // read the next line from the file
		
		// the following is by now also done by NextLine() but was left in for completeness
		if (Line.find ("#") == 0) continue;   // skip lines starting with # (comments)
		if (Line.length()   == 0) continue;   // skip empty lines
		
		if (Line.substr (0, 14) == "$CONFIGURATION")
			ErrorCode = Dichro::ReadInputSection ( &InFile, "$CONFIGURATION" );
		
		else if (Line.substr (0, 11) == "$PARAMETERS")
			ErrorCode = Dichro::ReadInputSection ( &InFile, "$PARAMETERS" );
		
		else if (Line.substr (0, 13) == "$CHROMOPHORES")
			ErrorCode = Dichro::ReadInputSection ( &InFile, "$CHROMOPHORES" );
		
		else if (Line.substr (0, 12) == "$COORDINATES")
			ErrorCode = Dichro::ReadInputSection ( &InFile, "$COORDINATES" );
		
		else {
			cerr << "\nERROR: Could not interpret this line in the input file "
			     << DC_InFile << ":\n" << Line << "\n\n";
			DC_Error = "Unknown line in input file";
			DC_ErrorCode = 103;
			return 103;
		}
		
		if (ErrorCode != 0) return ErrorCode;
	}
	
	InFile.close();
	
	if (DC_Debug > 3) Dichro::OutputInputClass();
	
	return 0;
} // of Dichro::ReadInput


// ================================================================================


int Dichro::ReadInputSection ( ifstream *InFile, string CurrentBlock )
// reads an parses a block ($BLOCK ... $END) in the input file
{
	string Line, Label;
	vector<string> Fields;
	unsigned int i, AtomIndex;
	
	while ( ( Line.compare(0, 4, "$END") != 0) && ( not InFile->eof() ) )
	{
		SplitNextLine (InFile, Line, Fields, " ");
		
		if ( Line.compare(0, 4, "$END") == 0 ) break;   // check for $END flag
		if ( Line.compare(0, 1, "$") == 0 )             // catch beginning of next block
		{
			cerr << "\nERROR: " << CurrentBlock << " was not terminated with $END flag.\n\n";
			DC_Error = "Missing $END in input file";
			DC_ErrorCode = 110;
			return 110;
		}
		
		Label = Fields.at(0);
		
		if (CurrentBlock.compare("$CONFIGURATION") == 0)
		{
			if (Fields.size() < 2) { ColumnError (DC_InFile, Line, 2); return 113; }
			
			     if ( StringInsCompare (Label, "BBTrans"))
				DC_Input.Configuration.BBTrans = atoi ( Fields.at(1).c_str() );
			
			else if ( StringInsCompare (Label, "CTTrans"))
				DC_Input.Configuration.CTTrans = atoi ( Fields.at(1).c_str() );
			
			else if ( StringInsCompare (Label, "Factor"))
				DC_Input.Configuration.Factor = atoi ( Fields.at(1).c_str() );
			
			else if ( StringInsCompare (Label, "MinWL"))
				DC_Input.Configuration.MinWL = atoi ( Fields.at(1).c_str() );
			
			else if ( StringInsCompare (Label, "MaxWL"))
				DC_Input.Configuration.MaxWL = atoi ( Fields.at(1).c_str() );
			else
			{
				cerr << "\nERROR: Unknown option " << Fields.at(0)
					<< " found in input file " << DC_InFile << ".\n";
				DC_Error = "Unknown option in input file";
				DC_ErrorCode = 112;
				return 112;
			}
		} // of if (CurrentBlock.compare("$CONFIGURATION") == 0)
		
		else if (CurrentBlock.compare("$PARAMETERS") == 0)
		{
			// check that there are indeed two or more fields
			// (comments may follow after parameter set and number of transitions)
			if ( Fields.size() < 2 )
			{
				cerr << "\nERROR: In the " << CurrentBlock
				     << " block, each parameter set is defined by its name\n";
				cerr << "       and number of transitions (2 fields). Error in line\n";
				cerr << "       " << Line << "\n\n";
				DC_Error = "Error in " + CurrentBlock + " block.";
				DC_ErrorCode = 114;
				return 114;
			}
			
			DC_Input.Parameters.Name.push_back  (Fields[0]);
			DC_Input.Parameters.Trans.push_back ( atoi (Fields[1].c_str()) );
		} // of if (CurrentBlock.compare("$PARAMETERS") == 0)
		
		else if (CurrentBlock.compare("$CHROMOPHORES") == 0)
		{
			if (Fields.size() < 4) { ColumnError (DC_InFile, Line, 4); return 115; }
			
			// the first column is the chromophore type, convert it to an integer value
			// this integer is the index of the chromophore in DC_Input.Parameters (both .Name and .Trans)
			DC_Input.Chromophores.Type.push_back ( atoi (Fields.at(0).c_str()) );
			Fields.erase( Fields.begin() );  // delete first element
			
			// create a new array for the atom numbers
			vector<int> NewChrom;
			
			// add the atom number one by one to the created array
			for (i = 0; i < Fields.size(); i++) {
				AtomIndex = atoi (Fields.at(i).c_str() ); // convert the field content to an integer
				--AtomIndex;                              // decrease to use it as array index
				NewChrom.push_back (AtomIndex);           // add it to the chromophore array
			}
				
			// add the array to the Atoms vector (i.e. it is a vector of vectors!)
			DC_Input.Chromophores.Atoms.push_back ( NewChrom );
		} // of if (CurrentBlock.compare("$CHROMOPHORES") == 0)
		
		else if (CurrentBlock.compare("$COORDINATES") == 0)
		{
			if (Fields.size() < 3) { ColumnError (DC_InFile, Line, 3); return 116; }
			
			// create a new array for the coordinates
			vector<double>  Coordinates;
			
			// add the coordinates one by one to the created array
			for (i = 0; i < 3; i++)
				Coordinates.push_back ( atof (Fields.at(i).c_str() ) );
			
			// add the array to the Groups vector (i.e. it is a vector of vectors!)
			DC_Input.Coordinates.Groups.push_back ( Coordinates );
			DC_Input.Coordinates.Atoms.push_back ( atoi(Fields.at(4).c_str())-1 ); // decreased by 1!
			DC_Input.Coordinates.Labels.push_back ( Fields.at(5).c_str() );
		} // of if (CurrentBlock.compare("$COORDINATES") == 0)
		
	} // of while ( ( Line.substr(0, 4) != "$END") && ( not InFile->eof() ) )
	
	if ( (not StringInsCompare (Line, "$END")) && InFile->eof() ) {
		cerr << "\nERROR: " << CurrentBlock << " was not terminated with $END flag.\n\n";
		DC_Error = "Missing $END in input file";
		DC_ErrorCode = 110;
		return 110;
	}
	
	return 0;
} // of Dichro::ReadInputSection


// ================================================================================


int Dichro::CheckInputData ( void )
// checks the input data for inconsistencies
{
	unsigned int Group, CurAtom, Atom, Index, Type;
	
	if (DC_Verbose) printf ("   Checking input data\n");
	
	for (Group = 0; Group < DC_Input.Chromophores.Type.size(); Group++) {
		Index = DC_Input.Chromophores.Type.at(Group);
		
		if ( Index > DC_Input.Parameters.Name.size() ) {
			cerr << "\nERROR: Chromophore type " << Index
			     << " referenced in $CHROMOPHORES block, but only "
			     << DC_Input.Parameters.Name.size()
			     << " types defined in $PARAMETERS block.\n\n";
			DC_Error = "Missing parameters in input file";
			DC_ErrorCode = 120;
			return 120;
		}
	} // of for (Group = 0; Group < DC_Input.Chromophores.Type.size(); Group++)
	
	for (Group = 0; Group <  DC_Input.Chromophores.Atoms.size(); Group++) {
		for (Atom = 0; Atom < DC_Input.Chromophores.Atoms.at(Group).size(); Atom++) {
			CurAtom = DC_Input.Chromophores.Atoms.at(Group).at(Atom);
			
			if (CurAtom > DC_Input.Coordinates.Groups.size()) {
				cerr << "\nERROR: Atom number " << CurAtom
				     << " is referenced in $CHROMOPHORES block, but only\n       "
				     << DC_Input.Coordinates.Groups.size()
				     << " atoms are present in $COORDINATES block.\n\n";
				DC_Error = "Missing atom coordinates in input file";
				DC_ErrorCode = 123;
				return 123;
			}
		}
	} // of for (Group = 0; Group < DC_Input.Chromophores.Type.size(); Group++)
	
	// loop over all parameter sets
	for (Type = 0; Type < DC_Input.Parameters.Trans.size(); Type++) {
		// if a specific backbone transition was requested and this is the first parameter set
		// type (i.e. the backbone parameters)
		if (DC_Input.Configuration.BBTrans > -1 and Type == 0) {
			if (DC_Input.Parameters.Trans.at(Type) != 1) {
				cerr << "\nERROR: A specific backbone transition was requested (BBTrans != -1)\n";
				cerr <<   "       but the number of transitions is != 1.\n\n";
				DC_Error = "Specific backbone transition requested and multiple transitions.";
				DC_ErrorCode = 126;
				return 126;
			}
		}
		
		// if a specific charge-transfer transition was requested and this is a CT chromophore
		if (DC_Input.Configuration.CTTrans > -1 and
		           DC_Input.Parameters.Name.at(Type).substr (0, 2) == "CT") {
			if (DC_Input.Parameters.Trans.at(Type) != 1) {
				cerr << "\nERROR: "
					  << "A specific charge-transfer transition was requested (CTTrans != -1)\n";
				cerr <<   "       but the number of transitions is != 1.\n\n";
				DC_Error = "Specific CT transition requested and multiple transitions.";
				DC_ErrorCode = (127);
				return (127);
			}
		}
	} // of for (Type = 0; Type < DC_Parameters.Name.size(); Type++) {
	
	// DEBUG OUTPUT
	// the individual blocks do not necessarily have to be in read in this order, but
	// need to be read in first as, for example, DC_Input.Chromophores needs information
	// from DC_Input.Parameters
	// Dichro::OutputInputConfigurationClass();
	// Dichro::OutputInputParametersClass();
	// Dichro::OutputInputChromophoresClass();
	// Dichro::OutputInputCoordinatesClass();
	
	return 0;
} // of Dichro::CheckInputData


// ================================================================================


int Dichro::ReadParameters ( void )
// reads and parses all parameter sets specified in the input file
{
	unsigned int i, CurFile, State, Trans, FilePos;
	bool Permanent;
	string Filename, Line;
	vector<string> Fields;
	vector<string> FileLines;
	vector< vector<string> > FileFields;
	
	if (DC_Verbose)   Dichro::NewTask ( "Reading Parameter Files" );
	if (DC_Debug > 3) Dichro::NewFileTask ( DC_DbgFile, "Reading Parameter Files" );
	
	// if no directory with the parameters is given, use the current directory
	// (should have been changed by the constructor, just to keep this routine general)
	if (DC_Params.size() == 0) DC_Params = ".";
	
	// if a slash (/) is found at the end of the directory name
	if (DC_Params.rfind("/") == DC_Params.size()-1)
		// delete the last character of the string
		DC_Params.erase(DC_Params.size()-1, 1);
	
	vector<string> ParFiles;
	
	if (not ReadDir (DC_Params, ".par", &ParFiles)) {
		cerr << "\nERROR: Could not read directory " << DC_Params << "\n\n";
		DC_Error = "Error reading directory with parameter files";
		DC_ErrorCode = 130;
		return 130;
	}
	
	if (DC_Debug > 3) {
		fprintf (DC_DbgFile, "   Reading parameter files in directory %s\n", DC_Params.c_str() );
		
		if (DC_Debug > 4) {
			fprintf (DC_DbgFile, "   %lu files found:\n", ParFiles.size() );
			
			for (i = 0; i < ParFiles.size(); i++)
				fprintf (DC_DbgFile, "   %s\n", ParFiles.at(i).c_str() );
		}
		
		fprintf (DC_DbgFile, "\n");
		
		if (DC_Input.Parameters.Name.size() == 1)
			fprintf (DC_DbgFile,
			         "   %lu parameter set is read in.\n\n", DC_Input.Parameters.Name.size() );
		else
			fprintf (DC_DbgFile,
			         "   %lu parameter sets are read in.\n\n", DC_Input.Parameters.Name.size() );
	}
	
	// run over each parameter set given in the input file, read and parse it into
	// the DC_ParSets vector
	for (CurFile = 0; CurFile < DC_Input.Parameters.Name.size(); CurFile++) {
		Filename = DC_Params + "/" + DC_Input.Parameters.Name.at(CurFile) + ".par";
		
		if ( not FileExists (Filename) ) {
			cerr << "\nERROR: Parameter file " << Filename << " not found.\n\n";
			DC_Error = "Parameter file " + Filename + " not found";
			DC_ErrorCode = 131;
			return 131;
		}
		
		if (DC_Verbose) printf ("   Reading %s\n", Filename.c_str() );
		if (DC_Debug > 3)
			fprintf (DC_DbgFile, "   Reading %s", Filename.c_str() );
		
		// create a new instance of the ParSet class
		ParSet CurParSet;
		
		CurParSet.Name = DC_Input.Parameters.Name.at(CurFile);
		
		// Charge-transfer parameters contain the four local transitions at first, which have
		// to be skipped later on. Therefore, a parameter set is marked as charge-transfer, if
		// the first two letters of the name are 'CT', this makes it easier later.
		if (CurParSet.Name.substr (0, 2) == "CT")
			CurParSet.ChargeTransfer = true;
		else
			CurParSet.ChargeTransfer = false;
		
		ifstream File;          // reinitialize each time, otherwise File.eof() remains true on Linux
		File.open (Filename.c_str(), ios::in);
		FileLines.clear();      // all lines of the current file as strings
		FileFields.clear();     // all lines of the current file as columns
		
		// The matrix method parameter sets were created for the use with a FORTRAN
		// program and are designed to be read on a line-by-line bases instead of
		// a block-wise fashion.
		
		// read and parse the complete file
		while ( not File.eof() ) {
			// split each line into columns and add to FileFields and FileLines
			SplitNextLine (&File, Line, Fields, " ");
			FileLines.push_back (Line);
			FileFields.push_back (Fields);
		}
		
		File.close();   // only FileFields and FileLines are used from now on
		
		if (FileFields.size() == 0) {
			cerr << "\nERROR: Could not read file " << Filename.c_str() << "\n\n";
			DC_Error = "Error reading in parameter set file";
			DC_ErrorCode = 132;
			return 132;
		}
		
		FilePos = 0;    // generally used to access a specific position in the file content
		Fields = FileFields.at(FilePos);
		
		// the first line should contain the filename (case-sensitive!)
		if (Fields.at(0).find (CurParSet.Name) == string::npos) {
			cerr << "\nERROR: In file " << Filename << " the first line does not contain\n"
			     <<   "       the parameter set name " << CurParSet.Name << ".\n\n";
			DC_Error = "Format error in parameter set file.";
			DC_ErrorCode = 133;
			return 133;
		}
		
		// --------------------------------------------------------------------------------
		
		// the second line contains the number of atoms to be read in the following lines
		Fields = FileFields.at(++FilePos); // FIRST increase FilePos and THEN get the line
		
		CurParSet.NumberOfAtoms = atoi ( Fields.at(0).c_str() );
		
		if (CurParSet.NumberOfAtoms == 0) {
			cerr << "\nERROR: In file " << Filename
			     << " the number of atoms could not be interpreted in line\n"
			     <<   FileLines.at(FilePos) << "\n\n";
			DC_Error = "Format error in parameter set file.";
			DC_ErrorCode = 134;
			return 134;
		}
		
		// initialize the coordinates of the reference point
		for (i = 0; i < 3; i++) CurParSet.Reference.push_back (0);
		double Weighting = 0;
		int Atom;
		
		// now read as many atoms as stated in the line before
		for (Atom = 0; Atom < CurParSet.NumberOfAtoms; Atom++) {
			Fields = FileFields.at(++FilePos); // FIRST increase FilePos and THEN get the line
			if (Fields.size() < 6) {
				ColumnError (CurParSet.Name, FileLines.at(FilePos), 6);
				return 135;
			}
			
			ParSetAtom Atom;
			Atom.Coord.push_back ( atof (Fields.at(0).c_str() ) );
			Atom.Coord.push_back ( atof (Fields.at(1).c_str() ) );
			Atom.Coord.push_back ( atof (Fields.at(2).c_str() ) );
			Atom.Weighting =       atof (Fields.at(3).c_str() );
			Atom.Label     = Fields.at(5);
			CurParSet.Atoms.push_back (Atom);
			
			CurParSet.Reference.at(0) += Atom.Coord.at(0) * Atom.Weighting;
			CurParSet.Reference.at(1) += Atom.Coord.at(1) * Atom.Weighting;
			CurParSet.Reference.at(2) += Atom.Coord.at(2) * Atom.Weighting;
			Weighting = Weighting + Atom.Weighting;
		}
		
		// divide each coordinate of the reference point by the weighting factor
		for (i = 0; i < 3; i++) CurParSet.Reference.at(i) /= Weighting;
		
		// --------------------------------------------------------------------------------
		
		Fields = FileFields.at(++FilePos); // FIRST increase FilePos and THEN get the line
		
		// after the atoms, a state should follow, starting with &TRANSITION
		if (Fields.at(0).find ("&TRANSITION") == string::npos) {
			cerr << "\nERROR: Label &TRANSITION expected in file "
			     << Filename << " in line\n       " << FileLines.at(FilePos) << "\n\n";
			DC_Error = "Format error in parameter set";
			DC_ErrorCode = 137;
			return 137;
		}
		
		State = 0;
		
		// this loop takes care of the states (started with a &TRANSITION label)
		while ( (FilePos+1 < FileFields.size() )     // 'EOF' not reached and
		         and         // and current Line contains &TRANSITION or &PERMANENT
			     ( (FileLines.at(FilePos).find ("&TRANSITION") != string::npos)
			       or (FileLines.at(FilePos).find ("&PERMANENT") != string::npos) ) ) {
			
			if (FileLines.at(FilePos).find ("&PERMANENT") != string::npos)
				Permanent = true;
			else
				Permanent = false;
			
			++FilePos; // advance to the next line, the start of the first transition of the state
			
			vector<ParSetTrans> NewState;
			CurParSet.States.push_back (NewState);
			Trans = 0;
			
			// this loop handles single transitions within a state
			while ( (FilePos+1 < FileFields.size() )  // 'EOF' not reached and
			         and         // and current line contains neither &TRANSITION nor &PERMANENT
			        ( (FileLines.at(FilePos).find ("&TRANSITION") == string::npos)
			          and (FileLines.at(FilePos).find ("&PERMANENT") == string::npos) ) ) {
			
				++Trans;
				
				// DEBUG OUTPUT
				// if (DC_Debug > 5) {
				// 	if (Permanent)
				// 		printf ("      - permanent moments of state %d\n", Trans);
				// 	else
				// 		printf ("      - state %d, transition %d\n", State, Trans);
				// }
				
				// create a new instance of the transition class and read the next transition
				ParSetTrans CurTrans;
				Dichro::ReadTransition ( FileLines, FileFields, CurParSet.Name,
				                         &FilePos, &CurTrans, Trans, Permanent );
				CurParSet.States.at(State).push_back ( CurTrans );
				
				++FilePos;  // go to the next line for the while loop checks
			} // of inner while loop
			
			++State;
		} // of outer while loop
		
		// DEBUG OUPUT: print the complete data
		if (DC_Debug > 4) Dichro::OutputParSetClass ( &CurParSet );
		
		// save the current set with the same index as in the input file
		DC_ParSets.push_back (CurParSet);
	} // of for (i = 0; i < DC_Input.Parameters.Name.size(); i++)
	
	return 0;
} // of Dichro::ReadParameters


// ================================================================================


int Dichro::ReadTransition ( vector<string> FileLines,  vector< vector<string> > FileFields,
                             string ParSetName,  unsigned int *FilePos,  ParSetTrans *CurTrans,
                             int Trans,  bool Permanent )
// reads a complete transition from a parameter set file
{
	int i;
	string Line;
	vector<string> Fields;
	
	// the first line of a transition contains the number of monopoles and the energy
	Fields = FileFields.at(*FilePos);  // no increase necessary here, it is already the first line
	if (Fields.size()<2) {Dichro::ColumnError (ParSetName,FileLines.at(*FilePos),2); return 135;}
	
	CurTrans->NumberOfMonopoles = atoi ( Fields.at(0).c_str() );
	CurTrans->Energy         = atof ( Fields.at(1).c_str() );
	if (CurTrans->Energy != 0)
		CurTrans->Wavelength     = 1E7 / CurTrans->Energy;
	else
		CurTrans->Wavelength     = 0;
	
	// second line is the electric transition dipole moment and a scale factor
	Fields = FileFields.at(++*FilePos);  // FIRST increase FilePos and THEN get the line
	if (Fields.size()<3) {Dichro::ColumnError (ParSetName,FileLines.at(*FilePos),3); return 135;}
	
	CurTrans->EDM.push_back ( atof (Fields.at(0).c_str()) );
	CurTrans->EDM.push_back ( atof (Fields.at(1).c_str()) );
	CurTrans->EDM.push_back ( atof (Fields.at(2).c_str()) );
	
	// the permanent moment do not have scale factors and magnetic dipole moments
	if (Permanent)
		CurTrans->Permanent = true;
	else {
		CurTrans->Permanent = false;
		
		if (Fields.size()<4) {Dichro::ColumnError (ParSetName,FileLines.at(*FilePos),4); return 135;}
		
		CurTrans->ScaleFactor = atof (Fields.at(3).c_str());
		
		// second line is the magnetic transition dipole moment
		Fields = FileFields.at(++*FilePos);  // FIRST increase FilePos and THEN get the line
		if (Fields.size()<3) {Dichro::ColumnError (ParSetName,FileLines.at(*FilePos),3); return 135;}
		CurTrans->MDM.push_back ( atof (Fields.at(0).c_str()) );
		CurTrans->MDM.push_back ( atof (Fields.at(1).c_str()) );
		CurTrans->MDM.push_back ( atof (Fields.at(2).c_str()) );
	}
	
	// now a monopole on each line, as many as specified before
	for (i = 0; i < CurTrans->NumberOfMonopoles; i++) {
		ParSetMonopole Mono;
		Fields = FileFields.at(++*FilePos);  // FIRST increase FilePos and THEN get the line
		if (Fields.size()<4) {Dichro::ColumnError (ParSetName,FileLines.at(*FilePos),4); return 135;}
		Mono.Coord.push_back ( atof (Fields.at(0).c_str()) );
		Mono.Coord.push_back ( atof (Fields.at(1).c_str()) );
		Mono.Coord.push_back ( atof (Fields.at(2).c_str()) );
		Mono.Charge =          atof (Fields.at(3).c_str());
		
		CurTrans->Monopoles.push_back (Mono);
	}
	
	// DEBUG OUTPUT: print the whole transition data
	// not necessarily needed, complete parameter set is printed out at the end of parsing
	// if (DC_Debug > 5) Dichro::OutputParSetTransClass (CurTrans, Trans);
	
	return 0;
} // of Dichro::ReadTransition


// ================================================================================


int Dichro::ColumnError (string File, string Line, int Columns)
{
	cerr << "\nERROR: In file " << File << " at least "
	     << Columns << " columns were expected in line\n       " << Line << "\n\n";
	DC_Error = "Format error in file " + File;
	DC_ErrorCode = 135;
	return 135;
} // of Dichro::ColumnError


// ================================================================================


