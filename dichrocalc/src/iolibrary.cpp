// #################################################################################################
//
//  Program:      iolibrary.cpp
//
//  Function:     Part of DichroCalc:
//                Input/output routines used by several parts of dichrocalc
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


bool FileExists ( string FileName )
// returns true if a file exists and false if it does not
{
	FILE* FilePointer = NULL;
	
	FilePointer = fopen ( FileName.c_str(), "rb" );
	
	if( FilePointer != NULL ) {
		fclose ( FilePointer );
		return true;
	}
	
	return false;
} // of FileExists


// ================================================================================


bool FileExtension ( string Filename, string Extension )
// checks a filename for a given file extension
{
	int Length = Filename.size();
	unsigned int ExtLength = Extension.size();
	
	if ( Filename.rfind(".par") == Length-ExtLength )
		return true;
	else
		return false;
} // of FileExtension


// ================================================================================


bool ReadDir ( string Dir, string Extension, vector<string>* Files )
// returns a vector with all files of a certain extension contained in a directory
{
	DIR *DirPointer;
	struct dirent *DirP;
	string File;
	
	if ( (DirPointer = opendir (Dir.c_str() ) ) == NULL ) { return false; }
	
	while ( (DirP = readdir (DirPointer) ) != NULL) {
		File = string(DirP->d_name);
		if (FileExtension (File, Extension))  // if the file extension is correct
			Files->push_back ( File );         // add the file to the array of parameter files
	}
	
	closedir (DirPointer);
	// PrintVector ( Files );
	return true;
} // of ReadDir


// ================================================================================


void VectorDiff (vector<double>* Vect1, vector<double>* Vect2, vector<double>* DiffVector)
// calculates the difference between two vectors
{
	for (int Coord = 0; Coord < 3; Coord++)
		DiffVector->at(Coord) = Vect1->at(Coord) - Vect2->at(Coord);
} // of VectorDiff


// ================================================================================


double VectorNorm (vector<double>* Vector)
// calculates the norm of a vector
{
	double Norm = 0.0;
	
	for (int Coord = 0; Coord < 3; Coord++)
		Norm = Norm + pow (Vector->at(Coord), 2);
	
	return sqrt(Norm);
} // of VectorNorm


// ================================================================================


void CrossProduct (vector<double>* Vect1, vector<double>* Vect2, vector<double>* Cross)
// calculates the difference between two vectors
{
	Cross->at(0) = (Vect1->at(1) * Vect2->at(2)) - (Vect2->at(1) * Vect1->at(2));	
	Cross->at(1) = (Vect1->at(2) * Vect2->at(0)) - (Vect2->at(2) * Vect1->at(0));	
	Cross->at(2) = (Vect1->at(0) * Vect2->at(1)) - (Vect2->at(0) * Vect1->at(1));	
} // of CrossProduct


// ================================================================================


double PointDistance (vector<double>* Vect1, vector<double>* Vect2)
// calculates the distance between two points
{
	return sqrt( pow((Vect1->at(0) - Vect2->at(0)), 2) +
	             pow((Vect1->at(1) - Vect2->at(1)), 2) +
	             pow((Vect1->at(2) - Vect2->at(2)), 2) );
} // of PointDistance


// ================================================================================


string tostring ( unsigned int Integer )
// cast to a string
{
	stringstream ss;
	ss << Integer;
	return ss.str();
}


string tostring ( int Integer )
// cast to a string
{
	stringstream ss;
	ss << Integer;
	return ss.str();
}


string tostring ( double Double )
// cast to a string
{
	stringstream ss;
	ss << Double;
	return ss.str();
}

// ================================================================================


void Dichro::NewTask ( string Message )
{
	cout << "\n-> " << Message << endl;
	// underline the message with a double line
	// for (unsigned int i = 0; i < Message.length(); i++) { cout << "="; }
	// cout << endl;
	
	return;
} // of Dichro::NewTask


void Dichro::NewFileTask ( FILE* File, string Message )
{
	if (DC_Debug < 1) {
		cerr << "\nERROR: The debug option is required to use Dichro::NewFileTask.\n\n";
		return;
	}
	
	fprintf (File,
	"\n\n################################################################################\n\n");
	
	fprintf (File, "\n%s\n", Message.c_str());
	// underline the message with a double line
	for (unsigned int i = 0; i < Message.length(); i++) { fprintf (File, "="); }
	fprintf (File, "\n\n");
	
	return;
} // of Dichro::NewFileTask


// ================================================================================


void Dichro::OutputHeadline ( string Headline )
// prints an underlined headline with some vertical space
{
	printf ("\n\n%s\n", Headline.c_str());
	// underline the message with a single line
	for (unsigned int i = 0; i < Headline.length(); i++) { cout << "-"; }
	cout << endl;
	
	return;
} // of Dichro::OutputHeadline


void Dichro::OutputFileHeadline ( FILE* File, string Headline, bool Lines )
// prints an underlined headline with some vertical space
{
	if (Lines) fprintf (File, "\n\n");
	
	fprintf (File, "%s\n", Headline.c_str());
	
	// check for leading spaces
	int i = Headline.find_first_not_of (" ");          // count leading blanks
	for (int j = 0; j < i; j++) fprintf (File, " ");   // print as many as found
	TrimSpaces (&Headline);                            // remove all leading blanks
	
	// underline the message with a double line
	for (unsigned int i = 0; i < Headline.length(); i++) { fprintf (File, "-"); }
	fprintf (File, "\n\n");
	
	return;
} // of Dichro::OutputFileHeadline


void Dichro::OutputFileSeparator ( FILE* File, int Indent )
// prints a separator line into a file
{
	if (DC_Debug < 1) {
		cerr << "\nERROR: The debug option is required to use Dichro::OutputFileSeparator.\n\n";
		return;
	}
	
	fprintf (File, "\n");
	
	// indent the line with some blanks
	for (int i = 0; i < Indent; i++) fprintf (File, " ");
	
	fprintf (File,
	         "--------------------------------------------------------------------------------\n");
} // of Dichro::OutputFileSeparator


// ================================================================================


void CString2String ( char* cLine, string* Line )
// converts a C-style string into a string
{
	*Line = "";                   // clear the string
	*Line = *Line + cLine;        // convert to a real string
} // of CString2String


// ================================================================================


void TrimSpaces ( string* Line )
// removes leading and trailing blanks and tabs
{
	// find the first character position after excluding leading blanks and tabs
	size_t FirstChar = Line->find_first_not_of (" \t");
	
	// find the first character position from the end of the string
	size_t LastChar  = Line->find_last_not_of (" \t");
	
	if ( ( string::npos == FirstChar ) || ( string::npos == LastChar) )
		// return an empty string if it was all blanks
		*Line = "";
	else
		// return the substring between the first and the last character
		*Line = Line->substr ( FirstChar, LastChar-FirstChar+1 );
} // of TrimSpaces


// ================================================================================


void TrimLeadingSpaces ( string* Line )
// removes leading blanks and tabs
{
	// find the first character position after excluding leading blanks and tabs
	size_t FirstChar = Line->find_first_not_of (" \t");
	if( string::npos != FirstChar ) *Line = Line->substr ( FirstChar );
} // of TrimLeadingSpaces


// ================================================================================


void TrimTrailingSpaces ( string* Line )
// removes trailing blanks and tabs
{
	// find the first character position at the end of the string
	size_t LastChar = Line->find_last_not_of (" \t");
	if( string::npos != LastChar )
	*Line = Line->substr( 0, LastChar+1 );
} // of TrimTrailingSpaces


// ================================================================================


void SplitString ( const string& Line, vector<string>& Fields, const string& Delimiter = " ")
// splits a string into its single 'fields' or columns, based on a given delimiter
{
	// clear the vector before filling it
	Fields.clear();
	// Skip delimiters at beginning.
	string::size_type lastPos = Line.find_first_not_of (Delimiter, 0);
	// Find first 'non-delimiter'
	string::size_type pos     = Line.find_first_of (Delimiter, lastPos);
	
	while (string::npos != pos || string::npos != lastPos) {
		// found a token, add it to the vector
		Fields.push_back (Line.substr (lastPos, pos - lastPos));
		// skip delimiters
		lastPos = Line.find_first_not_of (Delimiter, pos);
		// find next 'non-delimiter'
		pos = Line.find_first_of (Delimiter, lastPos);
	}
	return;
} // of SplitString


// ================================================================================


string NextLine ( ifstream* File )
// reads the next line from a file stream and returns it as whitespace-trimmed string
{
	char   cLine[255];
	string Line = "";
	
	// while it's not EOF and it's not an empty or a commented line
	while ( (not File->eof() ) && ( (Line.compare ("") == 0 ) || (Line.compare (0, 1, "#") == 0) ) )
	{
		File->getline (cLine, 255);       // get the current line as C-string
		CString2String (cLine, &Line);    // convert the c-string to a real string
		TrimSpaces ( &Line );             // remove leading and trailing blanks
	}
	
	return Line;
} // of NextLine


// ================================================================================


void SplitNextLine ( ifstream* File, string& Line, vector<string>& Fields, const string& Delimiter )
// reads the next line from a file stream and returns it as whitespace-trimmed string
{
	Line = NextLine (File);
	SplitString (Line, Fields, Delimiter);
	return;
} // of SplitNextLine


// ================================================================================


bool caseInsCharCompare (char a, char b) { return (toupper(a) == toupper(b)); }

bool StringInsCompare (const string& String1, const string& String2)
// compares two strings case-insensitively
{
	
	return ( (String1.size() == String2.size()) &&
	       equal (String1.begin(), String1.end(), String2.begin(), caseInsCharCompare) );
} // of StringInsCompare


// ================================================================================


void Dichro::OutputInputConfigurationClass ( void )
{
	if (DC_Debug < 1) {
		cerr << "\nERROR: The debug option is required to use Dichro::OutputInputConfigurationClass.\n\n";
		return;
	}
	
	Dichro::OutputFileHeadline (DC_DbgFile,
	        "   $DC_Input.Configuration: Data read from the $CONFIGURATION section");
	
	fprintf (DC_DbgFile, "      %-14s %4d\n", "BBTrans:",  DC_Input.Configuration.BBTrans );
	fprintf (DC_DbgFile, "      %-14s %4d\n", "CTTrans:",  DC_Input.Configuration.CTTrans );
	fprintf (DC_DbgFile, "      %-14s %4d\n", "Factor:",   DC_Input.Configuration.Factor  );
	fprintf (DC_DbgFile, "      %-14s %4d\n", "MinWL:",    DC_Input.Configuration.MinWL   );
	fprintf (DC_DbgFile, "      %-14s %4d\n", "MaxWL:",    DC_Input.Configuration.MaxWL   );
	return;
} // of Dichro::OutputInputConfigurationClass


// ================================================================================


void Dichro::OutputInputParametersClass ( void )
{
	unsigned int i;
	
	if (DC_Debug < 1) {
		cerr << "\nERROR: The debug option is required to use Dichro::OutputInputParametersClass.\n\n";
		return;
	}
	
	Dichro::OutputFileHeadline (DC_DbgFile,
	        "   $DC_Input.Parameters: Data read from the $PARAMETERS section");
	
	fprintf (DC_DbgFile, "      Name       Transitions   Index\n");
	
	for (i = 0; i < DC_Input.Parameters.Name.size(); i++) {
		fprintf (DC_DbgFile,"      %8s       %2d          %2d\n",
		         DC_Input.Parameters.Name[i].c_str(), DC_Input.Parameters.Trans[i], i);
	}
	
	return;
} // of Dichro::OutputInputParametersClass


// ================================================================================


void Dichro::OutputInputChromophoresClass ( void )
{
	unsigned int Groups = DC_Input.Chromophores.Type.size();
	unsigned int i, j, Chrom;
	
	if (DC_Debug < 1) {
		cerr << "\nERROR: The debug option is required to use Dichro::OutputInputChromophoresClass.\n\n";
		return;
	}
	
	Dichro::OutputFileHeadline (DC_DbgFile,
	        "   $DC_Input.Parameters: Data read from the $CHROMOPHORES section");
	
	fprintf (DC_DbgFile, "    index  name         atom indeces in $COORDINATES section\n");
	
	for (i = 0; i < Groups; i++) {
		// save the index of the chromophore type
		Chrom = DC_Input.Chromophores.Type.at(i);
		
		// the current chromophore type as index and string (-1 because access starts counting at 0)
		fprintf (DC_DbgFile, "    %5d (%s):  ", Chrom, DC_Input.Parameters.Name.at(Chrom).c_str());
		
		// each atom number
		for (j = 0; j < DC_Input.Chromophores.Atoms.at(i).size(); j++)
			fprintf (DC_DbgFile, "%5d  ", DC_Input.Chromophores.Atoms.at(i).at(j) );
		
		fprintf (DC_DbgFile, "\n");
	}
	
	return;
} // of Dichro::OutputInputChromophoresClass


// ================================================================================


void Dichro::OutputInputCoordinatesClass ( void )
{
	unsigned int Groups = DC_Input.Coordinates.Groups.size();
	unsigned int Group, j;
	
	if (DC_Debug < 1) {
		cerr << "\nERROR: The debug option is required to use Dichro::OutputInputCoordinatesClass.\n\n";
		return;
	}
	
	Dichro::OutputFileHeadline (DC_DbgFile,
	        "   $DC_Input.Coordinates: Data read from the $COORDINATES section");
	
	fprintf (DC_DbgFile, "           x              y              z         atom   label\n");
	
	for (Group = 0; Group < Groups; Group++) {
		// print each coodinate
		for (j = 0; j < DC_Input.Coordinates.Groups.at(Group).size(); j++)
			fprintf (DC_DbgFile, "   %12f", DC_Input.Coordinates.Groups.at(Group).at(j));
			
			fprintf (DC_DbgFile, "     %4d     %3s\n",
			         DC_Input.Coordinates.Atoms.at(Group),
			         DC_Input.Coordinates.Labels.at(Group).c_str() );
	}
	
	return;
} // of Dichro::OutputInputCoordinatesClass


// ================================================================================


void Dichro::OutputInputClass ( void )
{
	
	if (DC_Debug < 1) {
		cerr << "\nERROR: The debug option is required to use Dichro::OutputInputCoordinatesClass.\n\n";
		return;
	}
	
	Dichro::OutputInputConfigurationClass ();
	Dichro::OutputInputParametersClass ();
	Dichro::OutputInputChromophoresClass ();
	Dichro::OutputInputCoordinatesClass ();
	
	return;
} // of Dichro::OutputInputClass

// ================================================================================


void Dichro::OutputParSetMonopoleClass ( ParSetMonopole* Monopole )
{
		fprintf (DC_DbgFile, "   %14.8f  %14.8f  %14.8f        %14.8f\n",
		   Monopole->Coord.at(0), Monopole->Coord.at(1), Monopole->Coord.at(2),
		   Monopole->Charge );
} // of Dichro::OutputParSetMonopoleClass


// ================================================================================


void Dichro::OutputParSetTransClass ( ParSetTrans* CurTrans, int Trans, int State, string Name )
{
	int Mono;
	
	if (DC_Debug < 1) {
		cerr << "\nERROR: The debug option is required to use Dichro::OutputParSetTransClass.\n\n";
		return;
	}
	
	// build the headline step by step
	string Header = "   ";
	if (Name != "#") Header = Header + Name + "   -   ";
	if (State > -1)  Header = Header + "State " + tostring(State) + "   -   ";
	
	if (CurTrans->Permanent) Header = Header + "Permanent Moment ";
	                    else Header = Header + "Transition ";
	
	if (Trans > -1) Header = Header + tostring(Trans);
	Dichro::OutputFileHeadline (DC_DbgFile, Header);
	
	fprintf (DC_DbgFile, "   Monopoles:       %3d\n",      CurTrans->NumberOfMonopoles);
	fprintf (DC_DbgFile, "   Energy:       %10.3f cm-1\n", CurTrans->Energy);
	fprintf (DC_DbgFile, "   Wavelength:     %8.3f nm\n",   CurTrans->Wavelength);
	fprintf (DC_DbgFile, "   Scale factor:   %8.3f\n",      CurTrans->ScaleFactor);
	
	fprintf (DC_DbgFile, "\n");
	fprintf (DC_DbgFile, "   Elec. trans. dip. mom.:   %12.8f   %12.8f   %12.8f\n",
	      CurTrans->EDM.at(0), CurTrans->EDM.at(1), CurTrans->EDM.at(2) );
	
	// check whether a magnetic dipole moment exists, permanent moments don't have any
	if (CurTrans->MDM.size() >= 1)
		fprintf (DC_DbgFile, "   Magn. trans. dip. mom.:   %12.8f   %12.8f   %12.8f\n",
		   CurTrans->MDM.at(0), CurTrans->MDM.at(1), CurTrans->MDM.at(2) );
	
	if (DC_Debug > 4) {
		fprintf (DC_DbgFile,
		         "\n            x               y               z                 Charge\n");
		
		for (Mono = 0; Mono < CurTrans->NumberOfMonopoles; Mono++)
			Dichro::OutputParSetMonopoleClass ( &CurTrans->Monopoles.at(Mono) );
		
		fprintf (DC_DbgFile, "\n");
	}
	
	return;
} // of Dichro::OutputParSetTransClass


// ================================================================================


void Dichro::OutputParSetClass ( ParSet* CurParSet )
{
	unsigned int State, Trans;
	int i;
	
	if (DC_Debug < 1) {
		cerr << "\nERROR: The debug option is required to use Dichro::OutputParSetClass.\n\n";
		return;
	}
	
	Dichro::OutputFileHeadline (DC_DbgFile, "$DC_ParSets: Parameter set");
	
	fprintf (DC_DbgFile, "   Name:   %s", CurParSet->Name.c_str() );
	if (CurParSet->ChargeTransfer)
		fprintf (DC_DbgFile, " (charge-transfer)\n");
	else
		fprintf (DC_DbgFile, "\n");
	
	fprintf (DC_DbgFile, "   Atoms:  %d\n", CurParSet->NumberOfAtoms);
	fprintf (DC_DbgFile, "\n");
	
	for (i = 0; i < CurParSet->NumberOfAtoms; i++) {
		fprintf (DC_DbgFile,
		         "   Atom %2d: %8.4f %8.4f %8.4f     Weighting: %6.2f     Label: %3s\n",
		         i, CurParSet->Atoms.at(i).Coord.at(0), CurParSet->Atoms.at(i).Coord.at(1),
		            CurParSet->Atoms.at(i).Coord.at(2), CurParSet->Atoms.at(i).Weighting,
		            CurParSet->Atoms.at(i).Label.c_str() );
	}
	
	fprintf (DC_DbgFile, "\n");
	
	fprintf (DC_DbgFile, "   Reference point: %8.4f %8.4f %8.4f\n",
	                       CurParSet->Reference.at(0),
	                       CurParSet->Reference.at(1),
	                       CurParSet->Reference.at(2));
	
	Dichro::OutputFileSeparator (DC_DbgFile, 3);
	
	for (State = 0; State < CurParSet->States.size(); State++) {
		
		for (Trans = 0; Trans < CurParSet->States.at(State).size(); Trans++) {
			Dichro::OutputParSetTransClass (&CurParSet->States.at(State).at(Trans),
			                                Trans, State, CurParSet->Name);
		}
		
		// don't add the separator after the last transition
		if (State < CurParSet->States.size() -1)
			Dichro::OutputFileSeparator (DC_DbgFile, 3);
	}
	
	return;
} // of Dichro::OutputParSetClass


// ================================================================================


void Dichro::OutputSystemTransitionClass ( SystemTransition* CurTrans, int Group, int Trans )
{
	unsigned int Mono;
	
	if (DC_Debug < 1) {
		cerr << "\nERROR: The debug option is required to use Dichro::OutputSystemTransitionClass.\n\n";
		return;
	}
	
	if (not CurTrans->Permanent)
		fprintf (DC_DbgFile, "\n   Transition %d  (Group %d)\n\n", Trans, Group);
	else
		fprintf (DC_DbgFile, "\n   Permanent Moment %d  (Group %d)\n\n", Trans, Group);
	
	fprintf (DC_DbgFile, "      Parameter set:   %s\n", CurTrans->Origin.c_str());
	fprintf (DC_DbgFile, "      Energy:          %8.2f cm-1\n", CurTrans->Energy);
	fprintf (DC_DbgFile, "      Wavelength:      %8.2f nm\n", CurTrans->Wavelength);
	fprintf (DC_DbgFile, "      Monopoles:       %8d\n", CurTrans->NumberOfMonopoles);
	fprintf (DC_DbgFile, "      EDM:         %12.6f   %12.6f   %12.6f\n",
	         CurTrans->EDM.at(0), CurTrans->EDM.at(1), CurTrans->EDM.at(2));
	
	if (not CurTrans->Permanent)
		fprintf (DC_DbgFile, "      MDM:         %12.6f   %12.6f   %12.6f\n",
		      CurTrans->MDM.at(0), CurTrans->MDM.at(1), CurTrans->MDM.at(2));
	
	if (DC_Debug > 4) {
		fprintf (DC_DbgFile,
		         "\n   Monopoles:        x              y              z              q\n");
		
		for (Mono = 0; Mono < CurTrans->Monopoles.size(); Mono++)
			fprintf (DC_DbgFile, "             %12.6f   %12.6f   %12.6f   %12.6f\n",
			         CurTrans->Monopoles.at(Mono).Coord.at(0),
			         CurTrans->Monopoles.at(Mono).Coord.at(1),
			         CurTrans->Monopoles.at(Mono).Coord.at(2),
			         CurTrans->Monopoles.at(Mono).Charge);
	}
	
	return;
} // of Dichro::OutputSystemTransitionClass


// ================================================================================


void Dichro::OutputSystemGroupClass ( SystemGroup* CurGroup, int Group )
{
	unsigned int Atom, Trans, Perm, NumberOfTransitions;
	SystemTransition* CurTrans;
	SystemTransition* CurPerm;
	
	if (DC_Debug < 1) {
		cerr << "\nERROR: The debug option is required to use Dichro::OutputSystemGroupClass.\n\n";
		return;
	}
	
	Dichro::OutputFileSeparator (DC_DbgFile);
	
	if (Group > -1) { // -1 is the default value
		if (CurGroup->ChargeTransfer)
			Dichro::OutputFileHeadline (DC_DbgFile,
			        "$DC_System.Groups.[" + tostring(Group) + "]: Chromophore " + tostring(Group) + " (Charge-transfer)");
		else
			Dichro::OutputFileHeadline (DC_DbgFile,
			        "$DC_System.Groups.[" + tostring(Group) + "]: Chromophore " + tostring(Group));
	}
	else {
		if (CurGroup->ChargeTransfer)
			Dichro::OutputFileHeadline (DC_DbgFile, "$DC_System.Groups: Chromophore (Charge-transfer)");
		else
			Dichro::OutputFileHeadline (DC_DbgFile, "$DC_System.Groups: Chromophore");
	}
	
	fprintf (DC_DbgFile, "   Parameter set: %s\n\n", CurGroup->ParameterSet.c_str());
	
	// print the atoms
	for (Atom = 0; Atom < CurGroup->Atoms.size(); Atom++)
		fprintf (DC_DbgFile, "      Atom %2d (Index %2d): %12.4f %12.4f %12.4f\n",
		         Atom, CurGroup->AtomIndices.at(Atom),
		         CurGroup->Atoms.at(Atom).at(0),
		         CurGroup->Atoms.at(Atom).at(1),
		         CurGroup->Atoms.at(Atom).at(2));
	
	fprintf (DC_DbgFile, "\n      Reference point:    %12.4f %12.4f %12.4f\n",
	         CurGroup->Reference.at(0), CurGroup->Reference.at(1), CurGroup->Reference.at(2));
	
	NumberOfTransitions = CurGroup->NumberOfTransitions;
	
	// print the transitions
	for (Trans = 0; Trans < CurGroup->Trans.size(); Trans++) {
		CurTrans = &CurGroup->Trans.at(Trans);
		Dichro::OutputSystemTransitionClass (CurTrans, Group, Trans);
		Dichro::OutputFileSeparator (DC_DbgFile, 3);
	} // of for (Trans = 0; Trans < CurGroup.NumberOfTransitions; Trans++)
	
	// print the permanent moment for each ground state transition
	for (Perm = 0; Perm < CurGroup->Perm.size(); Perm++) {
		CurPerm = &CurGroup->Perm.at(Perm);
		Dichro::OutputSystemTransitionClass (CurPerm, Group, Perm);
		
		if (Perm < CurGroup->Perm.size()-1)
		Dichro::OutputFileSeparator (DC_DbgFile, 3);
	} // of for (Perm = 0; Perm < CurGroup.NumberOfTransitions; Perm++)
	
	return;
} // of Dichro::OutputSystemGroupClass


// ================================================================================


void Dichro::OutputSystemClass ( void )
// prints out the DC_System class (there is only one, hence the void)
{
	unsigned int Group, Atom;
	
	if (DC_Debug < 1) {
		cerr << "\nERROR: The debug option is required to use Dichro::OutputSystemClass.\n\n";
		return;
	}
	
	Dichro::NewFileTask (DC_DbgFile, "$DC_System: System Data");
	
	fprintf (DC_DbgFile, "Number of atoms:          %5d\n", DC_System.NumberOfAtoms);
	fprintf (DC_DbgFile, "Number of groups:         %5d\n", DC_System.NumberOfGroups);
	fprintf (DC_DbgFile, "Number of transitions:    %5d\n", DC_System.NumberOfTransitions);
	fprintf (DC_DbgFile, "Matrix dimension:         %5d\n", DC_System.MatrixDimension);
	
	fprintf (DC_DbgFile, "\nAtoms representing the system:\n\n");
	for (Atom = 0; Atom < DC_System.Atoms.size(); Atom++) {
		fprintf (DC_DbgFile, "Atom %2d: %12.4f %12.4f %12.4f\n", Atom,
		         DC_System.Atoms.at(Atom).at(0),
		         DC_System.Atoms.at(Atom).at(1),
		         DC_System.Atoms.at(Atom).at(2));
	}
	
	for (Group = 0; Group < DC_System.Groups.size(); Group++) {
		Dichro::OutputSystemGroupClass (&DC_System.Groups.at(Group), Group);
	}
	
	fprintf (DC_DbgFile, "\n\n");
	
	return;
} // of Dichro::OutputSystemClass


// ================================================================================


void Dichro::OutputResultsTransClass ( void )
{
	int Trans;
	
	if (DC_Debug < 1) {
		cerr << "\nERROR: The debug option is required to use Dichro::OutputResultsTransClass.\n\n";
		return;
	}
	
	Dichro::OutputFileSeparator (DC_DbgFile, 3);
	fprintf (DC_DbgFile,
	         "   $DC_Results.Trans: Data Sorted According to the Transition Sequence");
	Dichro::OutputFileSeparator (DC_DbgFile, 3);
	
	Dichro::OutputFileHeadline (DC_DbgFile, "   Dipole Moments");
	
	fprintf (DC_DbgFile, "                                          electric");
	fprintf (DC_DbgFile, "                                         magnetic\n");
	fprintf (DC_DbgFile, "   Index  group trans          x              y              z");
	fprintf (DC_DbgFile, "                  x              y              z\n");
	
	for (Trans = 0; Trans < DC_Results.NumberOfTransitions; Trans++) {
		fprintf (DC_DbgFile,
		     "   %5d   %3d , %-3d   %12.6f   %12.6f   %12.6f       %12.6f   %12.6f   %12.6f\n",
			  Trans,
			  DC_Results.Trans.GroupSequence.at(Trans),
			  DC_Results.Trans.TransSequence.at(Trans),
			
			  DC_Results.Trans.EDM.at(Trans).at(0),
			  DC_Results.Trans.EDM.at(Trans).at(1),
			  DC_Results.Trans.EDM.at(Trans).at(2),
			
			  DC_Results.Trans.MDM.at(Trans).at(0),
			  DC_Results.Trans.MDM.at(Trans).at(1),
			  DC_Results.Trans.MDM.at(Trans).at(2));
	}
	

	Dichro::OutputFileHeadline (DC_DbgFile, "   Energies, Wavelengths & Reference Vectors");
	
	fprintf (DC_DbgFile, "   Index  group trans    par. set      energy     wavelength");
	fprintf (DC_DbgFile, "   dip.str.    rot.str.  ");
	fprintf (DC_DbgFile, "        x              y              z\n");
	
	for (Trans = 0; Trans < DC_Results.NumberOfTransitions; Trans++) {
		fprintf (DC_DbgFile,
		     "   %5d   %3d , %-3d   %10s    %9.3f   %9.3f   %9.3f   %9.3f   %12.6f   %12.6f   %12.6f\n",
		     Trans,
			  DC_Results.Trans.GroupSequence.at(Trans),
			  DC_Results.Trans.TransSequence.at(Trans),
			  DC_Results.Trans.ParSetSequence.at(Trans).c_str(),
			  
		     DC_Results.Trans.Energy.at(Trans), DC_Results.Trans.Wavelength.at(Trans),
			  DC_Results.Trans.DipoleStrength.at(Trans),
		     DC_Results.Trans.RotationalStrength.at(Trans),
			  DC_Results.Trans.Reference.at(Trans).at(0),
			  DC_Results.Trans.Reference.at(Trans).at(1),
		     DC_Results.Trans.Reference.at(Trans).at(2));
	}
	
	Dichro::OutputFileHeadline (DC_DbgFile, "   Polarization Vectors");
	
	fprintf (DC_DbgFile, "   Index  group trans          x              y              z");
	fprintf (DC_DbgFile, "            Oscill. str.\n");
	
	for (Trans = 0; Trans < DC_Results.NumberOfTransitions; Trans++) {
		fprintf (DC_DbgFile,
		     "   %5d   %3d , %-3d   %12.6f   %12.6f   %12.6f       %12.6f\n",
			  Trans,
		     DC_Results.Trans.GroupSequence.at(Trans),
			  DC_Results.Trans.TransSequence.at(Trans),
			  
		     DC_Results.Trans.PolarizationVector.at(Trans).at(0),
			  DC_Results.Trans.PolarizationVector.at(Trans).at(1),
		     DC_Results.Trans.PolarizationVector.at(Trans).at(2),
		     
			  DC_Results.Trans.OscillatorStrength.at(Trans));
	}
} // of Dichro::OutputResultsTransClass ( void )


// ================================================================================


void Dichro::OutputResultsGroupClass ( void )
{
	unsigned int Group;
	int Trans;
	
	if (DC_Debug < 1) {
		cerr << "\nERROR: The debug option is required to use Dichro::OutputResultsGroupClass.\n\n";
		return;
	}
	
	Dichro::OutputFileSeparator (DC_DbgFile, 3);
	fprintf (DC_DbgFile, "   $DC_Results.Groups: Data Sorted by Groups");
	Dichro::OutputFileSeparator (DC_DbgFile, 3);
	
	for (Group = 0; Group < DC_Results.Groups.size(); Group++) {
		ResultsGroup* CurGroup = &DC_Results.Groups.at(Group);
		
		if (not CurGroup->ChargeTransfer)
			Dichro::OutputFileHeadline (DC_DbgFile, "   Group " + tostring(Group));
		else
			Dichro::OutputFileHeadline (DC_DbgFile,
			                            "   Group " + tostring(Group) + "(Charge-transfer)");
		
		fprintf (DC_DbgFile, "      Parameter set:        %s\n", CurGroup->ParameterSet.c_str());
		fprintf (DC_DbgFile, "      Transitions:          %d\n", CurGroup->NumberOfTransitions);
		fprintf (DC_DbgFile, "      Reference vector: %12.6f %12.6f %12.6f\n\n",
		      CurGroup->Reference.at(0), CurGroup->Reference.at(1), CurGroup->Reference.at(2));
		
		fprintf (DC_DbgFile, "                                   electric");
		fprintf (DC_DbgFile, "                                         magnetic\n");
		fprintf (DC_DbgFile, "   group trans          x              y              z");
		fprintf (DC_DbgFile, "                  x              y              z\n");
		
		for (Trans = 0; Trans < CurGroup->NumberOfTransitions; Trans++) {
			fprintf (DC_DbgFile,
			  "    %3d , %-3d   %12.6f   %12.6f   %12.6f       %12.6f   %12.6f   %12.6f\n",
			  Group, Trans,
			
			  CurGroup->EDM.at(Trans).at(0),
			  CurGroup->EDM.at(Trans).at(1),
			  CurGroup->EDM.at(Trans).at(2),
			
			  CurGroup->MDM.at(Trans).at(0),
			  CurGroup->MDM.at(Trans).at(1),
		     CurGroup->MDM.at(Trans).at(2));
		}
		
		fprintf (DC_DbgFile, "\n   group trans    energy     wavelength");
		fprintf (DC_DbgFile, "    dip.str.    rot.str.\n");
		
		for (Trans = 0; Trans < CurGroup->NumberOfTransitions; Trans++) {
			fprintf (DC_DbgFile,
			  "    %3d , %-3d    %9.3f   %9.3f   %9.3f   %9.3f\n",
			  Group, Trans,
			  CurGroup->Energy.at(Trans), DC_Results.Trans.Wavelength.at(Trans),
			  CurGroup->DipoleStrength.at(Trans), CurGroup->RotationalStrength.at(Trans));
		}
		
		fprintf (DC_DbgFile,"\n                                           polarization vectors\n");
		fprintf (DC_DbgFile,"   group trans    wavelength          x              y              z");
		fprintf (DC_DbgFile,"            Oscill. str.\n");
		
		for (Trans = 0; Trans < CurGroup->NumberOfTransitions; Trans++)
			fprintf (DC_DbgFile,
				"    %3d , %-3d  %12.3f   %12.6f   %12.6f   %12.6f       %12.6f\n",
				Group, Trans,
				DC_Results.Trans.Wavelength.at(Trans),
				
				CurGroup->PolarizationVector.at(Trans).at(0),
				CurGroup->PolarizationVector.at(Trans).at(1),
				CurGroup->PolarizationVector.at(Trans).at(2),
				
				CurGroup->OscillatorStrength.at(Trans));
		
		fprintf (DC_DbgFile, "\n");
		
		fprintf (DC_DbgFile, "   Submatrix of the group\n\n");
		FilePrintMatrix (DC_DbgFile, &CurGroup->Submatrix);
		fprintf (DC_DbgFile, "\n");
		fprintf (DC_DbgFile, "\n");
	}
	
	return;
} // of Dichro::OutputResultsGroupClass ( void )


// ================================================================================


void Dichro::OutputResultsClass ( void )
{
	if (DC_Debug < 1) {
		cerr << "\nERROR: The debug option is required to use Dichro::OutputResultsClass.\n\n";
		return;
	}
	
	Dichro::NewFileTask (DC_DbgFile, "$DC_Results: Calculated Results");
	
	Dichro::OutputFileHeadline (DC_DbgFile, "   Statistics");
	
	fprintf (DC_DbgFile, "   Number of atoms:          %5d\n", DC_System.NumberOfAtoms);
	fprintf (DC_DbgFile, "   Number of groups:         %5d\n", DC_System.NumberOfGroups);
	fprintf (DC_DbgFile, "   Number of transitions:    %5d\n", DC_System.NumberOfTransitions);
	fprintf (DC_DbgFile, "   Matrix dimension:         %5d\n", DC_System.MatrixDimension);
	
	Dichro::OutputResultsTransClass();  // access by transitions
	Dichro::OutputResultsGroupClass();  // access by groups
	
	if (DC_Debug > 4) {
		// the Hamiltonian matrix
		Dichro::OutputFileHeadline (DC_DbgFile, "$DC_Results.Hamiltonian: Hamiltonian Matrix");
		FilePrintMatrix (DC_DbgFile, &DC_Results.Hamiltonian);
		
		// the eigenvectors
		Dichro::OutputFileHeadline (DC_DbgFile, "$DC_Results.Eigenvectors: Eigenvectors");
		FilePrintMatrix (DC_DbgFile, &DC_Results.Eigenvectors);
		
		// the eigenvalues
		Dichro::OutputFileHeadline (DC_DbgFile, "$DC_Results.Eigenvalues: Eigenvalues");
		FilePrintMatrix (DC_DbgFile, &DC_Results.Eigenvalues, false );
	}
	
	fprintf (DC_DbgFile, "\n\n");
	
	return;
} // of Dichro::OutputResultsClass


// ================================================================================


void Dichro::OutputSystemData ( void )
// prints out all kinds of data to different files for plotting
{
	int Atom;
	string FileName;
	FILE* AtomsFile;
	
	if (DC_Debug < 1) {
		cerr << "\nERROR: The debug option is required to use Dichro::OutputSystemData.\n\n";
		return;
	}
	
	FileName = DC_InFileBaseName + ".atoms";
	
	AtomsFile = fopen (FileName.c_str(), "w");
	for (Atom = 0; Atom < DC_System.NumberOfAtoms; Atom++) {
		fprintf (AtomsFile, "%8.4f %8.4f %8.4f\n", DC_System.Atoms.at(Atom).at(0),
		           DC_System.Atoms.at(Atom).at(1), DC_System.Atoms.at(Atom).at(2));
	}
	fclose (AtomsFile);
	
	return;
} // of Dichro::OutputSystemData


// ================================================================================


void PrintCoord ( vector<double> *Vector, bool Norm )
// prints cartesian coordinates
{
	for ( unsigned int i = 0; i < Vector->size(); i++ )
		printf ("%12.6f", Vector->at(i));

	if (Norm) printf ("     Norm:  %12.6f", VectorNorm (Vector));
	printf ("\n");
	return;
} // of PrintCoord


void FilePrintCoord ( FILE* File, vector<double> *Vector, bool Norm )
// prints cartesian coordinates
{
	for ( unsigned int i = 0; i < Vector->size(); i++ )
		fprintf (File, "%12.6f", Vector->at(i));
	
	if (Norm) fprintf (File, "     Norm:  %12.6f", VectorNorm (Vector));
	fprintf (File, "\n");
	return;
} // of FilePrintCoord


// ================================================================================


void PrintVector ( vector<int>* Vector )
// prints all elements of a given vector
{
	cout << "\n-> Vector of int\n\n";
	
	for ( unsigned int i = 0; i < Vector->size(); i++ )
		cout << "  Vector[" << i << "]: " << Vector->at(i) << endl;
	
	cout << "\n------------------------------------\n";
	return;
} // of PrintVector


void PrintVector ( vector<double>* Vector )
// prints all elements of a given vector
{
	cout << "\n-> Vector of double\n\n";
	
	for ( unsigned int i = 0; i < Vector->size(); i++ )
		cout << "  Vector[" << i << "]: " << Vector->at(i) << endl;
	
	cout << "\n------------------------------------\n\n";
	return;
} // of PrintVector


void PrintVector ( vector<string>* Vector )
// prints all elements of a given vector
{
	cout << "\n-> Vector of string\n\n";
	
	for ( unsigned int i = 0; i < Vector->size(); i++ )
		cout << "  Vector[" << i << "]: " << Vector->at(i) << endl;
	
	cout << "\n------------------------------------\n\n";
	return;
} // of PrintVector


void PrintVector ( vector< vector<int> >* Vector )
// prints all elements of a given vector
{
	cout << "\n-> Vector of vector of int\n\n";
	
	for ( unsigned int i = 0; i < Vector->size(); i++ ) {
		cout << endl;
		
		for ( unsigned int k = 0; k < Vector->at(i).size(); k++ )
			cout << "  Vector[" << i << "][" << k << "] : "
			     << Vector->at(i).at(k) << endl;
	}
	
	cout << "\n------------------------------------\n\n";
	return;
} // of PrintVector


void PrintVector ( vector< vector<double> >* Vector )
// prints all elements of a given vector
{
	cout << "\n-> Vector of vector of double\n\n";
	
	for ( unsigned int i = 0; i < Vector->size(); i++ ) {
		cout << endl;
		
		for ( unsigned int k = 0; k < Vector->at(i).size(); k++ )
			cout << "  Vector[" << i << "][" << k << "] : "
			     << Vector->at(i).at(k) << endl;
	}
	
	cout << "\n------------------------------------\n\n";
	return;
} // of PrintVector


void PrintVector ( vector< vector<string> >* Vector )
// prints all elements of a given vector
{
	cout << "\n-> Vector of vector of string\n\n";
	
	for ( unsigned int i = 0; i < Vector->size(); i++ ) {
		cout << endl;
		
		for ( unsigned int k = 0; k < Vector->at(i).size(); k++ )
			cout << "  Vector[" << i << "][" << k << "] : "
			     << Vector->at(i).at(k) << endl;
	}
	
	cout << "\n------------------------------------\n\n";
	return;
} // of PrintVector


void Dichro::DebugOutput ( vector< vector<double> >* Vector )
// prints all elements of a given vector
{
	if (DC_Debug < 1) {
		cerr << "\nERROR: The debug option is required to use Dichro::DebugOutput.\n\n";
		return;
	}
	
	for ( unsigned int i = 0; i < Vector->size(); i++ ) {
		fprintf (DC_DbgFile, "\n");
		
		for ( unsigned int k = 0; k < Vector->at(i).size(); k++ )
			fprintf (DC_DbgFile, "     %f", Vector->at(i).at(k));
	}
	
	return;
} // of Dichro::DebugOutput


// ================================================================================


void PrintMatrix ( Matrix* InMatrix )
{
	int row, col;
	
	for (row = 0; row < InMatrix->nrows(); row++) {
		for (col = 0; col < InMatrix->ncols(); col++)
			printf ("%12.4f", InMatrix->element(row, col));
		
		printf ("\n");
	}
} // of PrintMatrix


void PrintMatrix ( SymmetricMatrix* InMatrix )
{
	int row, col;
	
	for (row = 0; row < InMatrix->nrows(); row++) {
		for (col = 0; col < InMatrix->ncols(); col++)
			printf ("%12.4f", InMatrix->element(row, col));
		
		printf ("\n");
	}
} // of PrintMatrix


void PrintMatrix ( DiagonalMatrix* InMatrix, bool Indent )
{
	int row, col;
	
	for (row = 0; row < InMatrix->nrows(); row++) {
		if (Indent) {
			// indent the diagonal element about as many columns as rows have been printed
			for (col = 0; col < row; col++) printf ("                 ");
		}
		
		// print the actual element
		printf ("%17f\n", InMatrix->element(row));
	}
} // of PrintMatrix


void FilePrintMatrix ( FILE* File, Matrix* InMatrix )
{
	int row, col;
	
	for (row = 0; row < InMatrix->nrows(); row++) {
		for (col = 0; col < InMatrix->ncols(); col++)
			fprintf (File, "%17f", InMatrix->element(row, col));
		
		fprintf (File, "\n");
	}
} // of FilePrintMatrix


void FilePrintMatrix ( FILE* File, SymmetricMatrix* InMatrix )
{
	int row, col;
	
	for (row = 0; row < InMatrix->nrows(); row++) {
		for (col = 0; col < InMatrix->ncols(); col++)
			fprintf (File, "%17f", InMatrix->element(row, col));
		
		fprintf (File, "\n");
	}
} // of FilePrintMatrix


void FilePrintMatrix ( FILE* File, DiagonalMatrix* InMatrix, bool Indent )
{
	int col, i;
	
	for (col = 0; col < InMatrix->ncols(); col++) {
		if (Indent) {
			// indent the diagonal element about as many columns as rows have been printed
			for (i = 0; i < col; i++) fprintf (File, "                 ");
		}
		
		// print the actual element
		fprintf (File, "%17f\n", InMatrix->element(col));
	}
} // of FilePrintMatrix


// ================================================================================


void dp ( string String )
{
	cout << "\n\nHappily died a controlled death in DebugPrint:\n\n" << String << " (string)\n\n";
	exit (1);
} // of dp

void dp ( int Integer )
{
	cout << "\n\nHappily died a controlled death in DebugPrint:\n\n" << Integer << " (int)\n\n";
	exit (1);
} // of dp

void dp ( double Double )
{
	cout << "\n\nHappily died a controlled death in DebugPrint:\n\n" << Double << " (double)\n\n";
	exit (1);
} // of dp

void dp ( size_t Size )
{
	cout << "\n\nHappily died a controlled death in DebugPrint:\n\n" << Size << " (size_t)\n\n";
	exit (1);
} // of dp

void dp ( vector<int> Vector )
{
	cout << "\n\nHappily died a controlled death in DebugPrint:\n\n";
	PrintVector ( &Vector );
	cout << "\n\n";
	exit (1);
} // of dp

void dp ( vector<double> Vector )
{
	cout << "\n\nHappily died a controlled death in DebugPrint:\n\n";
	PrintVector ( &Vector );
	cout << "\n\n";
	exit (1);
} // of dp

void dp ( vector<string> Vector )
{
	cout << "\n\nHappily died a controlled death in DebugPrint:\n\n";
	PrintVector ( &Vector );
	cout << "\n\n";
	exit (1);
} // of dp

void dp ( vector< vector<int> > Vector )
{
	cout << "\n\nHappily died a controlled death in DebugPrint:\n\n";
	PrintVector ( &Vector );
	cout << "\n\n";
	exit (1);
} // of dp

void dp ( vector< vector<double> > Vector )
{
	cout << "\n\nHappily died a controlled death in DebugPrint:\n\n";
	PrintVector ( &Vector );
	cout << "\n\n";
	exit (1);
} // of dp

void dp ( vector< vector<string> > Vector )
{
	cout << "\n\nHappily died a controlled death in DebugPrint:\n\n";
	PrintVector ( &Vector );
	cout << "\n\n";
	exit (1);
} // of dp


// ================================================================================


