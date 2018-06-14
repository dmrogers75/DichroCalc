// #################################################################################################
//
//  Program:      dichrocalc
//
//  Function:     Front-end for DichroCalc routines to perform matrix method calculations
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

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string.h>
#include <vector>
#include <getopt.h>

using namespace std;

#include "../include/dichrocalc.h"

// ================================================================================


void PrintArguments ( void );
int  Usage ( void );
int  ProcessCommandLineOptions ( int argc, char **argv );

class CommandLineArguments {
	public:
		int    Debug;
		bool   Verbose;
		bool   PrintVec;
		bool   PrintPol;
		bool   PrintMat;
		string InFile;
		string Params;
};

CommandLineArguments GlobalArgs;


// ================================================================================


void PrintArguments ( void )
// mainly for debugging purposes
{
	cout << "InFile   = " << GlobalArgs.InFile   << endl
	     << "Params   = " << GlobalArgs.Params   << endl
	     << "Verbose  = " << GlobalArgs.Verbose  << endl
	     << "PrintVec = " << GlobalArgs.PrintVec << endl
	     << "PrintPol = " << GlobalArgs.PrintPol << endl
	     << "PrintMat = " << GlobalArgs.PrintMat << endl
	     << "Debug    = " << GlobalArgs.Debug    << endl
	     << "\n\n";
	return;
} // of PrintArguments


// ================================================================================


int Usage ( void )
{
	cout << "\n";
	cout << "Usage: dichrocalc [options]\n\n";
	cout << "       -i , --input inputfile  filename of the input file to process (mandatory)\n";
	cout << "       -p , --params           directory with the parameter files (*.par)\n";
	cout << "       -v , --verbose          verbose output\n";
	cout << "       -d , --debug            set level of debug output (0-5)\n";
	cout << "            --vec              create .vec file (for absorbance/LD)\n";
	cout << "            --pol              create .pol file (transition polarizations)\n";
	cout << "            --mat              create .mat file (matrix, eigenvectors, eigenvalues)\n";
	cout << "       -h , --help, -?         usage output\n";
	cout << "\n";
	return 0;
} // of Usage


// ================================================================================


int ProcessCommandLineOptions ( int argc, char **argv )
{
	int NextOption;
	vector<string> FileNames;
	
	const char *const ShortOptions = "h?vd:i:p:";
	const struct option LongOptions[] = {
		// if NULL is given in the 3rd column, the value in the 4th column is returned if the
		// long option is found
		{ "input",   required_argument, NULL, 'i' },
		{ "params",  required_argument, NULL, 'p' },
		{ "help",    no_argument,       NULL, 'h' },
		{ "debug",   required_argument, NULL, 'd' },
		{ "verbose", no_argument,       NULL, 'v' },
		{ "vec",     no_argument,       NULL,  1  },
		{ "pol",     no_argument,       NULL,  2  },
		{ "mat",     no_argument,       NULL,  3  },
		{ NULL,      no_argument,       NULL,  0  },
	};
	
	do {
		NextOption = getopt_long (argc, argv, ShortOptions, LongOptions, NULL);
		
		switch (NextOption) {
			case 1:
				GlobalArgs.PrintVec = true;
				break;
			case 2:
				GlobalArgs.PrintPol = true;
				break;
			case 3:
				GlobalArgs.PrintMat = true;
				break;
			case 'i':
				GlobalArgs.InFile = string (optarg);
				
				if ( ! FileExists ( GlobalArgs.InFile ) ) {
					// try adding the .inp extension
					GlobalArgs.InFile += ".inp";
					
					if ( ! FileExists ( GlobalArgs.InFile ) ) {
						cerr << "\nERROR: File " << GlobalArgs.InFile << " not found.\n\n";
						return 5;
					}
				}
				
				break;
			case 'p':
				GlobalArgs.Params = string (optarg);
				
				ReadDir ( GlobalArgs.Params.c_str(), ".par", &FileNames );
				break;
			case 'd':
				GlobalArgs.Debug = atoi (optarg); // convert string to integer
				break;
			case 'v':
				GlobalArgs.Verbose = true;
				break;
			case 'h':
			case '?':
				Usage ();
				return 10;
				break;
		}
	} while (NextOption != -1);
	
	// check for mandatory parameters
	if (GlobalArgs.InFile == "")  {
		cerr << "\nERROR: No input file given via -i or --input.\n\n";
		return 15;
	}
	
	return 0;
} // of ProcessCommandLineOptions


// ================================================================================


int main ( int argc, char **argv )
{
	// if no parameter is given, display usage and quit
	if (argc == 1) {
		Usage ();
		return 1;
	}
	
	// initialize the variables
	GlobalArgs.InFile   = "";
	GlobalArgs.Params   = "";
	GlobalArgs.Debug    = 0;
	GlobalArgs.Verbose  = false;
	GlobalArgs.PrintVec = false;
	GlobalArgs.PrintPol = false;
	GlobalArgs.PrintMat = false;

	if (ProcessCommandLineOptions (argc, argv) != 0) { return 25; }
	
	// PrintArguments ();
	
	if (GlobalArgs.Verbose) {
		cout << "\nMatrix Method Calculations\n";
		cout <<   "==========================\n";
	}
	
	// -----------------------------------------------------------------------------
	
	// the object variable
	Dichro *DichroCalc;
	// construct the object with the required parameters
	DichroCalc = new Dichro ( GlobalArgs.InFile,  GlobalArgs.Params,    GlobalArgs.Verbose,
	                          GlobalArgs.Debug,   GlobalArgs.PrintVec,  GlobalArgs.PrintPol,
	                          GlobalArgs.PrintMat );
	
	if (GlobalArgs.Verbose) cout << "\n\n";
		
	cout << "Hamiltonian\n\n"  << setw(15) << DichroCalc->DC_Results.Hamiltonian  << "\n\n";
	cout << "Eigenvectors\n\n" << setw(15) << DichroCalc->DC_Results.Eigenvectors << "\n\n";
	cout << "Eigenvalues\n\n"  << setw(15) << DichroCalc->DC_Results.Eigenvalues  << "\n\n";
	
	return 0;
} // of main


// ================================================================================

