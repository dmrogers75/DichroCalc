// #################################################################################################
//
//  Header:       iolibrary.h
//
//  Author:       Benjamin M. Bulheller
//
//  Version:      $Revision: 4602 $, $Date: 2009-06-27 05:16:06 -0700 (Sat, 27 Jun 2009) $
//
// #################################################################################################

bool   FileExists    ( string FileName );
bool   FileExtension ( string Filename, string Extension );
bool   ReadDir ( string Dir, string Extension, vector<string>* Files );

void   VectorDiff    (vector<double>* Vect1, vector<double>* Vect2, vector<double>* Diff);
double VectorNorm (vector<double>* Vector);
void   CrossProduct  (vector<double>* Vect1, vector<double>* Vect2, vector<double>* Cross);
double PointDistance (vector<double>* Vect1, vector<double>* Vect2);

string tostring ( unsigned int Integer );
string tostring ( int Integer );
string tostring ( double Float );

void   TrimSpaces ( string* Line );
void   TrimLeadingSpaces  ( string* Line );
void   TrimTrailingSpaces ( string* Line );

void   SplitString ( const string &Line, vector<string> &Fields, const string &Delimiter );
void   CString2String ( char* cLine, string* Line );

string NextLine ( ifstream* File );
void   SplitNextLine ( ifstream* File, string& Line, vector<string>& Fields,
                       const string& Delimiter );

bool   StringInsCompare (const string& String1, const string& String2);

void   PrintCoord  ( vector<double>* Vector, bool Norm = false );
void   FilePrintCoord ( FILE* File, vector<double> *Vector, bool Norm = false );
void   PrintVector ( vector<string>* Vector  );
void   PrintVector ( vector<int>*    Vector  );
void   PrintVector ( vector<double>* Vector );
void   PrintVector ( vector< vector<int> >* Vector );
void   PrintVector ( vector< vector<double> >* Vector );
void   PrintVector ( vector< vector<string> >* Vector );

void   PrintMatrix ( Matrix* InMatrix );
void   PrintMatrix ( SymmetricMatrix* InMatrix );
void   PrintMatrix ( DiagonalMatrix* InMatrix, bool Indent = true );
void   FilePrintMatrix ( FILE* File, Matrix* InMatrix );
void   FilePrintMatrix ( FILE* File, SymmetricMatrix* InMatrix );
void   FilePrintMatrix ( FILE* File, DiagonalMatrix* InMatrix, bool Indent = true );

void   dp ( string String  );
void   dp ( int    Integer );
void   dp ( double Double  );
void   dp ( size_t Size    );

void   dp ( vector<int>    Vector );
void   dp ( vector<double> Vector );
void   dp ( vector<string> Vector );
void   dp ( vector< vector<int> > Vector );
void   dp ( vector< vector<double> > Vector );
void   dp ( vector< vector<string> > Vector );

