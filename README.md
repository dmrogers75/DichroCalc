# DichroCalc
Repository for the C++ version of DichroCalc.

DichroCalc C++ Version 1.0

Dependencies:          C++ compiler (to compile the DichroCalc .cpp source code). Linux OS or Mac OS (for the third-party library NewMat). Perl (to execute the DichroCalc perl scripts).


Description of DichroCalc files and directories present in git repo.


chromophores.dat       a list of parameter sets. The file can be placed in ~/bin and is required by, e.g., the perl script "./scripts/dcinput".


dichrocalc/            a directory containing sub-directories and files: 
dichrocalc/include/    DichroCalc header files;
dichrocalc/lib/        third-party header files and libraries (NewMat);
dichrocalc/makefile    the makefile;
dichrocalc/obj/        empty;
dichrocalc/src/        DichroCalc .cpp source code files.


documentation/         the DichroCalc documentation.


LICENSE                the GNU General Public License v3.0.


params/                a directory containing the parameter sets. The directory can be placed in ~/bin.


perllib/               a directory containing libraries required for the perl scripts (in ./scripts/).


REFERENCES	           a file listing the references you should cite if you present any data in a publication.


scripts/               a directory containing various perl scripts.


diabatisation/         a directory containing code for diabatisation of the diamide chromophore.
