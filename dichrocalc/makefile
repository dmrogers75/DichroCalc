####################################################################################################
#
# Program:    makefile
#
# Function:   Makefile for dichrocalc
#
# Author:     Benjamin M. Bulheller
#
# Version:    $Revision: 4453 $, $Date: 2009-05-17 01:01:54 -0700 (Sun, 17 May 2009) $
#
# Date:       May 2009
#
####################################################################################################

# filename of the binary to compile
BINS = dichrocalc
# filename of the library to produce
LIBS = libdichrocalc.a

# the used compiler
CC   = g++
# the used archiver to create a library
AR   = ar
# directory for dichrocalc header files and the library
INC  = include
# directory for required third-party header files and library
INC3 = lib
# directory for object files
OBJ  = obj
# directory for .cpp source files
SRC  = src

# flags for the compiler, -O0: do not use optimization, -O3: use full optimization
CPPFLAGS = -Wall  -O0

# all .cpp files that have to be compiled for the library
LIBOBJS = $(OBJ)/iolibrary.o     \
          $(OBJ)/readinput.o     \
          $(OBJ)/fitparameters.o \
          $(OBJ)/matrix.o        \
          $(OBJ)/dichroism.o

# all .cpp files that have to be compiled for the main program
BINOBJS = $(OBJ)/dichrocalc.o $(LIBOBJS)

# linker flags (only for the main program)
LDFLAGS  = -lnewmat  -ldichrocalc  -lm

# directories containing libraries and header files
LIBDIRS  = -I./lib/ -L./lib   -I./include/ -L./include

# all NewMat header files used in some of the programs
NEWMAT = ${INC}/dichrocalc.h  \
         ${INC3}/newmat.h    \
         ${INC3}/newmatio.h  \
         ${INC3}/newmatap.h

# link the main program
# CAUTION: On linux, the sequence matters: first objects, then CPP and linker flags, then output file!!!
$(BINS): $(INC)/$(LIBS)  ${BINOBJS}
	@echo
	@echo "=> Building dichrocalc binary"
	@echo "   --------------------------"
	@echo
	$(CC)  $(BINOBJS) \
	$(CPPFLAGS)  $(LIBDIRS)  $(LDFLAGS)  -o ${BINS}
	@echo
	@echo

$(INC)/$(LIBS): $(LIBOBJS)
	@echo
	@echo "=> Building $(INC)/libdichrocalc.a"
	@echo "   ---------------------------------"
	@echo
	$(AR) ruv $(INC)/$(LIBS)  $(LIBOBJS)
	@echo

$(OBJ)/dichrocalc.o: $(SRC)/dichrocalc.cpp $(SRC)/readinput.cpp
	$(CC)  $(CPPFLAGS)  $(LIBDIRS)  -c $(SRC)/dichrocalc.cpp      -o $(OBJ)/dichrocalc.o

$(OBJ)/iolibrary.o: $(SRC)/iolibrary.cpp ${INC}/dichrocalc.h  $(NEWMAT)
	$(CC)  $(CPPFLAGS)  $(LIBDIRS)  -c  $(SRC)/iolibrary.cpp      -o $(OBJ)/iolibrary.o

$(OBJ)/readinput.o: $(SRC)/readinput.cpp $(INC)/dichrocalc.h  $(SRC)/iolibrary.cpp  $(NEWMAT)
	$(CC)  $(CPPFLAGS)  $(LIBDIRS)  -c  $(SRC)/readinput.cpp      -o $(OBJ)/readinput.o

$(OBJ)/fitparameters.o: $(SRC)/fitparameters.cpp ${INC}/dichrocalc.h  $(SRC)/iolibrary.cpp  $(NEWMAT)
	$(CC)  $(CPPFLAGS)  $(LIBDIRS)  -c  $(SRC)/fitparameters.cpp  -o $(OBJ)/fitparameters.o

$(OBJ)/matrix.o: $(SRC)/matrix.cpp $(INC)/dichrocalc.h  $(SRC)/iolibrary.cpp  $(NEWMAT)
	$(CC)  $(CPPFLAGS)  $(LIBDIRS)  -c  $(SRC)/matrix.cpp         -o $(OBJ)/matrix.o

$(OBJ)/dichroism.o: $(SRC)/dichroism.cpp ${INC}/dichrocalc.h  $(SRC)/iolibrary.cpp  $(NEWMAT)
	$(CC)  $(CPPFLAGS)  $(LIBDIRS)  -c  $(SRC)/dichroism.cpp      -o $(OBJ)/dichroism.o


clean:
	@echo
	@echo "=> Cleaning directories"
	@echo "   --------------------"
	@echo
	rm -rf $(BINS)  $(INC)/$(LIBS)
	rm -rf $(OBJ)/*.o
	@echo


new:
	@make clean
	@make $(BINS)


