package ReadMatmac;

#######################################################################################################################
#
# Package:         ReadMatmac
#
# Function:        Reads the contents of the matmac output files (.kor, .dbg)
#
# Author:          Benjamin Bulheller
#
# Research Group:  Prof. Jonathan Hirst
#                  School of Physical Chemistry
#                  University of Nottingham
#
# Funded by:       EPSRC
#
# Version:         $Revision: 4661 $, $Date: 2009-08-01 18:13:42 +0100 (Sat, 01 Aug 2009) $
#
# Date:            March 2006 - October 2008
#
# CAUTION:         The atom label determination basically works for the Hirst peptide set only. The Woody set has a
#                  different order of the atoms, so the C, N and O labels would be wrong in that case.
#
# Usage:
#
#   use ReadMatmac;
#
#   ReadKorFile ($BaseName, \@KorData); # $BaseName is the filename without extension
#   ReadDbgFile ($BaseName, \%DbgData);
#   ReadAbiFile ($BaseName, \@AbiData);
#   ReadMatFile ($BaseName, \%MatData); # NOTE: The matrix is sorted $Matrix->[Row][Column]
#   ReadCdlFile ($BaseName, \@CdlData);
#   ReadVecFile ($BaseName, \@VecData);
#   ReadPolFile ($BaseName, \@PolData);
#   ReadStfFile ($BaseName, \%StfData);
#   ReadFmtFile ($BaseName, \@FmtData);
#   ReadChromophoresDat (\@Chromophores);
#   WriteChromophoresDat (\@Chromophores);
#   TransitionSequence ($BaseName, \%StfData, $Sequence);
#
# @KorData
#    |
#    |- 0 (first Group)
#    |  |
#    |  |- Atoms
#    |  |    |
#    |  |    |- 0 (x, y, z, Atom)
#    |  |    |- 1 (x, y, z, Atom)
#    |  |    `- 2 (x, y, z, Atom)
#    |  |
#    |  |- XS (x, y, z)
#    |  |
#    |  |- Extrema
#    |  |    |
#    |  |    |- MinX, MinY, MinZ
#    |  |    `- MaxX, MaxY, MaxZ
#    |  |
#    |  |- AbInitio
#    |  |    |
#    |  |    |- 0 (nPi*)
#    |  |    |  |
#    |  |    |  |- Energy
#    |  |    |  |- VM  (x, y, z)
#    |  |    |  `- VMU (X, y, z)
#    |  |    |
#    |  |    `- 1 (PiPi*)
#    |  |       |
#    |  |       |- Energy
#    |  |       |- VM  (x, y, z)
#    |  |       `- VMU (X, y, z)
#    |  |
#    |  |
#    |  `- Transitions
#    |        |
#    |        |- 0 (first Transition)
#    |        |  |
#    |        |  |- VM (x, y, z)  (magnetic moment before calculation)
#    |        |  |- VMU (x, y, z) (electric moment before calculation)
#    |        |  |
#    |        |  |- VMCoD (x, y, z) (VM Center of Dipole)
#    |        |  |- VMUCoD (x, y, z) (VM Center of Dipole)
#    |        |  |
#    |        |  |- Monopoles
#    |        |  |      |
#    |        |  |      |- 0 (x, y, z, q, Atom, Pos, Distance) # x, y, z are the coordinates; q ist the charge, Atom (with group number) is the atom label
#    |        |  |      |- 1 (x, y, z, q, Atom, Pos, Distance) # Pos is the relative position of the monopole to the atom it belongs to, Distance the distance to it
#    |        |  |      |
#    |        |  |     etc
#    |        |  |
#    |        |  `- BalancePoints
#    |        |         |
#    |        |         |- 0 (x, y, z)
#    |        |         |- 1 (x, y, z)
#    |        |         |
#    |        |        etc
#    |        |
#    |        |- 1 (second Transition)
#    |        |  |
#    |        | etc
#    |       etc
#    |
#    |- 1 (second Group)
#    |  |
#    | etc
#   etc
#
#
# ---------------------------------------------------------------------------------------------------------------------
#
# The EM, EMU and ab initio vectors do not belong to a specific  group and therefore don't fit into @KorData
# The EM and EMU vectors are dumped from matmac after the diagonaliztion and derive from the VM and VMU, which
# in turn come from the single groups
#
# %DbgData                          %AbiData
#    |                                 |
#    |- DipStrength                    |- DipStrength
#    |                                 |
#    |- EM                             |- VM
#    |   |                             |   |
#    |   |- 0 (x,y,z)                  |   |- 0 (x,y,z)
#    |   |- 1 (x,y,z)                  |   |- 1 (x,y,z)
#    |   |                             |   |
#    |  etc                            |  etc
#    |                                 |
#    |- EMU                            |- VMU
#    |   |                             |   |
#    |   |- 0 (x,y,z)                  |   |- 0 (x,y,z)
#    |   |- 1 (x,y,z)                  |   |- 1 (x,y,z)
#    |   |                             |   |
#    |  etc                            |  etc
#    |                                 |
#    |- Energy (in cm^-1)              |- Energy (in cm^-1)
#    |   |                             |   |
#    |   |- 0                          |   |- 0
#    |   |- 1                          |   |- 1
#    |   |                             |   |
#    |   etc                           |   etc
#    |                                 |
#    |- RotStrength                    |- RotStrength
#    |                                 |
#    `- Wavelength (in nm)             `- Wavelength (in nm)
#        |                                 |
#        |- 0                              |- 0
#        |- 1                              |- 1
#        |                                 |
#       etc                               etc
#
#
#
# %StfData
#    |
#    |- Assignments
#    |      |
#    |      |- atom1[1, 3, 4]
#    |      |- atom2[5, 6, 7]
#    |      |
#    |     etc
#    |       
#    |- Groups         (the number of groups)
#    |- Interval       (the interval)
#    |- TypeNumber     (the number of types)
#    |- Factor         (the scale factor)
#    |- AtomNumber     (the number of atoms)
#    |- Max            (maximum wavelength)
#    |- Min            (minimum wavelength)
#    |- BBTRANS        (backbone transition to start with)
#    |- CTTRANS        (charge transfer transition to use)
#    |- ConvFactor     (conversion factor to Angstrom)
#    |         
#    |- Names
#    |    |
#    |    `-[NMA4FIT2, PHECASI3, TYRCASY4, ...]
#    |         
#    |- Types
#    |    |
#    |    `-[1, 1, 1, 2, 4, 1, 2, ...]
#    |
#    `- Transitions
#         |
#         `-[2, 2, 2, 4, 4, 4, 2, ...]



#######################################################################################################################


use strict;                         # alway use this!!!
use FindBin qw/$Bin/;               # sets $Bin the directory of the script
use lib "$ENV{HOME}/bin/perllib";   # add ~/bin/perllib to the library path
use lib $Bin;                       # add the script's directory to the library path
use Data::Dumper;	                  # for easy printout of arrays and hashes
use VectorMath;                     # for several vector calculations
use DebugPrint;                     # handy during debugging
use ParsePDB;                       # to parse PDB files

$Data::Dumper::Sortkeys = 1; # sort the hash keys
# $Data::Dumper::Indent = 3;

require Exporter;
our @ISA    = qw( Exporter );
our @EXPORT = qw( ReadKorFile ReadDbgFile ReadPolFile ReadAbiFile ReadMatFile ReadCdlFile ReadVecFile ReadStfFile 
                  ReadFmtFile ReadChromophoresDat WriteChromophoresDat TransitionSequence );

my $ERROR = "\nERROR (ReadMatmac.pm)";

#######################################################################################################################


sub ReadKorFile { # reads the .kor file and splits the content in its parts and sorts everything into the @KorData-Array
	my $BaseName = shift;
	my $KorData  = shift;
	my $GroupCount = 0;
	my $LineCount  = 0;
	my $Transition = 0;
	my (@Content, $AtomCount, $Extrema);
	
	if (not $BaseName) {
		print STDERR "$ERROR: No file name given to ReadKorFile!\n\n";
		exit 120;
	}
	
	if (not -f "$BaseName.kor") {
		if (-f $BaseName) {
			$BaseName =~ s/\.kor$//;
		}
		else {
			return undef;
		}
	}
	
	open  KOR, "<$BaseName.kor" or die "$ERROR: Could not open $BaseName.kor: $!";
	@Content = <KOR>;
	close KOR;
	chomp @Content;
	
	if ($KorData !~ m/^ARRAY\(0x[\dA-zA-Z]+\)/) {
		print STDERR "$ERROR: \$KorData needs to be an array reference for ReadKorFile!\n\n";
		exit 122;
	}
	
	@{$KorData} = (); # make sure that the array is empty
	
	$Extrema = { MinX =>  1E20, MinY =>  1E20, MinZ =>  1E20,
	             MaxX => -1E20, MaxY => -1E20, MaxZ => -1E20 };
	
	# read the information from the input file
	until ($LineCount > $#Content) {
		if ($Content[$LineCount] =~ m/Group:/) { # the following lines belong to the atoms of the group
         # skip line till the "ATOM" line
			until ( ($Content[$LineCount] =~ m/ATOM/) or ($LineCount > $#Content) ) { ++$LineCount }
			
			$AtomCount = 0;
			
			until ( ($Content[$LineCount] !~ m/ATOM/) or ($LineCount > $#Content) ) {
				my $CoordHash = &GetXYZ ($Content[$LineCount], $Extrema);
				
				push @{$KorData->[$GroupCount]{Atoms}}, $CoordHash;
				
				++$AtomCount;
				++$LineCount;
			}

			# # if three atoms were found, it is a peptide bond and the atom labels can be assigned
			# if ($AtomCount == 3) {
			# 	for $AtomCount ( 0 .. 2 ) {
			# 		if		($AtomCount == 0) { $KorData->[$GroupCount]{Atoms}[$AtomCount]{Atom} = "C" }
			# 		elsif ($AtomCount == 1) { $KorData->[$GroupCount]{Atoms}[$AtomCount]{Atom} = "O" }
			# 		elsif ($AtomCount == 2) { $KorData->[$GroupCount]{Atoms}[$AtomCount]{Atom} = "N" }
			# 	}
			# }
			# else {
			# 	for $AtomCount ( 0 .. $#{$KorData->[$GroupCount]{Atoms}} ) {
			# 		$KorData->[$GroupCount]{Atoms}[$AtomCount]{Atom} = " ";
			# 	}
			# }
				
         # skip lines till the "VECT  XS" line
			until ( ($Content[$LineCount] =~ m/VECT  XS/) or ($LineCount > $#Content) ) { ++$LineCount }
			my $CoordHashXS = &GetXYZ ($Content[$LineCount], $Extrema);
			$KorData->[$GroupCount]{XS} = $CoordHashXS;
			
         # skip lines till the "REF   XS" line
			until ( ($Content[$LineCount] =~ m/REF   XS/) or ($LineCount > $#Content) ) { ++$LineCount }
			my $CoordHashRefXS = &GetXYZ ($Content[$LineCount], $Extrema);
			$KorData->[$GroupCount]{RefXS} = $CoordHashRefXS;
			
         # skip lines till the "VECT  VM" line, this is then the first transition
			until ( ($Content[$LineCount] =~ m/VECT  VM/) or ($LineCount > $#Content) ) { ++$LineCount }
			$Transition = 0;
			
			until ( ($Content[$LineCount] =~ m/---------------------/) or ($LineCount > $#Content) ) {
				if		($Content[$LineCount] =~ m/VECT  VM/) {
					my $CoordHash = &GetXYZ ($Content[$LineCount], $Extrema);
					
					$KorData->[$GroupCount]{Transitions}[$Transition]{VM} = $CoordHash;
				}
				elsif ($Content[$LineCount] =~ m/VECT VMU/) {
					my $CoordHash = &GetXYZ ($Content[$LineCount], $Extrema);
					
					$KorData->[$GroupCount]{Transitions}[$Transition]{VMU} = $CoordHash;
				}
				elsif ($Content[$LineCount] =~ m/MONOPOLE/) {
					my $CoordHash = &GetXYZ ($Content[$LineCount], $Extrema);
					push @{$KorData->[$GroupCount]{Transitions}[$Transition]{Monopoles}}, $CoordHash;
				}
				
				++$LineCount;
				if ($LineCount > $#Content) { last }
				if ($Content[$LineCount] =~ m/VECT  VM/) { ++$Transition }
			}
			
			++$GroupCount;
		} # of if ($Content[$LineCount] =~ m/Group:/)
		else { ++$LineCount }
	} # of until ($LineCount > $#Content)
	
	&DetermineBalancePoints ($KorData);
	&DetermineDipoleCenters ($KorData);
	&FindHydrogen ($KorData);
	&AssociateMonopoles ($KorData);
	
	&CleanExtrema ($Extrema);
	
	return $Extrema;
} # of sub ReadKorFile


sub ReadAbiFile { # reads in data from an .abi file, if it exists (tab separated table with ab inito results, sample at end of the routine) 
	my $BaseName = shift;
	my $AbiData  = shift;
	my (@Content, %AbiData,  $GroupCount, $Transition, $Line, @Fields, $Position, $i, $Dot);
	
	if (not $BaseName) {
		print "No file name given to ReadAbiFile!\n\n";
		exit 123;
	}
	
	if (not -f "$BaseName.abi") {
		if (-f $BaseName) {
			$BaseName =~ s/\.abi$//;
		}
		else {
			return;
		}
	}
	
	if ($AbiData !~ m/^HASH\(0x[\dA-zA-Z]+\)/) {
		print STDERR "$ERROR: \$AbiData needs to be a hash reference for ReadAbiFile!\n\n";
		exit 125;
	}
	
	%{$AbiData} = (); # make sure that the hash is empty
	
	open (FILE, "<$BaseName.abi");
	@Content = <FILE>;
	close FILE;
	chomp @Content; # remove all line feeds
	
	foreach $Line (@Content) {
		@Fields = split /\t/, $Line;
		$Position = 0;
		
		# the matmac eigenvalues are sorted by energy (in cm^-1, ascending). In order to make the values
		# compareable, the ab initio values also have to be in ascending order
		if ($Fields[3] =~ m/(pi\*)/) {
			for $i (0 .. $#{$AbiData->{Energy}}) { # if first value => $#array = -1 => loop is skipped => $Position = 0
				
				if ($AbiData->{Energy}[$i] > $Fields[5]) { # if the energy $i is bigger, insert the current energy at this position
					$Position = $i;
					last;
				}
				else { ++$Position } # increase $Position in case that no existing item is bigger than the new energy 
			}
			
			if (not @{$AbiData->{Energy}}) {
					push @{$AbiData->{Energy}},     $Fields[5];
					push @{$AbiData->{Wavelength}}, $Fields[4];
					push @{$AbiData->{EMU}},        { x => $Fields[6], y => $Fields[7],  z => $Fields[8]  };
					push @{$AbiData->{EM}},         { x => $Fields[9], y => $Fields[10], z => $Fields[11] };
					
					# to convert units from MOLCASS (conversion factor 2.54)
					# $AbiData->{EM} [$#{$AbiData{EM}}] = &Product (2.54, $AbiData{EM} [$#{$AbiData{EM}}]);
					# $AbiData->{EMU} [$#{$AbiData{EMU}}] = &Product (2.54, $AbiData{EMU}[$#{$AbiData{EMU}}]);
					
					$AbiData->{Phi} = $Fields[1];
					$AbiData->{Psi} = $Fields[2];
			}
			else {
				# add the item to the arrays at the determined position (without overwriting anything => "0")
				splice (@{$AbiData->{Energy}},     $Position, 0, $Fields[5]);
				splice (@{$AbiData->{Wavelength}}, $Position, 0, $Fields[4]);
				splice (@{$AbiData->{EMU}},        $Position, 0, { x => $Fields[6], y => $Fields[7],  z => $Fields[8]  });
				splice (@{$AbiData->{EM}},         $Position, 0, { x => $Fields[9], y => $Fields[10], z => $Fields[11] });
				
				# to convert units from MOLCAS (conversion factor 2.54)
				# $AbiData->{EM} [$Position] = &Product (1/2.54, $AbiData{EM} [$Position]);
				# $AbiData->{EMU} [$Position] = &Product (1/2.54, $AbiData{EMU}[$Position]);
				
				$AbiData->{Phi} = $Fields[1];
				$AbiData->{Psi} = $Fields[2];
			}
		}
	} # of foreach $Line (@Content)
	
	for $i (0 .. $#{$AbiData->{EMU}}) {
		# calculate the rotational strength
		$Dot = &DotProduct ($AbiData->{EM}[$i], $AbiData->{EMU}[$i]);
		push @{$AbiData->{ROTS}}, sprintf ("%.8f", $Dot);
		
		# calculate the dipole strenght
		$Dot = &DotProduct ($AbiData->{EMU}[$i], $AbiData->{EMU}[$i]);
		push @{$AbiData->{DIPS}}, sprintf ("%.8f", $Dot);
	}

# Sample .abi Input File:
# ~~~~~~~~~~~~~~~~~~~~~~~
#  Geom.	Phi	Psi	Trans.	Wavelength	Energy	VMUx	VMUy	VMUz	Vx	Vy	Vz
#  01	180	180	n1>pi*1	221.7906389	45087.56569	0.004338324	0.018306266	0.161073244	-0.78488841	0.035848442	0.003755469
#  01	180	180	n2>pi*2	216.5906558	46170.0435	0.006838205	-0.014341866	-0.169189801	-0.930077279	0.061665519	-0.001288004
#  01	180	180	pi2>pi*2	189.7015206	52714.39031	-1.173269129	3.582033526	-0.239185241	-0.01184243	-0.240296061	-0.011108707
#  01	180	180	pi1>pi*1	185.0664503	54034.6453	3.623236812	1.471595788	-0.090903162	-0.11467423	-1.494670756	0.021011155

	return 1;
} # of sub ReadAbiFile


sub ReadDbgFile { # reads the .dbg file if it exists
	my $BaseName = shift;
	my $DbgData  = shift;
	my (@Content, %DbgData, $Line, $xCoord, $yCoord, $zCoord, $ROTS, $DIPS, $Energy, $Wavelength);
	
	if (not $BaseName) {
		print "No file name given to ReadDbgFile!\n\n";
		exit 125;
	}
	
	if (not -f "$BaseName.dbg") {
		if (-f $BaseName) {
			$BaseName =~ s/\.dbg$//;
		}
		else {
			return;
		}
	}
	
	open (FILE, "<$BaseName.dbg");
	@Content = <FILE>;
	close (FILE);
	chomp @Content; # remove all line feeds
	
	if ($DbgData !~ m/^HASH\(0x[\dA-zA-Z]+\)/) {
		print STDERR "$ERROR: \$DbgData needs to be a hash reference for ReadDbgFile!\n\n";
		exit 130;
	}
	
	%{$DbgData} = (); # make sure that the hash is empty
	
	until ( (not @Content) or ($Content[$#Content] =~ m/CDCALC/) ) {
		$Line = pop @Content; # start at the end, .dbg files can be VERY large
		$Line =~ s/^\s+//; # remove all leading blanks
		
		if ($Line =~ m/^(EM  \(|EMU \()/) {
			$xCoord = substr ($Line, 16, 12);
			$xCoord =~ s/\s//g;
			
			$yCoord = substr ($Line, 28, 12);
			$yCoord =~ s/\s//g;
			
			$zCoord = substr ($Line, 40, length $Line);
			$zCoord =~ s/\s//g;
		}
		elsif ($Line =~ m/^ROTS/) {
			$ROTS = substr ($Line, 16, 12);
			$ROTS =~ s/\s//g;
			
			$DIPS = substr ($Line, 41, length $Line);
			$DIPS =~ s/\s//g;
			
			unshift @{$DbgData->{ROTS}}, $ROTS;
			unshift @{$DbgData->{DIPS}}, $DIPS;
			
			next;
		}
		elsif ($Line =~ m/^VEC1/) {
			$Energy = substr ($Line, 16, 12);
			$Energy =~ s/\s//g;
			
			$Wavelength = substr ($Line, 41, length $Line);
			$Wavelength =~ s/\s//g;
			
			unshift @{$DbgData->{Energy}}, $Energy;
			unshift @{$DbgData->{Wavelength}}, $Wavelength;
			
			next;
		}
		else { next }
		
		if ($Line =~ m/^EM  \(/) {
			unshift @{$DbgData->{EM}}, { x => $xCoord, y => $yCoord, z => $zCoord };
		}
		elsif ($Line =~ m/^EMU \(/) {
			unshift @{$DbgData->{EMU}}, { x => $xCoord, y => $yCoord, z => $zCoord };
		}
	} # of until ($Content[0] =~ m/^CDCALC/)

	return 1;
} # of sub ReadDbgFile


sub ReadPolFile { # reads the .pol file if it exists
	my $BaseName = shift;
	my $PolData  = shift;
	my (@Content, $Line, @Fields);
	
	if (not $BaseName) {
		print "No file name given to ReadPolFile!\n\n";
		die;
	}
	
	if (not -f "$BaseName.pol") {
		if (-f $BaseName) {
			$BaseName =~ s/\.pol$//;
		}
		else {
			return;
		}
	}
	
	open (FILE, "<$BaseName.pol");
	@Content = <FILE>;
	close (FILE);
	chomp @Content; # remove all line feeds
	
	if ($PolData !~ m/^ARRAY\(0x[\dA-Za-z]+\)/) {
		print STDERR "$ERROR: \$PolData needs to be an array reference for ReadPolFile!\n\n";
		exit 135;
	}
	
	@{$PolData} = (); # make sure that the array is empty
	
	$Line = "";
	until ( (not @Content) or ($Line =~ m/^Trans/) ) { $Line = shift @Content }
	
	while (@Content) {
		$Line = shift @Content;
		&StripLine (\$Line);
		
		# if the line does not begin with a number, skip it
		if ($Line !~ m/^\d/) { next }
		
		@Fields = split /[\t\s]+/, $Line;
		
		# if less then 9 columns are found, skip the line
		if (scalar @Fields < 8) { next }
		
		push @{$PolData},
		{ Wavelength => $Fields[1],
		  RotStrength => $Fields[8],
		  VMU   => { x => $Fields[2], y => $Fields[3], z => $Fields[4] },
		  nPi   => $Fields[5],
		  PiPi  => $Fields[6],
		  Total => $Fields[7]
	   };
	}
	
	return 1;
} # of sub ReadPolFile


sub ReadCdlFile { # reads the .cdl file if it exists
	my $BaseName = shift;
	my $CdlData  = shift;
	my (@Content, $Line, @Fields);
	
	if (not $BaseName) {
		print "No file name given to ReadCdlFile!\n\n";
		exit 127;
	}
	
	if ($CdlData !~ m/^ARRAY\(0x[\dA-zA-Z]+\)/) {
		print STDERR "$ERROR: \$CdlData needs to be an array reference for ReadCdlFile!\n\n";
		exit 140;
	}
	
	if (not -f "$BaseName.cdl") { return }
	
	@{$CdlData} = (); # make sure that the array is empty
	
	open (FILE, "<$BaseName.cdl");
	@Content = <FILE>;
	close (FILE);
	chomp @Content; # remove all line feeds
	
	while (@Content) {
		$Line = shift @Content;
		&StripLine (\$Line);
		if ($Line eq "") { next } # skip empty lines
		
		@Fields = split /\t/, $Line;
		push @{$CdlData}, [$Fields[0], $Fields[1]];
	}

	return 1;
} # of sub ReadCdlFile


sub ReadVecFile { # reads the .vec file if it exists
	my $BaseName = shift;
	my $VecData  = shift;
	my (@Content, $Line, @Fields);
	
	if (not $BaseName) {
		print "No file name given to ReadVecFile!\n\n";
		exit 129;
	}
	
	if (not -f "$BaseName.vec") { return }
	
	if ($VecData !~ m/^ARRAY\(0x[\dA-zA-Z]+\)/) {
		print STDERR "$ERROR: \$VecData needs to be an array reference for ReadVecFile!\n\n";
		exit 145;
	}
	
	@{$VecData} = (); # make sure that the array is empty
	
	open (FILE, "<$BaseName.vec");
	@Content = <FILE>;
	close (FILE);
	chomp @Content; # remove all line feeds
	
	while (@Content) {
		$Line = shift @Content;
		&StripLine (\$Line);
		if ($Line eq "") { next } # skip empty lines
		
		@Fields = split /[\t\s]+/, $Line;
		
		if (scalar @Fields < 4) {
			print STDERR "$ERROR: ReadVecFile expected 4 coloums in .vec (wavelength, x, y, z), but only $#Fields were found.\n\n";
			exit 150;
		}
		
		push @{$VecData}, { Wavelength => $Fields[0], x => $Fields[1], y => $Fields[2], z => $Fields[3]};
	}
	
	return 1;
} # of sub ReadVecFile


sub ReadMatFile { # reads the .mat file if it exists
	my $BaseName = shift;
	my $MatData = shift;
	my (@Matrix, @Vectors, @Values);
	my (@Content, $i);
	
	if (not $BaseName) {
		print "No file name given to ReadMatFile!\n\n";
		exit 131;
	}
	
	if (not -f "$BaseName.mat") { return }
	
	if ($MatData !~ m/^HASH\(0x[\dA-zA-Z]+\)/) {
		print STDERR "$ERROR: \$MatData needs to be a hash reference for ReadMatFile!\n\n";
		exit 155;
	}
	
	open (FILE, "<$BaseName.mat");
	@Content = <FILE>;
	close (FILE);
	chomp @Content; # remove all line feeds
	
	# remove all leading, trailing and multiple spaces from the array
	for ($i = 0; $i <= $#Content; $i++) {
		$Content[$i] =~ s/\s+/ /g;         # substitute a double space with a single space
		$Content[$i] =~ s/^\s+|\s+$//g;    # remove all leading and trailing spaces to avoid problem with split
		$Content[$i] =~ s/[\t\n]+//g;      # remove all tabs and line feeds just to be sure
	}
	
	&ReadMatrix (\@Content, \@Matrix);
	&ReadMatrix (\@Content, \@Vectors);
	&ReadValues (\@Content, \@Values);
	
	$MatData->{Matrix}  = \@Matrix;
	$MatData->{Vectors} = \@Vectors;
	$MatData->{Values}  = \@Values;

	return 1;
} # of sub ReadMatFile


sub ReadStfFile { # reads the stf file if it is found
	my $BaseName = shift;
	my $StfData  = shift;
	my (@Content, $LineCount, $Line, @Fields);
	
	if (not -f "$BaseName.stf") { return }
	
	if ($StfData !~ m/^HASH\(0x[\dA-zA-Z]+\)/) {
		print STDERR "$ERROR: \$StfData needs to be a hash reference for ReadStfFile!\n\n";
		exit 160;
	}
	
	open STF, "<$BaseName.stf";
	@Content = <STF>;        # read the whole file
	chomp @Content;          # remove all line feeds
	close STF;
	
	$LineCount = 0;
	
	while (@Content) {
		$Line = shift @Content;
		$Line =~ s/^\s+|\s+$//g; # remove all leading and trailing blanks
		
		if ( ($Line =~ m/^#/) or ($Line eq "") ) { next }         # if line is coment or blank, skip it
		                                    else { ++$LineCount } # count each non-comment line
		
		@Fields = split /\s+/, $Line;
		
		if ($LineCount == 1) { # line after the first comment
			$StfData->{Groups}      = $Fields[0];
			$StfData->{TypeNumber}  = $Fields[1];
			$StfData->{BBTRANS}     = $Fields[2];
			$StfData->{CTTRANS}     = $Fields[3];
			$StfData->{Factor}      = $Fields[4];
		}
		elsif ($LineCount == 2) { # line after the second comment
			$StfData->{ConvFactor} = $Fields[0];
		}
		elsif ($LineCount == 3) { # line after the third comment
			$StfData->{Min}      = $Fields[0];
			$StfData->{Max}      = $Fields[1];
			$StfData->{Interval} = $Fields[2];
		}
		elsif ($LineCount == 4) { # names of the monopole sets
			# there's only one name per line, thus $Line is used instead of @Fields
			$Line =~ s/^-|-$//g;             # remove all leading and trailing dashes
			push @{$StfData->{Names}}, $Line;
			
			until ($Content[0] =~ m/^#/) {   # read the lines until the next comment
				$Line = shift @Content;       # read the next line
				$Line =~ s/^\s+|\s+$//g;      # remove all leading and trailing blanks
				$Line =~ s/^-|-$//g;          # remove all leading and trailing dashes
				push @{$StfData->{Names}}, $Line;
			}
		}
		elsif ($LineCount == 5) { # numbers of transitions per group
			push @{$StfData->{Transitions}}, @Fields;
			
			until ($Content[0] =~ m/^#/) {   # read the lines until the next comment
				$Line = shift @Content;       # read the next line
				$Line =~ s/^\s+|\s+$//g;      # remove all leading and trailing blanks
				@Fields = split /\s+/, $Line;
				push @{$StfData->{Transitions}}, @Fields;
			}
		}
		elsif ($LineCount == 6) { # types of groups
			push @{$StfData->{Types}}, @Fields;
			
			until ($Content[0] =~ m/^#/) {   # read the lines until the next comment
				$Line = shift @Content;       # read the next line
				$Line =~ s/^\s+|\s+$//g;      # remove all leading and trailing blanks
				@Fields = split /\s+/, $Line;
				push @{$StfData->{Types}}, @Fields;
			}
		}
		elsif ($LineCount == 7) { # number of atoms
			# push @{$StfData->{AtomNumber}}, $Fields[0];
			$StfData->{AtomNumber} = $Fields[0];
		}
		elsif ($LineCount == 8) { # atomic assignments
			my @Group = @Fields;   # create a copy of the array (since a reference need to be used!)
			push @{$StfData->{Assignments}}, \@Group; # add the reference to the new array
			
			until (not @Content or $Content[0] =~ m/^#/) {   # read the lines until the next comment
				$Line = shift @Content;          # read the next line
				$Line =~ s/^\s+|\s+$//g;         # remove all leading and trailing blanks
				my @Group = split /\s+/, $Line;  # create a new array for the split line
				push @{$StfData->{Assignments}}, \@Group;
			}
		}
	}

	return 1;
} # of sub ReadStfFile


sub ReadFmtFile { # reads the .fmt file if it exists
	my $BaseName = shift;
	my $FmtData  = shift;
	my (@Content, $Line, @Fields);
	
	if (not $BaseName) {
		print "No file name given to ReadFmtFile!\n\n";
		exit 127;
	}
	
	if (not -f "$BaseName.fmt") { return }
	
	if ($FmtData !~ m/^ARRAY\(0x[\dA-zA-Z]+\)/) {
		print STDERR "$ERROR: \$FmtData needs to be an array reference for ReadFmtFile!\n\n";
		exit 165;
	}
	
	@{$FmtData} = (); # make sure that the array is empty
	
	open (FILE, "<$BaseName.fmt");
	@Content = <FILE>;
	close (FILE);
	chomp @Content; # remove all line feeds
	
	while (@Content) {
		$Line = shift @Content;
		&StripLine (\$Line);
		if ($Line eq "") { next } # skip empty lines
		
		@Fields = split /[\s\t]+/, $Line;
		push @{$FmtData}, [$Fields[0], $Fields[1], $Fields[2]];
	}
	
	return 1;
} # of sub ReadFmtFile


sub ReadChromophoresDat { # reads chromophores.dat
	my $Chromophores = shift;
	my ($File, @Content, $Line, @Fields);
	
	if ($Chromophores !~ m/^ARRAY\(0x[\dA-zA-Z]+\)/) {
		print STDERR "$ERROR: \$Chromophores needs to be an array reference for ReadChromophores!\n\n";
		exit 170;
	}
	
	@{$Chromophores} = (); # make sure that the array is empty
	
	if (-f "chromophores.dat") { # try at first in the current directory
		$File = "chromophores.dat";
	}
	else { # otherwise check at ~/bin
		if (-f "$ENV{HOME}/bin/chromophores.dat") {
			$File = "$ENV{HOME}/bin/chromophores.dat";
		}
		else {
			print STDERR "$ERROR: chromophores.dat not found in current directory\n";
			print STDERR "                             or $ENV{HOME}/bin!\n\n";
			return;
		}
	}
	
	open DATFILE, "<$File" or die "Error opening $File: $!";
	@Content = <DATFILE>;
	chomp @Content;
	close DATFILE;
	
	while (@Content) {
		$Line = shift @Content;
		$Line =~ s/^\s+|\s+$//g;   # remove all leading and trailing blanks
		
		if ($Line eq "")    { next }  # skip emtpy lines
		if ($Line =~ m/^#/) { next }  # skip comment lines
		
		@Fields = split /\s+/, $Line;
		if (scalar @Fields != 5) {
			print STDERR "$ERROR: Entry $Line in chromophores.dat could not be interpreted.\n\n";
		}
		
		my $CurChrom             = {}; # create an anonymous hash
		$CurChrom->{Name}        = $Fields[0];
		$CurChrom->{Type}        = $Fields[1];
		$CurChrom->{Transitions} = $Fields[2];
		$CurChrom->{Phi}         = $Fields[3];
		$CurChrom->{Psi}         = $Fields[4];
		
		push @{$Chromophores}, $CurChrom;
	}

	return 1;
} # of sub ReadChromophoresDat


sub WriteChromophoresDat { # create the chromophores.dat for dcinput
	my $Chromophores = shift;
	my $DATFORMAT = "%8s    %3s     %d   %6.1f   %6.1f\n"; # printf format for chromophores.dat
	my $Chrom;
	
	if ($Chromophores !~ m/^ARRAY\(0x[\dA-zA-Z]+\)/) {
		print STDERR "$ERROR: \$Chromophores needs to be an array reference for WriteChromophores!\n\n";
		exit 173;
	}
	
	open DAT, ">chromophores.dat";
	
	print DAT "# Name     Type   Trans   Phi      Psi\n";
	
	foreach $Chrom ( @{$Chromophores} ) {
		printf DAT $DATFORMAT,
			$Chrom->{Name}, $Chrom->{Type}, $Chrom->{Transitions},
			$Chrom->{Phi},  $Chrom->{Psi};
	}
	close DAT;
} # of sub WriteChromophoresDat


sub TransitionSequence { # determines the sequence of transitions in the Hamiltonian Matrix
	my $BaseName = shift;
	my $StfData  = shift;
	my $Sequence = shift;
	my $Groups   = shift;
	
	my ($Trans, $Type, $Name, $Transitions, $TransNumber, $PDB, $Index, $Group, $Atom, $ResidueNumber, $Label);
	my ($GroupIndex);
	
	if ($StfData !~ m/^HASH\(0x[\dA-zA-Z]+\)/) {
		print STDERR "$ERROR: \$StfData needs to be a hash reference containing the result from ReadStfFile!\n\n";
		exit 174;
	}
	
	if ($Sequence !~ m/^ARRAY\(0x[\dA-zA-Z]+\)/) {
		print STDERR "$ERROR: \$Sequence needs to be a hash reference containing the result from ReadStfFile!\n\n";
		exit 174;
	}
	
	
	if (not defined $StfData->{Types}) {
		print STDERR "$ERROR: \$StfData must already contain the result from ReadStfFile!\n\n";
		exit 175;
	}
	
	# check whether the PDB file is available
	if (not -f "$BaseName.pdb") {
		print STDERR "$ERROR: $BaseName.pdb was not found but is required for the analysis of the Hamiltonian matrix.\n\n";
		exit 177;
	}
	
	# the following routine was taken from dcinput and will produce the same PDB
	$PDB = ParsePDB->new (FileName => $BaseName, Verbose => 0, NoHETATM => 1, NoANISIG => 1); # create a new object
	
	$PDB->Parse;
	$PDB->RemoveInsertedResidues (Model => 0, Intensive => 1);
	$PDB->RemoveAtomLocations (Model => 0, AtomLocations => "First");
	
	# In the parser it is ensured that the atom number in the assignment section corresponds with the actual
	# line number in the PDB. To get the original residue number, the residues are not renumbered in this case.
	# $PDB->RenumberChains (Model => 0);                          # correct or add the ChainIDs
	# $PDB->RenumberResidues (Model => 0, KeepInsertions => 0);   # renumber the residues
	# $PDB->RenumberAtoms (Model => 0, IgnoreTER => 1);           # renumber the atoms
	
	$Transitions = {
		"NMA4FIT2" => [qw/nPi*      PiPi*      Pib_Pi*      n'Pi*/],
		"NMA99WDY" => [qw/nPi*      PiPi*/],
		"CT01009A" => [qw/CT-P1P2   CT-n1P2    CT-n2P1      CT-P2P1/],
		"CTBE009A" => [qw/CT-P1P2   CT-n1P2    CT-n2P1      CT-P2P1/],
		"CT-BE09A" => [qw/CT-P1P2   CT-n1P2    CT-n2P1      CT-P2P1/],
		"CT05009A" => [qw/CT-P1P2   CT-n1P2    CT-P2P1      CT-n2P1/],
		"CT-0509A" => [qw/CT-P1P2   CT-n1P2    CT-P2P1      CT-n2P1/],
		"CT02009B" => [qw/CT-P1P2   CT-n1P2    CT-n2P1      CT-P2P1/],
		"CT-0209B" => [qw/CT-P1P2   CT-n1P2    CT-n2P1      CT-P2P1/],
		"CT08009A" => [qw/CT-P1P2   CT-n1P1    CT-n2P1      CT-P2P1/],
		"CT-0809A" => [qw/CT-P1P2   CT-n1P2    CT-n2P1      CT-P2P1/],
		"CT03009A" => [qw/CT-P1P2   CT-n1P2    CT-n2P2      CT-P2P1/],
		"CT-0309A" => [qw/CT-P1P2   CT-n1P2    CT-n2P1      CT-P2P1/],
		"CT31009A" => [qw/CT-n1P2   CT-P1P2    CT-P2P1      CT-n2P1/],
		"CT-3109A" => [qw/CT-n1P2   CT-P1P2    CT-P2P1      CT-n2P1/],
		"CTBT009B" => [qw/CT-P1P2   CT-n2P1    CT-n1P2      CT-P2P1/],
		"CT-BT09B" => [qw/CT-P1P2   CT-n2P1    CT-n1P2      CT-P2P1/],
		"CTAL009A" => [qw/CT-n2P1   CT-P1P2    CT-P2P1      CT-n1P2/],
		"CT-AL09A" => [qw/CT-n2P1   CT-P1P2    CT-P2P1      CT-n1P2/],
		"CT15009A" => [qw/CT-P1P2   CT-n2P1    CT-P2P1      CT-n1P2/],
		"CT-1509A" => [qw/CT-P1P2   CT-n2P1    CT-P2P1      CT-n1P2/],
		"PHECASI3" => [qw/PHE-1Lb   PHE-1La    PHE-1Ba      PHE-1Bb/],
                "PHEtest1" => [qw/hah-1   hah-2    PHE-1Ba      PHE-1Bb     hah-3    hah-4    hah-5/],
		"TYRCASY4" => [qw/TYR-1Lb   TYR-1La    TYR-1Ba      TYR-1Bb/],
		"TRPCASH1" => [qw/TRP-1Lb   TRP-1La    TRP-1Bb      TRP-1Ba/],
		"ASNHIRST" => [qw/ASN-nPi*  ASN-PiPi*  ASN-Pib_Pi*  ASN-n'Pi*/],
		"GLNHIRST" => [qw/GLN-nPi*  GLN-PiPi*  GLN-Pib_Pi*  GLN-n'Pi*/],
		"ASP_CNVE" => [qw/ASP-nPi*  ASP-PiPi*/],
		"GLU_CNVE" => [qw/GLU-nPi*  GLU-PiPi*/],
		"NAP2S2FC" => [qw/B2uFC1 B2uFC2 B3uFC1 B3uFC2/],
		"NAP2TRAN" => [qw/B2u B3u/],
		"NAPH_10X" => [qw/1B1u 1B3u 1B2u 2B1u 2B3u 2B2u 3B3u/],
		"FCWEIGH7" => [qw/B3u-1 B3u-2 B3u-3 B2u-2 B2u-1/],
		"ADE-FUL3" => [qw/2A' 4A'/],
		"GUAWFUL3" => [qw/2A' 6A'/],
		"CYT-FUL3" => [qw/2A' 4A'/],
		"THY-FUL3" => [qw/2A' 4A'/],
		"URA-FUL3" => [qw/2A' 5A'/],
	};
	
	for $Index ( 0 .. $#{$StfData->{Types}} ) {
		# used in the label to make the output better readable, starting at 1
		$GroupIndex = $Index + 1;
		
		$Type = $StfData->{Types}[$Index];
		
		# all type numbers have to be decreased by 1 to be used as the array index to access the name
		--$Type;
		
		$Name = $StfData->{Names}[$Type];
		$TransNumber = $StfData->{Transitions}[$Type];
		
		# if the exact transitions/states were not defined above, simply number all transitions
		if (not defined $Transitions->{$Name}) {
			# print STDERR "WARNING (ReadMatmac.pm): The chromophore type $Name was not recognized in TransitionSequence.\n";
			
			for $Trans ( 0 .. $TransNumber ) {
				push @{$Transitions->{$Name}}, "Trans" . ($Trans+1);
			}
		}
		
		# $Type is the current index in the {Types} array. {Types}, {Names} and {Transitions} all contain
		# the data of one group with this index.
		
		# read the array with the atom numbers of the group (e.g. 3 4 7 for C, O, N)
		$Group = $StfData->{Assignments}[$Index];
		
		# read just the atom number of the C atom, this corresponds with the line number which is hence the 
		# Atom (NOT the AtomNumber) in the parser
		$Atom = $Group->[0] - 1; # -1 because the parser starts counting at 0
		
		# retrieve the $Atom-th atom from the parser, AtomIndex returns the readily parsed atom hash
		# (it's an array nevertheless since whole residues could be retrieved like that, with the atoms as array elements)
		my @AtomHash = $PDB->Get (Model => 0, Atom => $Atom, AtomIndex => 1);
		$ResidueNumber = $AtomHash[0]->{ResidueNumber};
		
		# verify some residues
		for $Label (qw/PHE TYR TRP ASP ASN/) {
			if ($StfData->{Names}[$Type] =~ m/$Label/ and $AtomHash[0]->{ResidueLabel} !~ m/$Label/) {
				print "ERROR: Atom number $Atom was searched, it should be $StfData->{Names}[$Type], but returned was\n";
				print Dumper @AtomHash;
				exit 177;
			}
		}
		
		if ($StfData->{Names}[$Type] =~ m/NMA4FIT2|NMA99WDY/) {
			if ($StfData->{BBTRANS} != 0) {
				print STDERR "$ERROR: BBTRANS not equal 0, this is not yet implemented in TransitionSequence!\n";
				exit 178;
			}
		}
		
		if ($StfData->{Names}[$Type] =~ m/^CT/) {
			if ($StfData->{CTTRANS} != 0) {
				print STDERR "$ERROR: CTTRANS not equal 0, this is not yet implemented in TransitionSequence!\n";
				exit 179;
			}
		}
		
		for $Trans ( 0 .. $TransNumber-1) { # -1: Fortran starts counting at 1, Perl at 0
			# mind that the matrix-analysis script splits up the interaction label to determine the interaction residues
			# so changes to the label format have to be considered there
			push @{$Sequence}, "R$ResidueNumber-G$GroupIndex-$Transitions->{$Name}[$Trans]";
			push @{$Groups}, $StfData->{Assignments}[$Index];
		}
	}
} # of sub TransitionSequence


sub ReadMatrix { # reads the hamiltonian matrix or eigen vectors from the .mat file (also handles the splitting of if over several blocks)
	my $Content = shift;
	my $Matrix = shift;
	my ($Line, @Fields, $i, @ColValues);
	
	my $Debug = 0;
	
	if ($Matrix !~ m/^ARRAY\(0x[\dA-zA-Z]+\)/) {
		print STDERR "$ERROR: \$Matrix needs to be an array reference for ReadMatrix!\n\n";
		exit 175;
	}
	
	# all leading and trailing blanks have already been removed, the beginning of the matrix
	# is thus marked by the column indices:
	# 1 2 3 4 5 6
	until ( (not @{$Content}) or ($Content->[0] =~ m/^1 /) ) { shift @{$Content} }
	
	if (not @{$Content}) { return }
	
	if ($Debug) { print Dumper "Content:", $Content; }
	
	while (@{$Content} and $Content->[0] !~ m/^Eigen/i) {    # read until "Eigenvectors" or "Eigenvalues" is reached
		$Line = shift @{$Content};
		if ($Line eq "") { next }             # e.g. an empty line between two blocks
		
		if ($Debug) { print Dumper " Current line:  $Line" }
		
		# read the first line of the block with the column indices, e.g. "1 2 3 4 5 6"
		@ColValues = split (/\s+/, $Line);      # these values are taken as array indices
		foreach (@ColValues) { --$_ }           # decrease all ColValues by 1 (Fortran starts counting at 1, Perl at 0)
		
		# allow an empty line after the column values, remove it if present
		if ($Content->[0] eq "") { shift @{$Content} }
		
		if ($Debug) { print Dumper "Column values:", \@ColValues; }
		
		until (not @{$Content} or $Content->[0] eq "") {         # until a new block begins
			$Line = shift @{$Content};
			@Fields = split /\s+/, $Line;      # $Fields[0] is the yValue
			--$Fields[0];                      # decrease the RowValue by 1 to speak Perl and not Fortran
			
			if ($Debug) { print Dumper $Line, "Fields", @Fields, $#Fields; }
			
			for $i (1 .. $#Fields) {           # add the read values to the matrix
			
				if ($Debug) {
					print "Fields[0] $Fields[0]  ColValues[i-1] $ColValues[$i-1]\n";
				}
				
				# NOTE: The matrix is sorted $Matrix->[Row][Column]
				$Matrix->[ $Fields[0] ][ $ColValues[$i-1] ] = $Fields[$i];
			}
		}
	}
	
	if ($Debug) { dpf ($Matrix) }
	
	return 1;
} # of sub ReadMatrix


sub ReadValues {
	my $Content = shift;
	my $Values = shift;
	my ($Line, @Fields);
	
	if ($Values !~ m/^ARRAY\(0x[\dA-zA-Z]+\)/) {
		print STDERR "$ERROR: \$Values needs to be an array reference for ReadValues!\n\n";
		exit 180;
	}
	
	until ( (not @{$Content}) or ($Content->[0] =~ m/^1 /) ) { shift @{$Content} }
	
	if (not @{$Content}) { return }
	
	while ( (@{$Content}) and ($Content->[0] ne "") ) {
		$Line = shift @{$Content};
		@Fields = split /\s+/, $Line;
		push @{$Values}, $Fields[1];
	}

	return 1;
} # of sub ReadValues


sub DetermineBalancePoints { # calculates the balance points of each "monopole set" and the center of the dipoles
	my $KorData = shift;
	my ($GroupCount, $Transition, $Monopole, $CurMono, $OldMono, $BalancePoint, $Atom);
	my ($MonopoleAdded);
	
	if ($KorData !~ m/^ARRAY\(0x[\dA-zA-Z]+\)/) {
		print STDERR "$ERROR: \$KorData needs to be an array reference for DetermineBalancePoints!\n\n";
		exit 185;
	}
	
	for $GroupCount (0 .. $#{$KorData}) { # loop over all groups
		
		for $Transition (0 .. $#{$KorData->[$GroupCount]{Transitions}} ) { # loop over all transitions
			my @Atoms; # actually the array contains the balance points of the monopoles
			
			for $Monopole (0 .. $#{$KorData->[$GroupCount]{Transitions}[$Transition]{Monopoles} } ) { # loop over all monopoles 
				$CurMono = $KorData->[$GroupCount]{Transitions}[$Transition]{Monopoles}[$Monopole];
				$MonopoleAdded = 0;
				
				if (not @Atoms) { # if it is the first monopole
					push @Atoms, [$CurMono]; # add the monopole to the 0th array in the array @Atoms
				}
				else { # if monopoles have already been added to @Atoms
					
					foreach $BalancePoint (@Atoms) { # $BalancePoint is actually a reference to an array!
						foreach $OldMono ( @{$BalancePoint} ) { # $OldMono is a monopole that has been already saved in the Atoms array
							if (&Distance ($CurMono, $OldMono) < 0.5) {
								push @{$BalancePoint}, $CurMono;
								$MonopoleAdded = 1;
								last;
							}
						}
						if ($MonopoleAdded) { last }
					}
					
					if (not $MonopoleAdded) { push @Atoms, [$CurMono] } # if the monopole has not been added yet, push it into a new item
				}
			} # of for $Monopole (0 .. $#{$KorData[$GroupCount]{Transitions}[$Transition]{Monopoles} } ) { # loop over all monopoles 
			
			foreach $BalancePoint ( @Atoms ) {
				my %ThisBalancePoint = &DetermineBalancePoint ($KorData->[$GroupCount]{Atoms}, $GroupCount, $BalancePoint);
				push @{$KorData->[$GroupCount]{Transitions}[$Transition]{BalancePoints}}, \%ThisBalancePoint;
			}
			
		} # of for $Transition (0 .. $#{$KorData->[$GroupCount]{Transitions}} ) { # loop over all transitions
		
	} # of for $GroupCount (0 .. $#{$KorData})
	
	return 1;
} # of sub DetermineBalancePoints


sub FindHydrogen { # determines the coordinates of the hydrogen atom
	my $KorData = shift;
	my ($GroupCount, $BalancePoint, $Atom, $Distance, $ShortestDistance, $Number);
	
	if ($KorData !~ m/^ARRAY\(0x[\dA-zA-Z]+\)/) {
		print STDERR "$ERROR: \$KorData needs to be an array reference for FindHydrogen!\n\n";
		exit 190;
	}
	
	$Number = $#{$KorData->[0]{Transitions}[0]{Monopoles}};
	if ($Number < 18) { return } # if less than 20 monopoles are present, it's the woody set, without hydrogen monopoles
	
	for $GroupCount (0 .. $#{$KorData}) { # loop over all groups
		foreach $BalancePoint ( @{$KorData->[$GroupCount]{Transitions}[0]{BalancePoints}} ) {
				
				$ShortestDistance = 1E20; # something big
				
				# Now check the distance to the atoms of the group. Currently there are only C,N and O in the list, so if
				# the shortest distance of the point is larger than 0.2 Angstroms, it's the balance point at the hydrogen	
				foreach $Atom ( @{$KorData->[$GroupCount]{Atoms}} ) {
					$Distance = &Distance ($BalancePoint, $Atom);
					if ($Distance < $ShortestDistance) {
							$ShortestDistance = $Distance;
					}
				}
				
				# if the shortest Distance to the next atom is larger than the threshold, the balance point belongs to the H (which is not included in the .kor file)
				if ($ShortestDistance > 0.2) {
					$BalancePoint->{Atom} = "H";
					push @{$KorData->[$GroupCount]{Atoms}}, $BalancePoint;
					last;
				}
				
			} # of for $Monopole (0 .. $#{$KorData[$GroupCount]{Transitions}[$Transition]{Monopoles} } ) { # loop over all monopoles 
		
	} # of for $GroupCount (0 .. $#KorData)

	return 1;
} # of sub FindHydrogen


sub DetermineDipoleCenters { # determines the center of the VMU vectors (the electrical dipole moments)
	my $KorData = shift;
	my ($GroupCount, $Transition);
	
	if ($KorData !~ m/^ARRAY\(0x[\dA-zA-Z]+\)/) {
		print STDERR "$ERROR: \$KorData needs to be an array reference for DetermineDipoleCenters!\n\n";
		exit 200;
	}
	
	for $GroupCount (0 .. $#{$KorData}) { # loop over all groups
		for $Transition (0 .. $#{$KorData->[$GroupCount]->{Transitions}} ) { # loop over all transitions
			my $CoordVMU = &GetXYZ ($KorData->[$GroupCount]{Transitions}[$Transition]{VMU});
			&MoveOrigin ($CoordVMU, $KorData->[$GroupCount]{XS});
			
			# The center of the dipole is at the half length of the vector, the starting point is defined in XS.
			
			$CoordVMU->{x} = ( $KorData->[$GroupCount]{XS}{x} + $CoordVMU->{x} ) / 2;
			$CoordVMU->{y} = ( $KorData->[$GroupCount]{XS}{y} + $CoordVMU->{y} ) / 2;
			$CoordVMU->{z} = ( $KorData->[$GroupCount]{XS}{z} + $CoordVMU->{z} ) / 2;
			
			$KorData->[$GroupCount]{Transitions}[$Transition]{VMUCoD} = $CoordVMU;
			
			my $CoordVM	= &GetXYZ ($KorData->[$GroupCount]{Transitions}[$Transition]{VM});
			&MoveOrigin ($CoordVM, $KorData->[$GroupCount]{XS});
			
			$CoordVM->{x} = ( $KorData->[$GroupCount]{XS}{x} + $CoordVM->{x} ) / 2;
			$CoordVM->{y} = ( $KorData->[$GroupCount]{XS}{y} + $CoordVM->{y} ) / 2;
			$CoordVM->{z} = ( $KorData->[$GroupCount]{XS}{z} + $CoordVM->{z} ) / 2;
			
			$KorData->[$GroupCount]{Transitions}[$Transition]{VMCoD} = $CoordVM;
		}
	}
	
	return 1
} # of sub DetermineDipoleCenters


sub AssociateMonopoles { # associates each monopole with an atom 
	my $KorData = shift;
	my ($GroupCount, $Transition, $Monopole, $CurMono, $Atom, $Distance, $ShortestDistance, $AtomLabel, $AtomCoord, $Point, $PointCoord, $Number);
	
	if ($KorData !~ m/^ARRAY\(0x[\dA-zA-Z]+\)/) {
		print STDERR "$ERROR: \$KorData needs to be an array reference for AssociateMonopoles!\n\n";
		exit 210;
	}
	
	for $GroupCount (0 .. $#{$KorData}) { # loop over all groups
		
		for $Transition (0 .. $#{$KorData->[$GroupCount]{Transitions}} ) { # loop over all transitions
			# determine the number of monopoles
			$Number = $#{$KorData->[$GroupCount]{Transitions}[$Transition]{Monopoles}};
			
			for $Monopole (0 .. $#{$KorData->[$GroupCount]{Transitions}[$Transition]{Monopoles} } ) { # loop over all monopoles 
				$CurMono = $KorData->[$GroupCount]{Transitions}[$Transition]{Monopoles}[$Monopole];
				$ShortestDistance = 1E20; # something big
				
				# now check which atom is next to the monopole
				foreach $Atom ( @{$KorData->[$GroupCount]{Atoms}} ) {
					$Distance = &Distance ($CurMono, $Atom);
				
					if ($Distance < $ShortestDistance) {
							$ShortestDistance = $Distance;
							if (not $Atom->{Atom}) {
								print "Atom label of group $GroupCount not defined!\n";
								next;
							}
							
							$AtomLabel = $Atom->{Atom} . $GroupCount;
							$AtomCoord = $Atom;
					}
				}
				
				# if the shortest Distance to the next atom is larger than the threshold, the monopole belongs to the H (which is not included in the .kor file)
#				if ($ShortestDistance > 0.5) {
#					$AtomLabel = "H$GroupCount";
#				}
				
				# now check which balance point is next to the monopole
				if ($Number > 18) { # if it's the NMA4FIT2 set, the Woody one doesn't really work with the balance point thing
					foreach $Point ( @{$KorData->[$GroupCount]{Transitions}[$Transition]{BalancePoints}} ) {
						$Distance = &Distance ($CurMono, $Point);
					
						if ($Distance < $ShortestDistance) {
							$ShortestDistance = $Distance;
							$PointCoord = $Point;
						}
					}
				}
				
				# now determine the relative position of the the monopole to the atom, using the balance point calculated before
				# the use of the balance points comes in handy to cope with the lack of the coordinates for the hydrogen
				if ($PointCoord) {
					if ($ShortestDistance < 0.01) { # then it is assumed to be equal to an atom
						$CurMono->{xPos} = "=";
						$CurMono->{yPos} = "=";
						$CurMono->{zPos} = "=";
					}
					else { # if the distance is bigger than 0.05
						if    ( $CurMono->{x} < $PointCoord->{x} ) { $CurMono->{xPos} = "-" }
						elsif ( $CurMono->{x} > $PointCoord->{x} ) { $CurMono->{xPos} = "+" }
#						                                      else { $CurMono->{xPos} = "=" }
						
						
						if    ( $CurMono->{y} < $PointCoord->{y} ) { $CurMono->{yPos} = "-" }
						elsif ( $CurMono->{y} > $PointCoord->{y} ) { $CurMono->{yPos} = "+" }
#						                                      else { $CurMono->{yPos} = "=" }
						
						
						if    ( $CurMono->{z} < $PointCoord->{z} ) { $CurMono->{zPos} = "-" }
						elsif ( $CurMono->{z} > $PointCoord->{z} ) { $CurMono->{zPos} = "+" }
#						                                      else { $CurMono->{zPos} = "=" }
					}
				}
				
				$CurMono->{Distance} = sprintf ("%.6f", $ShortestDistance); # add the distance to the next atom
				$CurMono->{Atom} = $AtomLabel;
			} # of for $Monopole (0 .. $#{$KorData[$GroupCount]{Transitions}[$Transition]{Monopoles} } ) { # loop over all monopoles 
		} # of for $Transition (0 .. $#{$KorData[$GroupCount]{Transitions}} ) { # loop over all transitions
		
	} # of for $GroupCount (0 .. $#KorData)

	return 1
} # of sub AssociateMonopoles


sub GetXYZ { # sets the global xyz variables
	my (@Fields, $Hash, $Type);
	my $Thingy  = shift;
	my $Extrema = shift;
		
	if ($Thingy =~ m/^HASH\(0x[\dA-zA-Z]+\)/) { # if the variable is a reference to a hash
		$Hash = { x => $Thingy->{x},
		          y => $Thingy->{y},
		          z => $Thingy->{z},
		        };
	}
	else { # OK, then it's a string...
		$Type   = substr ($Thingy, 0, 11); # save the first 12 characters
		$Thingy = substr ($Thingy, 12);    # cut away the first 12 characters
		
		$Type   =~ s/^\s+|\s+$//g;
		$Thingy =~ s/^\s+|\s+$//g;
		$Thingy =~ s/\s+/ /g;
		
		@Fields = split /\s+/, $Thingy;
		
		$Hash = { x => $Fields[0],
		          y => $Fields[1],
		          z => $Fields[2],
		        };
		
		if ($Type =~ m/MONO/) { $Hash->{q}    = $Fields[3] }
		if ($Type =~ m/ATOM/) { $Hash->{Atom} = $Fields[3] }
	}
	
	if ($Extrema) {
		if ($Hash->{x} < $Extrema->{MinX}) { $Extrema->{MinX} = $Hash->{x} }
		if ($Hash->{x} > $Extrema->{MaxX}) { $Extrema->{MaxX} = $Hash->{x} }
		
		if ($Hash->{y} < $Extrema->{MinY}) { $Extrema->{MinY} = $Hash->{y} }
		if ($Hash->{y} > $Extrema->{MaxY}) { $Extrema->{MaxY} = $Hash->{y} }
		
		if ($Hash->{z} < $Extrema->{MinZ}) { $Extrema->{MinZ} = $Hash->{z} }
		if ($Hash->{z} > $Extrema->{MaxZ}) { $Extrema->{MaxZ} = $Hash->{z} }
	}
	
	return $Hash;
} # of sub GetXYZ


sub MoveOrigin { # moves the first vector about the coordinates of the second one 
	my $Vector = shift;
	my $Origin = shift;
	
	$Vector->{x} = $Vector->{x} + $Origin->{x};
	$Vector->{y} = $Vector->{y} + $Origin->{y};
	$Vector->{z} = $Vector->{z} + $Origin->{z};
} # of sub MoveOrigin


sub DetermineBalancePoint { # calculates the balance point from a given set of points
	my $Atoms = shift;
	my $Group = shift;
	my $Points  = shift;
	my ($Point, $AtomLabel, $Atom, $Distance);
	my $ShortestDistance = 1E20;
	my %BalancePoint = ( x => 0, y => 0, z => 0, q => 0 );
	
	foreach $Point ( @{$Points} ) {
		%BalancePoint = (
			x => $BalancePoint{x} + $Point->{x},
			y => $BalancePoint{y} + $Point->{y},
			z => $BalancePoint{z} + $Point->{z},
			q => $BalancePoint{q} + $Point->{q},
		);
	}
	
	%BalancePoint = (
		x => $BalancePoint{x} / scalar @{$Points},
		y => $BalancePoint{y} / scalar @{$Points},
		z => $BalancePoint{z} / scalar @{$Points},
		q => $BalancePoint{q},
	);
	
	# now check which atom is next to the monopole
	foreach $Atom ( @{$Atoms} ) {
		$Distance = &Distance (\%BalancePoint, $Atom);
		
		if ($Distance < $ShortestDistance) {
			$ShortestDistance = $Distance;
			$AtomLabel = $Atom->{Atom};
		}
	}
	
	if (not $AtomLabel) {
		print "$AtomLabel of group $Group not defined!\n";
	}
	else {
		$BalancePoint{Atom} = $AtomLabel . $Group;
	}
	
	return %BalancePoint;
} # of sub DetermineBalancePoint


sub StripLine { # prepares a line for split
	my $Line = shift;
	# remove all leading and trailing tabs and blanks, carriage returns and line feeds
	${$Line} =~ s/^[\t\s]+|[\t\s\r\n]+$//g;
} # of sub StripLine


sub CleanExtrema { # sets the determined extrema to integers
	my $Extrema = shift;
	my ($xRange, $yRange, $zRange);
	
	$xRange = $Extrema->{MaxX} - $Extrema->{MinX};
	$yRange = $Extrema->{MaxY} - $Extrema->{MinY};
	$zRange = $Extrema->{MaxZ} - $Extrema->{MinZ};
	
	# decrease the minimum about 10% of the range
	$Extrema->{MinX} = $Extrema->{MinX} - ($xRange * 0.1);
	$Extrema->{MinY} = $Extrema->{MinY} - ($yRange * 0.1);
	$Extrema->{MinZ} = $Extrema->{MinZ} - ($zRange * 0.1);
	
	# increase the maximum about 10% of the range
	$Extrema->{MaxX} = $Extrema->{MaxX} + ($xRange * 0.1);
	$Extrema->{MaxY} = $Extrema->{MaxY} + ($yRange * 0.1);
	$Extrema->{MaxZ} = $Extrema->{MaxZ} + ($zRange * 0.1);
	
	# remove the floating point part
	foreach (keys %{$Extrema}) { $Extrema->{$_} = int $Extrema->{$_} }
} # of sub CleanExtrema

1;
