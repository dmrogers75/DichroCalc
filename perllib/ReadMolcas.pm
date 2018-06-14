package ReadMolcas;

####################################################################################################
#
# Program:    ReadMolcas
#
# Function:   Provides routines to read in molcas in- and output files
#
# Author:     Benjamin Bulheller
#
# Version:    $Revision: 3774 $, $Date: 2009-01-18 23:44:44 +0000 (Sun, 18 Jan 2009) $
#
# Date:       January 2007
#
# Usage:      ReadSewardInp  ($HashRef, "inputfile");
#             ReadSewardOut  ($HashRef, "outfile");
#             WriteSewardInp ("ouputfile");
#             ReadScfOut     ($HashRef);
#             ReadRasReadOut ($HashRef);
#             ReadCaspt2Out  ($HashRef);
#             ReadRassiOut   ($HashRef);
#             ReadRasScfOut  ($HashRed);
#
# Data Structures:
#
# $Seward
#    |
#    |-- Atoms
#    |     |
#    |     |- Angstrom
#    |     |    |-- 0
#    |     |    |   |- Label
#    |     |    |   |- x
#    |     |    |   |- y
#    |     |    |   `- z
#    |     |   ...
#    |     |
#    |     `- Bohr
#    |          |-- 0
#    |          |   |- Label
#    |          |   |- x
#    |          |   |- y
#    |          |   `- z
#    |         ...
#    |
#    |-- NuclearPotentialEnergy
#    |
#    |
#    `-- BasisFunctions
#          |
#          |-- 0
#          |   |- Basis
#          |   |- BasisFunctions
#          |   |    |-- 0
#          |   |   ...
#          |   |
#          |   |- Center
#          |   |- IrredRep
#          |   |- Label
#          |   |- Phase
#          |   `- Type
#         ...
#
#
# $Rassi
#   |
#   |-- Energies
#   |      |
#   |      |- 0 = Energy root 1
#   |      `- 1 = Energy root 2
#   |
#   |-- MLTPL = Electric dipole moments
#   |      |
#   |      |-- 1 = root 1
#   |      |   |
#   |      |   |- 1 => {x, y, z}
#   |      |   `- 2 => {x, y, z}
#   |      |
#   |      |-- 2 = root 2
#   |      |   |
#   |      |   |- 1 => {x, y, z}
#   |      |   `- 2 => {x, y, z}
#   |      |
#   |      `- Origin => {x, y, z}
#   |
#   `-- VELOCITY
#          |
#          |-- 1 = root 1
#          |   |
#          |   |- 1 => {x, y, z}
#          |   `- 2 => {x, y, z}
#          |
#          |-- 2 = root 2
#          |   |
#          |   |- 1 => {x, y, z}
#          |   `- 2 => {x, y, z}
#          |
#          `- Origin => {x, y, z}
#
#
####################################################################################################

use strict;                         # alway use this!!!
use FindBin qw/$Bin/;               # sets $Bin the directory of the script
use lib "$ENV{HOME}/bin/perllib";   # add ~/bin/perllib to the library path
use lib $Bin;                       # add the script's directory to the library path
use Data::Dumper;	                  # for easy printout of arrays and hashes
use DebugPrint;                     # handy during debugging

$Data::Dumper::Sortkeys = 1; # sort the hash keys

require Exporter;
our @ISA    = qw( Exporter );
our @EXPORT = qw( ReadSewardInp ReadSewardOut WriteSewardInp
                  ReadScfOut ReadRasReadOut ReadCaspt2Out ReadRassiOut ReadRasScfOut );

our $ERROR  = "\nERROR (ReadMolcas.pm) ";

#######################################################################################################################


sub ReadSewardInp {
	my $Seward = shift;
	my $InFile = shift;
	my (@Content, $Line, @Fields, @Atoms, @Charges);
	
	if (not -f "$InFile") {
		print STDERR "$ERROR in ReadSewardInp: The file $InFile has not been found!\n";
		return 0;
	}
	
	if ($Seward !~ m/^HASH\(0x[\dA-zA-Z]+\)/) {
		print STDERR "$ERROR: A hash reference is required for ReadSewardInp!\n\n";
		exit 125;
	}
	
	open FILE, "<$InFile";
	@Content = <FILE>;
	close FILE;
	
	while (@Content) {
		$Line = shift @Content;
		
		if ($Line =~ m/TITLE/) {
			$Line = shift @Content;
			chomp $Line;
			$Seward->{Title} = $Line;
		}
		elsif ($Line =~ m/MULTIPOLES/) {
			$Line = shift @Content;
			chomp $Line;
			$Seward->{Multipoles} = $Line;
		}
		elsif ($Line =~ m/^BASIS SET/) {
			my $Set = {};  # create a new hash for the current basis set
			my @CurAtoms;
			
			$Line = shift @Content;   # the basis set definition
			$Set->{Basis} = $Line;
			
			while ($Content[0] !~ m/^END OF BASIS/) {
				$Line = shift @Content;
				@Fields = &SplitLine ($Line);
				
				push @CurAtoms, {Label => $Fields[0], x => $Fields[1], y => $Fields[2], z => $Fields[3],
				                 Unit => $Fields[4], Line => $Line};
				
				push @Atoms, {Label => $Fields[0], x => $Fields[1], y => $Fields[2], z => $Fields[3],
				              Unit => $Fields[4], Line => $Line};
			}
			
			$Set->{Atoms} = \@CurAtoms;
			push @{$Seward->{BasisSets}}, $Set;
			
			# read the END OF INPUT LINE
			$Line = shift @Content;
		}
		elsif ($Line =~ m/^XFIELD/) {
			push @{$Seward->{Body}}, $Line;
			$Line = shift @Content;
			push @{$Seward->{Body}}, $Line;
			
			while ($Content[0] !~ m/^END OF INPUT/) {
				$Line = shift @Content;
				@Fields = &SplitLine ($Line);
				push @Charges, {x => $Fields[0], y => $Fields[1], z => $Fields[2], q => $Fields[3],
				                Unit => $Fields[4], Line => $Line};
				
				# add the reference to the current charge hash to the body for the use with WriteSewardInp
				push@{$Seward->{Body}}, $Charges[$#Charges];
			}
			
			# read the END OF INPUT LINE
			$Line = shift @Content;
			push @{$Seward->{Body}}, $Line;
		}
		else {
			# if it is not a BASIS SET or an XFIELD
#			push @{$Seward->{Body}}, $Line;
		}
	}
	
	$Seward->{Atoms}   = \@Atoms;
	$Seward->{Charges} = \@Charges;
	
	return 1;
} # of sub ReadSewardInp


sub ReadSewardOut {
	my $Seward = shift;
	my $InFile = shift;
	my (@Content, $Line, @Fields, @Atoms, @Chargesf, $IrredRep, $BasisFunctionsIrrep);
	
	if (not -f "$InFile") {
		print STDERR "$ERROR in ReadSewardOut: The file $InFile has not been found!\n";
		return 0;
	}
	
	if ($Seward !~ m/^HASH\(0x[\dA-zA-Z]+\)/) {
		print STDERR "$ERROR: A hash reference is required for ReadSewardInp!\n\n";
		exit 125;
	}
	
	open FILE, "<$InFile";
	@Content = <FILE>;
	close FILE;
	
	while ( @Content ) {
		$Line = shift @Content;
		
		if ($Line =~ m!Cartesian Coordinates / Bohr, Angstrom!) {
			until (not @Content or $Content[0] =~ m/Center  Label/) { shift @Content }
			$Line = shift @Content; # shift away the "Center  Label" line
			
			until ($Content[0] =~ m/^(\s+)?$/) { # until and empty line is found
				$Line = shift @Content;
				@Fields = &SplitLine ($Line);
				push @{$Seward->{Atoms}{Bohr}}, { Label => $Fields[1], x => $Fields[2], y => $Fields[3], z => $Fields[4] };
				push @{$Seward->{Atoms}{Angstrom}}, { Label => $Fields[1], x => $Fields[5], y => $Fields[6], z => $Fields[7] };
			}
		}
		elsif ($Line =~ m/Symmetry adapted Basis Functions/ or
		       $Line =~ m/Petite list Basis Functions/) {
			until ($Content[0] =~ m/Nuclear Potential Energy/) {
				$Line = shift @Content;
				
				if ($Line =~ m/Basis Label/) {
					until ($Content[0] =~ m/^(\s+)?$/) { # until an empty line is found
						$Line = shift @Content;
						@Fields = &SplitLine ($Line);
						
						# if there are two different centers it's actually two different atoms
						# related by symmetry
						if ($Fields[3]) {
							push @{$Seward->{BasisFunctions}},
								{ Basis    => $Fields[0], Label  => $Fields[1],
									Type     => $Fields[2], Center => $Fields[3], Phase => $Fields[4],
									IrredRep => $IrredRep,  BasisFunctions => $BasisFunctionsIrrep,
								};
						}
						
						if ($Fields[5]) {
							push @{$Seward->{BasisFunctions}},
								{ Basis    => $Fields[0], Label  => $Fields[1],
									Type     => $Fields[2], Center => $Fields[5], Phase => $Fields[6],
									IrredRep => $IrredRep,  BasisFunctions => $BasisFunctionsIrrep,
								};
						}
						
						if ($Fields[7]) {
							push @{$Seward->{BasisFunctions}},
								{ Basis    => $Fields[0], Label  => $Fields[1],
									Type     => $Fields[2], Center => $Fields[7], Phase => $Fields[8],
									IrredRep => $IrredRep,  BasisFunctions => $BasisFunctionsIrrep,
								};
						}
						
						if ($Fields[9]) {
							push @{$Seward->{BasisFunctions}},
								{ Basis    => $Fields[0], Label  => $Fields[1],
									Type     => $Fields[2], Center => $Fields[9], Phase => $Fields[10],
									IrredRep => $IrredRep,  BasisFunctions => $BasisFunctionsIrrep,
								};
						}
						
						if ($Fields[11]) {
							print STDERR "\n$ERROR in ReadSewardOut:\n";
							print STDERR "In the list with the basis functions of the atoms there are more columns\n";
							print STDERR "than currently supported. This was triggered by the line\n$Line\n";
							exit 126;
						}
					}
				}
				elsif ($Line =~ m/Irreducible representation/) {
					$Line = substr $Line, index ($Line, ":") + 1;  # copy everything from the char after the colon
					$Line =~ s/\s+//g;  # remove all blanks
					$IrredRep = $Line;
				}
				elsif ($Line =~ m/Basis function\(s\) of irrep/) {
					$Line = substr $Line, index ($Line, ":") + 1;  # copy everything from the char after the colon
					$Line =~ s/\s+//g;  # remove all blanks
					$BasisFunctionsIrrep = [split /,/, $Line];
				}
			}
		}
		elsif ($Line =~ m/Nuclear Potential Energy/) {
			chomp $Line;
			# cut out the number between "Nuclear Potential Energy" and "au"
			$Line =~ s/^(\s+)?Nuclear Potential Energy(\s+)([\d\.-]+)(\s+)au/$3/;
			$Seward->{NuclearPotentialEnergy} = $Line;
		}
	}
}  # of sub ReadSewardOut


sub WriteSewardInp {
	my $Seward = shift;
	my ($ChargesAdded, $Charge, $Line);
	
	print STDERR "\n$ERROR: The function ReadSewardInp had to be significantly changed to\n";
	print STDERR "adapt it to the format changes in MOLCAS 7. WriteSewardInp is not yet\n";
	print STDERR "adapted to those alterations in the hash structure and not functional.\n\n";
	exit 100;
	
	if (not $Seward->{OutFile}) {
		print STDERR "\nERROR: No outfile given!\n";
		return 0;
	}
	
	if (defined $Seward->{Charges}) { $ChargesAdded = 0 }
	
	open FILE, ">$Seward->{OutFile}";
	my $Content = $Seward->{Body};
	
	foreach $Line (@{$Seward->{Body}}) {
		# if it is a hash reference, take the values from it
		if ($Line =~ m/^HASH\(0x[^)]+\)/) {
			printf FILE "%-6s%10.6f   %10.6f   %10.6f %s\n",
			$Line->{Label},
			$Line->{x}, $Line->{y}, $Line->{z},
			$Line->{Unit};
		}
		else {
			if ($Line =~ m/END OF INPUT/ and defined $Seward->{Charges} and not $ChargesAdded) {
				print FILE "XFIELD\n";
				print FILE scalar @{$Seward->{Charges}}, " ANGSTROM\n";
				
				foreach $Charge ( @{$Seward->{Charges}} ) {
					printf FILE "%12.6f  %12.6f  %12.6f  %12.6f  0 0 0\n",
					       $Charge->{x}, $Charge->{y}, $Charge->{z}, $Charge->{q};
				}
				
				print FILE "END OF INPUT\n";
			}
			else {
				print FILE $Line;
			}
		}
	}
	close FILE;
	
	return 1;
} # of sub WriteSewardInp


sub ReadScfOut {
	my $OutFile = "scf.out";
	my $Data = shift;
	
	my ($Line, $Value);
	
	if (not -f $OutFile) {
		print STDERR "\nERROR: $OutFile not found!\n";
		exit 2;
	}
	
	&CheckLanding ($OutFile);
	
	open FILE, "<$OutFile";
	while ($Line = <FILE>) {
		if ($Line =~ m/Total SCF energy/) {
			$Value = substr $Line, 34;
			$Value =~ s/\s//g;
			&RoundValue (\$Value, $Data->{RoundNumbers});
		}
	}
	close FILE;
	
	&CheckValue ($OutFile, $Value);
	
	$Data->{ScfEnergy} = $Value;
	
	return 1;
} # of sub ReadScfOut


sub ReadRasReadOut {
	my $OutFile = "rasread.out";
	my $Data = shift;
	
	my ($Line, $Value, @Fields);
	
	if (not -f $OutFile) {
		print STDERR "\nERROR: $OutFile not found!\n";
		exit 2;
	}
	
	&CheckLanding ($OutFile);
	
	open FILE, "<$OutFile";
	while ($Line = <FILE>) {
		if ($Line =~ m/Final energy of root number/) {
			@Fields = &SplitLine ($Line);
			
			&RoundValue (\$Fields[7], $Data->{RoundNumbers});
			$Data->{RasRead}{$Fields[5]} = $Fields[7];
			$Value = 1;
		}
	}
	close FILE;
	
	&CheckValue ($OutFile, $Value);
	
	return 1;
} # of sub ReadRasReadOut


sub ReadRassiOut {
	my $OutFile = "rassi.out";
	my $Data = shift;
	
	my ($Line, $Value, @Fields, $Origin, $States);
	
	if (not -f $OutFile) {
		print STDERR "\nERROR: $OutFile not found!\n";
		exit 2;
	}
	
	&CheckLanding ($OutFile);
	
	open FILE, "<$OutFile" or die "\n$ERROR: File $OutFile could not be opened: $!";
	while ($Line = <FILE>) {
		if ($Line =~ m/SPIN-FREE ENERGIES/) {
			# skip to the first line with an energy
			do { $Line = <FILE> } until ( $Line =~ m/^(\s+)?1\s+/ );
			
			do {
				@Fields = &SplitLine ($Line);
				
				&RoundValue (\$Fields[1], $Data->{RoundNumbers});
				$Data->{Rassi}{Energies}{$Fields[0]}{Total}      = $Fields[1];
				$Data->{Rassi}{Energies}{$Fields[0]}{eV}         = $Fields[2];
				$Data->{Rassi}{Energies}{$Fields[0]}{WaveNumber} = $Fields[3];
				$Value = 1;
				$Line = <FILE>;
			} until ( $Line !~ m/^(\s+)?\d+\s+/ );
		}
		
		if ($Line =~ m/MATRIX ELEMENTS OF 1-ELECTRON OPERATORS/) {
			# this section heading changed from MOLCAS 5.2 to 7.0, to maintain the
			# compartibility to both versions the check is split in two comparisons
			$Line = <FILE>;
			
			# it is not the spin-free section, continue with sifting through the file
			if ($Line !~ m/SPIN-FREE EIGENSTATES/) { next }
			
			# now search for MLTPL or VELOCITY blocks until ****** is found which end that section
			do {
				$Line = <FILE>;
				
				if ($Line =~ m/PROPERTY:(\s+)MLTPL(\s+).(\s+)COMPONENT:(\s+)\d/) {
					# the asterisk is needed to pass a file handle to a sub routine
					&ReadMoments ($Data, $Line, *FILE);
				}
				
				if ($Line =~ m/PROPERTY:(\s+)VELOCITY(\s+).(\s+)COMPONENT:(\s+)\d/) {
					# the asterisk is needed to pass a file handle to a sub routine
					&ReadMoments ($Data, $Line, *FILE);
				}
				
				if ($Line =~ m/PROPERTY:(\s+)ANGMOM(\s+).(\s+)COMPONENT:(\s+)\d/) {
					# the asterisk is needed to pass a file handle to a sub routine
					&ReadMoments ($Data, $Line, *FILE);
				}
				
			} until ($Line =~ m/\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*/);
			
		}
		
	}
	close FILE;
	
	return 1;
} # of sub ReadRassiOut


sub ReadMoments {
	my $Data   = shift;
	my $Line   = shift;
	my $Handle = shift;
	
	my ($Type, @Fields, @States, $State, $Origin, $Component, $Coord, $i);
	
	@Fields = &SplitLine ($Line);
	$Type = $Fields[1];
	
	   if ($Type eq "MLTPL"   ) { $Component = $Fields[4] }
	elsif ($Type eq "VELOCITY") { $Component = $Fields[3] }
	elsif ($Type eq "ANGMOM"  ) { $Component = $Fields[3] }
	else {
		print STDERR "$ERROR (ReadMoments): Could not interpret type \"$Type\".\n";
		print STDERR "Only MLTPL, ANGMOM or VELOCITY components are currently supported.\n\n";
		exit 112;
	}
	
	# read the origin
	$Line = <$Handle>;
	@Fields = &SplitLine ($Line);
	$Origin = { x => $Fields[2], y => $Fields[3], z => $Fields[4] };
	
	   if ($Component eq "1") { $Coord = "x" }
	elsif ($Component eq "2") { $Coord = "y" }
	elsif ($Component eq "3") { $Coord = "z" }
	else {
		print STDERR "$ERROR (ReadMoments): Could not interpret \n$Component\" as a component of a dipole moment.\n";
		print STDERR "(1, 2 or 3 is expected for x, y or z component).\n\n";
		exit 114;
	}
	
	# The origin may be overwritten over and over again for multiple blocks
	# but it's assumed that it is always the same
	$Data->{Rassi}{$Type}{Origin} = $Origin;
	
	# read the state indices
	$Line = <$Handle>;
	@States = &SplitLine ($Line);
	shift @States;  # STATE
	shift @States;  # :
	
	until ($Line =~ m/^(\s+)?1\s+/) { $Line = <$Handle>; }
	
	do {
		@Fields = &SplitLine ($Line);
		
		$i = 0;  # this is the index of the state in the $Data hash
		
		foreach $State ( @States ) { # these are the state number, e.g. 1, 2, 3, 4
			++$i;
			$Data->{Rassi}{$Type}{$State}{$Fields[0]}{$Coord} = $Fields[$i];
		}
		
		$Line = <$Handle>;
	} until ($Line !~ m/^(\s+)?\d\d?\s+/);
	
} # of sub ReadMoments


sub ReadComponent {
	my $Data      = shift;
	my $Component = shift;
	my $Handle    = shift;
	
	my ($Line, @Fields);
	$Line = ""; # to avoid an undefined value in the first pattern match
	
	# read lines until $Line starts with "    1    "
	($Line = <$Handle>) until ($Line =~ m/^\s+1\s+/);
	
	@Fields = &SplitLine ($Line);
	&RoundValue (\$Fields[1], $Data->{RoundNumbers});
	$Data->{1.1}{Matrix}{$Component} = $Fields[1];
	
	$Line = <$Handle>;
	@Fields = &SplitLine ($Line);
	&RoundValue (\$Fields[1], $Data->{RoundNumbers});
	$Data->{1.2}{Matrix}{$Component} = $Fields[1];
	
	$Line = <$Handle>;
	@Fields = &SplitLine ($Line);
	&RoundValue (\$Fields[1], $Data->{RoundNumbers});
	$Data->{2.1}{Matrix}{$Component} = $Fields[1];
	
	return 1;
} # of sub ReadComponent


sub ReadCaspt2Out {
	my $Data = shift;
	
	my ($OutFile, @Content, @Values, @Fields);
	
	$OutFile = "caspt2.1.2.out";
	
	if (not -f $OutFile) {
		print STDERR "\nERROR: $OutFile not found!\n";
		exit 2;
	}
	
	&CheckLanding ($OutFile);
	
	open (FILE, "<$OutFile");
	@Content = <FILE>;
	@Values = grep /E2 \(Variational\):/, @Content;
	
	if (scalar @Values != 2) {
		print STDERR "\n\nERROR: $OutFile did not contain exactly two correlation energies!\n\n";
		print STDERR "Found energies:\n", join ("\n", @Values), "\n\n";
		exit 6;
	}
	
	@Fields = &SplitLine ($Values[0]);
	&RoundValue (\$Fields[2], $Data->{RoundNumbers});   # round the values, if requested
	$Data->{"1.1"}{Caspt2} = $Fields[2];                # save the value
	
	@Fields = &SplitLine ($Values[1]);
	&RoundValue (\$Fields[2], $Data->{RoundNumbers});   # round the values, if requested
	$Data->{"1.2"}{Caspt2} = $Fields[2];                # save the value
	close FILE;
	
	$OutFile = "caspt2.2.1.out";
	
	if (not -f $OutFile) {
		print STDERR "\nERROR: $OutFile not found!\n";
		exit 2;
	}
	
	&CheckLanding ($OutFile);
	
	open (FILE, "<$OutFile");
	@Content = <FILE>;
	@Values = grep /E2 \(Variational\):/, @Content;
	
	if (scalar @Values != 1) {
		print STDERR "\n\nERROR: $OutFile did not contain exactly one correlation energy!\n\n";
		print STDERR "Found energies:\n", join ("\n", @Values), "\n\n";
		exit 7;
	}
	
	@Fields = &SplitLine ($Values[0]);
	&RoundValue (\$Fields[2], $Data->{RoundNumbers});   # round the values, if requested
	$Data->{"2.1"}{Caspt2} = $Fields[2];                # save the value
	close FILE;
	
	return 1;
} # of sub ReadCaspt2Out


sub ReadRasScfOut {
	my $Data = shift;
	my $OutFile = "rasscf.1.2.out";
	
	my (@Content, @Values, @Fields);
	
	if (not -f $OutFile) {
		print STDERR "\nERROR: $OutFile not found!\n";
		exit 2;
	}
	
	open (FILE, "<$OutFile");
	@Content = <FILE>;
	@Values = grep /^\s+root number  [12] E = /, @Content;
	
	if (scalar @Values != 2) {
		print STDERR "\n\nERROR: $OutFile did not contain exactly 2 energies!\n\n";
		print STDERR "Found energies:\n", join ("\n", @Values), "\n\n";
		exit 6;
	}
	
	@Fields = &SplitLine ($Values[0]);
	$Data->{1.1}{RasScf} = $Fields[5];
	
	@Fields = &SplitLine ($Values[1]);
	$Data->{1.2}{RasScf} = $Fields[5];
	close FILE;
	
	
	$OutFile = "rasscf.2.1.out";
	
	if (not -f $OutFile) {
		print STDERR "\nERROR: $OutFile not found!\n";
		exit 2;
	}
	
	open (FILE, "<$OutFile");
	@Content = <FILE>;
	@Values = grep /^\s+root number  [12] E = /, @Content;
	
	if (scalar @Values != 1) {
		print "\n\nERROR: $OutFile did not contain exactly one energy!\n\n";
		print "Found energies:\n", join ("\n", @Values), "\n\n";
		exit 7;
	}
	
	$Values[0] =~ s/^\s+//;  # remove leading blanks
	@Fields = &SplitLine ($Values[0]);
	$Data->{2.1}{RasScf} = $Fields[5];
	close FILE;
	
	return 1;
} # of sub ReadRasScfOut


sub RoundValue {
	my $Value = shift;
	my $RoundNumbers = shift;
	
	# the reference of $Value is passes, thus it needs to be dereferenced to change it
	if ($RoundNumbers) { ${$Value} = sprintf "%.".$RoundNumbers."f", ${$Value}  }
} # of RoundValue


sub CheckLanding {
	my $File = shift;
	my (@Content, @Landing);
	
	open FILE, "<$File" or die "\nERROR: Could not open $File!:$!";
	@Content = <FILE>;
	close FILE;
	
	@Landing = grep /Happy Landing!/i, @Content;
	
	if (@Landing) { return 1 }
	else {
		print "\nERROR: No Happy Landing found in file $File!\n\n";
		exit 10;
	}
} # of sub CheckLanding


sub CheckValue {
	my $OutFile = shift;
	my $Value = shift;
	
	if (not defined $Value) {
		print "\nERROR in $OutFile, required value not found!\n\n";
		exit 3;
	}
} # of sub CheckValue


sub SplitLine {
	my $Line = shift;
	
	$Line =~ s/^\s+//g;  # remove leading blanks
	chomp $Line;
	return split (/\s+/, $Line);
} # of sub SplitLine

1; # return true (mandatory for all packages)
