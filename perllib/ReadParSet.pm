package ReadParSet;

#######################################################################################################################
#
# Package:         ReadParSet
#
# Function:        Reads in a parameter set from a given file into a hash
#
# Author:          Benjamin Bulheller
#
# Research Group:  Prof. Jonathan Hirst
#                  School of Physical Chemistry
#                  University of Nottingham
#
# Funded by:       EPSRC
#
# Version:         $Revision: 3773 $, $Date: 2009-01-18 23:43:12 +0000 (Sun, 18 Jan 2009) $
#
# Date:            March 2006
#
# Usage:           use ReadParSet;
#                  ReadParSet ($InFile, $ParameterSet, \%ParData);
#                  WriteParSet (\%ParData, $OutFile);
#
#                  $InFile and $OutFile contain the full file name (e.g. parsets.944)
#                  $ParSetName contains the name of the requested parameter set (e.g. NMA4FIT2)
#                  \%ParData is the reference to a hash in which the information will be stored
#
#
# %ParData
#   |
#   |- {Header} => [ array with header lines of the parameter set ]
#   |
#   |- {Atoms}
#   |     |
#   |     |- 0 => { x, y, z, q, Comment }
#   |     |
#   |    etc
#   |
#   |- {TransHeader} => [ array with header lines for each transition ]
#   |
#   |- {Trans} (transition moments)
#   |     |
#   |     |- 0 (transition state, e.g. 1 => ... )
#   |     |  |
#   |     |  |- 0 (transition, e.g. 1 => 2 )
#   |     |  |  |
#   |     |  |  |- Energy
#   |     |  |  |- Header => [ array with header lines ]
#   |     |  |  |- VM     => { x, y, z, q, Comment }
#   |     |  |  |- VMU    => { x, y, z, q, Comment }
#   |     |  |  |
#   |     |  |  |- Monopoles
#   |     |  |  |    |
#   |     |  |  |    |- 0 => { x, y, z, q }
#   |     |  |  |    |
#   |     |  |  |   etc
#   |     |  |  |
#   |     |  | etc
#   |     |  |
#   |     | etc
#   |     |
#   |    etc
#   |
#   `- {Perm}
#         |
#         |- 0 (there is only transition state 0 in the Permanent Moments, to be "compatible" to the {Trans} Hash
#         |  |
#         |  `- 0 (transition, e.g. 1 => 2 )
#         |     |
#         |     |- VMU  => { x, y, z, q, Comment }
#         |     |
#         |     |- Monopoles
#         |     |    |
#         |     |    |- 0 => { x, y, z, q }
#         |     |    |
#         |     |   etc
#         |     |
#         |    etc
#         |
#         |
#        etc
#
#
#######################################################################################################################


use strict;                           # always use this!!!
use FindBin qw/$Bin/;                 # sets $Bin to the directory of the script
use lib $Bin;                         # adds the script's directory to the library path
use lib "$ENV{HOME}/bin/perllib";     # adds ~/bin/perllib to the library path
use Data::Dumper;	                    # for easy printout of arrays and hashes
use DebugPrint;                       # handy during debugging

$Data::Dumper::Sortkeys = 1; # sort the hash keys

require Exporter;
our @ISA    = qw( Exporter );
our @EXPORT = qw( ReadParSet WriteParSet );


#######################################################################################################################

sub ReadParSet {
	my $InFile       = shift;
	my $ParameterSet = shift;
	my $ParData      = shift;
	my ($Line, @Content, @ParContent);
	
	if (not $InFile or not $ParameterSet or not $ParData) {
		print "Error: parameter missing for ReadParSet!\n\n";
		print "Usage: ReadParSet (\$InFile, \$ParSetName, \\\%ParData);\n";
		print "       \$InFile contains the full file name (e.g. parsets.par)\n";
		print "       \$ParSetName contains the name of the requested parameter set (e.g. NMA4FIT2)\n";
		print "       \\\%ParData is the reference to a hash in which the information will be stored\n";
		die;
	}
	
	if (not -f $InFile) {
		print "\n$InFile could not be found!\n\n";
		die;
	}
	
	open FILE, "<$InFile";
	@Content = <FILE>;
	chomp @Content;
	close FILE;
	
	# reject the leading lines until the requested parameter set is found
	until ( (not @Content) or ($Content[0] =~ m/-$ParameterSet-/) ) { shift @Content }
	
	if (not @Content) {
		print "\n-$ParameterSet- has not been found in $InFile!\n\n";
		die;
	}
	
	$Line = shift @Content; # the first line containing -ParameterSet-
	$Line =~ s///g;
	push @ParContent, $Line;
	
	# Push each line from Content into ParContent, until Content is empty or the next parset begins.
	# CAUTION:
	#    - The RegEx contains all characters allowed in the parameter set name and might need
	#      to be amended if new parameter set names use additional character.
	#    - The length of the parameter set name is fixed to 8 characters. This should be safe
	#      as matmac requires this (defined in subroutine FINDTY)
	#    - A problem could be two negative numbers whose minus signs happen to be 8 characters apart.
	until ( (not @Content) or ($Content[0] =~ m/^(\s+)?-[-_A-Z0-9]{8}-(\s+)?/) ) {
		$Line = shift @Content;
		$Line =~ s///g;
		push @ParContent, $Line;
	}
	
	&ParseData ($ParData, $ParameterSet, \@ParContent);
} # of sub ReadParSet


sub ParseData { # parses the content and sorts all information into $ParData
	my $ParData = shift;
	my $ParameterSet = shift;
	my $ParContent = shift;
	my ($i, $Line, $Number, @Fields, @TransHeader, $PermHeader);
	my $StateCount = 0;
	
	while ( @{$ParContent} ) {
		$Line = &NextLine ($ParContent);
		
		if ($Line =~ m/-$ParameterSet-/) {
			my @Header;
			push @Header, $Line;
			
			while ($ParContent->[0] =~ m/^(\s+)?#/) { # if a description of the dataset is given
				$Line = shift @{$ParContent};
				push @Header, $Line;
			}
			
			# read the line with the number of atoms
			$Line = &NextLine ($ParContent);
			push @Header, $Line;
			$ParData->{Header} = \@Header;
			
			$Number = $Line;
			
			# cut away leading blanks if present and everything after the number (a comment usually)
			$Number =~ s/(\s+)?(\d+)[^\d]+$/$2/;
			
			for $i (0 .. $Number-1) {
				@Fields = &NextFields ($ParContent);
				push @{$ParData->{Atoms}}, { x => $Fields[0], y => $Fields[1], z => $Fields[2], q => $Fields[3], Comment => $Fields[4] };
			}
		}
		elsif ($Line =~ m/&TRANSITION/i) {
			$TransHeader[$StateCount] = $Line;
			
			my @Trans = &ReadState ($ParData, $ParContent, 1);
			push @{$ParData->{Trans}[$StateCount]}, @Trans;
			++$StateCount;
		}
		elsif ($Line =~ m/&PERMANENT/i) {
			$PermHeader = $Line;
			
			my @Perm  = &ReadState ($ParData, $ParContent, 0);
			push @{$ParData->{Perm}[0]}, @Perm; # item [0] is taken to be consistent with the {Trans} hash an be able to use the same write procedure later on
		}
		else {
			print "\n\nERROR: Line not recognized:\n\n$Line\n\n";
			die;
		}
		
	} # until (not @{$ParContent});
	
	$ParData->{PermHeader}  = $PermHeader;
	$ParData->{TransHeader} = \@TransHeader;
	
	return 1;
} # of sub ParseData


sub ReadState { # reads in a complete state with all its transitions
	my $ParData    = shift;
	my $ParContent = shift;
	my $IsTrans    = shift;
	my (@AllTrans, $Line, $Number, $Energy, @Fields, $i);
	
	while ( defined $ParContent->[0] and ($ParContent->[0] !~ m/^&/) and @{$ParContent} ) {
		my (@Header, %Trans);
		
		$Line = &NextLine ($ParContent);
		
		$Line =~ s/^\s+//g;  # remove leading blanks
		@Fields = split /\s+/, $Line;
		
		$Number = $Fields[0];
		$Energy = $Fields[1];
		
		push @Header, $Line;
		$Trans{Header} = \@Header;
		
		@Fields = &NextFields ($ParContent);
		$Trans{VMU} = { x => $Fields[0],
		                y => $Fields[1],
		                z => $Fields[2],
		                q => $Fields[3],
		                Comment => $Fields[4],
		              };
		
		if ($IsTrans) { # if it is not a permanent moment, read the magnetic transition dipole moment
			$Trans{Ewn} = $Energy;
			
			if ($Energy != 0) { $Trans{Enm} = 1E7 / $Energy }
			             else { $Trans{Enm} = 0             }
			
			@Fields = &NextFields ($ParContent);
			$Trans{VM}  = { x => $Fields[0],
			                y => $Fields[1],
			                z => $Fields[2],
			                q => $Fields[3],
			              };
		}
		
		for $i (0 .. $Number-1) {
			@Fields = &NextFields ($ParContent);
			my $Monopole = { x => $Fields[0], y => $Fields[1], z => $Fields[2], q => $Fields[3] };
			
			&DetermineAtom ($ParData->{Atoms}, $Monopole);
			
			push @{$Trans{Monopoles}}, $Monopole;
		}
		
		
		if (%Trans and defined $Trans{Header}) {
			push @AllTrans, \%Trans
		}
	}
	
	return @AllTrans;
} # of sub ReadState


sub DetermineAtom { # determines the closest atom to the monopole
	my $Atoms  = shift;
	my $Monopole = shift;
	my ($Atom, $Dist, $ShortestDist, $ClosestAtom);
	
	$ShortestDist = 1000;
	
	for $Atom (0 .. $#{$Atoms}) {
		if (defined $Atoms->[$Atom]{x} and defined $Monopole->{x}) {
			$Dist = sqrt (   ( $Atoms->[$Atom]{x} - $Monopole->{x} )**2
								+ ( $Atoms->[$Atom]{y} - $Monopole->{y} )**2
								+ ( $Atoms->[$Atom]{z} - $Monopole->{z} )**2 );
			
			if (defined $Dist and ($Dist < $ShortestDist) ) {
				$ClosestAtom  = $Atom;
				$ShortestDist = $Dist;
			}
		}
	}
	
	$Monopole->{Atom}     = $ClosestAtom;
	$Monopole->{AtomDist} = sprintf "%.8f", $ShortestDist;
	
	return 1;
} # of sub DetermineAtom


sub NextLine { # returns the next line from @Content which is not a comment line
	my $ParContent = shift;
	my $Line;
	
	do {
		if (not @{$ParContent}) { return undef }
		$Line = shift @{$ParContent};
	} until ($Line !~ m/^(\s+)?#/);
	
	return $Line;
} # of sub NextLine


sub NextFields { # returns the next line from @Content splitted into its fields
	my $ParContent = shift;
	my ($Line, @Fields, $Position, $Comment);
		
	$Line = &NextLine ($ParContent);
	
	if (defined $Line) {
		if ( rindex ($Line, "#") > -1 ) { # if there is a comment at the end of the line
			$Line =~ m/\s+#/; # just match the "#" including all spaces before it
			$Position = rindex ($Line, $&); # $& contains the complete last regex match
			
			$Comment = substr ($Line, $Position, length $Line);
			$Line = substr ($Line, 0, $Position-1); # delete the comment
		}
		
		$Line =~ s/\s+/ /g;
		
		@Fields = split (" ", $Line);
		if ( $Comment ) { push @Fields, $Comment }
		
		return @Fields;
	}
	else {
		return undef;
	}
} # of sub NextFields


sub WriteParSet {
	my $ParData = shift;
	my $OutFile = shift;
	
	my ($State, $i, $j, $k);
	
	if (not $ParData or not $OutFile) {
		print "\nError: \%ParData and \$OutFile required as parameters for WriteParSet!\n\n";
		die;
	}
	
	open (FILE, ">$OutFile") or die "Cannot write to file $OutFile: $!";
	
	print FILE join ("\n", @{$ParData->{Header}}), "\n";
	
	for $i (0 .. $#{$ParData->{Atoms}}) {
		printf FILE "  %11.4f %11.4f %11.4f %11.4f%s\n",
		       $ParData->{Atoms}[$i]{x}, $ParData->{Atoms}[$i]{y}, $ParData->{Atoms}[$i]{z},
		       $ParData->{Atoms}[$i]{q}, $ParData->{Atoms}[$i]{Comment};
	}
	
	for $State (0 .. $#{$ParData->{TransHeader}}) {
		print FILE $ParData->{TransHeader}[$State], "\n";
		&WriteState ($ParData, "Trans", $State);
	}
	
	print FILE $ParData->{PermHeader}, "\n";
	&WriteState ($ParData, "Perm", 0);
	
	close FILE;

	return 1;
} # of sub WriteData


sub WriteState {
	my $ParData = shift;
	my $Type = shift;
	my $State = shift;
	my ($i, $k);
	
	for $i (0 .. $#{$ParData->{$Type}[$State]}) {
	
		if ($Type eq "Trans") { # if it is a transition
			# print the "main" transition header (&TRANSITION x->)
			print FILE join ("\n", @{$ParData->{$Type}[$State][$i]{Header}}), "\n";
			
			# print the electrical moment
			printf FILE "  %11.4f %11.4f %11.4f %7.3f%s\n",
			       $ParData->{$Type}[$State][$i]{VMU}{x},
			       $ParData->{$Type}[$State][$i]{VMU}{y},
			       $ParData->{$Type}[$State][$i]{VMU}{z},
			       $ParData->{$Type}[$State][$i]{VMU}{q},
			       $ParData->{$Type}[$State][$i]{VMU}{Comment};
			
			# print the magnetic moment
			printf FILE "  %11.4f %11.4f %11.4f%s\n",
			       $ParData->{$Type}[$State][$i]{VM}{x},
			       $ParData->{$Type}[$State][$i]{VM}{y},
			       $ParData->{$Type}[$State][$i]{VM}{z},
			       $ParData->{$Type}[$State][$i]{VM}{q};
		}
		else { # if it is a permanent moment
			# print the electrical moment
			print FILE join ("\n", @{$ParData->{$Type}[$State][$i]{Header}}), "\n";
			printf FILE "  %11.4f %11.4f %11.4f%s\n",
			       $ParData->{$Type}[$State][$i]{VMU}{x},
			       $ParData->{$Type}[$State][$i]{VMU}{y},
			       $ParData->{$Type}[$State][$i]{VMU}{z},
			       $ParData->{$Type}[$State][$i]{VMU}{q};
		}
		
		# print the monopoles
		for $k (0 .. $#{$ParData->{$Type}[$State][$i]{Monopoles}}) {
			printf FILE "  %15.8f %15.8f %15.8f %15.8f\n",
			       $ParData->{$Type}[$State][$i]{Monopoles}[$k]{x},
			       $ParData->{$Type}[$State][$i]{Monopoles}[$k]{y},
			       $ParData->{$Type}[$State][$i]{Monopoles}[$k]{z},
			       $ParData->{$Type}[$State][$i]{Monopoles}[$k]{q};
		}
	}
	
	return 1;
} # of sub WriteTransition


1;

