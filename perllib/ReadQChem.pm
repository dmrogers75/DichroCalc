package ReadQChem;

####################################################################################################
#
# Program:    ReadQChem
#
# Function:   To read QChem input and output files
#
# Author:     Benjamin Bulheller
#
# Version:    $Revision: 3183 $, $Date: 2008-08-12 23:01:41 +0100 (Tue, 12 Aug 2008) $
#
# Date:       July 2008
#
# Usage:      ReadQChemInp ($File, $HashRef);
#             ReadQChemOut ($File, $HashRef);
#
####################################################################################################

use strict;			                  # alway use this!!!
use FindBin qw/$Bin/;               # sets $Bin the directory of the script
use lib "$ENV{HOME}/bin/perllib";   # add ~/bin/perllib to the library path
use lib $Bin;                       # add the script's directory to the library path
use Data::Dumper;	                  # for easy printout of arrays and hashes
use DebugPrint;                     # handy during debugging

$Data::Dumper::Sortkeys = 1; # sort the hash keys

require Exporter;
our @ISA    = qw( Exporter );
our @EXPORT = qw( ReadQChemInp ReadQChemOut);

our $ERROR    = "\nERROR (ReadQChem.pm) ";
our $WARNING  = "\nWARNING (ReadQChem.pm) ";

#######################################################################################################################


sub ReadQChemInp {
	my $File = shift;
	my $Data = shift;

	my ($Line, $i, @Fields);
	
	my $Content = &ReadFile ($File, $Data, "inp");

	&ReadRem      ($Data, $Content);
	&ReadMolecule ($Data, $Content);
	&ReadSolute   ($Data, $Content);
} # of sub ReadQChemInp


sub ReadQChemOut {
	my $File = shift;
	my $Data = shift;

	my $Content = &ReadFile ($File, $Data, "out");

	&ReadRem      ($Data, $Content);
	&ReadMolecule ($Data, $Content);
	&ReadSolute   ($Data, $Content);
	&ReadResults  ($Data, $Content);	
} # of sub ReadQChemOut


####################################################################################################
# Internal Routines
####################################################################################################


sub ReadFile {
	my $File = shift;
	my $Data = shift;
	my $Ext  = shift;

	if (not -f $File) {
		if (-f "$File.$Ext") {
			$File = $File . ".$Ext";
		}
		else {
			print STDERR "$ERROR: File $File not found.\n\n";
			exit 130;
		}
	}
	
	if ($Data !~ m/^HASH\(0x[\dA-zA-Z]+\)/) {
		print STDERR "$ERROR: A hash reference is required for ReadQChemInp!\n\n";
		exit 125;
	}
	
	open OUT, "<$File" or die "$ERROR: Could not open file $File: $!";
	my @Content = <OUT>;
	close OUT;

	if (not @Content) {
		print STDERR "$ERROR: File $File appears to be empty.\n\n";
		exit 155;
	}
	
	$Data->{File} = $File;
	$Data->{BaseName} = $File;
	$Data->{BaseName} =~ s/\.$Ext$//;
	
	return \@Content;
}

####################################################################################################

sub ReadRem {
	my $Data    = shift;
	my $Content = shift;
	
	my ($Line, $i, @Fields);
	

	################################################################################
	# Read $rem block
	################################################################################
	
	$i = 0; # start at line 0
	
	until ($i == $#{$Content} or $Content->[$i] =~ m/^\$rem/i) { ++$i }
	
	++$i; # $rem
	
	if ($i > $#{$Content}) {
		print "$WARNING: End of $Data->{File} reached while searching for \$rem.\n\n";
	}
	else {	
		my @rem;

		until ($i == $#{$Content} or $Content->[$i] =~ m/\$end/) {
			$Line = $Content->[$i];
			++$i;
			
			$Line =~ s/^\s+//;      # remove leading blanks
			chomp $Line;            # remove the line feed;
			
			# split into 3 fields using blanks or tabs as delimiter
			@Fields = split /\s+|\t+/, $Line, 3;
			
			my $Item = {
				Keyword => $Fields[0],
				Value   => $Fields[1],
				Comment => $Fields[2],
			};
			
			if (not $Item->{Comment}) { $Item->{Comment} = ""; }
			
			push @rem, $Item;
		}
		
		$Data->{rem} = \@rem;
	}
} # of sub ReadRem

####################################################################################################
	
sub ReadMolecule {
	my $Data    = shift;
	my $Content = shift;
	
	my ($Line, $i, @Fields);

	
	################################################################################
	# Read $molecule block
	################################################################################
	
	$i = 0; # start at line 0
	
	until ($i == $#{$Content} or $Content->[$i] =~ m/^\$molecule/i) { ++$i }
	
	if ($i > $#{$Content}) {
		print "$WARNING: End of $Data->{File} reached while searching for \$molecule.\n\n";
	}
	else {	
		++$i; # $molecule
		++$i; # 0 1

		my @Atoms;
		
		until ($i >= $#{$Content} or $Content->[$i] =~ m/\$end/) {
			$Line = $Content->[$i];
			++$i;
			
			@Fields = split /\s+/, $Line;
			my $Atom = {
				Label => $Fields[0],
				x => $Fields[1],
				y => $Fields[2],
				z => $Fields[3],
			};

			push @Atoms, $Atom;
		}
		
		$Data->{molecule} = \@Atoms;
	}
} # of sub ReadMolecule

####################################################################################################

sub ReadSolute {
	my $Data    = shift;
	my $Content = shift;
	
	my ($Line, $i, @Fields);

	
	################################################################################
	# Read $solute block
	################################################################################
	
	$i = 0; # start at line 0
	
	until ($i == $#{$Content} or $Content->[$i] =~ m/^\$solute/i) { ++$i }
	
	if ($i > $#{$Content}) {
		print "$WARNING: End of $Data->{File} reached while searching for \$solute.\n\n";
	}
	else {	
		++$i; # $solute

		my @Solute;
		
		until ($i >= $#{$Content} or $Content->[$i] =~ m/\$end/) {
			$Line = $Content->[$i];
			++$i;
			
			$Line =~ s/[\s\t]//g;
			chomp $Line;

			push @Solute, $Line;
		}
		
		$Data->{solute} = \@Solute;
	}
} # of sub ReadSolute

####################################################################################################

sub ReadResults { # reads the SCF and DFT results
	my $Data    = shift;
	my $Content = shift;
	
	my ($Line, $i, @Fields);

	
	################################################################################
	# Read SCF and DFT results
	################################################################################
	
	$i = 0; # start at line 0
	
	until ($i >= $#{$Content} or $Content->[$i-1] =~ m/Roots Converged/) {
		$Line = $Content->[$i];
		++$i;
		
		chomp $Line;
		$Line =~ s/^\s+//;

		if ($Line =~ m/Convergence criterion met/) {
			@Fields = split /\s+/, $Line;
			$Data->{SCF}{SCFiter}    = $Fields[0];
			$Data->{SCF}{SCFenergy}  = $Fields[1];
			$Data->{SCF}{SCF1}       = $Fields[2];
			$Data->{SCF}{SCF2}       = $Fields[3];
			# print $Line;
		}
		
		if ($Line =~ m/Requested basis set is/) {
			$Line =~ s/^(\s+)?Requested basis set is //;
			chomp $Line;
			$Data->{BasisSet} = $Line;
		}
		
		if ($Line =~ m/Roots Converged/) {
			@Fields = split /\s+/, $Line;
			$Data->{DFT}{DFTiter}  = $Fields[0];
			$Data->{DFT}{DFT1}     = $Fields[1];
			$Data->{DFT}{DFT2}     = $Fields[2];
			$Data->{DFT}{DFT3}     = $Fields[3];
			$Data->{DFT}{DFT4}     = $Fields[4];
			# print $Line;
		}
	}

	################################################################################
	# Read TDDFT results
	################################################################################
	
	# continue reading the file after "Roots Converged"
	until ($i == $#{$Content} or $Content->[$i] =~ m!TDDFT/TDA Excitation Energies!i) { ++$i }
	until ($i == $#{$Content} or $Content->[$i] =~ m!Excited state!i) { ++$i }
	
	my $Trans = 0;
	my $Value;
	
	if ($i < $#{$Content}) {
		until ($Line =~ m/---------------/) {
			$Line = $Content->[$i];
			++$i;

			if ($Line =~ m/Excited state/) { 
				until ($Line =~ m/^\n$/) {
					if ($Line =~ m/Excited state/) {
						$Value = &GetNumber ($Line);
						$Data->{Trans}[$Trans]{ExEnergy} = $Value;
					}
					elsif ($Line =~ m/Total energy for state/) {
						$Value = &GetNumber ($Line);
						$Data->{Trans}[$Trans]{TotalEnergy} = $Value;
					}
					elsif ($Line =~ m/Trans\. Mom\./) {
						$Value = $Line;
						$Value = substr $Line, index ($Line, ":")+1;
						
						my $Vec = {};
							
						$Vec->{x} = substr $Value, 0, index ($Value, "X");
						$Vec->{x} =~ s/\s//g;
						$Value = substr $Line, index ($Line, "X")+1;

						$Vec->{y} = substr $Value, 0, index ($Value, "Y");
						$Vec->{y} =~ s/\s//g;
						$Value = substr $Line, index ($Line, "Y")+1;

						$Vec->{z} = substr $Value, 0, index ($Value, "Z");
						$Vec->{z} =~ s/\s//g;

						$Data->{Trans}[$Trans]{TransMom} = $Vec;
					}
					elsif ($Line =~ m/Strength/) {
						$Value = &GetNumber ($Line);
						$Data->{Trans}[$Trans]{Strength} = $Value;
					}
					elsif ($Line =~ m/Multiplicity/) {
						$Value = $Line;
						$Value =~ s/[\s\t\n\r]+$//;
						$Value = substr $Value, rindex ($Value, " ")+1, length ($Value)-rindex ($Value, " ");
						$Data->{Trans}[$Trans]{Multiplicity} = $Value;
					}
					
					push @{$Data->{Trans}[$Trans]{Result}}, $Line;
					++$i;
					$Line = $Content->[$i];
				}

				++$Trans;
			} # of if ($Line =~ m/Excited state/) 
		} # of until ($Line =~ m/---------------/) 
	} # of if ($i < $#{$Content})
} # of sub ReadResults

####################################################################################################

sub GetNumber { # reads the last number from a line
	my $Value = shift;
	
	$Value =~ s/[\s\t\n\r]+$//;
	$Value = substr $Value, rindex ($Value, " "), length ($Value)-rindex ($Value, " ");
	return $Value;
}

1; # return true (mandatory for all packages)
