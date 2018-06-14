package ReadSpectrum;

####################################################################################################
#
# Package:    ReadSpectrum
#
# Function:   Reads xy data like spectra from a file and returns it as array. The order of the
#             the first column is checked to be ascending and reversed, if not.
#
# Usage:      @Spectrum = ReadSpectrum ("filename");
#             ($xMin, $xMax) = xExtrema (\@Spectrum);  # add ",1" to paramters for symmetric range
#             ($yMin, $yMax) = yExtrema (\@Spectrum);  # add ",1" to paramters for symmetric range
#
# Author:     Benjamin Bulheller
#
# Version:    $Revision: 3089 $, $Date: 2008-07-13 21:46:51 +0200 (Sun, 13 Jul 2008) $
#
# Date:       March 2006
#
####################################################################################################

use strict;                   # always use this!!!
use Data::Dumper;             # to print out arrays and hashes

require Exporter;
our @ISA    = qw( Exporter );
our @EXPORT = qw( ReadSpectrum xExtrema yExtrema );

sub ReadSpectrum {
	my $FileName = shift;
	my (@Content, @Fields, $Line, @Data);
	
	if (not $FileName) {
		print "\nERROR (ReadSpectrum.pm): No file name given!\n\n";
		return undef;
	}
	
	if (not -f $FileName) {
		print "\nERROR (ReadSpectrum.pm): File $FileName not found!\n\n";
		return undef;
	}
	
	open FILE, "<$FileName";
	@Content = <FILE>;       # read the file
	chomp @Content;          # remove line feeds
	close FILE;
	
	# if only one line and DOS-line feeds are found
	if (scalar @Content == 1 and $Content[0] =~ m//) {
		# replace all DOS-line feeds with proper unix ones
		$Content[0] =~ s//\n/g;
		@Content = split "\n", $Content[0];
	}
	
	while (@Content) {
		$Line = shift @Content;
		
		$Line =~ s/^\s+//;             # remove leading spaces
		$Line =~ s/\s+$//;             # remove trailing spaces
		
		if ($Line eq "")    { next }   # skip empty lines
		if ($Line =~ m/^#/) { next }   # skip comment lines	

		@Fields = split /[\s\t]+/, $Line; # split at multiple blanks or tabs
		
		push @Data, [@Fields];
	}
	
	# if the first and last field are both numbers and the very first wavelength
	# is bigger than the very last one, reverse the sorting
	if ($Data[0]->[0] =~ m/\d+/g and $Data[0]->[0] =~ m/\d+/ and
		$Data[0]->[0] > $Data[$#Data]->[0]) {
		@Data = reverse @Data;
	}
	
	return @Data;
} # of sub ReadSpectrum


sub xExtrema { # returns the absolute minumum and maximum of the x axis
	my $Spectrum  = shift;
	my $Symmetric = shift;

	return &Extrema ($Spectrum, 0, $Symmetric);
} # of sub xExtrema


sub yExtrema { # returns the absolute minumum and maximum of the y axis
	my $Spectrum = shift;
	my $Symmetric = shift;

	return &Extrema ($Spectrum, 1, $Symmetric);
} # of sub yExtrema


sub Extrema { # returns the absolute minimum and maximum of a given colum
	my $Spectrum = shift;
	my $Column   = shift;
	my $Symmetric = shift;
	
	my $Min = +1E20;
	my $Max = -1E20;
	my $Point;

	foreach $Point ( @{$Spectrum} ) {
		if ($Point->[$Column] < $Min) { $Min = $Point->[$Column] }
		if ($Point->[$Column] > $Max) { $Max = $Point->[$Column] }
	}
	
	if ($Symmetric) {
		if (abs $Max > abs $Min) {
			$Max = abs $Max;
			$Min = -$Max;
		}
		else {
			$Max = abs $Min;
			$Min = $Min;
		}
	}

	return ($Min, $Max);
} # of sub Extrema


1;
