package GetBaseName;

####################################################################################################
#
# Package:    GetBaseName
#
# Function:   To split a file name into its base name and the extension (similar to the fileparse
#             command).
#
# Author:     Benjamin Bulheller
#
# Version:    $Revision: 3382 $, $Date: 2008-10-12 14:59:58 +0200 (Sun, 12 Oct 2008) $
#
# Date:       December 2005
#
# Usage:
#
# $BaseName  = GetBaseName  ($FileName);    # returns everything before the last dot
#                                             (or $Filename, if not dot is found)
#
# $Extension = GetExtension ($FileName);    # returns everything after the last dot
#                                             (or undef, if no dot is found)
#
# ($BaseName, $Extension) = SplitFileName ($FileName); # both in one command
#
####################################################################################################

use strict;

require Exporter;
our @ISA    = qw( Exporter );
our @EXPORT = qw( GetBaseName GetExtension SplitFileName );

####################################################################################################

sub GetBaseName { # returns the base name of a given file name
	my $FileName = shift;
	my ($DotPosition, $BaseName);
	
	if (not $FileName) { return undef }
	
	# separate path, file name and extension
	$DotPosition = rindex ($FileName, "."); # determine the position of the last dot
	
	if ($DotPosition != -1) { # if a dot has been found
		$BaseName = substr ($FileName, 0, $DotPosition); # take everything till the last dot
		return $BaseName;
	}
	else {
		return $FileName;
	}
} # of sub GetBaseName


sub GetExtension { # returns the extension (everything after the last dot, without the dot itself
	my $FileName = shift;
	my ($DotPosition, $Extension);
	
	if (not $FileName) { return undef }
	
	# separate path, file name and extension
	$DotPosition = rindex ($FileName, "."); # determine the position of the last dot
	
	if ($DotPosition != -1) { # if a dot has been found
		$Extension = substr ($FileName, $DotPosition+1); # take everything after the last dot
		return $Extension;
	}
	else {
		return undef;
	}

} # of sub GetExtension


sub SplitFileName { # returns the base name and the extension
	my $FileName = shift;

	my $BaseName  = &GetBaseName  ($FileName);
	my $Extension = &GetExtension ($FileName);

	return ($BaseName, $Extension);
} # of sub SplitFileName


1;

