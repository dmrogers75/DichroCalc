package DebugPrint;

####################################################################################################
#
# Package:         DebugPrint
#
# Function:        Provides the command dp, which simply dumps the given variables and lets the
#                  program die afterwards. The command dpf prints the output to the file debub.txt
#                  in case of large variable contents. Both commands are handy during debugging as
#                  they can be faster used than print statements.
#
# Author:          Benjamin Bulheller
#
# Website:         www.bulheller.com
#
# Mail address:    webmaster.-at-.bulheller.com
#
# Version:         $Revision: 4262 $, $Date: 2009-04-27 03:10:17 +0200 (Mon, 27 Apr 2009) $
#
# Date:            December 2005
#
# Usage:           use DebugPrint;
#
#                  dp ($Scalar);
#                  dp (\@Array);   # arrays and hashes have to be passed as references
#                  dp (\%Hash);
#                  dp ($Scalar, \@Array1, \@Array2);
#
#                  dpf (\@Array); # dumps the output into the file debug.txt
#
# Licence:         This program is free software: you can redistribute it and/or modify
#                  it under the terms of the GNU General Public License as published by
#                  the Free Software Foundation, either version 3 of the License, or
#                  (at your option) any later version.
#
#                  This program is distributed in the hope that it will be useful,
#                  but WITHOUT ANY WARRANTY; without even the implied warranty of
#                  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#                  GNU General Public License for more details.
#
#                  You should have received a copy of the GNU General Public License
#                  along with this program.  If not, see http://www.gnu.org/licenses/.
#
####################################################################################################

use strict;                             # always use this!!!
use Data::Dumper;                       # for easy printout of arrays and hashes
use FindBin qw/$Bin/;                   # sets $Bin to the directory of the script
use lib $Bin;                           # add the script's directory to the library path

$Data::Dumper::Sortkeys = 1;            # sort the hash keys
# $Data::Dumper::Indent = 3;

require Exporter;
our @ISA    = qw( Exporter );
our @EXPORT = qw( dp dpf );

our $OutputFile = "debug.txt";

####################################################################################################


sub dp { # DebugPrint
	my @args = @_;
	my ($Var, $File);
	
	if ($args[$#args] eq "DumpIntoFile") {   # if it was called by dpf
		open (OUTPUT, ">$OutputFile");        # dump variable into file
		pop @args;                            # remove last element
		$File = 1;
	}
	else {
		open (OUTPUT, ">&STDOUT");            # dump variable to STDOUT
		$File = 0;
	}
	
	while (@args) {
		$Var = shift @args;
		
		print OUTPUT "\n\n";
		
		if (ref $Var) {
			if		($Var =~ m/ARRAY/)  {
				chomp @{$Var};
				print OUTPUT Dumper @{$Var};
			}
			elsif ($Var =~ m/HASH/)   { print OUTPUT Dumper %{$Var} }
			elsif ($Var =~ m/SCALAR/) { print OUTPUT Dumper ${$Var} }
			else  { print "Error determining reference type!\n"     }
		}
		else {
			print OUTPUT Dumper $Var;
		}
		
		print OUTPUT "\n";
	
	} # of while (@args)
	
	close OUTPUT;
	
	if ($File) { die "\nDied happily a controlled death in debug print. Output saved to $OutputFile\n\n" }
	      else { die "\nDied happily a controlled death in debug print.\n\n"                             }
} # of sub dp


sub dpf { # DebugPrintFile, prints to debug.txt
	my @args = @_;
	
	push @args, "DumpIntoFile";
	dp (@args);
} # of sub dpf

1;
