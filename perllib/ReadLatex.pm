package ReadLatex;

####################################################################################################
#
# Package:         ReadLatex
#
# Function:        Provides methods for dealing with LaTeX files
#
# Usage:           ReadLatex (file.tex);
#
#                  Returns an array with the latex source code, included files via input are considered
#
# Author:          Benjamin Bulheller
#
# Website:         www.bulheller.com
#
# Mail address:    webmaster.-at-.bulheller.com
#
# Version:         $Revision: 4553 $, $Date: 2009-06-06 07:56:44 +0100 (Sat, 06 Jun 2009) $
#
# Date:            May 2007
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

use strict;                          # always use this!!!
use Data::Dumper;                    # for easy printout of arrays and hashes
use FindBin qw/$Bin/;                # sets $Bin to the directory of the script
use lib $Bin;                        # add the script's directory to the library path
use lib "$ENV{HOME}/bin/perllib/";   # adds ~/bin/perllib to the library path
use DebugPrint;                      # handy during debugging

$Data::Dumper::Sortkeys = 1; # sort the hash keys
# $Data::Dumper::Indent = 3;

require Exporter;
our @ISA    = qw( Exporter );
our @EXPORT = qw( ReadLatex );

####################################################################################################
	
sub ReadLatex { # reads in a LaTeX source code and deals with include directives
	my $File   = shift;
	my $Silent = shift;
	my (@Content, @Source, $InputFile, @Input, $Line);

	if (not -f $File) {
		if (-f "$File.tex") {
			$File = $File . ".tex";
		}
		else {
			print STDERR "ERROR (ReadLatex): File $File could not be found!\n\n";
			exit 100;
		}
	}
	
	open FILE, "<$File";
	@Content = <FILE>;
	close FILE;
	
	while (@Content) {
		$Line = shift @Content;
		
		if ($Line =~ m/\\input(\s+)?\{[^}]+\}/) {
			$InputFile = $Line;
			
			# cut away everything till after "input"
			$InputFile = substr $InputFile, index ($InputFile, "input") + 5;
			# cut away everything till after the next "{"
			$InputFile = substr $InputFile, index ($InputFile, "{") + 1;
			# take only everything between char 0 and the next "}"
			$InputFile = substr $InputFile, 0, index ($InputFile, "}");

			if (not $Silent) { print "Including $InputFile...\n"; }
			
			# if there is an \input{filename} command, cut out the filename
			# $InputFile =~ s/^.*\\input(\s+)?\{([\w\.\/]+)\}.*$/$2/;
			
			if ($InputFile) {
				# check whether the extension has to be added
				if (-f "$InputFile.tex") { $InputFile = $InputFile . ".tex" }
				
				if (not -f $InputFile) {
					print "WARNING (ReadLatex): input file $InputFile not found!\n";
				}
				else {
					open FILE, "<$InputFile";
					@Input = <FILE>;
					close FILE;
					
					unshift @Content, @Input;
					next;
				}
			}
		}
		
		push @Source, $Line;
	}
	
	return @Source;
} # of sub ReadLatex


1;
