package GetParameters;

####################################################################################################
#
# Package:         GetParameters
#
# Author:          Benjamin Bulheller
#
# Website:         www.bulheller.com
#
# Mail address:    webmaster.-at-.bulheller.com
#
# Version:         $Revision: 4552 $, $Date: 2009-06-06 07:55:42 +0100 (Sat, 06 Jun 2009) $
#
# Date:            November 2005
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

require Exporter;
our @ISA    = qw( Exporter );
our @EXPORT = qw( GetParameters );
our $ERROR  = "\nERROR (GetParameters.pm):";

####################################################################################################


sub GetParameters {
	my @Parameters;  # hold the keys of the given Parameter Hash, i.e. all possible defined parameters
	my ($CurPar, $Boolean, $String, $File);
	my ($ParameterHash, $Help, $OptionsHash);
	my (@FileList, $WasGiven);
	
	
	if (scalar @_ < 3) { # if less than 3 parameters are passed to the sub routine
		print STDERR "$ERROR The command line parameters or the help message is not defined.\n\n";
		exit 100;
	} # of if (scalar(@_) < 2)
	
	$ParameterHash = shift;  # the definitions of the parameters (that is, which are possible, and what types they are
	$OptionsHash   = shift;  # the hash that is filled with the parsed arguments
	$Help          = shift;  # the help message
	
	# if -h is not defined as a program parameter, use is as help request and check for it 
	if (not $ParameterHash->{h}) {  
		foreach (@ARGV) { # if -h parameters is given, display help screen
			if ($_ =~ m/^-h/) {
				print STDERR $Help;
				exit 101;
			}
		} # of foreach
	}

	# check whether the "rest" array is defined (in which everything without a preceeding
	# -option is collected). This has to be a list.
	if (defined $ParameterHash->{rest} and $ParameterHash->{rest} !~ m/list/) {
		print STDERR "$ERROR If the 'rest' parameter is defined it has to be a list.\n\n";
		exit 113;
	}
	
	# display help if no parameter is given at all
	if (not @ARGV) {
		print STDERR $Help;
		exit 102;
	}
	
	# parse the parameters
	until (not @ARGV) {
		$CurPar = shift @ARGV;
		
		# if the parameters end with =0 or =1
		if ($CurPar =~ m/=[01]$/) {
			# set $Boolean to the last number after the "="
			$Boolean = substr $CurPar, (length $CurPar)-1;
			# remove the =[01]
			$CurPar =~ s/=[01]//;
		}
		else {
			$Boolean = undef;
		}
		
		if ( ($CurPar =~ m/^-/) and ($CurPar !~ m/^-\d/)  # if the argument begins with a dash but no number follows
			  or ($ParameterHash->{rest}) )  # or the {rest} key is defined and this parameter is going to go in there
		{
			
			# If the read parameter in $CurPar does not begin with "-" or is a negative number (-1.3456)
			# it is going to be added to the {rest} array. If {rest} was predefined as e.g. a integerlist
			# or a filelist then the type of this parameter has to be checked for being an integer or a
			# file. In order to go through the following checks, it has to be faked that the user actually
			# gave this parameter after a -rest option and $CurPar is added to the beginning of @ARGV again
			# where it will be read further down
			if ($ParameterHash->{rest} and ($CurPar !~ m/^-/ or $CurPar =~ m/^-\d/) ) {
				unshift @ARGV, $CurPar;    # add it again to the beginning of @ARGV
				$CurPar = "-rest";         # "fake" the -rest option
			}
			
			$CurPar =~ s/^-//;   # remove the leading dash
			
			# the "--" parameter is only used to end list input
			if ($CurPar eq "-") { next }
			
			# if parameter has not been defined in the options hash
			if (not $ParameterHash->{$CurPar}) {
				print STDERR "$ERROR -$CurPar is not a defined Parameter!\n\n";
				print STDERR $Help;
				exit 103;
			}
			
			if ($ParameterHash->{$CurPar} =~ m/^switch\*?$/) {
				# if =0 or =1 was given, set the parameter to the chosen value, otherwise the switch is "true"
				if (defined $Boolean) {
					$OptionsHash->{$CurPar} = $Boolean;
					$OptionsHash->{GivenParams}{$CurPar} = $Boolean;
				}
				else {
					$OptionsHash->{$CurPar} = 1;
					$OptionsHash->{GivenParams}{$CurPar} = 1;
				}
			}
			
			elsif ($ParameterHash->{"$CurPar"} =~ m/^(string|integer|real|file)?list/) {
				# Parameters given via the command line overwrite default parameters. It needs to be checked, whether
				# default array elements had been declared, which will be erased if this parameter is given. On the 
				# other hand, list parameters may also be split, e.g. "prog -l as er ty -f file -l fg cx" would result
				# in a list "l" containing the five given elements.
				# The parameter's hash entry in $WasGiven will be true, if it had been passed via the command line before,
				# i.e. then existing elements of this array will be kept. If the entry is false, then existing values were
				# given as default parameters and will be erased.
				
				# delete default values
				if (not $WasGiven->{$CurPar}) { $OptionsHash->{$CurPar} = [] }
				
				# 'mark' this parameter, so that more elements can be added later without deleting these ones
				$WasGiven->{$CurPar} = 1;  # set this parameter true
				
				@FileList = ();
				
				# add the arguments following the list parameter to @FileList
				# the -- parameter to end the list input is caught by the two regexes
				until ( (not @ARGV) or ( ($ARGV[0] =~ m/^-/) and ($ARGV[0] !~ m/^-\d/) ) ) {
					$String = shift (@ARGV);                    # until the next parameter appears or @ARGV is empty
					
					if ($ParameterHash->{"$CurPar"} =~ m/^integerlist/) { &CheckInteger ($String, $CurPar) }
					if ($ParameterHash->{"$CurPar"} =~ m/^reallist/)    { &CheckReal    ($String, $CurPar) }
					if ($ParameterHash->{"$CurPar"} =~ m/^filelist/)    { $String = &CheckFile    ($String, $ParameterHash->{$CurPar}) }
					
					push @FileList, $String;
				} # of until (($ARGV[0] =~ m/^-/) or (not @ARGV))
				
				push @{$OptionsHash->{$CurPar}}, @FileList;
				push @{$OptionsHash->{GivenParams}{$CurPar}}, @FileList;
			} # of elsif (${$ParameterHash}{$CurPar} eq "list")
			
			elsif ($ParameterHash->{$CurPar} =~ m/^string\*?$/) {
				$String = shift @ARGV; # read the parameter that followed -$CurPar
				$OptionsHash->{$CurPar} = $String;
				$OptionsHash->{GivenParams}{$CurPar} = $String;
			}
			
			elsif ($ParameterHash->{$CurPar} =~ m/^integer\*?$/) {
				$String = shift @ARGV; # read the parameter that followed -$CurPar
				&CheckInteger ($String, $CurPar);
				$OptionsHash->{$CurPar} = $String;
				$OptionsHash->{GivenParams}{$CurPar} = $String;
			}
			
			elsif ($ParameterHash->{$CurPar} =~ m/^real\*?$/) {
				$String = shift @ARGV; # read the parameter that followed -$CurPar
				&CheckReal ($String, $CurPar);
				$OptionsHash->{$CurPar} = $String;
				$OptionsHash->{GivenParams}{$CurPar} = $String;
			}
			
			elsif ($ParameterHash->{$CurPar} =~ m/^file(\[[^\]]+\])?\*?/) {
				$File = shift @ARGV; # read the parameter that followed -$CurPar
				
				$File = &CheckFile ($File, $ParameterHash->{$CurPar});
				print Dumper $File;
			}
			
			else { # if the defintion of the parameter type is wrong (that is, not switch, string or list
				print STDERR "\nParameter type \"", $ParameterHash->{"$CurPar"}, "\" illegal.\n";
				print STDERR "Either \"switch\", \"string\" or \"list\" is expected.\n\n";
				exit 104;
			}
		} # if ($CurPar =~ m/^-/)
		
		else { # if the argument does not begin with a dash
			push @{$OptionsHash->{rest}}, $CurPar;
			push @{$OptionsHash->{GivenParams}{rest}}, $CurPar;
		} # of else
		
	} # of until (not @ARGV)
	
	# check for missing mandatory parameters
	foreach $CurPar (keys %{$ParameterHash}) {
		# if the parameter type declaration contains an asterisk
		if ($ParameterHash->{$CurPar} =~ m/\*/) {
			# but the parameter was not given
			if (not defined $OptionsHash->{$CurPar}) {
				if ($CurPar eq "rest") {
					print STDERR "\nAdditional values (without any parameter or switch) are required.\n$Help";
				}
				else {
					print STDERR "\nMandatory parameter -$CurPar not given!\n$Help";
				}
				
				exit 105;
			}
		}
	} # of foreach $CurPar (keys %{$ParameterHash})
	
	# check for min, max numbers of elements in lists
	foreach $CurPar (keys %{$ParameterHash}) {
		# if the parameter type declaration defines a list
		if ($ParameterHash->{$CurPar} =~ m/(integer|real|string|file)?list/) {
			# if the definition contains a min/max-value or range
			if ($ParameterHash->{$CurPar} =~ m/\[(\d+)?(,)?(\d+)?\]/) {
				# save the backreferences of the RegEx
				my $Min   = $1; if (not defined $Min) { $Min = 0 }
				my $Range = $2; # true if a comma was given, false otherwise
				my $Max   = $3;
				
				if (defined $Max and $Max == 0) {
					print STDERR "\nERROR: For the parameter -$CurPar a maximum amount of 0 options was defined.\n";
					print STDERR "       This has to be greater than zero.\n\n";
					exit 140;
				}
				
				if (not defined $Max) { $Max = 0 }	
				my $Size;
				
				# determine the size of the array in question
				if (defined $OptionsHash->{$CurPar}) {
					$Size = scalar @{$OptionsHash->{$CurPar}};
				}
				else {
					next;
				}
				
				if (not $Range and $Size != $Min) {
					if ($Min == 1) { $Min = "1 option"     }
					          else { $Min = "$Min options" }
					
					if ($CurPar eq "rest") { 
						print STDERR "\nERROR: Only $Min (values without any parameter or switch) are allowed!\n$Help!";
					}
					else {
						print STDERR "\nERROR: Parameter -$CurPar requires $Min!\n$Help!";
					}
					
					exit 108;
				}
				
				if ($Size < $Min) {
					if ($Min == 1) { $Min = "1 option"     }
					          else { $Min = "$Min options" }
					
					if ($CurPar eq "rest") { 
						print STDERR "\nERROR: At least $Min (values without any parameter or switch) are required!\n$Help!";
					}
					else {
						print STDERR "\nERROR: Parameter -$CurPar needs at least $Min!\n$Help!";
					}
					
					exit 106;
				}
				
				if ($Max > 0 and $Size > $Max) {
					if ($Max == 1) { $Max = "1 option"     }
					          else { $Max = "$Min options" }
					
					if ($CurPar eq "rest") { 
						print STDERR "\nERROR: At most $Max (values without any parameter or switch) are allowed!\n$Help!";
					}
					else {
						print STDERR "\nERROR: Parameter -$CurPar can have at most $Max!\n$Help!";
					}
					
					exit 107;
				}
			}
		}
	} # of foreach $CurPar (keys %{$ParameterHash})
	
	return 1;
} # of sub GetParameters


sub CheckInteger { # checks whether a given string is an integer number
	my $String = shift;
	my $CurPar = shift;
	
	if (not defined $String or $String !~ m/^-{0,1}\d+$/) {
		print STDERR "$ERROR -$CurPar must be an integer!\n\n";
		exit 109;
	}
} # of sub CheckInteger


sub CheckReal { # checks whether a given string is a real number
	my $String = shift;
	my $CurPar = shift;
	
	if (not defined $String or $String !~ m/^-{0,1}\d*\.{0,1}\d+$/) {
		print STDERR "$ERROR -$CurPar must be a real number!\n\n";
		exit 110;
	}
} # of sub CheckReal

sub CheckFile { # checks for a given file using multiple extensions if given
	my $File   = shift;
	my $CurPar = shift;
	my ($Extension, @Extensions);
	
	if (-f $File) { # if the file was found
		return $File;
	}
	# if  the file was not found but extensions were given,
	# containing "{" and "}" with multiple characters (not "}") in between
	elsif ($CurPar =~ m/\{[^\}]+\}/) {
		# cut out everything between { and }, syntax is substr (string, startpos, length)
		$Extension = substr $CurPar, # the parameter definition
		                    index ($CurPar, "{")+1, # the position of the first "{" (including it)
		                    # the position of the last "}"               again the position of the first "{"
		                    rindex ($CurPar, "}") - (index ($CurPar, "{")+1);
		                    # the last line is the number of chars in between { and } (the "lenght" of the string")
		
		# split the string into multiple extensions, separates by "|" (this also works with only one extension)
		@Extensions = split /\|/, $Extension;
		
		foreach $Extension ( @Extensions ) {
			# if the file is found with the current extension
			if (-f "$File.$Extension") {
				return "$File.$Extension";
			}
		}
	}
	
	# if the file was not found and no extensions were defined
	print STDERR "$ERROR File $File not found.\n";
	
	if (@Extensions) {
		print STDERR "       The following extensions were tested: ", join (", ", @Extensions), "\n";
	}
	
	print "\n";
	exit 114;
} # of sub CheckFile

1;
