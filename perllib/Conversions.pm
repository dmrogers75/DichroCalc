package Conversions;

####################################################################################################
#
# Package:    Conversions
#
# Function:   Provides methods for converting units (mostly energy/wavelengths etc.)
#
# Usage:      $Wavelength = eV2nm ( $ElectronVolts );
#
# Author:     Benjamin Bulheller
#
# Version:    $Revision: 3181 $, $Date$
#
# Date:       July 2008
#
####################################################################################################

use strict;                          # always use this!!!
use Data::Dumper;                    # for easy printout of arrays and hashes
use FindBin qw/$Bin/;                # sets $Bin to the directory of the script
use lib $Bin;                        # add the script's directory to the library path
use lib "$ENV{HOME}/bin/perllib/";   # adds ~/bin/perllib to the library path
use DebugPrint;                      # handy during debugging

require Exporter;
our @ISA    = qw( Exporter );
our @EXPORT = qw( eV2nm nm2eV nm2wn eV2hartree hartree2eV wn2nm wn2eV eV2wn );

####################################################################################################
	

sub eV2nm {	1239.84 / $_[0] }

sub nm2eV {	1239.84 / $_[0] }


sub nm2wn { 1E7 / $_[0] }

sub wn2nm { 1E7 / $_[0] }


sub wn2eV { (1239.84 * $_[0]) / 1E7 }

sub eV2wn { 1E7/(1239.84/$_[0]) }


sub eV2hartree { $_[0] / 27.211383 }

sub hartree2eV { $_[0] * 27.211383 }

