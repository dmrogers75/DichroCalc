package VectorMath;

#######################################################################################################################
#
# Package:         VectorMath
#
# Function:        Provides several functions for vector calculations
#
# Author:          Benjamin Bulheller
#
# Research Group:  Prof. Jonathan Hirst
#                  School of Physical Chemistry
#                  University of Nottingham
#
# Version:         $Revision: 3775 $, $Date: 2009-01-18 23:45:29 +0000 (Sun, 18 Jan 2009) $
#
# Date:            March 2006
#
# Usage:           $Scalar = Distance ($Point1, $Point2);
#                  $Scalar = DotProduct ($Vector1, $Vector2);
#                  $Vector = CrossProduct ($Vector1, $Vector2);
#                  $Vector = Product ($Scalar, $Vector);
#                  $Scalar = Norm ($Vector);
#                  $Scalar = Angle ($Vector1, $Vector2);
#                  $Scalar = DihedralAngle ($Point1, Point2, Point3, $Point4);
#                  $Vector = VecDiff ($Vector1, $Vector2);
#                  $Vector = VecSum  ($Vector1, $Vector2);
#                  $Vector = Normalize ($Vector);
#                  $Scalar = Deg2Rad ($Angle);
#                  $Scalar = Rad2Deg ($Angle);
#                  $Scalar = xComponent ($Vector);
#                  $Scalar = yComponent ($Vector);
#                  $Scalar = zComponent ($Vector);
#
#                  Rotation  ($Vector, $Angle, $AxisVector); # rotation about arbitrary axis
#                  xRotation ($Vector, $Value);
#                  yRotation ($Vector, $Value);
#                  zRotation ($Vector, $Value);
#                  xMatRotation ($Vector, $Value);
#                  yMatRotation ($Vector, $Value);
#                  zMatRotation ($Vector, $Value);
#                  Translation  ($Vector1, $Vector2);  # translation about arbitrary x,y,z values
#                  xTranslation ($Vector, $Value);
#                  yTranslation ($Vector, $Value);
#                  zTranslation ($Vector, $Value);
#
#######################################################################################################################

use strict;        # always use this!!!
use POSIX;         # for asin and acos functions
use Data::Dumper;  # for printing arrays and hashes

$Data::Dumper::Sortkeys = 1; # sort the hash keys

require Exporter;
our @ISA    = qw( Exporter );
our @EXPORT = qw( Distance DotProduct CrossProduct Product Norm Angle DihedralAngle VecDiff VecSum Normalize
                  xComponent yComponent zComponent
                  Rotation xRotation yRotation zRotation xMatRotation yMatRotation zMatRotation
                  Translation xTranslation yTranslation zTranslation Deg2Rad Rad2Deg );

our $ERROR  = "\nERROR (VectorMath.pm) ";

#######################################################################################################################


sub Distance { # calculates the distance of two points
	my $Point1 = shift;
	my $Point2 = shift;
	my $Distance;
	
	if ( (not $Point1) or (not $Point2) ) {
		print STDERR "$ERROR in Distance: Less than two points have been given!\n";
		return undef;
	}
	
	if ( ($Point1 !~ m/^HASH/) or ($Point2 !~ m/^HASH/) ) {
		print STDERR "$ERROR in Distance: Hash references are required!\n";
		return undef;
	}
	
	for ( qw/x y z/ ) {
		if (not defined $Point1->{$_} or not defined $Point2->{$_}) { return undef }
	}
	
	$Distance = sqrt (   ( $Point1->{x} - $Point2->{x} )**2
						    + ( $Point1->{y} - $Point2->{y} )**2
						    + ( $Point1->{z} - $Point2->{z} )**2 );
						
	return $Distance;
} # of sub Distance


sub DotProduct { # calculates the dot product of two vectors
	my $Vector1 = shift;
	my $Vector2 = shift;
	my ($Theta, $DotProduct);
	
	if ( (not $Vector1) or (not $Vector2) ) {
		print STDERR "$ERROR in DotProduct: Less than two vectors have been given!\n";
		return undef;
	}
	
	if ( ($Vector1 !~ m/^HASH/) or ($Vector2 !~ m/^HASH/) ) {
		print STDERR "$ERROR in DotProduct: Hash references are required!\n";
		return undef;
	}
	
	&Normalize ($Vector1);
	&Normalize ($Vector2);
	
	$DotProduct = (   $Vector1->{x} * $Vector2->{x}
					    + $Vector1->{y} * $Vector2->{y}
					    + $Vector1->{z} * $Vector2->{z} );
					
	return $DotProduct;
} # of sub DotProduct


sub CrossProduct { # calculates the cross product of two vectors
	my $Vector1 = shift;
	my $Vector2 = shift;
	my $CrossProduct;
	
	if ( (not $Vector1) or (not $Vector2) ) {
		print STDERR "$ERROR in CrossProduct: Less than two vectors have been given!\n";
		return undef;
	}
	
	if ( ($Vector1 !~ m/^HASH/) or ($Vector2 !~ m/^HASH/) ) {
		print STDERR "$ERROR in CrossProduct: Hash references are required!\n";
		return undef;
	}
	
	$CrossProduct = { x => ($Vector1->{y} * $Vector2->{z} - $Vector2->{y} * $Vector1->{z}),
	                  y => ($Vector1->{z} * $Vector2->{x} - $Vector2->{z} * $Vector1->{x}),
	                  z => ($Vector1->{x} * $Vector2->{y} - $Vector2->{x} * $Vector1->{y}) };
	
	return $CrossProduct;
} # of sub CrossProduct


sub Product { # multiplies a vector with a scalar value
	my $Scalar = shift;
	my $Vector = shift;
	my $Product;
	
	if ( (not $Scalar) or (not $Vector) ) {
		print STDERR "$ERROR in Product: No scalar and/or no vector have been given!\n";
		return undef;
	}
	
	if ($Vector !~ m/^HASH/) {
		print STDERR "$ERROR in Product: A hash reference is required!\n";
		return undef;
	}
	
	$Product = {
		x => $Vector->{x} * $Scalar,
		y => $Vector->{y} * $Scalar,
		z => $Vector->{z} * $Scalar
	};
	
	return $Product;
} # sub Product


sub Norm { # calculates the length of a vector
	my $Vector = shift;
	my $Norm;
	
	if (not defined $Vector->{x} or not defined $Vector->{y} or not defined $Vector->{z}) {
		print STDERR "$ERROR in Norm: Not all coordinates of vector are defined!\n";
		return undef;
	}
	
	if ($Vector !~ m/^HASH/) {
		print STDERR "$ERROR in Norm: A hash reference is required!\n";
		return undef;
	}
	
	$Norm = sqrt ( $Vector->{x}**2 + $Vector->{y}**2 + $Vector->{z}**2 );
	
	return $Norm;
} # of sub &Norm


sub Normalize { # normalizes a vector
	my $Vector = shift;
	my $Norm;
	
	if (not defined $Vector->{x} or not defined $Vector->{y} or not defined $Vector->{z}) {
		print STDERR "$ERROR in Normalize: Not all coordinates of vector are defined!\n";
		return undef;
	}
	
	if ($Vector !~ m/^HASH/) {
		print STDERR "$ERROR in Normalize: A hash reference is required!\n";
		return undef;
	}
	
	$Norm = &Norm ($Vector);
	
	if ($Norm != 0) {
		$Vector->{x} = $Vector->{x} / $Norm;
		$Vector->{y} = $Vector->{y} / $Norm;
		$Vector->{z} = $Vector->{z} / $Norm;
	}
} # of sub Normalize


sub Angle { # calculates the angle between two vectors
	my $Vector1 = shift;
	my $Vector2 = shift;
	my ($Theta, $DotProduct);
	
	if (not $Vector1 or not $Vector2) {
		print STDERR "$ERROR in Angle: Less than two vectors have been given!\n";
		return undef;
	}
	
	if ( ($Vector1 !~ m/^HASH/) or ($Vector2 !~ m/^HASH/) ) {
		print STDERR "$ERROR in Angle: Hash references are required!\n";
		return undef;
	}
	
	if (&Norm($Vector1) == 0) {
		print STDERR "$ERROR in Angle: norm of vector 1 == 0, division by zero while calculating angle!\n\n";
		return 360;
	}
	elsif (&Norm($Vector2) == 0) {
		print STDERR "$ERROR in Angle: norm of vector 2 == 0, division by zero while calculating angle!\n\n";
		return 360;
	}
	
	$Theta = &DotProduct ($Vector1, $Vector2) / ( &Norm ($Vector1) * &Norm ($Vector2) );
	$Theta = acos ($Theta) * 57.29578; # 180/Pi = 57.29578
	
	return $Theta;
} # of sub Angle


sub DihedralAngle { # calculates the dihedral angle between four given angles
	my $Atom1 = shift;
	my $Atom2 = shift;
	my $Atom3 = shift;
	my $Atom4 = shift;
	
	my $Result = 360.0;
	
	my ($Vector12, $Vector23, $Vector43, $CrossP, $CrossX, $CrossY, $DotX, $DotY, $DotPX, $DotPY);
	
	my $Pi = 3.141593;
	my $Radian = 57.29578;  # = 180 / Pi
	
	if ( (not $Atom1) or (not $Atom2) or (not $Atom3) or (not $Atom4) ) {
		print STDERR "$ERROR in DihedralAngle: Less than four atoms defined!\n";
		return 360;
	}
	
	$Vector12 = &VecDiff ($Atom1, $Atom2);
	$Vector23 = &VecDiff ($Atom2, $Atom3);
	$Vector43 = &VecDiff ($Atom4, $Atom3);
		
	&Normalize ($Vector12);
	&Normalize ($Vector23);
	&Normalize ($Vector43);
	
	$CrossP = &CrossProduct ($Vector23, $Vector12);
	$CrossX = &CrossProduct ($Vector23, $Vector43);
	$CrossY = &CrossProduct ($Vector23, $CrossX);
	
	$DotX = &DotProduct ($CrossX, $CrossX);
	$DotY = &DotProduct ($CrossY, $CrossY);
	
	if ( ($DotX <= 0.0) or ($DotY <= 0.0) ) { return $Result }
	
	$DotPX = &DotProduct ($CrossP, $CrossX) / sqrt ($DotX);
	$DotPY = &DotProduct ($CrossP, $CrossY) / sqrt ($DotY);
	
	if ( ($DotPX == 0.0) or ($DotPY == 0.0) ) { return $Result }
	
	$Result = atan2 ($DotPY, $DotPX) * $Radian;
	
	return $Result;
} # of sub DihedralAngle


sub VecDiff { # subtracts a given vector from another one
	my $Vector1 = shift;
	my $Vector2 = shift;
	my $ResultVector;
	
	if ( (not $Vector1) or (not $Vector2) ) {
		print STDERR "$ERROR in VecDiff: Less than two vectors have been given!\n";
		return undef;
	}
	
	$ResultVector = { x => $Vector1->{x} - $Vector2->{x},
	                  y => $Vector1->{y} - $Vector2->{y},
	                  z => $Vector1->{z} - $Vector2->{z},
	                };
	
	return $ResultVector;
} # of sub VecDiff

sub VecSum { # subtracts a given vector from another one
	my $Vector1 = shift;
	my $Vector2 = shift;
	my $ResultVector;
	
	if ( (not $Vector1) or (not $Vector2) ) {
		print STDERR "$ERROR in VecSum: Less than two vectors have been given!\n";
		return undef;
	}
	
	$ResultVector = { x => $Vector1->{x} + $Vector2->{x},
	                  y => $Vector1->{y} + $Vector2->{y},
	                  z => $Vector1->{z} + $Vector2->{z},
	                };
	
	return $ResultVector;
} # of sub VecDiff


sub xComponent { # calculates the x component of the unit vector along the given vector
	my $Vector = shift;
	my $xAxis  = { x => 1, y => 0, z => 0 };

	my $xAngle = Angle ($Vector, $xAxis);
	return cos ( Deg2Rad ($xAngle));
} # of sub xComponent


sub yComponent { # calculates the y component of the unit vector along the given vector
	my $Vector = shift;
	my $yAxis  = { x => 0, y => 1, z => 0 };

	my $yAngle = Angle ($Vector, $yAxis);
	return cos ( Deg2Rad ($yAngle));
} # of sub yComponent


sub zComponent { # calculates the z component of the unit vector along the given vector
	my $Vector = shift;
	my $zAxis  = { x => 0, y => 0, z => 1 };

	my $zAngle = Angle ($Vector, $zAxis);
	return cos ( Deg2Rad ($zAngle));
} # of sub zComponent


sub xTranslation { # translation of the x coordinate
	my $Coord = shift;
	my $Value = shift;
	
	$Coord->{x} = $Coord->{x} + $Value;
} # of sub xTranslation


sub yTranslation { # translation of the y coordinate
	my $Coord = shift;
	my $Value = shift;
	
	$Coord->{y} = $Coord->{y} + $Value;
} # of sub yTranslation


sub zTranslation { # translation of the z coordinate
	my $Coord = shift;
	my $Value = shift;
	
	$Coord->{z} = $Coord->{z} + $Value;
} # of sub zTranslation


sub Translation { # translation about an arbitrary vector
	my $Vector = shift;
	my $Trans  = shift;
	
	if (not $Vector) {
		print STDERR "$ERROR in Translation: vector not defined!\n";
		return undef;
	}
	
	if (not $Trans) {
		print STDERR "$ERROR in Translation: scalar not defined!\n";
		return undef;
	}
	
	$Vector->{x} = $Vector->{x} + $Trans->{x};
	$Vector->{y} = $Vector->{y} + $Trans->{y};
	$Vector->{z} = $Vector->{z} + $Trans->{z};
} # of sub Translation


sub xRotation { # rotation about x axis
	my $Coord = shift;
	my $Angle = shift;
	
	my ($x, $y, $z);
	
	if (not $Coord) {
		print STDERR "$ERROR in xRotation: vector not defined!\n";
		return undef;
	}
	
	if (not defined $Angle) {
		print STDERR "$ERROR in xRotation: angle not defined!\n";
		return undef;
	}
	
	$Angle = &Deg2Rad ($Angle);
	
	$x =  $Coord->{x};
	$y = ($Coord->{y} * cos ($Angle)) - ($Coord->{z} * sin ($Angle));
	$z = ($Coord->{y} * sin ($Angle)) + ($Coord->{z} * cos ($Angle));
	
	$Coord->{x} = $x;
	$Coord->{y} = $y;
	$Coord->{z} = $z;
} # of sub xRotation


sub yRotation { # rotation about y axis
	my $Coord = shift;
	my $Angle = shift;
	
	my ($x, $y, $z);
	
	if (not $Coord) {
		print STDERR "$ERROR in yRotation: vector not defined!\n";
		return undef;
	}
	
	if (not defined $Angle) {
		print STDERR "$ERROR in yRotation: angle not defined!\n";
		return undef;
	}
	
	$Angle = &Deg2Rad ($Angle);
	
	$x = ($Coord->{x} * cos ($Angle)) + ($Coord->{z} * sin ($Angle));
	$y =  $Coord->{y};
	$z = ($Coord->{z} * cos ($Angle)) - ($Coord->{x} * sin ($Angle));
	
	$Coord->{x} = $x;
	$Coord->{y} = $y;
	$Coord->{z} = $z;
} # of sub yRotation


sub zRotation { # rotation about z axis
	my $Coord = shift;
	my $Angle = shift;
	
	my ($x, $y, $z);
	
	if (not $Coord) {
		print STDERR "$ERROR in zRotation: vector not defined!\n";
		return undef;
	}
	
	if (not defined $Angle) {
		print STDERR "$ERROR in zRotation: angle not defined!\n";
		return undef;
	}
	
	$Angle = &Deg2Rad ($Angle);
	
	$x = ($Coord->{x} * cos ($Angle)) - ($Coord->{y} * sin ($Angle));
	$y = ($Coord->{x} * sin ($Angle)) + ($Coord->{y} * cos ($Angle));
	$z =  $Coord->{z};
	
	$Coord->{x} = $x;
	$Coord->{y} = $y;
	$Coord->{z} = $z;
} # of sub zRotation


sub xMatRotation { # rotation about x axis
	my $Coord = shift;
	my $Angle = shift;
	
	my ($x, $y, $z, @xMat, @yMat, @zMat);
	
	if (not $Coord) {
		print STDERR "$ERROR in xMatRotation: vector not defined!\n";
		return undef;
	}
	
	if (not defined $Angle) {
		print STDERR "$ERROR in xMatRotation: angle not defined!\n";
		return undef;
	}
	
	$Angle = &Deg2Rad ($Angle);
	
	# build the rotation matrix
	@xMat = (1,           0,           0);
	@yMat = (0,  cos $Angle,  sin $Angle);
	@zMat = (0, -sin $Angle,  cos $Angle);
	
	$x = $Coord->{x} * $xMat[0] + $Coord->{y} * $xMat[1] + $Coord->{z} * $xMat[2];
	$y = $Coord->{x} * $yMat[0] + $Coord->{y} * $yMat[1] + $Coord->{z} * $yMat[2];
	$z = $Coord->{x} * $xMat[0] + $Coord->{y} * $xMat[1] + $Coord->{z} * $xMat[2];
	
	$Coord->{x} = $x;
	$Coord->{y} = $y;
	$Coord->{z} = $z;
} # of sub xMatRotation


sub yMatRotation { # rotation about y axis
	my $Coord = shift;
	my $Angle = shift;
	
	my ($x, $y, $z, @xMat, @yMat, @zMat);
	
	if (not $Coord) {
		print STDERR "$ERROR in yMatRotation: vector not defined!\n";
		return undef;
	}
	
	if (not defined $Angle) {
		print STDERR "$ERROR in yMatRotation: angle not defined!\n";
		return undef;
	}
	
	$Angle = &Deg2Rad ($Angle);
	
	# build the rotation matrix
	@xMat = ( cos $Angle,      0, -sin $Angle);
	@yMat = (          0,      1,           0);
	@zMat = ( sin $Angle,      0,  cos $Angle);
	
	$x = $Coord->{x} * $xMat[0] + $Coord->{y} * $xMat[1] + $Coord->{z} * $xMat[2];
	$y = $Coord->{x} * $yMat[0] + $Coord->{y} * $yMat[1] + $Coord->{z} * $yMat[2];
	$z = $Coord->{x} * $xMat[0] + $Coord->{y} * $xMat[1] + $Coord->{z} * $xMat[2];
	
	$Coord->{x} = $x;
	$Coord->{y} = $y;
	$Coord->{z} = $z;
} # of sub yMatRotation


sub zMatRotation { # rotation about z axis
	my $Coord = shift;
	my $Angle = shift;
	
	my ($x, $y, $z, @xMat, @yMat, @zMat);
	
	if (not $Coord) {
		print STDERR "$ERROR in zMatRotation: vector not defined!\n";
		return undef;
	}
	
	if (not defined $Angle) {
		print STDERR "$ERROR in zMatRotation: angle not defined!\n";
		return undef;
	}
	
	$Angle = &Deg2Rad ($Angle);
	
	# build the rotation matrix
	@xMat = ( cos $Angle, sin $Angle, 0);
	@yMat = (-sin $Angle, cos $Angle, 0);
	@zMat = (          0,          0, 1);
	
	$x = $Coord->{x} * $xMat[0] + $Coord->{y} * $xMat[1] + $Coord->{z} * $xMat[2];
	$y = $Coord->{x} * $yMat[0] + $Coord->{y} * $yMat[1] + $Coord->{z} * $yMat[2];
	$z = $Coord->{x} * $xMat[0] + $Coord->{y} * $xMat[1] + $Coord->{z} * $xMat[2];
	
	$Coord->{x} = $x;
	$Coord->{y} = $y;
	$Coord->{z} = $z;
} # of sub zMatRotation


sub Rotation { # rotation about an arbitrary axis
	# Rotates a given point about an arbitrary axis u about a given angle. The axis starts at the origin
	# and points to the given end point.
	# See Anton & Rorres, Elementary Linear Algebra (Applications Version)
	# 8th Edition, page 181
	
	my $Coord = shift;
	my $Angle = shift;
	my $Axis  = shift;
	
	my ($x, $y, $z, @xMat, @yMat, @zMat, $ux, $uy, $uz, $Fac);
	
	if (not defined $Coord) {
		print STDERR "$ERROR in Rotation: \$Coord not defined!\n";
		return undef;
	}
	
	if (not defined $Angle) {
		print STDERR "$ERROR in Rotation: \$Angle not defined!\n";
		return undef;
	}
	
	if (not defined $Axis) {
		print STDERR "$ERROR in Rotation: \$Axis not defined!\n";
		return undef;
	}
	
	$Angle = &Deg2Rad ($Angle);
	&Normalize ($Axis);
	
	$ux = $Axis->{x};
	$uy = $Axis->{y};
	$uz = $Axis->{z};
	
	$Fac = 1 - cos $Angle;
	
	# build the rotation matrix
	@xMat = ( $ux*$ux * $Fac +       cos $Angle,   $uy*$ux * $Fac - $uz * sin $Angle,   $uz*$ux * $Fac + $uy * sin $Angle );
	@yMat = ( $ux*$uy * $Fac + $uz * sin $Angle,   $uy*$uy * $Fac +       cos $Angle,   $uz*$uy * $Fac - $ux * sin $Angle );
	@zMat = ( $ux*$uz * $Fac - $uy * sin $Angle,   $uy*$uz * $Fac + $ux * sin $Angle,   $uz*$uz * $Fac +       cos $Angle );
	
	# multiply the vector with the rotation matrix
	$x = $Coord->{x} * $xMat[0] + $Coord->{y} * $xMat[1] + $Coord->{z} * $xMat[2];
	$y = $Coord->{x} * $yMat[0] + $Coord->{y} * $yMat[1] + $Coord->{z} * $yMat[2];
	$z = $Coord->{x} * $zMat[0] + $Coord->{y} * $zMat[1] + $Coord->{z} * $zMat[2];
	
	$Coord->{x} = $x;
	$Coord->{y} = $y;
	$Coord->{z} = $z;
} # of sub Rotation


sub Deg2Rad { # converts degrees to readians
	my $Pi = atan2 (1,1) * 4;  # camel book, page 684
	my $Rad = $Pi/180;         # convert degrees to radians
	
	my $Angle = shift;
	return $Angle * $Rad;
} # of sub Deg2Rad


sub Rad2Deg { # converts radians to degrees
	my $Pi = atan2 (1,1) * 4;  # camel book, page 684
	my $Deg = 180/$Pi;         # convert radians to degrees
	
	my $Angle = shift;
	return $Angle * $Deg;
} # of sub Rad2Deg


1;
