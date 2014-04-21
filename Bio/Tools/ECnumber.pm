#
# BioPerl module for Bio::Tools::ECnumber
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Christian M. Zmasek <czmasek-at-burnham.org> or <cmzmasek@yahoo.com>
#
# (c) Christian M. Zmasek, czmasek-at-burnham.org, 2002.
# (c) GNF, Genomics Institute of the Novartis Research Foundation, 2002.
#
# You may distribute this module under the same terms as perl itself.
# Refer to the Perl Artistic License (see the license accompanying this
# software package, or see http://www.perl.com/language/misc/Artistic.html)
# for the terms under which you may use, modify, and redistribute this module.
#
# THIS PACKAGE IS PROVIDED "AS IS" AND WITHOUT ANY EXPRESS OR IMPLIED
# WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED WARRANTIES OF
# MERCHANTIBILITY AND FITNESS FOR A PARTICULAR PURPOSE.
#

# POD documentation - main docs before the code


=head1 NAME

Bio::Tools::ECnumber - representation of EC numbers (Enzyme Classification)

=head1 SYNOPSIS

  use Bio::Tools::ECnumber;

  # Creation of ECnumber objects
  my $EC1 = Bio::Tools::ECnumber->new( -ec_string => "4.3.2.1" );
  my $EC2 = Bio::Tools::ECnumber->new( -ec_string => "EC 1.1.1.1" );
  my $EC3 = Bio::Tools::ECnumber->new();

  # Copying
  my $EC4 = $EC1->copy();

  # Modification/canonicalization of ECnumber objects
  print $EC3->EC_string( "1.01.01.001" ); # Prints "1.1.1.1".

  # Stringify
  print $EC3->EC_string();
  # or
  print $EC3->to_string();

  # Test for equality
  # -- Against ECnumber object:
  if ( $EC3->is_equal( $EC2 ) ) { # Prints "equal".
      print "equal";
  }
  # -- Against string representation of EC number:
  if ( ! $EC3->is_equal( "1.1.1.-" ) ) { # Prints "not equal".
      print "not equal";
  }

  # Test for membership
  my $EC5 = Bio::Tools::ECnumber->new( -ec_string => "4.3.2.-" ); 
  # -- Against ECnumber object.
  if ( $EC1->is_member( $EC5 ) ) { # Prints "member".
      print "member"; 
  }
  # -- Against string representation of EC number.
  if ( ! $EC1->is_member( "4.3.1.-" ) ) { # Prints "not member".
      print "not member";
  }

=head1 DESCRIPTION

L<Bio::Tools::ECnumber> is a representation of EC numbers, 
the numerical heirarchy for Enzyme Classification.

See L<http://www.chem.qmul.ac.uk/iubmb/enzyme/> for more details.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists.  Your participation is much appreciated.

  bioperl-l@bioperl.org                  - General discussion
  http://bioperl.org/wiki/Mailing_lists  - About the mailing lists

=head2 Support 

Please direct usage questions or support issues to the mailing list:

I<bioperl-l@bioperl.org>

rather than to the module maintainer directly. Many experienced and 
reponsive experts will be able look at the problem and quickly 
address it. Please include a thorough description of the problem 
with code and data examples if at all possible.

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via the
web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR

Christian M. Zmasek

Email: czmasek-at-burnham.org  or  cmzmasek@yahoo.com

WWW:   http://monochrome-effect.net/

Address: 

  Genomics Institute of the Novartis Research Foundation
  10675 John Jay Hopkins Drive
  San Diego, CA 92121

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...

package Bio::Tools::ECnumber;
use strict;

use constant DEFAULT => "-";
use constant TRUE    => 1;
use constant FALSE   => 0;

use base qw(Bio::Root::Root);

=head2 new

 Title   : new
 Usage   : $EC1 = Bio::Tools::ECnumber->new( -ec_string => "4.3.2.1" );
           or
           $EC2 = Bio::Tools::ECnumber->new( -ec_string => "4.3.2.2",
                                             -comment   => "Is EC 4.3.2.2" );
           or                      
           $EC3 = Bio::Tools::ECnumber->new(); # EC3 is now "-.-.-.-"                      
 Function: Creates a new ECnumber object.
           Parses a EC number from "x.x.x.x", "EC x.x.x.x",
           "ECx.x.x.x", or "EC:x.x.x.x";
           x being either a positive integer or a "-".
 Returns : A new ECnumber object.
 Args    : A string representing a EC number, e.g. "4.3.2.1"
           or "EC 4.3.2.1" or "1.-.-.-".

=cut

sub new {
    my( $class, @args ) = @_;
    
    my $self = $class->SUPER::new( @args );

    my ( $EC_string, $comment )
    = $self->_rearrange( [ qw( EC_STRING COMMENT ) ], @args );

    $self->init(); 
 
    $EC_string && $self->EC_string( $EC_string );
    $comment   && $self->comment( $comment );
                        
    return $self;
    
} # new



=head2 init

 Title   : init()
 Usage   : $EC1->init(); # EC1 is now "-.-.-.-"    
 Function: Initializes this ECnumber to default values.
 Returns : 
 Args    :

=cut

sub init {
    my( $self ) = @_;

    $self->enzyme_class( DEFAULT );
    $self->sub_class( DEFAULT );
    $self->sub_sub_class( DEFAULT );
    $self->serial_number( DEFAULT );
    $self->comment( "" );
  
} # init



=head2 copy

 Title   : copy()
 Usage   : $EC2 = $EC1->copy();
 Function: Creates a new ECnumber object which is an exact copy
           of this ECnumber.
 Returns : A copy of this ECnumber.
 Args    :

=cut

sub copy {
    my( $self ) = @_;
    
    my $new_ec = $self->new();
    $new_ec->enzyme_class(  $self->enzyme_class() );
    $new_ec->sub_class(     $self->sub_class() );
    $new_ec->sub_sub_class( $self->sub_sub_class() );
    $new_ec->serial_number( $self->serial_number() );
    $new_ec->comment(       $self->comment() );
    return $new_ec; 

} # copy



=head2 EC_string

 Title   : EC_string
 Usage   : $EC3->EC_string( "1.1.1.-" );
           or
           print $EC3->EC_string();
 Function: Set/get for string representations of EC numbers.
           Parses a EC number from "x.x.x.x", "EC x.x.x.x",
           "ECx.x.x.x", or "EC:x.x.x.x";
           x being either a positive integer or a "-".
 Returns : A string representations of a EC number.
 Args    : A string representations of a EC number.

=cut

sub EC_string {
    my ( $self, $value ) = @_;

    if ( defined $value) {
        $value =~ s/\s+//g; # Removes white space.
        $value =~ s/^EC//i; # Removes "EC".
        $value =~ s/^://;   # Removes ":".

        if ( $value =~ /^([\d-]*)\.([\d-]*)\.([\d-]*)\.([\d-]*)$/ ) {
            $self->enzyme_class( $1 );
            $self->sub_class( $2 );
            $self->sub_sub_class( $3 );
            $self->serial_number( $4 );
        }
        else {
            $self->throw( "Illegal format error [$value]" );
        }
    }

    return $self->to_string();

} # EC_string



=head2 to_string

 Title   : to_string()
 Usage   : print $EC3->to_string();
 Function: To string method for EC numbers
           (equals the "get" functionality of "EC_string").
 Returns : A string representations of a EC number.
 Args    :

=cut

sub to_string {
    my ( $self ) = @_;

    my $s  = $self->enzyme_class() . ".";
    $s    .= $self->sub_class() . ".";
    $s    .= $self->sub_sub_class() . ".";   
    $s    .= $self->serial_number();
    return $s;
    
} # to_string



=head2 is_equal

 Title   : is_equal
 Usage   : if ( $EC3->is_equal( $EC2 ) )
           or
           if ( $EC3->is_equal( "1.1.1.-" ) )
 Function: Checks whether this ECnumber is equal to the argument
           EC number (please note: "1.1.1.1" != "1.1.1.-").
 Returns : True (1) or false (0).
 Args    : A ECnumber object or a string representation of a EC number.

=cut

sub is_equal {
    my ( $self, $value ) = @_;

    if ( $self->_is_not_reference( $value ) ) {
        $value = $self->new( -ec_string => $value );
    }
    else {
        $self->_is_ECnumber_object( $value );
    }
    
    unless ( $self->enzyme_class() eq $value->enzyme_class() ) {
        return FALSE;
    } 
    unless ( $self->sub_class() eq $value->sub_class() ) {
        return FALSE;
    } 
    unless ( $self->sub_sub_class() eq $value->sub_sub_class() ) {
        return FALSE;
    } 
    unless ( $self->serial_number() eq $value->serial_number() ) {
        return FALSE;
    } 
    return TRUE;

} # is_equal



=head2 is_member

 Title   : is_member
 Usage   : if ( $EC1->is_member( $EC5 ) )
           or
           if ( $EC1->is_member( "4.3.-.-" ) )
 Function: Checks whether this ECnumber is a member of the (incomplete)
           argument EC number (e.g. "1.1.1.1" is a member of "1.1.1.-"
           but not of "1.1.1.2").
 Returns : True (1) or false (0).
 Args    : A ECnumber object or a string representation of a EC number.

=cut

sub is_member {
    my ( $self, $value ) = @_;

    if ( $self->_is_not_reference( $value ) ) {
        $value = $self->new( -ec_string => $value );
    }
    else {
        $self->_is_ECnumber_object( $value );
    }
    $self->_check_for_illegal_defaults();
    $value->_check_for_illegal_defaults();

    unless ( $value->enzyme_class() eq DEFAULT
    ||       $self->enzyme_class() eq $value->enzyme_class() ) {
        return FALSE;
    } 
    unless (  $value->sub_class() eq DEFAULT 
    ||        $self->sub_class() eq $value->sub_class() ) {
        return FALSE;
    } 
    unless ( $value->sub_sub_class() eq DEFAULT
    ||       $self->sub_sub_class() eq $value->sub_sub_class() ) {
        return FALSE;
    } 
    unless ( $value->serial_number() eq DEFAULT
    ||       $self->serial_number() eq $value->serial_number() ) {
        return FALSE;
    } 
    return TRUE;

} # is_member 



=head2 enzyme_class

 Title   : enzyme_class
 Usage   : $EC1->enzyme_class( 1 );
           or 
           print $EC1->enzyme_class(); 
 Function: Set/get for the enzyme class number of ECnumbers.
 Returns : The enzyme class number of this ECnumber.
 Args    : A positive integer or "-".

=cut

sub enzyme_class {
    my ( $self, $value ) = @_;

    if ( defined $value) {
        $self->{ "_enzyme_class" } = $self->_check_number( $value );
    }
   
    return $self->{ "_enzyme_class" };
    
} # enzyme_class



=head2 sub_class

 Title   : sub_class
 Usage   : $EC1->sub_class( 4 );
           or 
           print $EC1->sub_class(); 
 Function: Set/get for the enzyme sub class number of ECnumbers.
 Returns : The enzyme sub class number of this ECnumber.
 Args    : A positive integer or "-".

=cut

sub sub_class {
    my ( $self, $value ) = @_;

    if ( defined $value) {
        $self->{ "_sub_class" } = $self->_check_number( $value );
    }
   
    return $self->{ "_sub_class" };
    
} # sub_class



=head2 sub_sub_class 

 Title   : sub_sub_class
 Usage   : $EC1->sub_sub_class( 12 );
           or 
           print $EC1->sub_sub_class(); 
 Function: Set/get for the enzyme sub sub class number of ECnumbers.
 Returns : The enzyme sub sub class number of this ECnumber.
 Args    : A positive integer or "-".

=cut

sub sub_sub_class {
    my ( $self, $value ) = @_;

    if ( defined $value) {
        $self->{ "_sub_sub_class" } = $self->_check_number( $value );
    }
   
    return $self->{ "_sub_sub_class" };
    
} # sub_sub_class



=head2 serial_number

 Title   : serial_number
 Usage   : $EC1->serial_number( 482 );
           or 
           print $EC1->serial_number(); 
 Function: Set/get for the serial number of ECnumbers.
 Returns : The serial number of this ECnumber.
 Args    : A positive integer or "-".

=cut

sub serial_number {
    my ( $self, $value ) = @_;

    if ( defined $value) {
        $self->{ "_serial_number" } = $self->_check_number( $value );
    }
   
    return $self->{ "_serial_number" };
    
} # serial_number



=head2 comment

 Title   : comment
 Usage   : $EC1->comment( "deprecated" );
           or 
           print $EC1->comment();
 Function: Set/get for a arbitrary comment.
 Returns : A comment [scalar].
 Args    : A comment [scalar].

=cut

sub comment {
    my ( $self, $value ) = @_;

    if ( defined $value) {
        $self->{ "_comment" } = $value;
    }
   
    return $self->{ "_comment" };
    
} # comment



# Title   : _check_number
# Function: Checks and standardizes the individual numbers of a EC number
#           (removes leading zeros, removes white spaces).
# Returns : A standardized number.
# Args    : A string representing a number in a EC number.
sub _check_number {
    my ( $self, $value ) = @_;
    
    my $original_value = $value;
    $value =~ s/\s+//g;   # Removes white space.
    if ( $value eq "" ) {
        $value = DEFAULT;  
    }
    $value =~ s/^0+//;    # Removes leading zeros.
    if ( $value eq "" ) { # If it was "0" (or "00"), it would be "" now.
        $value = "0";
    }
    elsif ( $value ne DEFAULT 
    &&      $value =~ /\D/ ) {
        $self->throw( "Illegal format error [$original_value]" );
    }
    return $value;

} # _check_number



# Title   : _check_for_illegal_defaults()
# Function: Checks for situations like "1.-.1.1", which
#           are illegal in membership tests.
# Returns :
# Args    :
sub _check_for_illegal_defaults {
    my ( $self ) = @_;
   
    if ( ( $self->sub_sub_class() eq DEFAULT 
    &&     $self->serial_number() ne DEFAULT ) ||
         ( $self->sub_class()     eq DEFAULT 
    &&     $self->sub_sub_class() ne DEFAULT ) ||
         ( $self->enzyme_class()  eq DEFAULT 
    &&     $self->sub_class()     ne DEFAULT ) ) {
        $self->throw( "Illegal format error for comparison ["
        . $self->to_string() . "]" );
    } 

} # _check_for_illegal_defaults



# Title   : _is_not_reference
# Function: Checks whether the argument is not a reference.
# Returns : True or false.
# Args    : A scalar.
sub _is_not_reference {
    my ( $self, $value ) = @_;

    return ( ! ref( $value ) );
    
} # _is_not_reference



# Title   : _is_ECnumber_object
# Function: Checks whether the arument is a ECnumber.
# Returns :
# Args    : A reference.
sub _is_ECnumber_object {
    my ( $self, $value ) = @_;

    unless( $value->isa( "Bio::Tools::ECnumber" ) ) {
        $self->throw( "Found [". ref( $value ) 
        ."] where [Bio::Tools::ECnumber] expected" );
    }   
    
} # _is_ECnumber_object



1;
