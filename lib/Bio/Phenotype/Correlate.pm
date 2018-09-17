#
# BioPerl module for Bio::Phenotype::Correlate
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
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Phenotype::Correlate - Representation of a correlating phenotype in a given species

=head1 SYNOPSIS

  use Bio::Phenotype::Correlate;

  $co = Bio::Phenotype::Correlate->new( -name        => "4(Tas1r3)",
                                        -description => "mouse correlate of human phenotype MIM 605865",
                                        -species     => $mouse,
                                        -type        => "homolog",
                                        -comment     => "type=homolog is putative" );

  print $co->name();
  print $co->description();
  print $co->species()->binomial();
  print $co->type();
  print $co->comment();

  print $co->to_string();

=head1 DESCRIPTION

This class models correlating phenotypes.
Its creation was inspired by the OMIM database where many human phenotypes
have a correlating mouse phenotype. Therefore, this class is intended
to be used together with a phenotype class. 


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
methods.

=cut


# Let the code begin...

package Bio::Phenotype::Correlate;
use strict;
use Bio::Species;

use base qw(Bio::Root::Root);


=head2 new

 Title   : new
 Usage   : $co = Bio::Phenotype::Correlate->new( -name        => "4(Tas1r3)",
                                                 -description => "mouse correlate of human phenotype MIM 605865",
                                                 -species     => $mouse,
                                                 -type        => "homolog",
                                                 -comment     => "type=homolog is putative" );                      
 Function: Creates a new Correlate object.
 Returns : A new Correlate object.
 Args    : -name        => a name or id
           -description => a description
           -species     => the species of this correlating phenotype [Bio::Species]
           -type        => the type of correlation
           -comment     => a comment

=cut

sub new {

    my( $class, @args ) = @_;
    
    my $self = $class->SUPER::new( @args );
 
    my ( $name, $desc, $species, $type, $comment )
    = $self->_rearrange( [ qw( NAME
                               DESCRIPTION
                               SPECIES
                               TYPE
                               COMMENT ) ], @args );
                         
    $self->init();                     
   
    $name    && $self->name( $name );
    $desc    && $self->description( $desc );
    $species && $self->species( $species );
    $type    && $self->type( $type );
    $comment && $self->comment( $comment );
   
    return $self;
    
} # new




=head2 init

 Title   : init()
 Usage   : $co->init();   
 Function: Initializes this Correlate to all "".
 Returns : 
 Args    :

=cut

sub init {

    my( $self ) = @_;

    $self->name( "" );
    $self->description( "" );
    my $species = Bio::Species->new();
    $species->classification( qw( species Undetermined ) );
    $self->species( $species );
    $self->type( "" );
    $self->comment( "" );
  
} # init




=head2 name

 Title   : name
 Usage   : $co->name( "4(Tas1r3)" );
           or
           print $co->name();
 Function: Set/get for the name or id of this Correlate.
 Returns : The name or id of this Correlate.
 Args    : The name or id of this Correlate (optional).

=cut

sub name {
    my ( $self, $value ) = @_;

    if ( defined $value ) {
        $self->{ "_name" } = $value;
    }

    return $self->{ "_name" };

} # name




=head2 description

 Title   : description
 Usage   : $co->description( "mouse correlate of human phenotype MIM 03923" );
           or
           print $co->description();
 Function: Set/get for the description of this Correlate.
 Returns : A description of this Correlate.
 Args    : A description of this Correlate (optional).

=cut

sub description {
    my ( $self, $value ) = @_;

    if ( defined $value ) {
        $self->{ "_description" } = $value;
    }

    return $self->{ "_description" };

} # description




=head2 species

 Title   : species
 Usage   : $co->species( $species );
           or
           $species = $co->species();
 Function: Set/get for the species of this Correlate.
 Returns : The Bio::Species of this Correlate [Bio::Species].
 Args    : The Bio::Species of this Correlate [Bio::Species] (optional).

=cut

sub species {

    my ( $self, $value )  = @_;

    if ( defined $value ) {
        $self->_check_ref_type( $value, "Bio::Species" );
        $self->{ "_species" } = $value;
    }
    
    return $self->{ "_species" };
    
} # species




=head2 type

 Title   : type
 Usage   : $co->type( "homolog" );
           or
           print $co->type();
 Function: Set/get for the type of this Correlate.
 Returns : The type of this Correlate.
 Args    : The type of this Correlate (optional).

=cut

sub type {
    my ( $self, $value ) = @_;

    if ( defined $value ) {
        $self->{ "_type" } = $value;
    }

    return $self->{ "_type" };

} # type




=head2 comment

 Title   : comment
 Usage   : $co->comment( "doubtful" );
           or 
           print $co->comment();
 Function: Set/get for an arbitrary comment about this Correlate.
 Returns : A comment.
 Args    : A comment (optional).

=cut

sub comment {
    my ( $self, $value ) = @_;

    if ( defined $value ) {
        $self->{ "_comment" } = $value;
    }
   
    return $self->{ "_comment" };
    
} # comment



=head2 to_string

 Title   : to_string()
 Usage   : print $co->to_string();
 Function: To string method for Correlate objects.
 Returns : A string representations of this Correlate.
 Args    :

=cut

sub to_string {

    my ( $self ) = @_;

    my $s = "";
    
    $s .= "-- Name:\n";
    $s .= $self->name()."\n";
    $s .= "-- Description:\n";
    $s .= $self->description()."\n";
    $s .= "-- Species:\n";
    $s .= $self->species()->binomial()."\n";
    $s .= "-- Type of correlation:\n";
    $s .= $self->type()."\n";
    $s .= "-- Comment:\n";
    $s .= $self->comment();
  
    return $s;
    
} # to_string




# Title   : _check_ref_type              
# Function: Checks for the correct type.
# Returns : 
# Args    : The value to be checked, the expected class.
sub _check_ref_type {
    my ( $self, $value, $expected_class ) = @_;

    if ( ! defined( $value ) ) {
        $self->throw( ( caller( 1 ) )[ 3 ] .": Found [undef" 
        ."] where [$expected_class] expected" );
    }
    elsif ( ! ref( $value ) ) {
        $self->throw( ( caller( 1 ) )[ 3 ] .": Found scalar"
        ." where [$expected_class] expected" );
    } 
    elsif ( ! $value->isa( $expected_class ) ) {
        $self->throw( ( caller( 1 ) )[ 3 ] .": Found [". ref( $value ) 
        ."] where [$expected_class] expected" );
    }    
} # _check_ref_type



1;
