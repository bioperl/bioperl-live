#
# BioPerl module for Bio::Phenotype::OMIM::OMIMentryAllelicVariant
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

Bio::Phenotype::OMIM::OMIMentryAllelicVariant - Representation of a allelic
variant of the OMIM database

=head1 SYNOPSIS

  use Bio::Phenotype::OMIM::OMIMentryAllelicVariant;

  $av = Bio::Phenotype::OMIM::OMIMentryAllelicVariant->new( -number               => ".0001",
                                                            -title                => "ALCOHOL INTOLERANCE",
                                                            -symbol               => "ALDH2*2",
                                                            -description          => "The ALDH2*2-encoded ...",
                                                            -aa_ori               => "GLU",
                                                            -aa_mut               => "LYS",
                                                            -position             => 487,
                                                            -additional_mutations => "IVS4DS, G-A, +1" );

=head1 DESCRIPTION

This class models the allelic variant of the OMIM database.
This class is intended to be used together with a OMIM entry class. 


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


package Bio::Phenotype::OMIM::OMIMentryAllelicVariant;
use strict;

use base qw(Bio::Root::Root);

=head2 new

 Title   : new
 Usage   : $av = Bio::Phenotype::OMIM::OMIMentryAllelicVariant->new( -number               => ".0001",
                                                                     -title                => "ALCOHOL INTOLERANCE",
                                                                     -symbol               => "ALDH2*2",
                                                                     -description          => "The ALDH2*2-encoded ...",
                                                                     -aa_ori               => "GLU",
                                                                     -aa_mut               => "LYS",
                                                                     -position             => 487,
                                                                     -additional_mutations => "IVS4DS, G-A, +1" );                      
 Function: Creates a new OMIMentryAllelicVariant object.
 Returns : A new OMIMentryAllelicVariant object.
 Args    : -number               => the OMIM allelic variant number
           -title                => the title
           -symbol               => a symbol
           -description          => a description
           -aa_ori               => the original amino acid
           -aa_mut               => the mutated amino acid
           -position             => the position of the mutation
           -additional_mutations => free form description of additional mutations

=cut

sub new {

    my( $class, @args ) = @_;
  
    my $self = $class->SUPER::new( @args );
   
    my ( $number, $title, $symbol, $desc, $ori, $mut, $pos, $am )
    = $self->_rearrange( [ qw( NUMBER
                               TITLE
                               SYMBOL
                               DESCRIPTION
                               AA_ORI
                               AA_MUT
                               POSITION
                               ADDITIONAL_MUTATIONS ) ], @args );

    $self->init(); 

    $number && $self->number( $number );
    $title  && $self->title( $title );
    $symbol && $self->symbol( $symbol );
    $desc   && $self->description( $desc );
    $ori    && $self->aa_ori( $ori );
    $mut    && $self->aa_mut( $mut );
    $pos    && $self->position( $pos );
    $am     && $self->additional_mutations( $am );
   
    return $self;

} # new 




=head2 init

 Title   : init()
 Usage   : $av->init();   
 Function: Initializes this OMIMentryAllelicVariant to all "".
 Returns : 
 Args    :

=cut

sub init {
    my( $self ) = @_;

    $self->number( "" );
    $self->title( "" );
    $self->symbol( "" );
    $self->description( "" );
    $self->aa_ori( "" );
    $self->aa_mut( "" );
    $self->position( "" );
    $self->additional_mutations( "" );
    
} # init




=head2 number

 Title   : number
 Usage   : $av->number( ".0001" );
           or
           print $av->number();
 Function: Set/get for the OMIM allelic variant number of this
           OMIMentryAllelicVariant.
 Returns : The OMIM allelic variant number.
 Args    : The OMIM allelic variant number (optional).

=cut

sub number {
    my ( $self, $value ) = @_;

    if ( defined $value ) {
        $self->{ "_number" } = $value;
    }

    return $self->{ "_number" };

} # number



=head2 title

 Title   : title
 Usage   : $av->title( "ALCOHOL INTOLERANCE" );
           or
           print $av->title();
 Function: Set/get for the title of this OMIMentryAllelicVariant.
 Returns : The title.
 Args    : The title (optional).

=cut

sub title {
    my ( $self, $value ) = @_;

    if ( defined $value ) {
        $self->{ "_title" } = $value;
    }

    return $self->{ "_title" };

} # title




=head2 symbol

 Title   : symbol
 Usage   : $av->symbol( "ALDH2*2" );
           or
           print $av->symbol();
 Function: Set/get for the symbol of this OMIMentryAllelicVariant.
 Returns : A symbol.
 Args    : A symbol (optional).

=cut

sub symbol {
    my ( $self, $value ) = @_;

    if ( defined $value ) {
        $self->{ "_symbol" } = $value;
    }

    return $self->{ "_symbol" };

} # symbol




=head2 description

 Title   : description
 Usage   : $av->description( "The ALDH2*2-encoded protein has a change ..." );
           or
           print $av->description();
 Function: Set/get for the description of this OMIMentryAllelicVariant.
 Returns : A description.
 Args    : A description (optional).

=cut

sub description {
    my ( $self, $value ) = @_;

    if ( defined $value ) {
        $self->{ "_description" } = $value;
    }

    return $self->{ "_description" };

} # description




=head2 aa_ori

 Title   : aa_ori
 Usage   : $av->aa_ori( "GLU" );
           or
           print $av->aa_ori();
 Function: Set/get for the original amino acid(s).
 Returns : The original amino acid(s).
 Args    : The original amino acid(s) (optional).

=cut

sub aa_ori {
    my ( $self, $value ) = @_;

    if ( defined $value ) {
        $self->{ "_aa_ori" } = $value;
    }

    return $self->{ "_aa_ori" };

} # aa_ori




=head2 aa_mut

 Title   : aa_mut
 Usage   : $av->aa_mut( "LYS" );
           or
           print $av->aa_mut();
 Function: Set/get for the mutated amino acid(s).
 Returns : The mutated amino acid(s).
 Args    : The mutated amino acid(s) (optional).

=cut

sub aa_mut {
    my ( $self, $value ) = @_;

    if ( defined $value ) {
        $self->{ "_aa_mut" } = $value;
    }

    return $self->{ "_aa_mut" };

} # aa_mut




=head2 position

 Title   : position
 Usage   : $av->position( 487 );
           or
           print $av->position();
 Function: Set/get for the position of the mutation.
 Returns : The position of the mutation.
 Args    : The position of the mutation (optional).

=cut

sub position {
    my ( $self, $value ) = @_;

    if ( defined $value ) {
        $self->{ "_position" } = $value;
    }

    return $self->{ "_position" };

} # position




=head2 additional_mutations

 Title   : additional_mutations
 Usage   : $av->additional_mutations( "1-BP DEL, 911T" );
           or
           print $av->additional_mutations();
 Function: Set/get for free form description of (additional) mutation(s).
 Returns : description of (additional) mutation(s).
 Args    : description of (additional) mutation(s) (optional).

=cut

sub additional_mutations {
    my ( $self, $value ) = @_;

    if ( defined $value ) {
        $self->{ "_additional_mutations" } = $value;
    }

    return $self->{ "_additional_mutations" };

} # additional_mutations



=head2 to_string

 Title   : to_string()
 Usage   : print $av->to_string();
 Function: To string method for OMIMentryAllelicVariant objects.
 Returns : A string representations of this OMIMentryAllelicVariant.
 Args    :

=cut

sub to_string {
    my( $self ) = @_;

    my $s = "";
    
    $s .= "-- Number:\n";
    $s .= $self->number()."\n";
    $s .= "-- Title:\n";
    $s .= $self->title()."\n";
    $s .= "-- Symbol:\n";
    $s .= $self->symbol()."\n";
    $s .= "-- Description:\n";
    $s .= $self->description()."\n";
    $s .= "-- Original AA(s):\n";
    $s .= $self->aa_ori()."\n";
    $s .= "-- Mutated AA(s):\n";
    $s .= $self->aa_mut()."\n";
    $s .= "-- Position:\n";
    $s .= $self->position()."\n";
    $s .= "-- Additional Mutation(s):\n";
    $s .= $self->additional_mutations();
  
    return $s;
 
} # to_string




1;
