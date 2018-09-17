#
# BioPerl module for Bio::Ontology::GOterm
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

Bio::Ontology::GOterm - representation of GO terms 

=head1 SYNOPSIS

  $term = Bio::Ontology::GOterm->new
    ( -go_id       => "GO:0016847",
      -name        => "1-aminocyclopropane-1-carboxylate synthase",
      -definition  => "Catalysis of ...",
      -is_obsolete => 0,
      -comment     => "" );

  $term->add_definition_references( @refs );
  $term->add_secondary_GO_ids( @ids );
  $term->add_aliases( @aliases );

  foreach my $dr ( $term->each_definition_reference() ) {
      print $dr, "\n";
  }

  # etc.

=head1 DESCRIPTION

This is "dumb" class for GO terms (it provides no functionality 
related to graphs). Implements Bio::Ontology::TermI.

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
the bugs and their resolution.  Bug reports can be submitted via the web:

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

package Bio::Ontology::GOterm;
use strict;

use constant GOID_DEFAULT => "GO:0000000";
use constant TRUE         => 1;
use constant FALSE        => 0;

use base qw(Bio::Ontology::Term);

=head2 new

 Title   : new
 Usage   : $term = Bio::Ontology::GOterm->new( 
       -go_id       => "GO:0016847",
       -name        => "1-aminocyclopropane-1-carboxylate synthase",
       -definition  => "Catalysis of ...",
       -is_obsolete => 0,
       -comment     => "" );                   
 Function: Creates a new Bio::Ontology::GOterm.
 Returns : A new Bio::Ontology::GOterm object.
 Args    : -go_id         => the goid of this GO term [GO:nnnnnnn] 
                             or [nnnnnnn] (nnnnnnn is a zero-padded
                             integer of seven digits)
           -name          => the name of this GO term [scalar]
           -definition    => the definition of this GO term [scalar]  
           -ontology      => the ontology for this term (a
                             Bio::Ontology::OntologyI compliant object)
           -version       => version information [scalar]
           -is_obsolete   => the obsoleteness of this GO term [0 or 1]   
           -comment       => a comment [scalar]

=cut

sub new {

    my( $class,@args ) = @_;
    
    my $self = $class->SUPER::new( @args );
   
    my ( $GO_id )
	= $self->_rearrange( [ qw( GO_ID ) ], @args );
   
    $GO_id && $self->GO_id( $GO_id );
  
                                                    
    return $self;
    
} # new


=head2 init

 Title   : init()
 Usage   : $term->init();   
 Function: Initializes this GOterm to all "" and empty lists.
 Returns : 
 Args    :

=cut

sub init {

    my $self = shift;

    # first call the inherited version to properly chain up the hierarchy
    $self->SUPER::init(@_);

    # then only initialize what we implement ourselves here
    #$self->GO_id( GOID_DEFAULT );
  
} # init


=head2 GO_id

 Title   : GO_id
 Usage   : $term->GO_id( "GO:0003947" );
           or
           print $term->GO_id();
 Function: Set/get for the goid of this GO term.

           This is essentially an alias to identifier(), with added
           format checking.

 Returns : The goid [GO:nnnnnnn].
 Args    : The goid [GO:nnnnnnn] or [nnnnnnn] (nnnnnnn is a
           zero-padded integer of seven digits) (optional).

=cut

sub GO_id {
    my $self = shift;
    my $value;

    if ( @_ ) {
        $value = $self->_check_go_id( shift );
	unshift(@_, $value);
    }

    return $self->identifier( @_ );

} # GO_id


=head2 get_secondary_GO_ids

 Title   : get_secondary_GO_ids
 Usage   : @ids = $term->get_secondary_GO_ids();
 Function: Returns a list of secondary goids of this Term.

           This is aliased to remove_secondary_ids().

 Returns : A list of secondary goids [array of [GO:nnnnnnn]]
           (nnnnnnn is a zero-padded integer of seven digits).
 Args    :

=cut

sub get_secondary_GO_ids {
    return shift->get_secondary_ids(@_);
} # get_secondary_GO_ids


=head2 add_secondary_GO_id

 Title   : add_secondary_GO_id
 Usage   : $term->add_secondary_GO_id( @ids );
           or
           $term->add_secondary_GO_id( $id );                  
 Function: Pushes one or more secondary goids into
           the list of secondary goids.

           This is aliased to remove_secondary_ids().

 Returns : 
 Args    : One secondary goid [GO:nnnnnnn or nnnnnnn] or a list
           of secondary goids [array of [GO:nnnnnnn or nnnnnnn]]
           (nnnnnnn is a zero-padded integer of seven digits).

=cut

sub add_secondary_GO_id {
    return shift->add_secondary_id(@_);
} # add_secondary_GO_id


=head2 remove_secondary_GO_ids

 Title   : remove_secondary_GO_ids()
 Usage   : $term->remove_secondary_GO_ids();
 Function: Deletes (and returns) the secondary goids of this Term.

           This is aliased to remove_secondary_ids().

 Returns : A list of secondary goids [array of [GO:nnnnnnn]]
           (nnnnnnn is a zero-padded integer of seven digits).
 Args    :

=cut

sub remove_secondary_GO_ids {
    return shift->remove_secondary_ids(@_);
} # remove_secondary_GO_ids


=head2 to_string

 Title   : to_string()
 Usage   : print $term->to_string();
 Function: to_string method for GO terms.
 Returns : A string representation of this GOterm.
 Args    :

=cut

sub to_string {
    my( $self ) = @_;

    my $s = "";

    $s .= "-- GO id:\n";
    $s .= ($self->GO_id() || '')."\n";
    $s .= "-- Name:\n";
    $s .= ($self->name() || '') ."\n";
    $s .= "-- Definition:\n";
    $s .= ($self->definition() || '') ."\n";
    $s .= "-- Category:\n";
    if ( defined( $self->ontology() ) ) {
        $s .= $self->ontology()->name()."\n";
    }
    else {
        $s .= "\n";
    }
    $s .= "-- Version:\n";
    $s .= ($self->version() || '') ."\n";
    $s .= "-- Is obsolete:\n";
    $s .= $self->is_obsolete()."\n";
    $s .= "-- Comment:\n";
    $s .= ($self->comment() || '') ."\n"; 
    $s .= "-- Definition references:\n";
    $s .= $self->_array_to_string( $self->get_dbxrefs() )."\n";
    $s .= "-- Secondary GO ids:\n";
    $s .= $self->_array_to_string( $self->get_secondary_GO_ids() )."\n";
    $s .= "-- Aliases:\n";
    $s .= $self->_array_to_string( $self->get_synonyms() );
    
    return $s;
    
} # to_string




# Title   : _check_go_id
# Function: Checks whether the argument is [GO:nnnnnnn].
#           If "GO:" is not present, it adds it.
# Returns : The canonical GO id.
# Args    : The value to be checked.
sub _check_go_id {
    my ( $self, $value ) = @_;
    unless ( $value =~ /^(GO:)?\d{7}$/ || $value eq GOID_DEFAULT ) {
        $self->throw( "Found [" . $value
        . "] where [GO:nnnnnnn] or [nnnnnnn] expected" );
    } 
    unless ( $value =~ /^GO:/ ) {
        $value = "GO:".$value;
    }
    return $value;
} # _check_go_id



# Title   : _array_to_string         
# Function:
# Returns : 
# Args    : 
sub _array_to_string {
    my( $self, @value ) = @_;

    my $s = "";
    
    for ( my $i = 0; $i < scalar( @value ); ++$i ) {
        if ( ! ref( $value[ $i ] ) ) {
            $s .= "#" . $i . "\n--  " . $value[ $i ] . "\n";
        }
    }
    
    return $s;
    
} # _array_to_string

#################################################################
# aliases or forwards to maintain backward compatibility
#################################################################

*each_secondary_GO_id = \&get_secondary_GO_ids;
*add_secondary_GO_ids = \&add_secondary_GO_id;

1;
