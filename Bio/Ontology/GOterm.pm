# $Id$
#
# BioPerl module for Bio::Ontology::GOterm
#
# Cared for by Christian M. Zmasek <czmasek@gnf.org> or <cmzmasek@yahoo.com>
#
# (c) Christian M. Zmasek, czmasek@gnf.org, 2002.
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

GOterm - representation of GO terms 

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

This is "dumb" class for GO terms (it provides no functionality related to graphs).
Implements Bio::Ontology::TermI.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists.  Your participation is much appreciated.

  bioperl-l@bioperl.org             - General discussion
  http://bio.perl.org/MailList.html - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via email
or the web:

  bioperl-bugs@bio.perl.org
  http://bugzilla.bioperl.org/

=head1 AUTHOR

Christian M. Zmasek

Email: czmasek@gnf.org  or  cmzmasek@yahoo.com

WWW:   http://www.genetics.wustl.edu/eddy/people/zmasek/

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
use vars qw( @ISA );
use strict;
use Bio::Ontology::Term;

use constant GOID_DEFAULT => "GO:-------";
use constant TRUE         => 1;
use constant FALSE        => 0;

@ISA = qw( Bio::Ontology::Term );




=head2 new

 Title   : new
 Usage   : $term = Bio::Ontology::GOterm->new( -go_id       => "GO:0016847",
                                               -name        => "1-aminocyclopropane-1-carboxylate synthase",
                                               -definition  => "Catalysis of ...",
                                               -is_obsolete => 0,
                                               -comment     => "" );                   
 Function: Creates a new Bio::Ontology::GOterm.
 Returns : A new Bio::Ontology::GOterm object.
 Args    : -go_id                 => the goid of this GO term [GO:nnnnnnn] 
                                     or [nnnnnnn] (nnnnnnn is a zero-padded
                                     integer of seven digits)
           -name                  => the name of this GO term [scalar]
           -definition            => the definition of this GO term [scalar]  
           -category              => a relationship between this Term and another Term [TermI or scalar]
           -version               => version information [scalar]
           -is_obsolete           => the obsoleteness of this GO term [0 or 1]   
           -comment               => a comment [scalar]

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

    $self->SUPER::init(@_);

    $self->GO_id( GOID_DEFAULT );
    $self->remove_dblinks();
    $self->remove_secondary_GO_ids();
    $self->remove_synonyms();
  
} # init




=head2 GOid

 Title   : GO_id
 Usage   : $term->GO_id( "GO:0003947" );
           or
           print $term->GO_id();
 Function: Set/get for the goid of this GO term.
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





=head2 each_dblink

 Title   : each_dblink()
 Usage   : @ds = $term->each_dblink();                 
 Function: Returns a list of each dblinks of this GO term.
 Returns : A list of dblinks [array of [scalars]].
 Args    :

=cut

sub each_dblink {
    my ( $self ) = @_;
    
    if ( $self->{ "_dblinks" } ) {
        return @{ $self->{ "_dblinks" } };
    }
    else {
        return my @a = (); 
    }
    
} # each_dblink


=head2 add_dblinks

 Title   : add_dblinks
 Usage   : $term->add_dblinks( @dbls );
           or
           $term->add_dblinks( $dbl );                  
 Function: Pushes one or more dblinks
           into the list of dblinks.
 Returns : 
 Args    : One  dblink [scalar] or a list of
            dblinks [array of [scalars]].

=cut

sub add_dblinks {
    my ( $self, @values ) = @_;
    
    return unless( @values );
        
    push( @{ $self->{ "_dblinks" } }, @values );
    
} # add_dblinks


=head2 remove_dblinks

 Title   : remove_dblinks()
 Usage   : $term->remove_dblinks();
 Function: Deletes (and returns) the definition references of this GO term.
 Returns : A list of definition references [array of [scalars]].
 Args    :

=cut

sub remove_dblinks {
    my ( $self ) = @_;
     
    my @a = $self->each_dblink();
    $self->{ "_dblinks" } = [];
    return @a;

} # remove_dblinks




=head2 each_secondary_GO_id

 Title   : each_secondary_GO_id()
 Usage   : @ids = $term->each_secondary_GO_id();                 
 Function: Returns a list of secondary goids of this Term.
 Returns : A list of secondary goids [array of [GO:nnnnnnn]]
           (nnnnnnn is a zero-padded integer of seven digits).
 Args    :

=cut

sub each_secondary_GO_id {
    my ( $self ) = @_;
    
    if ( $self->{ "_secondary_GO_ids" } ) {
        return @{ $self->{ "_secondary_GO_ids" } };
    }
    else {
        return my @a = (); 
    }
    
} # each_secondary_GO_id


=head2 add_secondary_GO_ids

 Title   : add_secondary_GO_ids
 Usage   : $term->add_secondary_GO_ids( @ids );
           or
           $term->add_secondary_GO_ids( $id );                  
 Function: Pushes one or more secondary goids into
           the list of secondary goids.
 Returns : 
 Args    : One secondary goid [GO:nnnnnnn or nnnnnnn] or a list
           of secondary goids [array of [GO:nnnnnnn or nnnnnnn]]
           (nnnnnnn is a zero-padded integer of seven digits).

=cut

sub add_secondary_GO_ids {
    my ( $self, @values ) = @_;
    
    return unless( @values );
    
    foreach my $value ( @values ) {  
        $value = $self->_check_go_id( $value );
    }
    
    push( @{ $self->{ "_secondary_GO_ids" } }, @values );
    
} # add_secondary_GO_ids


=head2 remove_secondary_GO_ids

 Title   : remove_secondary_GO_ids()
 Usage   : $term->remove_secondary_GO_ids();
 Function: Deletes (and returns) the secondary goids of this Term.
 Returns : A list of secondary goids [array of [GO:nnnnnnn]]
           (nnnnnnn is a zero-padded integer of seven digits).
 Args    :

=cut

sub remove_secondary_GO_ids {
    my ( $self ) = @_;
     
    my @a = $self->each_secondary_GO_id();
    $self->{ "_secondary_GO_ids" } = [];
    return @a;

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
    $s .= $self->GO_id()."\n";
    $s .= "-- Name:\n";
    $s .= $self->name() || '' ."\n";
    $s .= "-- Definition:\n";
    $s .= $self->definition() || '' ."\n";
    $s .= "-- Category:\n";
    if ( defined( $self->category() ) ) {
        $s .= $self->category()->name()."\n";
    }
    else {
        $s .= "\n";
    }
    $s .= "-- Version:\n";
    $s .= $self->version() || '' ."\n";
    $s .= "-- Is obsolete:\n";
    $s .= $self->is_obsolete()."\n";
    $s .= "-- Comment:\n";
    $s .= $self->comment() || '' ."\n"; 
    $s .= "-- Definition references:\n";
    $s .= $self->_array_to_string( $self->each_dblink() )."\n";
    $s .= "-- Secondary GO ids:\n";
    $s .= $self->_array_to_string( $self->each_secondary_GO_id() )."\n";
    $s .= "-- Aliases:\n";
    $s .= $self->_array_to_string( $self->each_synonym() );
    
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


1;
