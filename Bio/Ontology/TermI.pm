# $Id$
#
# BioPerl module for Bio::Ontology::TermI
#
# Cared for by Christian M. Zmasek <czmasek@gnf.org> or <zmasek@yahoo.com>
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

TermI - interface for ontology terms

=head1 SYNOPSIS

#get Bio::Ontology::TermI somehow.
  
  print $term->ontology_id(), "\n";
  print $term->name(), "\n";
  print $term->definition(), "\n";
  print $term->is_obsolete(), "\n";
  print $term->comment(), "\n";
  
  foreach my $alias ( $term->each_alias() ) {
       print $alias, "\n";
  }
    

=head1 DESCRIPTION

This is "dumb" interface for ontology terms providing basic methods
(it provides no functionality related to graphs).

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
  http://bio.perl.org/bioperl-bugs/

=head1 AUTHOR

Christian M. Zmasek

Email: czmasek@gnf.org  or  zmasek@yahoo.com

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

package Bio::Ontology::TermI;
use vars qw( @ISA );
use strict;
use Bio::Root::Object;

use constant TRUE    => 1;
use constant FALSE   => 0;

@ISA = qw( Bio::Root::Root );




=head2 ontology_id

 Title   : ontology_id
 Usage   : $term->ontology_id( "0003947" );
           or
           print $term->ontology_id();
 Function: Set/get for the ontology ID of this Term.
 Returns : The ontology ID [scalar].
 Args    : The ontology ID [scalar] (optional).

=cut

sub ontology_id {
    my ( $self, $value ) = @_;

    if ( defined $value ) {
        $self->{ "_ontology_id" } = $value;
    }

    return $self->{ "_ontology_id" };

} # ontology_id




=head2 name

 Title   : name
 Usage   : $term->name( "N-acetylgalactosaminyltransferase" );
           or
           print $term->name();
 Function: Set/get for the name of this Term.
 Returns : The name [scalar].
 Args    : The name [scalar] (optional).

=cut

sub name {
    my ( $self, $value ) = @_;

    if ( defined $value ) {
        $self->{ "_name" } = $value;
    }

    return $self->{ "_name" };

} # name





=head2 definition

 Title   : definition
 Usage   : $term->definition( "Catalysis of ..." );
           or
           print $term->definition();
 Function: Set/get for the definition of this Term.
 Returns : The definition [scalar].
 Args    : The definition [scalar] (optional).

=cut

sub definition {
    my ( $self, $value ) = @_;

    if ( defined $value ) {
        $self->{ "_definition" } = $value;
    }

    return $self->{ "_definition" };

} # definition




=head2 is_obsolete

 Title   : is_obsolete
 Usage   : $term->is_obsolete( 1 );
           or
           if ( $term->is_obsolete() )
 Function: Set/get for the obsoleteness of this Term.
 Returns : the obsoleteness [0 or 1].
 Args    : the obsoleteness [0 or 1] (optional).

=cut

sub is_obsolete {
    my ( $self, $value ) = @_;

    if ( defined $value ) {
        $self->_is_true_or_false( $value );
        $self->{ "_is_obsolete" } = $value;
    }

    return $self->{ "_is_obsolete" };

} # is_obsolete





=head2 comment

 Title   : comment
 Usage   : $term->comment( "Consider the term ..." );
           or 
           print $term->comment();
 Function: Set/get for an arbitrary comment about this Term.
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




=head2 each_alias

 Title   : each_alias()
 Usage   : @aliases = $term->each_alias();                 
 Function: Returns a list of aliases of this Term.
 Returns : A list of aliases [array of [scalar]].
 Args    :

=cut

sub each_alias {
    my ( $self ) = @_;
    
    if ( $self->{ "_aliases" } ) {
        return @{ $self->{ "_aliases" } };
    }
    else {
        return my @a = (); 
    }
    
} # each_alias


=head2 add_aliases

 Title   : add_aliases
 Usage   : $term->add_aliases( @aliases );
           or
           $term->add_aliases( $alias );                  
 Function: Pushes one or more aliases into the list of aliases.
 Returns : 
 Args    : One alias [scalar] or a list of aliases [array of [scalar]].

=cut

sub add_aliases {
    my ( $self, @values ) = @_;
    
    return unless( @values );
        
    push( @{ $self->{ "_aliases" } }, @values );
    
} # add_aliases


=head2 remove_aliases

 Title   : remove_aliases()
 Usage   : $term->remove_aliases();
 Function: Deletes (and returns) the aliases of this Term.
 Returns : A list of aliases [array of [scalar]].
 Args    :

=cut

sub remove_aliases {
    my ( $self ) = @_;
     
    my @a = $self->each_alias();
    $self->{ "_aliases" } = [];
    return @a;

} # remove_aliases




# Title   :_is_true_or_false
# Function: Checks whether the argument is TRUE or FALSE.
# Returns :
# Args    : The value to be checked.
sub _is_true_or_false {
    my ( $self, $value ) = @_;
    unless ( $value !~ /\D/ && ( $value == TRUE || $value == FALSE ) ) {
        $self->throw( "Found [" . $value
        . "] where " . TRUE . " or " . FALSE . " expected" );
    }
} # _is_true_or_false


1;
