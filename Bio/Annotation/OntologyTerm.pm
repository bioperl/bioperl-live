# $Id$
#
# BioPerl module for Bio::Annotation::OntologyTerm
#
# Cared for by bioperl <bioperl-l@bio.perl.org>
#
# Copyright bioperl
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Annotation::OntologyTerm - An ontology term adapted to AnnotationI

=head1 SYNOPSIS

   use Bio::Annotation::OntologyTerm;
   use Bio::Annotation::Collection;

   my $col = new Bio::Annotation::Collection;
   my $sv = new Bio::Annotation::OntologyTerm(-value => 'someval');   
   $col->add_Annotation('tagname', $sv);

=head1 DESCRIPTION

Ontology term annotation object 

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists. Your participation is much appreciated.

  bioperl-l@bioperl.org              - General discussion
  http://bio.perl.org/MailList.html  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via email
or the web:

  bioperl-bugs@bioperl.org
  http://bio.perl.org/bioperl-bugs/

=head1 AUTHOR - bioperl

Email bioperl-l@bio.perl.org

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Annotation::OntologyTerm;
use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::Root::Root

use Bio::AnnotationI;
use Bio::Ontology::TermI;
use Bio::Root::Root;

@ISA = qw(Bio::Root::Root Bio::AnnotationI Bio::Ontology::TermI);

=head2 new

 Title   : new
 Usage   : my $sv = new Bio::Annotation::OntologyTerm;
 Function: Instantiate a new OntologyTerm object
 Returns : Bio::Annotation::OntologyTerm object
 Args    : -term => $term to initialize the term data field [optional]

=cut

sub new{
   my ($class,@args) = @_;

   my $self = $class->SUPER::new(@args);

   my ($term) = $self->_rearrange([qw(TERM)], @args);
   $self->term($term) if $term;
   return $self;
}


=head1 AnnotationI implementing functions

=cut

=head2 as_text

 Title   : as_text
 Usage   : my $text = $obj->as_text
 Function: return the string "Name: $v" where $v is the name of the term
 Returns : string
 Args    : none


=cut

sub as_text{
   my ($self) = @_;

   return "Name: ".$self->name();
}

=head2 hash_tree

 Title   : hash_tree
 Usage   : my $hashtree = $value->hash_tree
 Function: For supporting the AnnotationI interface just returns the value
           as a hashref with the key 'value' pointing to the value
 Returns : hashrf
 Args    : none


=cut

sub hash_tree{
   my ($self) = @_;
   
   my $h = {};
   $h->{'name'} = $self->name();
   $h->{'identifier'} = $self->identifier();
   $h->{'definition'} = $self->definition();
   $h->{'synonyms'} = [$self->each_synonym()];
}


=head1 Methods for Bio::Ontology::TermI compliance

=cut

=head2 term

 Title   : term
 Usage   : $obj->term($newval)
 Function: Get/set the Bio::Ontology::TermI implementing object.

           We implement TermI by composition, and this method sets/gets the
           object we delegate to.
 Example : 
 Returns : value of term (a Bio::Ontology::TermI compliant object)
 Args    : new value (a Bio::Ontology::TermI compliant object, optional)


=cut

sub term{
    my ($self,$value) = @_;
    if( defined $value) {
	$self->{'term'} = $value;
    }
    if(! exists($self->{'term'})) {
	$self->{'term'} = Bio::Ontology::Term->new();
    }
    return $self->{'term'};
}

=head2 identifier

 Title   : identifier
 Usage   : $term->identifier( "0003947" );
           or
           print $term->identifier();
 Function: Set/get for the identifier of this Term.
 Returns : The identifier [scalar].
 Args    : The identifier [scalar] (optional).

=cut

sub identifier {
    return shift->term()->identifier(@_);
} # identifier

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
    return shift->term()->name(@_);
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
    return shift->term()->definition(@_);
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
    return shift->term()->is_obsolete(@_);
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
    return shift->term()->comment(@_);
} # comment

=head2 each_synonym

 Title   : each_synonym()
 Usage   : @aliases = $term->each_synonym();                 
 Function: Returns a list of aliases of this Term.
 Returns : A list of aliases [array of [scalar]].
 Args    :

=cut

sub each_synonym {
    return shift->term()->each_synonym(@_);
} # each_synonym

=head2 add_synonyms

 Title   : add_synonyms
 Usage   : $term->add_synonyms( @asynonyms );
           or
           $term->add_synonyms( $synonym );                  
 Function: Pushes one or more synonyms into the list of synonyms.
 Returns : 
 Args    : One synonym [scalar] or a list of synonyms [array of [scalar]].

=cut

sub add_synonyms {
    return shift->term()->add_synonyms(@_);
} # add_synonyms


=head2 remove_synonyms

 Title   : remove_synonyms()
 Usage   : $term->remove_synonyms();
 Function: Deletes (and returns) the synonyms of this Term.
 Returns : A list of synonyms [array of [scalar]].
 Args    :

=cut

sub remove_synonyms {
    return shift->term()->remove_synonyms(@_);
} # remove_synonyms


1;
