# $Id$
#
# BioPerl module for Bio::Annotation::OntologyTerm
#
# Cared for by Hilmar Lapp <hlapp at gmx.net>
#
# Copyright Hilmar Lapp
#
# You may distribute this module under the same terms as perl itself

# 
# (c) Hilmar Lapp, hlapp at gmx.net, 2002.
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

Bio::Annotation::OntologyTerm - An ontology term adapted to AnnotationI

=head1 SYNOPSIS

   use Bio::Annotation::OntologyTerm;
   use Bio::Annotation::Collection;
   use Bio::Ontology::Term;

   my $coll = new Bio::Annotation::Collection;

   # this also implements a tag/value pair, where tag _and_ value are treated
   # as ontology terms
   my $annterm = new Bio::Annotation::OntologyTerm(-label => 'ABC1',
                                                   -tagname => 'Gene Name');
   # ontology terms can be added directly - they implicitly have a tag
   $coll->add_Annotation($annterm);

   # implementation is by composition - you can get/set the term object
   # e.g.
   my $term = $annterm->term(); # term is-a Bio::Ontology::TermI
   print "ontology term ",$term->name()," (ID ",$term->identifier(),
         "), category ",$term->category()->name(),"\n";
   $term = Bio::Ontology::Term->new(-name => 'ABC2', -category => 'Gene Name');
   $annterm->term($term);

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
  http://bugzilla.bioperl.org/

=head1 AUTHOR - Hilmar Lapp

Email bioperl-l@bio.perl.org
Email hlapp at gmx.net


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
use Bio::Ontology::Term;
use Bio::Root::Root;

@ISA = qw(Bio::Root::Root Bio::AnnotationI Bio::Ontology::TermI);

=head2 new

 Title   : new
 Usage   : my $sv = new Bio::Annotation::OntologyTerm;
 Function: Instantiate a new OntologyTerm object
 Returns : Bio::Annotation::OntologyTerm object
 Args    : -term => $term to initialize the term data field [optional]
           Most named arguments that Bio::Ontology::Term accepts will work
           here too. -label is a synonym for -name, -tagname is a synonym for
           -category.

=cut

sub new{
    my ($class,@args) = @_;
    
    my $self = $class->SUPER::new(@args);
    
    my ($term,$name,$label,$identifier,$definition,$cat,$tag) =
	$self->_rearrange([qw(TERM
			      NAME
			      LABEL
			      IDENTIFIER
			      DEFINITION
			      CATEGORY
			      TAGNAME)],
			  @args);
    if($term) {
	$self->term($term);
    } else {
	$self->name($name || $label) if $name || $label;
	$self->identifier($identifier) if $identifier;
	$self->definition($definition) if $definition;
    }
    $self->category($cat || $tag) if $cat || $tag;

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

   return $self->tagname()."|".$self->name()."|".$self->identifier();
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


=head2 tagname

 Title   : tagname
 Usage   : $obj->tagname($newval)
 Function: Get/set the tagname for this annotation value.

           Setting this is optional. If set, it obviates the need to provide
           a tag to AnnotationCollection when adding this object.

           This is aliased to category() here.
 Example : 
 Returns : value of tagname (a scalar)
 Args    : new value (a scalar, optional)


=cut

sub tagname{
    my $self = shift;

    return $self->category(@_) if @_;
    # if in get mode we need to get the name from the category term
    my $cat = $self->category();
    return ref($cat) ? $cat->name() : $cat;
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

=head2 category

 Title   : category
 Usage   : $term->category( $top );
           or 
           $top = $term->category();
 Function: Set/get for a relationship between this Term and
           another Term (e.g. the top level of the ontology).
 Returns : The category of this Term [TermI].
 Args    : The category of this Term [TermI or scalar -- which
           becomes the name of the catagory term] (optional).

=cut

sub category {
    return shift->term()->category(@_);
}

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
