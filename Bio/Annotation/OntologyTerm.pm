#
# BioPerl module for Bio::Annotation::OntologyTerm
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
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

   my $coll = Bio::Annotation::Collection->new();

   # this also implements a tag/value pair, where tag _and_ value are treated
   # as ontology terms
   my $annterm = Bio::Annotation::OntologyTerm->new(-label => 'ABC1',
                                                   -tagname => 'Gene Name');
   # ontology terms can be added directly - they implicitly have a tag
   $coll->add_Annotation($annterm);

   # implementation is by composition - you can get/set the term object
   # e.g.
   my $term = $annterm->term(); # term is-a Bio::Ontology::TermI
   print "ontology term ",$term->name()," (ID ",$term->identifier(),
         "), ontology ",$term->ontology()->name(),"\n";
   $term = Bio::Ontology::Term->new(-name => 'ABC2',
                                    -ontology => 'Gene Name');
   $annterm->term($term);

=head1 DESCRIPTION

Ontology term annotation object

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists. Your participation is much appreciated.

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
the bugs and their resolution.  Bug reports can be submitted via
the web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Hilmar Lapp

Email hlapp at gmx.net


=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Annotation::OntologyTerm;
use strict;

# Object preamble - inherits from Bio::Root::Root

use Bio::Ontology::Term;

use base qw(Bio::Root::Root Bio::AnnotationI Bio::Ontology::TermI);

=head2 new

 Title   : new
 Usage   : my $sv = Bio::Annotation::OntologyTerm->new();
 Function: Instantiate a new OntologyTerm object
 Returns : Bio::Annotation::OntologyTerm object
 Args    : -term => $term to initialize the term data field [optional]
           Most named arguments that Bio::Ontology::Term accepts will work
           here too. -label is a synonym for -name, -tagname is a synonym for
           -ontology.

=cut

sub new{
    my ($class,@args) = @_;

    my $self = $class->SUPER::new(@args);
    my ($term,$name,$label,$identifier,$definition,$ont,$tag) =
	$self->_rearrange([qw(TERM
                          NAME
                          LABEL
                          IDENTIFIER
                          DEFINITION
                          ONTOLOGY
                          TAGNAME)],
                      @args);
    if($term) {
        $self->term($term);
    } else {
        $self->name($name || $label) if $name || $label;
        $self->identifier($identifier) if $identifier;
        $self->definition($definition) if $definition;
    }
    $self->ontology($ont || $tag) if $ont || $tag;
    return $self;
}


=head1 AnnotationI implementing functions

=cut

=head2 as_text

 Title   : as_text
 Usage   : my $text = $obj->as_text
 Function: Returns a textual representation of the annotation that
           this object holds. Presently, it is tag name, name of the
           term, and the is_obsolete attribute concatenated togather
           with a delimiter (|).

 Returns : string
 Args    : none


=cut

sub as_text{
   my ($self) = @_;

   return $self->tagname()."|".$self->name()."|".($self->is_obsolete()||'');
}

=head2 display_text

 Title   : display_text
 Usage   : my $str = $ann->display_text();
 Function: returns a string. Unlike as_text(), this method returns a string
           formatted as would be expected for te specific implementation.

           One can pass a callback as an argument which allows custom text
           generation; the callback is passed the current instance and any text
           returned
 Example :
 Returns : a string
 Args    : [optional] callback

=cut

{
  my $DEFAULT_CB = sub { $_[0]->identifier || ''};

  sub display_text {
    my ($self, $cb) = @_;
    $cb ||= $DEFAULT_CB;
    $self->throw("Callback must be a code reference") if ref $cb ne 'CODE';
    return $cb->($self);
  }

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
   $h->{'synonyms'} = [$self->get_synonyms()];
}


=head2 tagname

 Title   : tagname
 Usage   : $obj->tagname($newval)
 Function: Get/set the tagname for this annotation value.

           Setting this is optional. If set, it obviates the need to provide
           a tag to AnnotationCollection when adding this object.

           This is aliased to ontology() here.
 Example :
 Returns : value of tagname (a scalar)
 Args    : new value (a scalar, optional)


=cut

sub tagname{
    my $self = shift;

    return $self->ontology(@_) if @_;
    # if in get mode we need to get the name from the ontology
    my $ont = $self->ontology();
    return ref($ont) ? $ont->name() : $ont;
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

=head2 ontology

 Title   : ontology
 Usage   : $term->ontology( $top );
           or
           $top = $term->ontology();
 Function: Set/get for a relationship between this Term and
           another Term (e.g. the top level of the ontology).
 Returns : The ontology of this Term [TermI].
 Args    : The ontology of this Term [TermI or scalar -- which
           becomes the name of the catagory term] (optional).

=cut

sub ontology {
    return shift->term()->ontology(@_);
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

=head2 get_synonyms

 Title   : get_synonyms()
 Usage   : @aliases = $term->get_synonyms();
 Function: Returns a list of aliases of this Term.
 Returns : A list of aliases [array of [scalar]].
 Args    :

=cut

sub get_synonyms {
    return shift->term()->get_synonyms(@_);
} # get_synonyms

=head2 add_synonym

 Title   : add_synonym
 Usage   : $term->add_synonym( @asynonyms );
           or
           $term->add_synonym( $synonym );
 Function: Pushes one or more synonyms into the list of synonyms.
 Returns :
 Args    : One synonym [scalar] or a list of synonyms [array of [scalar]].

=cut

sub add_synonym {
    return shift->term()->add_synonym(@_);
} # add_synonym


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

=head2 get_dblinks

 Title   : get_dblinks()
 Usage   : @ds = $term->get_dblinks();
 Function: Returns a list of each dblinks of this GO term.
 Returns : A list of dblinks [array of [scalars]].
 Args    :
 Note    : this is deprecated in favor of get_dbxrefs(), which works with strings
           or L<Bio::Annotation::DBLink> instances

=cut

sub get_dblinks {
    my $self = shift;
    $self->deprecated('get_dblinks() is deprecated; use get_dbxrefs()');
    return $self->term->get_dbxrefs(@_);
} # get_dblinks

=head2 get_dbxrefs

 Title   : get_dbxrefs()
 Usage   : @ds = $term->get_dbxrefs();
 Function: Returns a list of each dblinks of this GO term.
 Returns : A list of dblinks [array of [scalars] or Bio::Annotation::DBLink instances].
 Args    :

=cut

sub get_dbxrefs {
    return shift->term->get_dbxrefs(@_);
} # get_dblinks

=head2 add_dblink

 Title   : add_dblink
 Usage   : $term->add_dblink( @dbls );
           or
           $term->add_dblink( $dbl );
 Function: Pushes one or more dblinks
           into the list of dblinks.
 Returns :
 Args    : One  dblink [scalar] or a list of
            dblinks [array of [scalars]].
 Note    : this is deprecated in favor of add_dbxref(), which works with strings
           or L<Bio::Annotation::DBLink> instances

=cut

sub add_dblink {
    my $self = shift;
    $self->deprecated('add_dblink() is deprecated; use add_dbxref()');
    return $self->term->add_dbxref(@_);
} # add_dblink

=head2 add_dbxref

 Title   : add_dbxref
 Usage   : $term->add_dbxref( @dbls );
           or
           $term->add_dbxref( $dbl );
 Function: Pushes one or more dblinks
           into the list of dblinks.
 Returns :
 Args    : 

=cut

sub add_dbxref {
    return shift->term->add_dbxref(@_);
} 

=head2 remove_dblinks

 Title   : remove_dblinks()
 Usage   : $term->remove_dblinks();
 Function: Deletes (and returns) the definition references of this GO term.
 Returns : A list of definition references [array of [scalars]].
 Args    :
 Note    : this is deprecated in favor of remove_dbxrefs(), which works with strings
           or L<Bio::Annotation::DBLink> instances

=cut

sub remove_dblinks {
    my $self = shift;
    $self->deprecated('remove_dblinks() is deprecated; use remove_dbxrefs()');
    return $self->term->remove_dbxrefs(@_);
} # remove_dblinks

=head2 remove_dbxrefs

 Title   : remove_dbxrefs()
 Usage   : $term->remove_dbxrefs();
 Function: Deletes (and returns) the definition references of this GO term.
 Returns : A list of definition references [array of [scalars]].
 Args    :

=cut

sub remove_dbxrefs {
    return shift->term->remove_dbxrefs(@_);
} 

=head2 get_secondary_ids

 Title   : get_secondary_ids
 Usage   : @ids = $term->get_secondary_ids();
 Function: Returns a list of secondary identifiers of this Term.

           Secondary identifiers mostly originate from merging terms,
           or possibly also from splitting terms.

 Returns : A list of secondary identifiers [array of [scalar]]
 Args    :

=cut

sub get_secondary_ids {
    return shift->term->get_secondary_ids(@_);
} # get_secondary_ids


=head2 add_secondary_id

 Title   : add_secondary_id
 Usage   : $term->add_secondary_id( @ids );
           or
           $term->add_secondary_id( $id );
 Function: Adds one or more secondary identifiers to this term.
 Returns :
 Args    : One or more secondary identifiers [scalars]

=cut

sub add_secondary_id {
    return shift->term->add_secondary_id(@_);
} # add_secondary_id


=head2 remove_secondary_ids

 Title   : remove_secondary_ids
 Usage   : $term->remove_secondary_ids();
 Function: Deletes (and returns) the secondary identifiers of this Term.
 Returns : The previous list of secondary identifiers [array of [scalars]]
 Args    :

=cut

sub remove_secondary_ids {
    return shift->term->remove_secondary_ids(@_);
} # remove_secondary_ids


1;
