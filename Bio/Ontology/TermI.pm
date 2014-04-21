#
# BioPerl module for Bio::Ontology::Term
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

Bio::Ontology::TermI - interface for ontology terms

=head1 SYNOPSIS

  #get Bio::Ontology::TermI somehow.

  print $term->identifier(), "\n";
  print $term->name(), "\n";
  print $term->definition(), "\n";
  print $term->is_obsolete(), "\n";
  print $term->comment(), "\n";

  foreach my $synonym ( $term->get_synonyms() ) {
      print $synonym, "\n";
  }


=head1 DESCRIPTION

This is "dumb" interface for ontology terms providing basic methods
(it provides no functionality related to graphs).

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
the bugs and their resolution.  Bug reports can be submitted via
the web:

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

package Bio::Ontology::TermI;
use strict;

use base qw(Bio::Root::RootI);



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
    shift->throw_not_implemented();
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
    shift->throw_not_implemented();
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
    shift->throw_not_implemented();
} # definition



=head2 ontology

 Title   : ontology
 Usage   : $ont = $term->ontology();
           or 
           $term->ontology( $ont );
 Function: Get the ontology this term is in.

           An implementation may not permit the value of this
           attribute to be changed once it is set, since that may have
           serious consequences (note that with the ontology in hand
           you can query for all related terms etc).

           Note for implementors: you will almost certainly have to
           take special precaution in order not to create cyclical
           references in memory.

 Returns : The ontology of this Term as a Bio::Ontology::OntologyI
           implementing object.
 Args    : On set, the  ontology of this Term as a Bio::Ontology::OntologyI
           implementing object or a string representing its name.

See L<Bio::Ontology::OntologyI>.

=cut

sub ontology {
    shift->throw_not_implemented();
} # ontology



=head2 version

 Title   : version
 Usage   : $term->version( "1.00" );
           or 
           print $term->version();
 Function: Set/get for version information.
 Returns : The version [scalar].
 Args    : The version [scalar] (optional).

=cut

sub version {
    shift->throw_not_implemented();
} # version




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
    shift->throw_not_implemented();
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
    shift->throw_not_implemented();
} # comment




=head2 get_synonyms

 Title   : get_synonyms
 Usage   : @aliases = $term->get_synonyms();
 Function: Returns a list of aliases of this Term.

           If an implementor of this interface permits modification of
           this array property, the class should define at least
           methods add_synonym() and remove_synonyms(), with obvious
           functionality.

 Returns : A list of aliases [array of [scalar]].
 Args    :

=cut

sub get_synonyms {
    shift->throw_not_implemented();
} # get_synonyms

=head2 get_dblinks

 Title   : get_dblinks()
 Usage   : @ds = $term->get_dblinks();
 Function: Returns a list of each dblink of this term.

           If an implementor of this interface permits modification of
           this array property, the class should define at least
           methods add_dblink() and remove_dblinks(), with obvious
           functionality.

 Returns : A list of dblinks [array of [scalars]].
 Args    : 
 Note    : This has been deprecated in favor of get_dbxrefs()

=cut

sub get_dblinks {
    shift->throw('get_dblinks() is deprecated, use get_dbxrefs() instead');
} # get_dblinks

=head2 get_dbxrefs

 Title   : get_dbxrefs()
 Usage   : @ds = $term->get_dbxrefs();
 Function: Returns a list of each link for this term.

           If an implementor of this interface permits modification of
           this array property, the class should define at least
           methods add_dbxref() and remove_dbxrefs(), with obvious
           functionality.

 Returns : A list of dblinks.  This can return a mixed 'bag' of scalars and
           L<Bio::Annotation::DBLink> instances, or specific subgroups
           can be returned based on passed arguments
 Args    : implementation-specific

=cut

sub get_dbxrefs {
    shift->throw_not_implemented();
} # get_dblinks

=head2 get_secondary_ids

 Title   : get_secondary_ids
 Usage   : @ids = $term->get_secondary_ids();
 Function: Returns a list of secondary identifiers of this Term.

           Secondary identifiers mostly originate from merging terms,
           or possibly also from splitting terms.

           If an implementor of this interface permits modification of
           this array property, the class should define at least
           methods add_secondary_id() and remove_secondary_ids(), with
           obvious functionality.

 Returns : A list of secondary identifiers [array of [scalar]]
 Args    :

=cut

sub get_secondary_ids {
    shift->throw_not_implemented();
} # get_secondary_ids


=head1  Deprecated methods

Used for looking up the methods that supercedes them.

=cut

=head2 category

 Title   : category
 Usage   :
 Function: This method is deprecated. Use ontology() instead. We provide
           an implementation here that preserves backwards compatibility,
           but if you do not have legacy code using it you should not be
           calling this method.
 Example :
 Returns :
 Args    :

=cut

sub category {
    my $self = shift;

    $self->warn("TermI::category is deprecated and being phased out. ".
		"Use TermI::ontology instead.");

    # called in set mode?
    if(@_) {
	# yes; what is incompatible with ontology() is if we were given
	# a TermI object
	my $arg = shift;
	$arg = $arg->name() if ref($arg) && $arg->isa("Bio::Ontology::TermI");
	return $self->ontology($arg,@_);
    } else {
	# No, called in get mode. This is always incompatible with ontology()
	# since category is supposed to return a TermI.
	my $ont = $self->ontology();
	my $term;
	if(defined($ont)) {
	    $term = Bio::Ontology::Term->new(-name => $ont->name(),
					     -identifier =>$ont->identifier());
	}
	return $term;
    }
} # category

1;
