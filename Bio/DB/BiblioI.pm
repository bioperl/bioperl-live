#
# BioPerl module for Bio::DB::BiblioI
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Martin Senger <senger@ebi.ac.uk>
# For copyright and disclaimer see below.
#

# POD documentation - main docs before the code

=head1 NAME

Bio::DB::BiblioI - An interface to a Bibliographic Query Service

=head1 SYNOPSIS

This is an interface module - you do not instantiate it.
Use I<Bio::Biblio> module:

  use Bio::Biblio;
  my $biblio = Bio::Biblio->new(@args);

=head1 DESCRIPTION

This interface describes the methods for accessing a bibliographic
repository, for querying it and for retrieving citations from it. The
retrieved citations are in XML format and can be converted to perl
objects using I<Bio::Biblio::IO>.

The interface complies (with some simplifications) with the
specification described in the B<OpenBQS> project. Its home page is at
http://www.ebi.ac.uk/~senger/openbqs/.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

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
of the bugs and their resolution. Bug reports can be submitted via the
web:

  https://redmine.open-bio.org/projects/bioperl/

=head1 AUTHOR

Martin Senger (martin.senger@gmail.com)

=head1 COPYRIGHT

Copyright (c) 2002 European Bioinformatics Institute. All Rights Reserved.

This module is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.

=head1 DISCLAIMER

This software is provided "as is" without warranty of any kind.

=head1 APPENDIX

This is actually the main documentation...

If you try to call any of these methods directly on this
Bio::DB::BiblioI object you will get a I<not implemented> error
message. You need to call them on a Bio::Biblio object.

=cut


# Let the code begin...

package Bio::DB::BiblioI;
use strict;

use base qw(Bio::Root::RootI);

# -----------------------------------------------------------------------------

=head2 get_collection_id

 Usage   : my $collection_id = $biblio->get_collection_id;
 Returns : string identifying a query collection
           represented by the $biblio object
 Args    : none

Every query collection is uniquely identify-able by its collection
ID. The returned value can be used to populate another $biblio object
and then to access that collection.

=cut

sub get_collection_id {
   my ($self,@args) = @_;
   $self->throw_not_implemented();
}


# -----------------------------------------------------------------------------

=head2 get_count

 Usage   : my $count = $biblio->get_count;
 Returns : integer
 Args    : none, or a string identifying a query collection

It returns a number of citations in the query collection represented
by the calling $biblio object, or in the collection whose ID is given
as an argument.

=cut

sub get_count { shift->throw_not_implemented(); }

# -----------------------------------------------------------------------------

=head2 find

 Usage   : my $new_biblio = $biblio->find ($keywords, $attrs);
           my $new_biblio = $biblio->find ('perl', 'abstract');
           my $new_biblio = $biblio->find ( [ 'perl', 'Java' ] );
 Returns : new Bio::Biblio object representing a new query
           collection
 Args    : $keywords - what to look for (mandatory)
            - a comma-delimited list of keywords, or
            - an array reference with keywords as elements
           $attrs - where to look in (optional)
            - a comma-delimited list of attribute names, or
            - an array reference with attribute names as elements

This is the main query method. It looks for the $keywords in a default
set of attributes, or - if $attrs given - only in the given
attributes.

Because it returns a new Bio::Biblio object which can be again queried
it is possible to chain together several invocations:

    $biblio->find ('Brazma')->find ('Robinson')->get_collection_id;

=cut

sub find { shift->throw_not_implemented; }

# -----------------------------------------------------------------------------

# TBD: AFAIK this method is not implemented on the server side.
#      Let's comment it out for the time being...
#sub query { shift->throw_not_implemented(); }

# -----------------------------------------------------------------------------

=head2 reset_retrieval

 Usage   : $biblio->reset_retrieval;
 Returns : nothing
 Args    : none

It sets an iterator stored in the $biblio object back to its
beginning. After this, the retrieval methods I<has_next>, I<get_next>
and I<get_more> start to iterate the underlying query collection
again from its start.

It throws an exception if this object does not represent any query
result (e.i. it does not contain a collection ID). Note that a
collection ID is created automatically when this object was returned
by a I<find> method, or it can be assigned in a constructor using
argument I<-collection_id>.

=cut

sub reset_retrieval { shift->throw_not_implemented; }

# -----------------------------------------------------------------------------

=head2 get_next

 Usage   : my $citation = $biblio->get_next;
 Returns : a citation in an XML format
 Args    : none

It returns the next available citation from the underlying query
collection. It throws an exception if there are no more citations. In
order to avoid this, use it together with the I<has_next> method:

  my $result = $biblio->find ('brazma', 'authors');
  while ( $result->has_next ) {
      print $result->get_next;
  }

It also throws an exception if this object does not represent any
query result - see explanation in the I<reset_retrieval> elsewhere in
this document.

=cut

sub get_next { shift->throw_not_implemented; }

# -----------------------------------------------------------------------------

=head2 get_more

 Usage   : my $r_citations = $biblio->get_more (5);
 Returns : an array reference - each element has a citation
           in an XML format
 Args    : an integer 'how_many' citations to return;
           default is 1 - but it is assigned with warning

It returns the next I<how_many> available citations from the
underlying query collection. It does not throw any exception if
'how_many' is more than currently available - it simply returns
less. However, it throws an exception if used again without calling
first I<reset_retrieval>.

It also throws an exception if this object does not represent any
query result - see explanation in method I<reset_retrieval> elsewhere
in this document.

=cut

sub get_more { shift->throw_not_implemented; }

# -----------------------------------------------------------------------------

=head2 has_next

 Usage   : my $is = $biblio->has_next;
 Returns : 1 or undef
 Args    : none

It returns 1 if there is a next citation available in the underlying
query collection. Otherwise it returns undef.

It throws an exception if this object does not represent any query
result - see explanation in method I<reset_retrieval> elsewhere in
this document.

=cut

sub has_next { shift->throw_not_implemented; }

# -----------------------------------------------------------------------------

=head2 get_all_ids

 Usage   : my $r_ids = $biblio->get_all_ids;
 Returns : an array reference - each element has
           a citation identifier
 Args    : none

The identifiers of all citations in the underlying query collection
are returned. A usual pattern is to use them then in the I<get_by_id>
method:

    my $biblio = $repository->find ('brazma')->find ('robinson');
    foreach my $id ( @{ $biblio->get_all_ids } ) {
        print $biblio->get_by_id ($id);
    }

It throws an exception if this object does not represent any query
result - see explanation in method I<reset_retrieval> elsewhere in
this document.

=cut

sub get_all_ids { shift->throw_not_implemented; }

# -----------------------------------------------------------------------------

=head2 get_by_id

 Usage   : my $citation = $biblio->get_by_id ('12368254');
 Returns : a citation in an XML format
 Args    : a citation identifier (PMID for Medline)

It returns a citation - disregarding if the citation is or is not in
the underlying query collection (of course, it must be in the
repository).

=cut

sub get_by_id { shift->throw_not_implemented; }

# -----------------------------------------------------------------------------

=head2 get_all

 Usage   : my $all = $biblio->get_all;
 Returns : a (big) string with all citations in an XML format
 Args    : none

It returns an XML valid string (which means that individual citations
are also surrounded by a "set" XML tag) representing all citations
from the underlying query collection.

Note that some servers may limit the number of citations which can be
returned by this method. In such case you need either to refine
further your query collection (using I<find> method) or to retrieve
results by iteration (methods I<has_next>, I<get_next>, I<get_more>).

It throws an exception if this object does not represent any query
result - see explanation in method I<reset_retrieval> elsewhere in
this document.

=cut

sub get_all { shift->throw_not_implemented; }

# -----------------------------------------------------------------------------

=head2 exists

 Usage   : my $exists = $biblio->exists;
 Returns : 1 or undef
 Args    : none

It returns 1 if the underlying query collection represented by the
$biblio object still exists (on the server side).

If you have a collection ID (e.g. stored or printed in a previous
session) but you do not have anymore a C<Bio::Biblio> object representing
it this is how you can check the collection existence:

    use Bio::Biblio;
    print
      Bio::Biblio->new(-collection_id => '1014324148861')->exists;

It throws an exception if this object does not represent any query
result - see explanation in method I<reset_retrieval> elsewhere in
this document.

=cut

sub exists { shift->throw_not_implemented; }

# -----------------------------------------------------------------------------

=head2 destroy

 Usage   : $biblio->destroy;
 Returns : nothing
 Args    : none

It sends a message to the remote server to forget (or free, or destroy
- whatever server choose to do) the query collection represented by
this object.

It throws an exception if this object does not represent any query
collection.

=cut

sub destroy { shift->throw_not_implemented; }

# -----------------------------------------------------------------------------

=head2 get_vocabulary_names

 Usage   : print join ("\n", @{ $biblio->get_vocabulary_names });
 Returns : an array reference - each element has a name
           of a controlled vocabulary
 Args    : none

The controlled vocabularies allow to introspect bibliographic
repositories and to find what citation resource types (such as journal
and book articles, patents or technical reports) are provided by the
repository, what attributes they have, eventually what attribute
values are allowed.

This method returns names of all available controlled
vocabularies. The names can than be used in other methods dealing with
vocabularies: I<contains>, I<get_entry_description>,
I<get_all_values>, and I<get_all_entries>.

=cut

sub get_vocabulary_names { shift->throw_not_implemented; }

# -----------------------------------------------------------------------------

=head2 contains

 Usage   : my $yes = $biblio->contains ($vocabulary_name, $value);
 Returns : 1 or undef
 Args    : $vocabulary_name defines a vocabulary where to look,
           and a $value defines what to look for

It returns 1 if the given controlled vocabulary contains the given
value.

For example, when you know, that a vocabulary
C<MEDLINE/JournalArticle/properties> contains value C<COUNTRY> you can
use it in the I<find> method:

    $biblio->find ('United States', 'COUNTRY');

=cut

sub contains { shift->throw_not_implemented; }

# -----------------------------------------------------------------------------

=head2 get_entry_description

 Usage   : $biblio->get_entry_description ($voc_name, $value);
 Returns : a string with a desciption
 Args    : $voc_name defines a vocabulary where to look,
           and a $value defines whose description to return

Each vocabulary entry has its value (mandatory attribute), and can
have a description (optional attribute). The description may be just a
human readable explanation of an attribute, or it can have more exact
meaning. For example, the server implementation of the bibliographic
query service provided by the EBI puts into attribute descriptions
words I<queryable> and/or I<retrievable> to distinguish the role of
the attributes.

It throws an exception if either vocabulary or value do not exist.

=cut

sub get_entry_description { shift->throw_not_implemented; }

# -----------------------------------------------------------------------------

=head2 get_all_values

 Usage   : $biblio->get_all_values ($vocabulary_name);
 Returns : an array reference - each element has a value (scalar)
           from the given controlled vocabulary
 Args    : $vocabulary_name defines a vocabulary whose values
           are being returned

It returns all values of the given vocabulary.  It throws an exception
if the vocabulary does not exist.

=cut

sub get_all_values { shift->throw_not_implemented; }

# -----------------------------------------------------------------------------

=head2 get_all_entries

 Usage   : $biblio->get_all_entries ($vocabulary_name);
 Returns : a hash reference - keys are vocabulary values
           and values are their descriptions
 Args    : $vocabulary_name defines a vocabulary whose entries
           are being returned

It returns pairs of values and their descriptions of the whole
vocabulary. It throws an exception if the vocabulary does not exist.

This is one way how to get it and print it:

    my $name = 'MEDLINE2005/JournalArticle/properties';
    use Data::Dumper;
    print Data::Dumper->Dump ( [$biblio->get_all_entries ($name)],
			       ['All entries']);

=cut

sub get_all_entries { shift->throw_not_implemented; }

# -----------------------------------------------------------------------------

=head2 VERSION and Revision

 Usage   : print $Bio::DB::BiblioI::VERSION;
           print $Bio::DB::BiblioI::Revision;

=cut

1;
__END__

