#
# BioPerl module Bio::DB::Biblio::eutils.pm
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Allen Day <allenday@ucla.edu>
# For copyright and disclaimer see below.

# POD documentation - main docs before the code

=head1 NAME

Bio::DB::Biblio::eutils - Access to PubMed's bibliographic query service

=head1 SYNOPSIS

Do not use this object directly, it is recommended to access it and use
it through the I<Bio::Biblio> module:

  use Bio::Biblio;
  use Bio::Biblio::IO;

  my $biblio = Bio::Biblio->new( -access => 'eutils' );
  $biblio->find("10336996");
  my $xml = $biblio->get_next;
  my $io = Bio::Biblio::IO->new( -data => $xml, 
                                 -format => 'medlinexml' );
  my $article = $io->next_bibref();

=head1 DESCRIPTION

This object contains the real implementation of a Bibliographic Query
Service as defined in L<Bio::DB::BiblioI>.

L<Bio::DB::BiblioI> is not implemented as documented in the interface,
particularly the find() method, which is not compatible with PubMed's
query language.

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
of the bugs and their resolution. Bug reports can be submitted via
the web:

  https://redmine.open-bio.org/projects/bioperl/

=head1 AUTHOR

Allen Day E<lt>allenday@ucla.eduE<gt>

=head1 COPYRIGHT

Copyright (c) 2004 Allen Day, University of California, Los Angeles.

This module is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.

=head1 DISCLAIMER

This software is provided "as is" without warranty of any kind.

=head1 BUGS AND LIMITATIONS

=over

=item *

More testing and debugging needed to ensure that returned citations
are properly transferred even if they contain foreign characters.

=item *

Maximum record count (MAX_RECORDS) returned currently hard coded to
100K.

=item *

Biblio retrieval methods should be more tightly integrated with
L<Bio::Biblio::Ref> and L<Bio::DB::MeSH>.

=back

=head1 SEE ALSO

 Pub Med Help:
 http://eutils.ncbi.nlm.nih.gov/entrez/query/static/help/pmhelp.html

 Entrez Utilities:
 http://eutils.ncbi.nlm.nih.gov/entrez/query/static/eutils_help.html

 Example code:
 examples/biblio/biblio-eutils-example.pl

=head1 APPENDIX

The main documentation details are to be found in
L<Bio::DB::BiblioI>.

Here is the rest of the object methods.  Interface methods first,
followed by internal methods.

=cut

# Let the code begin...


package Bio::DB::Biblio::eutils;
use vars qw($DEFAULT_URN);
use strict;

use LWP::Simple;
use XML::Twig;
use URI::Escape;
use base qw(Bio::Biblio Bio::DB::BiblioI);

our $EFETCH      = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi';
our $ESEARCH     = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi';
our $MAX_RECORDS = 100_000;

# -----------------------------------------------------------------------------

=head2 _initialize

 Usage   : my $obj = Bio::Biblio->new(-access => 'eutils' ...);
           (_initialize is internally called from this constructor)
 Returns : 1 on success
 Args    : none

This is an actual new() method (except for the real object creation
and its blessing which is done in the parent class Bio::Root::Root in
method _create_object).

Note that this method is called always as an I<object> method (never as
a I<class> method) - and that the object who calls this method may
already be partly initiated (from Bio::Biblio::new method); so if you
need to do some tricks with the 'class invocation' you need to change
Bio::Biblio::new method, not this one.

=cut

sub _initialize {
    my ($self, @args) = @_;

    #eutils doesn't need this code, but it doesn't hurt to leave it here... -ad

    # make a hashtable from @args
    my %param = @args;
    @param { map { lc $_ } keys %param } = values %param; # lowercase keys

    # copy all @args into this object (overwriting what may already be
    # there) - changing '-key' into '_key'
    my $new_key;
    foreach my $key (keys %param) {
	($new_key = $key) =~ s/^-/_/;
	$self->{ $new_key } = $param { $key };
    }


    # set up internal data
    $self->twig(XML::Twig->new());

    # finally add default values for those keys who have default value
    # and who are not yet in the object

    #AOK
    return 1;
}

=head2 db

 Title   : db
 Usage   : $obj->db($newval)
 Function: specifies the database to search.  valid values are:

           pubmed, pmc, journals

           it is also possible to add the following, and i will do
           so on request:

           genome, nucleotide, protein, popset, snp, sequence, taxonomy

           pubmed is default.

 Returns : value of db (a scalar)
 Args    : on set, new value (a scalar or undef, optional)


=cut

sub db{
    my($self,$arg) = @_;

    if($arg){
      my %ok = map {$_=>1} qw(pubmed pmc journals);
      if($ok{lc($arg)}){
        $self->{'db'} = lc($arg);
      } else {
        $self->warn("invalid db $arg, keeping value as ".$self->{'db'} || 'pubmed');
      }
    }
    return $self->{'db'};
}


=head1 Methods implementing Bio::DB::BiblioI interface

=head2 get_collection_id

  Title   : get_collection_id
  Usage   : $id = $biblio->get_collection_id();
  Function: returns WebEnv value from ESearch
  Returns : ESearch WebEnv value as a string
  Args    : none


=cut

sub get_collection_id {
   return shift->collection_id();
}

sub get_count {
  return shift->count();
}

sub get_by_id {
  my $self = shift;
  my $id = shift;
  my $db = $self->db || 'pubmed';
  $self->throw("must provide valid ID, not undef") unless defined($id);
  my $xml = get($EFETCH.'?rettype=abstract&retmode=xml&db='.$db.'&id='.$id);
  return $xml;
}

=head2 reset_retrieval

  Title   : reset_retrieval
  Usage   : $biblio->reset_retrieval();
  Function: reset cursor in id list, see cursor()
  Returns : 1
  Args    : none


=cut

sub reset_retrieval {
  shift->cursor(0);
  return 1;
}

=head2 get_next

  Title   : get_next
  Usage   : $xml = $biblio->get_next();
  Function: return next record as xml
  Returns : an xml string
  Args    : none


=cut

sub get_next {
  my $self = shift;

  return unless $self->has_next;

  my $xml = $self->get_by_id( @{ $self->ids }[$self->cursor] );
  $self->cursor( $self->cursor + 1 );

  return $xml;
}

=head2 get_more

  Title   : get_more
  Usage   : $xml = $biblio->get_more($more);
  Function: returns next $more records concatenated
  Returns : a string containing multiple xml documents
  Args    : an integer representing how many records to retrieve


=cut

sub get_more {
  my ($self,$more) = @_;

  my @return = ();

  for(1..$more){
    my $next = $self->get_next();
    last unless $next;
    push @return, $next;
  }

  return \@return;
}

=head2 has_next

  Title   : has_next
  Usage   : $has_next = $biblio->has_next();
  Function: check to see if there are more items to be retrieved
  Returns : 1 on true, undef on false
  Args    : none


=cut

sub has_next {
  my $self = shift;
  return ($self->cursor < $self->count) ? 1 : undef;
}



=head2 find

  Title   : find
  Usage   : $biblio = $biblio->find($pubmed_query_phrase);
  Function: perform a PubMed query using Entrez ESearch
  Returns : a reference to the object on which the method was called
  Args    : a PubMed query phrase.  See
            http://eutils.ncbi.nlm.nih.gov/entrez/query/static/help/pmhelp.html
            for help on how to construct a query.

=cut

sub find {
  my ($self,$query) = @_;

  $query = uri_escape($query);

  my $db = $self->db || 'pubmed';

  my $url = $ESEARCH."?usehistory=y&db=$db&retmax=$MAX_RECORDS&term=$query";

  my $xml = get($url) or $self->throw("couldn't retrieve results from $ESEARCH: $!");

  $self->twig->parse($xml);

  my @ids = map {$_->text} $self->twig->get_xpath('//IdList//Id');
  $self->ids(\@ids);

  ##
  #should we be using the ids, or the count tag?
  ##
  my($count_element)  = $self->twig->get_xpath('//Count');
  if (defined $count_element) {
    my $count = $count_element->text();
    $self->count(scalar(@ids));
  }

  my($retmax_element) = $self->twig->get_xpath('//RetMax');
  if (defined $retmax_element) {
    my $retmax = $retmax_element->text();
  }

  my($querykey_element) = $self->twig->get_xpath('//QueryKey');
  if (defined $querykey_element) {
    $self->query_key($querykey_element->text());
  }

  my($webenv_element) = $self->twig->get_xpath('//WebEnv');
  if (defined $webenv_element) {
    $self->collection_id($webenv_element->text());
  }

  #initialize/reset cursor
  $self->cursor(0);

  return $self;
}


=head2 get_all_ids

  Title   : get_all_ids
  Usage   : @ids = $biblio->get_all_ids();
  Function: return a list of PubMed ids resulting from call to find()
  Returns : a list of PubMed ids, or an empty list
  Args    : none


=cut

sub get_all_ids {
  my $self = shift;
  return $self->ids() if $self->ids();
  return ();
}

=head2 get_all

  Title   : get_all
  Usage   : $xml = $biblio->get_all();
  Function: retrieve all records from query
  Returns : return a large concatenated string of PubMed xml documents
  Args    : none


=cut

sub get_all {
  my ($self) = shift;

  my $db = $self->db || 'pubmed';

  my $xml = get($EFETCH.'?rettype=abstract&retmode=xml&db=pubmed&query_key='.
                $self->query_key.'&WebEnv='.$self->collection_id.
                '&retstart=1&retmax='.$MAX_RECORDS
               );

  return $xml;
}

=head2 exists

  Title   : exists
  Usage   : do not use
  Function: no-op.  this is here only for interface compatibility
  Returns : undef
  Args    : none


=cut

sub exists {
  return;

}

=head2 destroy

  Title   : destroy
  Usage   : do not use
  Function: no-op.  this is here only for interface compatibility
  Returns : undef
  Args    : none


=cut

sub destroy {
  return;

}

=head2 get_vocabulary_names

  Title   : get_vocabulary_names
  Usage   : do not use
  Function: no-op.  this is here only for interface compatibility
  Returns : empty arrayref
  Args    : none


=cut

sub get_vocabulary_names {
  return [];
}

=head2 contains

  Title   : contains
  Usage   : do not use
  Function: no-op.  this is here only for interface compatibility
  Returns : undef
  Args    : none


=cut

sub contains {
  return;
}

=head2 get_entry_description

  Title   : get_entry_description
  Usage   : do not use
  Function: no-op.  this is here only for interface compatibility
  Returns : undef
  Args    : none


=cut

sub get_entry_description {
  return;
}

=head2 get_all_values

  Title   : get_all_values
  Usage   : do not use
  Function: no-op.  this is here only for interface compatibility
  Returns : undef
  Args    : none


=cut

sub get_all_values {
  return;
}

=head2 get_all_entries

  Title   : get_all_entries
  Usage   : do not use
  Function: no-op.  this is here only for interface compatibility
  Returns : undef
  Args    : none


=cut

sub get_all_entries {
  return;
}

=head1 Internal methods unrelated to Bio::DB::BiblioI

=head2 cursor

  Title   : cursor
  Usage   : $obj->cursor($newval)
  Function: holds position in reference collection
  Returns : value of cursor (a scalar)
  Args    : on set, new value (a scalar or undef, optional)


=cut

sub cursor {
  my $self = shift;
  my $arg  = shift;

  return $self->{'cursor'} = $arg if defined($arg);
  return $self->{'cursor'};
}

=head2 twig

  Title   : twig
  Usage   : $obj->twig($newval)
  Function: holds an XML::Twig instance.
  Returns : value of twig (a scalar)
  Args    : on set, new value (a scalar or undef, optional)


=cut

sub twig {
  my $self = shift;

  return $self->{'twig'} = shift if @_;
  return $self->{'twig'};
}

=head2 ids

  Title   : ids
  Usage   : $obj->ids($newval)
  Function: store pubmed ids resulting from find() query
  Returns : value of ids (a scalar)
  Args    : on set, new value (a scalar or undef, optional)


=cut

sub ids {
  my $self = shift;

  return $self->{'ids'} = shift if @_;
  return $self->{'ids'};
}

=head2 collection_id

  Title   : collection_id
  Usage   : $obj->collection_id($newval)
  Function:
  Returns : value of collection_id (a scalar)
  Args    : on set, new value (a scalar or undef, optional)


=cut

sub collection_id {
  my $self = shift;

  return $self->{'collection_id'} = shift if @_;
  return $self->{'collection_id'};
}

=head2 count

  Title   : count
  Usage   : $obj->count($newval)
  Function:
  Returns : value of count (a scalar)
  Args    : on set, new value (a scalar or undef, optional)


=cut

sub count {
  my $self = shift;

  return $self->{'count'} = shift if @_;
  return $self->{'count'};
}

=head2 query_key

  Title   : query_key
  Usage   : $obj->query_key($newval)
  Function: holds query_key from ESearch document
  Returns : value of query_key (a scalar)
  Args    : on set, new value (a scalar or undef, optional)


=cut

sub query_key {
  my $self = shift;

  return $self->{'query_key'} = shift if @_;
  return $self->{'query_key'};
}


1;
