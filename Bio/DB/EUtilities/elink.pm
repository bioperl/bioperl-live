# $Id$
#
# BioPerl module for Bio::DB::EUtilities::elink
#
# Cared for by Chris Fields <cjfields at uiuc dot edu>
#
# Copyright Chris Fields
#
# You may distribute this module under the same terms as perl itself
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code
# 

=head1 NAME

Bio::DB::EUtilities::elink - check for and retrieve external or related ID's
from a list of one or more primary ID's, including relevancy scores.

=head1 SYNOPSIS

B<Do not use this module directly.>  Use it via the
L<Bio::DB::EUtilities|Bio::DB::EUtilities> class.

  # chain EUtilities for complex queries

  use Bio::DB::EUtilities;

  my $esearch = Bio::DB::EUtilities->new(-eutil      => 'esearch',
                                         -db         => 'pubmed',
                                         -term       => 'hutP',
                                         -usehistory => 'y');
  
  $esearch->get_response; # parse the response, fetch a cookie
  
  my $elink = Bio::DB::EUtilities->new(-eutil        => 'elink',
                                       -db           => 'protein,taxonomy',
                                       -dbfrom       => 'pubmed',
                                       -cookie       => $esearch->next_cookie,
                                       -cmd          => 'neighbor');
  
  # this retrieves the Bio::DB::EUtilities::ElinkData object
  
  my ($linkset) = $elink->next_linkset;
  my @ids;
  
  # step through IDs for each linked database in the ElinkData object
  
  for my $db ($linkset->get_databases) {   
    @ids = $linkset->get_LinkIds_by_db($db); #returns primary ID's
    # do something here
  }
  
  # multiple ID groups (for one-to-one-correspondence of IDs)

  my $elink = Bio::DB::EUtilities->new(-eutil        => 'elink',
                                       -db           => 'all',
                                       -dbfrom       => 'protein',
                                       -id           => [\@id1, @ids2],
                                       -multi_id     => 1,
                                       -cmd          => 'neighbor');

  for my $linkset ($elink->get_all_linksets) {
    for my $db ($linkset->get_databases) {
      my @ids = $linkset->get_LinkIds_by_db($db); #returns primary ID's
      # do something here
    }
  }
  
  # to retrieve scores for a linkset

  for my $linkset ($elink->get_all_linksets) {
    my @score_dbs = $linkset->has_scores; # retrieve databases with score values
    for my $db (@score_dbs) {
      my @ids = $linkset->get_LinkIds_by_db($db); #returns primary ID's
      $linkset->set_score_db($db);  # to current database containing scores
      for my $id (@ids) {
         my $score = get_score($id);  
         # do something here, like screen for IDs based on score
      }
    }
  }
  
  # or just receive a hash containing ID-score key-value pairs
  
  for my $linkset ($elink->get_all_linksets) {
    my @score_dbs = $linkset->has_scores; 
    for my $db (@score_dbs) {
      $linkset->set_score_db($db);
      %scores = $linkset->get_score_hash;
    }
  }
  
=head1 DESCRIPTION

B<WARNING>: Please do B<NOT> spam the Entrez web server with multiple requests.

The EUtility Elink is used to check for and retrieve external or related ID's
from a list of one or more primary ID's.  Using the C<cmd> parameter, one can
vary the returned data.  See the below command options for explanations on
returned XML output.  For certain command options one can retrieve one or more
L<Bio::DB::EUtilities::Cookie|Bio::DB::EUtilities::Cookie> objects to be used in
other EUtility searches or efetch primary IDs.  Other will return the ID
information and relevancy scores in one or more
L<Bio::DB::EUtilities::ElinkData|Bio::DB::EUtilities::ElinkData> objects.

=head2 Parameters

The following are a general list of parameters that can be used to take
advantage of ELink.  Up-to-date help for ELink is available at this URL
(the information below is a summary of the options found there):

  http://eutils.ncbi.nlm.nih.gov/entrez/query/static/elink_help.html

=over 3

=item C<db>

One or more database available through EUtilities. If set to 'all', will
retrieve all relelvant information from each database based on the C<cmd>
parameter (the default setting is to retrieve related primary ID's).  One
interesting behaviour is when C<db> and C<dbfrom> are set to the same database;
related IDs from database are retrieved along with a relevancy score.  This
score differs from database to database; if protein-protein elinks are sought,
the scores are generated from BLASTP

=item C<dbfrom>

originating database; useful only if using directly when querying with ID's

=item C<id>

a list of primary ID's

Below are a list of IDs which can be used with ELink:

B<PMID> (pubmed), B<MIM number> (omim), B<GI number> (nucleotide, protein),
B<Genome ID> (genome), B<Popset ID> (popset), B<SNP cluster ID> (snp),
B<UniSTS ID> (unists), B<UniGene cluster ID> (unigene), <MMDB-ID> (structure),
B<PSSM-ID> (cdd), B<3D SDI> (domains), B<TAXID> (taxonomy), B<GEO ID> (geo)

=item C<reldate>

limits results to the number of days preceding today's date

=item C<mindate>, C<maxdate>

limits results by dates (C<yyyy/mm/dd> format, or by year)

=item C<term>

limits results by Entrez query (only valid when C<cmd=neighbor> within a single
database)

=item C<retmode>

set to XML, but can be changed to ref when needed

=item C<cookie>

a Bio::DB::EUtilities::cookie object (see below)

=item C<cmd>

command values (see below)

=item C<holding>

list LinkOut URLs for specified holding provider; used with C<cmd=llinks>
or C<cmd=llinkslib> (rarely used)

=back

=head2 Command Values

Command values are set using the C<cmd> parameter. 

=over 3

=item C<prlinks>

List the hyperlink to the primary LinkOut provider for multiple IDs and
database. Each ID is processed separately.

=item C<prlinks&retmode=ref>

Create a hyperlink to the primary LinkOut provider for a single ID and database.

=item C<llinks>

List LinkOut URLs and Attributes, except PubMed libraries, for multiple IDs
and database. Each ID is processed separately.

=item C<llinkslib>

List LinkOut URLs and Attributes for multiple IDs and database.  Each ID is
processed separately.

=item C<lcheck>

Check for the existence (Y or N) of an external link in for multiple IDs and
database.

=item C<ncheck>

Check for the existence of a neighbor link for each ID within a database,
e.g., Related Articles in PubMed.

=item C<neighbor>

The default setting. Display neighbors and their scores within a database.
This module will parse XML output from an ELink query and will return a
L<Bio::DB::EUtilities::ElinkData> object, which contains IDs for every database
liked to using C<db> (see C<id> and C<db> for more details).  

=item C<neighbor_history>

Create history (WebEnv & query_key) for use in other EUtilities.

=item C<acheck>

Lists Entrez databases links for multiple IDs from a single database.

=back

=head2 Cookies

Some EUtilities (C<epost>, C<esearch>, or C<elink>) are able to retain information on
the NCBI server under certain settings.  This information can be retrieved by
using a B<cookie>.  Here, the idea of the 'cookie' is similar to the 'cookie' set
on a user's computer when browsing the Web.  XML data returned by these
EUtilities, when applicable, is parsed for the cookie information (the 'WebEnv'
and 'query_key' tags to be specific)  The information along with other identifying
data, such as the calling eutility, description of query, etc.) is stored as a
L<Bio::DB::EUtilities::cookie|Bio::DB::EUtilities::cookie> object in an internal
queue.  These can be retrieved one at a time by using the next_cookie method or
all at once in an array using get_all_cookies.  Each cookie can then be 'fed',
one at a time, to another EUtility object, thus enabling chained queries as
demonstrated in the synopsis.

By default, a EUtilities object will retrieve records using a cookie if the
cookie parameter is set.  Also, the object will use the database parameter C<db>
stored in the L<Bio::DB::EUtilities::cookie|Bio::DB::EUtilities::cookie> object
when that parameter isn't set upon instantiation:

  my $efetch = Bio::DB::EUtilities->new(-cookie       => $elink->next_cookie,
                                        -rettype      => 'fasta');

ELink, in particular, is capable of returning multiple cookies based on the
setting for the database; if C<db> is set to C<'all'>, you will retrieve a cookie for
each database with related records.  

=head1 CURRENT USES

=head2 complex queries

Chaining queries for retrieving related data using elink and other EUtilities is
now possible (see the L</"SYNOPSIS"> for an example).

=head2 Retrieving relevancy scores

Currently, this is supported for only one ID at a time!

When the C<db> and C<dbfrom> parameters are set to the same database, one can
retrieve relevancy scores for an ID.  These are based on several different
factors.  For proteins, they are precomputed TBLASTX scores, so this is actually
a quick way to get the best hits without having to run TBLASTX directly!

=head2 ID groups

These are not completely implemented yet but support is being added in the next
few weeks.

=head1 TODO

Using multiple ID lists for processing multiple groups and maintaining
record-to-record correspondence is not completely implemented.

http://www.ncbi.nlm.nih.gov/books/bv.fcgi?rid=coursework.section.elink-considerations

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the 
evolution of this and other Bioperl modules. Send
your comments and suggestions preferably to one
of the Bioperl mailing lists. Your participation
is much appreciated.

  bioperl-l@lists.open-bio.org               - General discussion
  http://www.bioperl.org/wiki/Mailing_lists  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to
help us keep track the bugs and their resolution.
Bug reports can be submitted via the web.

  http://bugzilla.open-bio.org/

=head1 AUTHOR 

Email cjfields at uiuc dot edu

=head1 APPENDIX

The rest of the documentation details each of the
object methods. Internal methods are usually
preceded with a _

=cut

# Let the code begin...

package Bio::DB::EUtilities::elink;

use strict;
use warnings;
use Bio::DB::EUtilities;
use Bio::DB::EUtilities::Cookie;
use Bio::DB::EUtilities::ElinkData;
use URI::Escape qw(uri_unescape);
use XML::Simple;
use Data::Dumper;

use vars qw(@ISA $EUTIL $VERSION %CMD);

@ISA = qw(Bio::DB::EUtilities);

BEGIN {
    #set as default
    $EUTIL = 'elink';
    $VERSION = '1';
    # cmd parameter options; these haven't been mapped yet
    %CMD = ('prlinks'   => 1,
            'llinks'    => 1,
            'llinkslib' => 1,
            'lcheck'    => 1,
            'ncheck'    => 1,
            'neighbor'  => 1,
            'neighbor_history'  => 1,
            'acheck'    => 1,
           );
}

sub _initialize {
    my ($self, @args ) = @_;
    $self->SUPER::_initialize(@args);
	my ($term, $field, $reldate, $mindate, $maxdate, $datetype, $multi_id, $retstart,
        $retmax, $report, $dbfrom, $cmd, $holding, $version, $retmode, $linkname) = 
	  $self->_rearrange([qw(TERM FIELD RELDATE MINDATE MAXDATE DATETYPE MULTI_ID
        RETSTART RETMAX REPORT DBFROM CMD HOLDING VERSION LINKNAME)], @args);
    # set by default
    $self->_eutil($EUTIL);
    # defaults which can be overridden
    # Note : retmode should be 'xml' for all elink queries except when cmd=prlinks
    $datetype ||= 'mdat';
    $self->datetype($datetype);
    $version ||= $VERSION; # DTD to use, should leave alone
    $self->version($version);
    # normal settings
    $term       && $self->term($term);
    $field      && $self->field($field);
    $reldate    && $self->reldate($reldate);
    $mindate    && $self->mindate($mindate);
    $maxdate    && $self->maxdate($maxdate);
    $retstart   && $self->retstart($retstart);
    $retmax     && $self->retmax($retmax);
    $report     && $self->report($report);
    $dbfrom     && $self->dbfrom($dbfrom);
    $cmd        && $self->cmd($cmd);
    $holding    && $self->holding($holding);
    $linkname   && $self->linkname($linkname);
	$multi_id	&& $self->multi_id($multi_id);
}

=head2 parse_response

 Title   : parse_response
 Usage   : $db->parse_response($content)
 Function: parse out response for cookie and/or id's
 Returns : none
 Args    : HTTP::Response object
 Throws  : 'NCBI elink nonrecoverable error'
           'No links' error

=cut

sub parse_response {
	# to add: parsing for dbfrom/dbto ids, tagging cookies with databases
    my $self    = shift;
    my $response = shift if @_;
    if (!$response || !$response->isa("HTTP::Response")) {
        $self->throw("Need HTTP::Response object");
    }
    my $xs = XML::Simple->new();
    my $simple = $xs->XMLin($response->content,
                            forcearray => [qw(LinkSet LinkSetDb LinkSetDbHistory)]);
    # check for errors
    if ($simple->{ERROR}) {
        $self->throw("NCBI elink nonrecoverable error: ".$simple->{ERROR});
    }
	$self->debug("Response dumper:\n".Dumper($simple));
    my $cmd = $self->cmd ? $self->cmd : 'neighbor'; # set default cmd
    # process possible cookies first
    if (defined($cmd) && $cmd eq 'neighbor_history') {
        # process each LinkSet hash, one at at time;  
        # No scores when using history (only ids)
        if (! $simple->{LinkSet} ) {
            $self->warn('No link history');
            return;
        }
        for my $linkset (@{ $simple->{LinkSet} }) {
            my $webenv = $linkset->{WebEnv};
            my $dbfrom =  $linkset->{DbFrom};
            my $from_ids = $linkset->{IdList}->{Id};
            for my $history (@{ $linkset->{LinkSetDbHistory} }) {
                my $query_key = $history->{QueryKey};
                next if (!$query_key || ($history->{Info} eq 'Empty result') );
                my $lname = $history->{LinkName};
                my $db = $history->{DbTo};
                my $cookie = Bio::DB::EUtilities::Cookie->new(
                                        -verbose   => $self->verbose,
                                        -webenv    => $webenv,
                                        -querykey  => $query_key,
                                        -eutil     => 'elink',
                                        -database  => $db,
                                        -dbfrom    => $dbfrom,
                                        -query_id  => $from_ids,
                                        -linkname  => $lname,
                                        );
                $self->add_cookie($cookie);
            }
        }
        return;
    }
    elsif ($cmd eq 'neighbor' || !$cmd) {
        if (!$simple->{LinkSet}) {
            $self->warn('No returned links.');
            return;
        }
        for my $linkset (@{ $simple->{LinkSet} }) {
            my $linkobj = Bio::DB::EUtilities::ElinkData->new
                                (-verbose => $self->verbose,
                                 -command =>$cmd);
            $linkobj->_add_set($linkset);
            $self->_add_linkset($linkobj);
        }
    } else {
        $self->warn("$cmd not yet supported; no parsing occurred");
        return;
        # need to add a few things for cmd=llinks
    }
}

=head2 multi_id

 Title   : multi_id
 Usage   : $db->multi_id(1);
 Function: gets/sets value (switch for using multiple ids)
 Returns : value evaluating to true or false
 Args    : value evaluating to true or false

=cut

sub multi_id {
	my $self = shift;
	return $self->{'_multi_id'} = shift if @_;
	return $self->{'_multi_id'};
}

sub next_linkset {
    my $self = shift;
    if ($self->{'_linksets'}) {
        return shift @{ $self->{'_linksets'} };
    }
    $self->{'_tot_linksets'}--;
}

sub get_all_linksets {
    my $self = shift;
    return @{ $self->{'_linksets'} };
}

sub reset_linksets{
    my $self = shift;
    $self->{'_linksets'} = [];
    $self->{'_tot_linksets'} = 0;
}

sub total_linksets {
    my $self = shift;
    return $self->{'_tot_linksets'};
}

sub _add_linkset {
    my $self = shift;
    if (@_) {
        my $data_links = shift;
        $self->throw("Expecting a Bio::DB::EUtilities::ElinkData, got $data_links.")
          unless $data_links->isa("Bio::DB::EUtilities::ElinkData");
        push @{ $self->{'_linksets'} }, $data_links;
        $self->{'_tot_linksets'}++;
    }
}

=head2 Methods inherited from L<Bio::DB::EUtilities|Bio::DB::EUtilities>

=head3 add_cookie

 Title   : cookie
 Usage   : $db->add_cookie($cookie)
 Function: adds an NCBI query cookie to the internal cookie queue
 Returns : none
 Args    : a Bio::DB::EUtilities::Cookie object

=cut

=head3 next_cookie

 Title   : next_cookie
 Usage   : $cookie = $db->next_cookie
 Function: return a cookie from the internal cookie queue
 Returns : a Bio::DB::EUtilities::Cookie object
 Args    : none

=cut

=head3 reset_cookies

 Title   : reset_cookie
 Usage   : $db->reset_cookie
 Function: resets the internal cookie queue
 Returns : none
 Args    : none

=cut

=head3 get_all_cookies

 Title   : get_all_cookies
 Usage   : @cookies = $db->get_all_cookies
 Function: retrieves all cookies from the internal cookie queue; this leaves
           the cookies in the queue intact 
 Returns : array of Bio::DB::EUtilities::cookie objects
 Args    : none

=cut

=head3 get_response

 Title   : get_response
 Usage   : $db->get_response($content)
 Function: main method to retrieve data stream; parses out response for cookie
 Returns : HTTP::Response object
 Args    : optional : Bio::DB::EUtilities::cookie from a previous search
 Throws  : 'not a cookie' exception, response errors (via HTTP::Response)

=cut

=head3 reset_parameters 

 Title   : reset_parameters
 Usage   : $db->reset_parameters(@args);
 Function: resets the parameters for a EUtility with args (in @args)
 Returns : none
 Args    : array of arguments (arg1 => value, arg2 => value)

=cut

=head3 get_db_ids

 Title   : get_db_ids
 Usage   : $count = $elink->get_db_ids($db); # gets array ref of IDs
           @count = $elink->get_db_ids($db); # gets array of IDs
           %hash  = $elink->get_db_ids(); # hash of databases (keys) and array_refs(value)
 Function: returns an array or array ref if a database is the argument,
           otherwise returns a hash of the database (keys) and id_refs (values)
 Returns : array or array ref of ids (arg=database) or hash of
           database-array_refs (no args); relies on wantarray
 Args    : database (string);

=cut

1;
__END__