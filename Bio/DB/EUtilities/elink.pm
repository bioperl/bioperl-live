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

  use Bio::DB::EUtilities;

  my $esearch = Bio::DB::EUtilities->new(-eutil      => 'esearch',
                                         -db         => 'pubmed',
                                         -term       => 'hutP',
                                         -usehistory => 'y');
  
  $esearch->get_response; # parse the response, fetch a cookie
  
  my $elink = Bio::DB::EUtilities->new(-eutil        => 'elink',
                                       -db           => 'protein',
                                       -dbfrom       => 'pubmed',
                                       -cookie       => $esearch->next_cookie,
                                       -cmd          => 'neighbor');
  
  my @prot_ids = $elink->get_db_ids('protein'); # retrieves protein UID's

=head1 DESCRIPTION

B<WARNING>: Please do B<NOT> spam the Entrez web server with multiple requests.

The EUtility Elink is used to check for and retrieve external or related ID's
from a list of one or more primary ID's.  There are some pretty variations on
what can be returned based on the parameters used.  

=head2 Parameters

The following are a general list of parameters that can be used to take
advantage of ELink.  Up-to-date help for ELink is available here (all information
below is from this location):

  http://eutils.ncbi.nlm.nih.gov/entrez/query/static/elink_help.html

=over 3

=item C<db>

one or more database available through EUtilities if set to 'all', will retrieve
all related ID's from each database (see method get_db_ids to retrieve these).

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
                            forcearray => [qw(LinkSetDb LinkSetDbHistory Link)]);
    # check for errors
    if ($simple->{ERROR}) {
        $self->throw("NCBI elink nonrecoverable error: ".$simple->{ERROR});
    }
	$self->debug("Response dumper:\n".Dumper($simple));
    my ($webenv, $dbfrom);
    my $cmd = $self->cmd;
    if ($self->multi_id) {
        $self->warn("Not implemented yet for multiple ID groups\n".
                    "No scores or IDs retained");
        return;
    }
    # go through to make sure this catches errors
    $dbfrom = $simple->{LinkSet}->{DbFrom};

    if ($cmd && $cmd eq 'neighbor_history') {
        $webenv   = $simple->{LinkSet}->{WebEnv};
        $self->throw("No links; make sure you are retrieving the cookie correctly")
		    unless $simple->{LinkSet}->{LinkSetDbHistory};
        for my $linkset  (@{ $simple->{LinkSet}->{LinkSetDbHistory} }) {
            my $query_key = $linkset->{QueryKey};
            my $db = $linkset->{DbTo};
            my $cookie = Bio::DB::EUtilities::Cookie->new(
												 -webenv    => $webenv,
												 -querykey  => $query_key,
												 -eutil     => 'elink',
												 -database  => $db,
												 -dbfrom    => $dbfrom,
												);
            $self->add_cookie($cookie);
        }
        
        return;
    }
    else { # this sets internal hash of databases and ids
		# need to rethink how to handle scores and so on with multiple IDs;
		# maybe do something like with get_db_ids (hash of id and id-score pairs? Oi!!)
		# will ned to rethink get_scores and get_ids_by_score based on same priciple
        for my $linkset  ( @{ $simple->{LinkSet}->{LinkSetDb} } ) {
            my $id_ref = [];
            my $db = $linkset->{DbTo};
            for my $id (@{ $linkset->{Link} }) {
                push @{ $id_ref }, $id->{Id};
                $self->_add_relevancy_ids($id->{Id},$id->{Score})
                    if ($id->{Score});
            }
            $self->throw("Missing database and/or id; something wrong with XML")
                if (!$db && !$id_ref);
                
            $self->_add_db_ids($db, $id_ref);
        }
        # toss the dbfrom ids (you should have these already!)
        return;
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

=head2 get_score

 Title   : get_score
 Usage   : $score = $db->get_score($id);
 Function: gets score for ID (if present)
 Returns : integer (score) 
 Args    : ID values

=cut

sub get_score {
    my $self = shift;
    my $id = shift if @_;
    $self->throw("No ID given") if !$id;
    return $self->{'_rel_ids'}->{$id} if $self->{'_rel_ids'}->{$id};
    $self->warn("No score for $id");
}

=head2 get_ids_by_score

 Title   : get_ids_by_score
 Usage   : @ids = $db->get_ids_by_score;  # returns IDs
           @ids = $db->get_ids_by_score($score); # get IDs by score
 Function: returns ref of array of ids based on relevancy score from elink;
           To return all ID's above a score, use the normal score value;
           to return all ID's below a score, append the score with '-';
 Returns : ref of array of ID's; if array, an array of IDs
 Args    : integer (score value); returns all if no arg provided

=cut

sub get_ids_by_score {
    my $self = shift;
    my $score = shift if @_;
    my @ids;
    if (!$score) {
        @ids = sort keys %{ $self->{'_rel_ids'} };
    }
    elsif ($score && $score > 0) {
        for my $id (keys %{ $self->{'_rel_ids'}}) {
            push @ids, $id if $self->{'_rel_ids'}->{$id} > $score;
        }
    }
    elsif ($score && $score < 0) {
        for my $id (keys %{ $self->{'_rel_ids'}}) {
            push @ids, $id if $self->{'_rel_ids'}->{$id} < abs($score);
        }
    }
    if (@ids) {
        @ids = sort {$self->get_score($b) <=> $self->get_score($a)} @ids;
        return @ids if wantarray;
        return \@ids;
    }
    # if we get here, there's trouble
    $self->warn("No returned IDs!");
}

=head1 Private methods

=head2 _add_relevancy_ids

 Title   : _add_relevancy_ids
 Usage   : $self->add_relancy_ids($db, $ids);
 Function: sets internal hash of ids with relevancy scores
 Returns : none
 Args    : two numbers: id (key) and relevancy score (value)

=cut

sub _add_relevancy_ids {
    my ($self, $id, $score) = @_;
    $self->throw ("Must have id-score pair for hash") unless ($id && $score);
    $self->throw ("ID, score must be scalar") if
         (ref($id) && ref($score));
    $self->{'_rel_ids'}->{$id} = $score;
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

B<Experimental method at this time>

=cut

=head3 get_db_ids

 Title   : get_db_ids
 Usage   : $count = $elink->get_db_ids($db); # gets array ref of IDs
           @count = $elink->get_db_ids($db); # gets array of IDs
           %hash  = $elink->get_db_ids(); # hash of databases (keys) and array_refs(value)
 Function: returns an array or array ref if a database is the argument,
           otherwise returns a hash of the database (keys) and id_refs (values)
 Returns : array or array ref of ids (arg=database) or hash of
           database-array_refs (no args)
 Args    : database string;

=cut

1;
__END__