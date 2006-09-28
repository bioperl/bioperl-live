# $Id$
#
# BioPerl module for Bio::DB::EUtilities::esummary
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

Bio::DB::EUtilities::esummary - retrieval of NCBI DocSum data from a list
of primary IDs or a Cookie

=head1 SYNOPSIS

B<Do not use this module directly.>
Use it via the L<Bio::DB::EUtilities|Bio::DB::EUtilities> class.

  use Bio::DB::EUtilities;

  my $esearch = Bio::DB::EUtilities->new(-eutil      => 'esearch',
                                         -db         => 'pubmed',
                                         -term       => 'hutP',
                                         -usehistory => 'y');

  $esearch->get_response; # parse the response, fetch a cookie

  my $esummary = Bio::DB::EUtilities->new(-eutil        => 'esummary',
                                       -cookie       => $esearch->next_cookie);

  print $esearch->get_response-content; # prints XML output

=head1 DESCRIPTION

B<WARNING>: Please do B<NOT> spam the Entrez web server with multiple requests.

The EUtility ESummary is used to retrieve ducument summaries from a list of
primary IDs or the user's history (stored on the remote server and accessible
using a L<Cookie|Bio::DB::EUtilities::Cookie>.  The returned data is processed
for errors, but no further processing is done at this time.

=over 3

=item C<db>

one or more database available through EUtilities if set to 'all', will retrieve
all related ID's from each database (see method get_db_ids to retrieve these)

=item C<id>

a list of primary ID's (see below)

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

a Bio::DB::EUtilities::Cookie object (see below)

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

package Bio::DB::EUtilities::esummary;
use strict;
use warnings;
use XML::Simple;
#use Data::Dumper;

use base qw(Bio::DB::EUtilities);

our $EUTIL = 'esummary';

sub _initialize {
    my ($self, @args ) = @_;
    $self->SUPER::_initialize(@args);
    # set by default
    $self->_eutil($EUTIL);
	my ($retstart, $retmax) =  $self->_rearrange([qw(RETSTART RETMAX)],@args);
    $retstart       && $self->retstart($retstart);
    $retmax         && $self->retmax($retmax);
}

=head2 parse_response

 Title   : parse_response
 Usage   : $db->parse_response($content)
 Function: parse out response for cookie and/or id's
 Returns : none
 Args    : HTTP::Response object
 Throws  : 'NCBI elink nonrecoverable error'

=cut

sub parse_response {
    my $self    = shift;
    my $response = shift if @_;
    if (!$response || !$response->isa("HTTP::Response")) {
        $self->throw("Need HTTP::Response object");
    }
    my $xs = XML::Simple->new();
    my $simple = $xs->XMLin($response->content);
    #$self->debug("Response dumper:\n".Dumper($simple));
    # check for errors
    if ($simple->{ERROR}) {
        $self->throw("NCBI esummary nonrecoverable error: ".$simple->{ERROR});
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
 Returns : none
 Args    : none

=cut

=head3 get_response

 Title   : get_response
 Usage   : $db->get_response($content)
 Function: main method to retrieve data stream; parses out response for cookie
 Returns : HTTP::Response object
 Args    : optional : Bio::DB::EUtilities::Cookie from a previous search
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

=head3 count

 Title   : count
 Usage   : $count = $db->count;
 Function: return count of number of entries retrieved by query
 Returns : integer
 Args    : none

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

=head3 get_score

 Title   : get_score
 Usage   : $score = $db->get_score($id);
 Function: gets score for ID (if present)
 Returns : integer (score) 
 Args    : ID values

=cut

=head3 get_ids_by_score

 Title   : get_ids_by_score
 Usage   : @ids = $db->get_ids_by_score;  # returns IDs
           @ids = $db->get_ids_by_score($score); # get IDs by score
 Function: returns ref of array of ids based on relevancy score from elink;
           To return all ID's above a score, use the normal score value;
           to return all ID's below a score, append the score with '-';
 Returns : ref of array of ID's; if array, an array of IDs
 Args    : integer (score value); returns all if no arg provided

=cut

1;
__END__
