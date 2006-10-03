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

  print $esummary->get_response->content; # prints XML output

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
use Bio::DB::EUtilities::DocSum;
use strict;
use warnings;
use XML::Simple;
use Data::Dumper;

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
    my $simple = $xs->XMLin($response->content,
                            forcearray => [qw(DocSum Item)]);
    #$self->debug("Response dumper:\n".Dumper($simple));
    # check for errors
    if ($simple->{ERROR}) {
        $self->throw("NCBI esummary nonrecoverable error: ".$simple->{ERROR});
    }
    if (!exists $simple->{DocSum}) {
        $self->warn('No returned docsums.');
        return;
    }
    print STDERR "reference:",ref($simple->{DocSum}),"\n";
    for my $data (@{ $simple->{DocSum} }) {
        my $ds = Bio::DB::EUtilities::DocSum->new(-verbose => $self->verbose);
        $ds->_add_data($data->{Item}) if exists $data->{Item};
        $ds->esummary_id($data->{Id}) if exists $data->{Id}; 
        $self->debug("DocSum Object:".Dumper($ds));
        $self->add_docsum($ds);
    }
}

sub add_docsum {
    my $self = shift;
    if (@_) {
        my $docsum = shift;
        $self->throw("Expecting a Bio::DB::EUtilities::DocSum, got $docsum.")
          unless $docsum->isa("Bio::DB::EUtilities::DocSum");
        push @{ $self->{'_docsums'} }, $docsum;
        $self->{'_docsum_ct'}++;
    }
}

sub next_docsum {
    my $self = shift;
    my $index = $self->_next_docsum_index;
    return if ($index > scalar($self->{'_docsums'}));
    return $self->{'_docsums'}->[$index] ;
}

=head2 get_all_linksets

 Title   : get_all_linksets
 Usage   : @ls = $elink->get_all_linksets;
 Function: returns array of Bio::DB::EUtilities::ElinkData objects
 Returns : array or array ref of Bio::DB::EUtilities::ElinkData objects
           based on wantarray
 Args    : None

=cut

sub get_all_docsums {
    my $self = shift;
    return @{ $self->{'_docsums'} } if wantarray;
    return $self->{'_docsums'};
}

=head2 reset_linksets

 Title   : reset_linksets
 Usage   : $elink->reset_linksets;
 Function: resets (empties) internal cache of Linkset objects
 Returns : None
 Args    : None

=cut

sub reset_docsums{
    my $self = shift;
    $self->{'_docsums'} = [];
    $self->rewind_docsums;
    $self->{'_docsum_ct'} = 0;
}

=head2 rewind_linksets

 Title   : rewind_linksets
 Usage   : $elink->rewind_linksets;
 Function: resets linkset index to 0 (starts over)
 Returns : None
 Args    : None

=cut

sub rewind_linksets{
    my $self = shift;
    $self->{'_docsumindex'} = 0;
}

sub _next_docsum {
    my $self = shift;
    return $self->{'_docsumindex'}++;
}


1;
__END__
