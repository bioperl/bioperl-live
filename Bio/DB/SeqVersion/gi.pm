#
# BioPerl module for Bio::DB::SeqVersion::gi
#
# Please direct questions and support issues to <bioperl-l@bioperl.org>
#
# Cared for by Brian Osborne
#
# Copyright Brian Osborne 2006
#
# You may distribute this module under the same terms as Perl itself
#
# POD documentation - main docs before the code

=head1 NAME

Bio::DB::SeqVersion::gi - interface to NCBI Sequence Revision History page

=head1 SYNOPSIS

Do not use this module directly, use Bio::DB::SeqVersion.

    use Bio::DB::SeqVersion;

    my $query = Bio::DB::SeqVersion->new(-type => 'gi');

    # all GIs, which will include the GI used to query
    my @all_gis = $query->get_all(2);

    # the most recent GI, which may or may not be the GI used to query
    my $live_gi = $query->get_recent(2);

    # get all the visible data on the Sequence Revision page
    my $array_ref = $query->get_history(11111111);

These methods can also take accession numbers as arguments, just like
the Sequence Revision page itself.

=head1 DESCRIPTION

All sequence entries at GenBank are identified by a pair of
identifiers, an accession and a numeric identifier, and this number is
frequently called a GI number (B<G>enInfo B<I>dentifier). The accession
is stable, but each new version of the sequence entry for the accession
receives a new GI number (see L<http://www.ncbi.nlm.nih.gov/Sitemap/sequenceIDs.html>
for more information on GenBank identifiers). One accession
can have one or more GI numbers and the highest of these is the most recent,
or "live", GI.

Information on an accession and its associated GI numbers is available at
the Sequence Revision History page at NCBI,
L<http://www.ncbi.nlm.nih.gov/entrez/sutils/girevhist.cgi>, this information is
not available in file format. This module queries the Web page and retrieves GI
numbers and related data given an accession (e.g. NP_111111, A11111, P12345) or
a GI number (e.g. 2, 11111111) as query.

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
the bugs and their resolution.  Bug reports can be submitted via the
web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Brian Osborne

Email E<lt> osborne at optonline dot net E<gt>

=head1 CONTRIBUTORS

Torsten Seemann - torsten.seemann AT infotech.monash.edu.au

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::DB::SeqVersion::gi;
use strict;
use Encode;
use HTML::TableExtract;
use base qw(Bio::DB::SeqVersion);

# Private class variables

# TODO: this may be an unstable setting (text is actually minimal XHTML)
my $URL = 'http://www.ncbi.nlm.nih.gov/nuccore/%s?report=girevhist&format=text';

=head2 new

 Title   : new
 Usage   : $gb = Bio::DB::SeqVersion::gi->new
 Function: Creates a new query object
 Returns : New query object

=cut

sub new {
    my ( $class, @args ) = @_;
    my $self = $class->SUPER::new(@args);
    $self->_initialize;
    return $self;
}

=head2 get_all

 Title   : get_all
 Usage   : my @gis = $q->get_all(2)
 Function: Get all GI numbers given a GI number
 Returns : An array of GI numbers, earliest GI number is the 0 element
 Args    : A single GI number (string)

=cut

sub get_all {
    my ( $self, $id ) = @_;
    my ( @arr, $ref );
    $id eq $self->{_last_id}
      ? $ref = $self->{_last_result}
      : $ref = $self->get_history($id);
    for my $row ( @{$ref} ) {
        push @arr, $$row[0];
    }
    @arr;
}

=head2 get_recent

 Title   : get_recent
 Usage   : my $newest_gi = $q->get_recent(2)
 Function: Get most recent GI given a single GI
 Returns : String
 Args    : A single GI number (string)

=cut

sub get_recent {
    my ( $self, $id ) = @_;
    my $ref;
    $id eq $self->{_last_id}
      ? $ref = $self->{_last_result}
      : $ref = $self->get_history($id);
    $ref->[0]->[0];
}

=head2 get_status

 Title   : get_status
 Usage   : my $newest_gi = $q->get_status(2)
 Function: Get most recent GI given a single GI
 Returns : String
 Args    : A single GI number (string)

=cut

sub get_status {
    my ( $self, $id ) = @_;
    $self->throw("Must pass an ID") if !defined $id;
    if ($id ne $self->{_last_id} ) {
        $self->get_history($id);
    }
    $self->{_last_status};
}

=head2 get_history

 Title   : get_history
 Usage   : my $ref = $query_obj->get_history()
 Function: Queries the NCBI Revision page, gets the data from the HTML table
 Returns : Reference to an array of arrays where element 0 refers to the most
           recent version and the last element refers to the oldest version.
           In the second dimension the elements are:

           0      GI number
           1      Version
           2      Update Date

           For example, to get the GI number of the first version:

           $ref->[$#{@$ref}]->[0]

           To get the Update Date of the latest version:

           $ref->[0]->[2]

 Args    : One identifier (string)
 Note    : Status of the GI was returned here previously as the last element in
           the row of elemnts above; however the status is currently only
           returned for the GI requested (e.g. a single value).  One can get
           the status for this using the get_status() method above

=cut

sub get_history {
    my ( $self, $id ) = @_;
    my $html = $self->_get_request($id);
    my ( $ref, $status ) = $self->_process_data($html);

    # store the very last result in case some other methods
    # are called using the same identifier
    $self->{_last_result} = $ref;
    $self->{_last_id}     = $id;
    $self->{_last_status} = $status;
    $ref;
}

=head2 _get_request

 Title   : _get_request
 Usage   : my $url = $self->_get_request
 Function: GET using NCBI Revision page URL, uses Root::HTTPget
 Returns : HTML
 Args    : One identifier (string)

=cut

sub _get_request {
    my ( $self, $id ) = @_;

    $self->throw("Must specify a single id to query") if ( !defined($id) || ref($id) );

    my $url = sprintf( $URL, $id );
    my $response = $self->get($url);
    if ( not $response->is_success ) {
        $self->throw( "Can't query $url: "
              . $response->status_line . "\n"
              . "ID likely does not exist" );
    }
    return $response->content;
}

=head2 _process_data

 Title   : _process_data
 Usage   : $self->_process_data($html)
 Function: extract data from HTML
 Args    : HTML from Revision History page
 Returns : reference to an array of arrays

=cut

sub _process_data {
    my ( $self, $html ) = @_;

    # Only one status is returned (not one per revision).  Setting once
    my $status;
    if ($html =~ /<div class="status">Current status:\s+(\S+)<\/div>/) {
        $status = $1;
    } else {
        $self->warn("No current status found, setting to 'unknown'");
        $status = 'unknown';
    }

    my $te = HTML::TableExtract->new(
        headers => ['Gi', 'Version', 'Update Date'] ,
        depth   => 0);
    $te->parse(decode_utf8($html));
    my $table = $te->first_table_found;
    $self->throw("No table found") unless defined $table;
    my $t = [$table->rows];
    ($t, $status);
}

1;

__END__
