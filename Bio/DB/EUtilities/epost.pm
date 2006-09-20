# $Id$
#
# BioPerl module for Bio::DB::EUtilities::epost
#
# Cared for by Chris Fields
#
# Copyright Chris Fields
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code
# 
# Part of the EUtilities BioPerl package

=head1 NAME

Bio::DB::EUtilities::epost - posting IDs on the remote NCBI server for batch
retrieval and chained queries

=head1 SYNOPSIS

    my $epost = Bio::DB::EUtilities->new(
                                          -eutil    => 'epost',
                                          -id       => \@ids,
                                          -db       => 'protein',
                                          );

    $epost->get_response;

=head1 DESCRIPTION

B<WARNING>: Please do B<NOT> spam the Entrez web server with multiple requests.

The EUtility EPost is used to post a list of primary IDs to the NCBI EUtilities
server for retrieval by L<EFetch|Bio::DB::EUtilities::efetch> or for using in
futher searches using L<ELink|Bio::DB::EUtilities::elink> or
L<ESearch|Bio::DB::EUtilities::esearch>.  The data is posted using:

    $epost->get_response;

When not used in void context, this will also return a
L<HTTP::Response|HTTP::Response> object for further processing.  This is not
necessary, as any posts made will automatically generate a
L<Cookie|Bio::DB::EUtilities::Cookie>,
which can be used to retrieve the posted information using
L<EFetch|Bio::DB::EUtilities::efetch>.

Using EPost is recommended for retrieving large lists of primary IDs and is
capable, when used repeatedly and in combination with EFetch, of retrieving
thousands of database entries.  

=head2 Parameters

The following are a general list of parameters that can be used to take
advantage of EPost.  Up-to-date help for EPost is available at this URL
(the information below is a summary of the options found there):

  http://eutils.ncbi.nlm.nih.gov/entrez/query/static/epost_help.html

=over 3

=item C<db>

The name of an Entrez database available through EUtilities. 

=item C<id>

a list of primary ID's

Below are a list of IDs which can be used with EPost:

B<PMID> (pubmed), B<MEDLINE UI> (NIH MedLine), B<MIM number> (omim),
B<GI number> (nucleotide, protein), B<MMDB-ID> (structure), B<TAXID> (taxonomy)

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

package Bio::DB::EUtilities::epost;
use strict;
use warnings;
use Bio::DB::EUtilities;
use Bio::DB::EUtilities::Cookie;
use XML::Simple;
#use Data::Dumper;

use vars qw(@ISA $EUTIL $RETMODE);

@ISA = qw(Bio::DB::EUtilities);

BEGIN {
    #set as default
    $EUTIL = 'epost';
    $RETMODE = 'xml';
}

sub _initialize {
    my ($self, @args) = @_;
    $self->SUPER::_initialize(@args);
    # set by default
    $self->_eutil($EUTIL);
    $self->retmode($RETMODE);
}

=head2 parse_response

 Title   : parse_response
 Usage   : $db->_parse_response($content)
 Function: parse out response for cookie
 Returns : empty
 Args    : none
 Throws  : 'unparseable output exception'

=cut

sub parse_response {
    my $self    = shift;
    my $response = shift if @_;
    if (!$response || !($response->isa("HTTP::Response"))) {
        $self->throw("Need HTTP::Response object");
    }
    my $xs = XML::Simple->new();
    my $simple = $xs->XMLin($response->content);
    #$self->debug("Response dumper:\n".Dumper($simple));
    # check for errors
    if ($simple->{ERROR}) {
        $self->throw("NCBI epost nonrecoverable error: ".$simple->{ERROR});
    }
    if ($simple->{InvalidIdList}) {
        $self->warn("NCBI epost error: Invalid ID List".$simple->{InvalidIdList});
    }
    my $db = $self->db;
    my $webenv    = $simple->{WebEnv};
    my $querykey  = $simple->{QueryKey};
    my $cookie = Bio::DB::EUtilities::Cookie->new(-webenv   => $webenv,
                                                  -querykey => $querykey,
                                                  -eutil    => 'epost',
                                                  -database => $db,
                                                  );
    $self->add_cookie($cookie);
    return $response;
}

1;
__END__