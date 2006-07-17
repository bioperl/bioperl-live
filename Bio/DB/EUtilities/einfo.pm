# $Id$
#
# BioPerl module for Bio::DB::EUtilities::einfo
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

Bio::DB::EUtilities::einfo - NCBI database information

=head1 SYNOPSIS

my $einfo = Bio::DB::EUtilities->new(
                                     -eutil    => 'einfo',
                                     -db       => 'pubmed'
                                      );

print $einfo->get_response->content;

=head1 DESCRIPTION

L<Bio::DB::EUtilities::einfo|Bio::DB::EUtilities::einfo> queries provide
information about NCBI databases.  At this time no data is postprocessed
for further information; this is currently implemented to complete access
to the NCBI EUtilities suite of web tools.  Using the base URL with no
parameters returns a list of all Entrez databases.

=head2 Parameters

The following are a general list of parameters that can be used to take
advantage of EInfo.  Up-to-date help for EInfo is available at this URL
(the information below is a summary of the options found there):

  http://eutils.ncbi.nlm.nih.gov/entrez/query/static/einfo_help.html

=over 3

=item C<db>

Database parameter.  This is optional; not setting this will return a list of
all the available Entrez databases.  If a database is specified, returned XML
data will provide the following information : field names, index term counts,
last update, and available links for each Entrez database.

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

package Bio::DB::EUtilities::einfo;
use strict;
use warnings;
use Bio::DB::EUtilities;
use XML::Simple;
use Data::Dumper;

use vars qw(@ISA $EUTIL);

@ISA = qw(Bio::DB::EUtilities);

BEGIN {
    $EUTIL = 'einfo';
}

sub _initialize {
    my ($self, @args ) = @_;
    $self->SUPER::_initialize(@args);
    # set by default
    $self->_eutil($EUTIL);
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
    if (!$response || !$response->isa("HTTP::Response")) {
        $self->throw("Need HTTP::Response object");
    }
    my $xs = XML::Simple->new();
    my $simple = $xs->XMLin($response->content);
    $self->debug("Response dumper:\n".Dumper($simple));
    # check for errors
    if ($simple->{ERROR}) {
        $self->throw("NCBI esummary nonrecoverable error: ".$simple->{ERROR});
    }
}

1;
__END__