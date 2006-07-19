# $Id$
#
# BioPerl module for Bio::DB::EUtilities::efetch
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

Bio::DB::EUtilities::efetch

=head1 SYNOPSIS

my $efetch = Bio::DB::EUtilities->new(
                                      -verbose => 1,
                                      -cookie   => $esearch->next_cookie,
                                      -retmax   => $retmax,
                                      -rettype  => 'fasta'
                                      );

print $efetch->get_response->content;

=head1 DESCRIPTION

L<EFetch|Bio::DB::EUtilities::efetch> retrieve data records from a list of ID's
from the user environment.  This can be accomplished directly (using C<id>) or
indirectly (by using a L<Cookie|Bio::DB::EUtilities::Cookie>.

=head2 Parameters

The following are a general list of parameters that can be used to take
advantage of EFetch.  Up-to-date help for EFetch is available at this URL
(the information below is a summary of the options found there):

  http://eutils.ncbi.nlm.nih.gov/entrez/query/static/efetch_help.html

=over 3

=item C<db>

Database parameter.  This should be set at all times; the only exception is
when setting a C<cookie>.

others to follow...

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

package Bio::DB::EUtilities::efetch;
use strict;
use warnings;
use Bio::DB::EUtilities;
use URI::Escape qw(uri_unescape);

use vars qw(@ISA $EUTIL);

@ISA = qw(Bio::DB::EUtilities);

BEGIN {
    $EUTIL = 'efetch';
}

sub _initialize {
    my ($self, @args ) = @_;
    $self->SUPER::_initialize(@args);
	my ($term, $field, $reldate, $mindate, $maxdate, $datetype, $rettype, $retstart, 
        $retmax, $report, $seq_start, $seq_stop, $strand, $complexity) = 
	  $self->_rearrange([qw(TERM FIELD RELDATE MINDATE MAXDATE DATETYPE RETTYPE
        RETSTART RETMAX REPORT SEQ_START SEQ_STOP STRAND COMPLEXITY)], @args);    
    # set by default
    $self->_eutil($EUTIL);
    $datetype ||= 'mdat';
    $self->datetype($datetype) if $datetype;
    defined($retstart)       && $self->retstart($retstart);
    $retmax         && $self->retmax($retmax);
    $rettype        && $self->rettype($rettype);
    $seq_start      && $self->seq_start($seq_start);
    $seq_stop       && $self->seq_stop($seq_stop);
    $strand         && $self->strand($strand);
    defined($complexity) && $self->complexity($complexity);
    $report         && $self->report($report);    
}

=head2 parse_response

 Title   : parse_response
 Usage   : $db->_parse_response($content)
 Function: parse out response for cookie
 Returns : empty
 Args    : none
 Throws  : 'unparseable output exception'

=cut

# this is NOOP b/c efetch returns raw data to be processed or saved;
# HTTP errors caught in get_response 

sub parse_response {
}

1;
__END__