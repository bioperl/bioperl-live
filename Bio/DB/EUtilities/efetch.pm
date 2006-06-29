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
    $self->eutil($EUTIL);
    $datetype ||= 'mdat';
    $self->datetype($datetype) if $datetype;
    $retstart       && $self->retstart($retstart);
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

# this is NOOP b/c efetch returns raw data to be processed or saved (i.e. no cookies)

sub parse_response {
}

1;
__END__