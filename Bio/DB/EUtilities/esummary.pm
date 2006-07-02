# $Id$
# POD to come...

# Let the code begin...

package Bio::DB::EUtilities::esummary;
use strict;
use warnings;
use Bio::DB::EUtilities;
use XML::Simple;
#use Data::Dumper;
use vars qw(@ISA $EUTIL);

@ISA = qw(Bio::DB::EUtilities);

BEGIN {
    #set as default
    $EUTIL = 'esummary';
}

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
    #$self->debug("Response dumper:\n".Dumper($simple));
    # check for errors
    if ($simple->{ERROR}) {
        $self->throw("NCBI esummary nonrecoverable error: ".$simple->{ERROR});
    }    
}

1;
__END__