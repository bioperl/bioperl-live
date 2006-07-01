
# Let the code begin...

package Bio::DB::EUtilities::einfo;
use strict;
use warnings;
use Bio::DB::EUtilities;

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
    my $content = $response->content;
    # go through to make sure this catches errors
    if (my ($warning) = $content =~ m!<ErrorList>(.+)</ErrorList>!s) {
        warn "Warning(s) from GenBank: $warning\n";
    }
    if (my ($error) = $content =~ /<OutputMessage>([^<]+)/) {
        $self->throw("Error from Genbank: $error");
    }
}

1;
__END__