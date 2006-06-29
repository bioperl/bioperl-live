package Bio::DB::EUtilities::epost;
use strict;
use warnings;
use Bio::DB::EUtilities;
use Bio::DB::EUtilities::Cookie;
use URI::Escape qw(uri_unescape);

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
    $self->eutil($EUTIL);
    $self->return_mode($RETMODE);
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
    my $content = $response->content;
    # go through to make sure this catches errors
    if (my ($warning) = $content =~ m!<ErrorList>(.+)</ErrorList>!s) {
        warn "Warning(s) from GenBank: $warning\n";
    }
    if (my ($error) = $content =~ /<OutputMessage>([^<]+)/) {
        $self->throw("Error from Genbank: $error");
    }
    my ($webenv)    = $content =~ m!<WebEnv>(\S+)</WebEnv>!;
    my ($querykey)  = $content =~ m!<QueryKey>(\d+)!;
    my $cookie = Bio::DB::EUtilities::Cookie->new(-webenv   => uri_unescape($webenv),
                                                  -querykey => $querykey,
                                                  -eutil    => 'epost');
    $self->add_cookie($cookie);
    return $response;
}

1;
__END__