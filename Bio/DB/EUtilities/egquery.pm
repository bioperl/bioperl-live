# $Id$
# POD to come...

# Let the code begin...

package Bio::DB::EUtilities::egquery;
use strict;
use warnings;
use Bio::DB::EUtilities;

use vars qw(@ISA $EUTIL);

@ISA = qw(Bio::DB::EUtilities);

BEGIN {
    #set as default
    $EUTIL = 'egquery';
}

sub _initialize {
    my ($self, @args ) = @_;
    $self->SUPER::_initialize(@args);
	my ($term) =  $self->_rearrange([qw(TERM)],@args);	
    # set by default
    $self->_eutil($EUTIL);
    $term	        && $self->term($term);
}

=head2 parse_response

 Title   : parse_response
 Usage   : $db->_parse_response($content)
 Function: parse out response for cookie
 Returns : empty
 Args    : none
 Throws  : 'unparseable output exception'

=cut

# EGQuery doesn't have error checking, so this is a no-op for now

sub parse_response {
}

1;
__END__