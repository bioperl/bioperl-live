
# simple object to hold NCBI cookie information and descriptions

package Bio::DB::EUtilities::Cookie;
use strict;
use warnings;
use URI::Escape qw(uri_unescape);
use Bio::Root::Root;
use vars '@ISA';

@ISA = 'Bio::Root::Root';

sub new {
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(@args);
    my ($webenv, $querykey, $desc, $eutil, $total) = $self->_rearrange
      ([qw(WEBENV QUERYKEY DESCRIPTION EUTIL TOTAL)], @args);
    unless ($webenv && $querykey) {
        my $missing;
        if (!$webenv) {
            $missing = 'WebEnv';
        } elsif (!$querykey) {
            $missing = 'query_key';
        } else {
            $self->throw("Abnormal cookie");
        }
        $self->throw("Missing ".$missing);
    }
    $self->cookie(uri_unescape($webenv), $querykey);
    # holds description of database, elink, etc if present
    $desc       && $self->description($desc);
    # holds total hits if present (i.e. esearch)
    $total      && $self->total($total);
    # holds originating eutil
    $eutil      && $self->eutil($eutil);
    return $self;
}

sub cookie {
    my $self = shift;
    if (@_) {
        my ($webenv, $querykey) = (shift, shift);
        $self->throw("Missing part of cookie!") if (!$webenv || !$querykey);
        return $self->{'_cookie'} = [$webenv, $querykey];
    } else {
        return $self->{'_cookie'};
    }
}

sub eutil {
    my $self = shift;
    return $self->{'_eutil'} = shift if @_;
    return $self->{'_eutil'};
}

sub description {
    my $self = shift;
    return $self->{'_description'} = shift if @_;
    return $self->{'_description'};
}

sub total {
    my $self = shift;
    return $self->{'_total'} = shift if @_;
    return $self->{'_total'};
}

1;
__END__