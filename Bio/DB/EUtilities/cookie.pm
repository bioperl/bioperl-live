# $Id$

# simple object to hold NCBI cookie information and descriptions
# POD to come...

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
    my ($webenv, $querykey, $database, $dbfrom, $eutil, $total, $desc) =
      $self->_rearrange
      ([qw(WEBENV QUERYKEY DATABASE DBFROM EUTIL TOTAL DESCRIPTION)], @args);
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
    # holds descriptions of database
    $database   && $self->database($database);
    # for elink only
    $dbfrom     && $self->dbfrom($dbfrom);
    # for esearch, to hold original search query
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

sub database {
    my $self = shift;
    return $self->{'_database'} = shift if @_;
    return $self->{'_database'};
}

sub dbfrom {
    my $self = shift;
    return $self->{'_dbfrom'} = shift if @_;
    return $self->{'_dbfrom'};
}

sub total {
    my $self = shift;
    return $self->{'_total'} = shift if @_;
    return $self->{'_total'};
}

sub description {
    my $self = shift;
    return $self->{'_description'} = shift if @_;
    return $self->{'_description'};
}

1;
__END__