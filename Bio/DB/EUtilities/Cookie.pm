# $Id$

# simple object to hold NCBI cookie information and descriptions
# POD to come...

=head1 NAME

Bio::DB::EUtilities::Cookie - simple object to hold NCBI cookie information and descriptions

=head1 DESCRIPTION

Some EUtilities (C<epost>, C<esearch>, or C<elink>) are able to retain information on
the NCBI server under certain settings.  This information can be retrieved by
using a B<cookie>.  Here, the idea of the 'cookie' is similar to the 'cookie' set
on a user's computer when browsing the Web.  XML data returned by these
EUtilities, when applicable, is parsed for the cookie information (the 'WebEnv'
and 'query_key' tags to be specific)  The information along with other identifying
data, such as the calling eutility, description of query, etc.) is stored as a
L<Bio::DB::EUtilities::cookie|Bio::DB::EUtilities::cookie> object in an internal queue.  These can be retrieved
one at a time by using the next_cookie method or all at once in an array using
get_all_cookies.  Each cookie can then be 'fed', one at a time, to another
EUtility object, thus enabling chained queries as demonstrated in the synopsis.

By default, a EUtilities object will retrieve records using a cookie if the
cookie parameter is set.  Also, the object will use the database parameter
stored in the L<Bio::DB::EUtilities::cookie|Bio::DB::EUtilities::cookie> object when the parameter isn't set
upon instantiation:

  my $efetch = Bio::DB::EUtilities->new(-cookie       => $elink->next_cookie,
                                        -rettype      => 'fasta');

ELink, in particular, is capable of returning multiple cookies based on the
setting for the database; if C<db> is set to C<'all'>, you will retrieve a cookie for
each database with related records.

=cut

package Bio::DB::EUtilities::Cookie;
use strict;
use warnings;
use URI::Escape qw(uri_unescape);

use base qw(Bio::Root::Root);

sub new {
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(@args);
    my ($webenv, $querykey, $database, $dbfrom, $query_id, $eutil,
      $total, $term, $linkname) = $self->_rearrange ([qw(WEBENV QUERYKEY
      DATABASE DBFROM QUERY_ID EUTIL TOTAL TERM LINKNAME)], @args);
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
    # holds originating eutil
    $eutil      && $self->eutil($eutil);
    # holds descriptions of database being queried
    $database   && $self->database($database);
    
    # for elink only, originating database
    $dbfrom     && $self->elink_dbfrom($dbfrom);
    # holds elink dbfrom ID's used for querys
    $query_id   && $self->elink_queryids($query_id);
    # holds elink linkname; information can be found using einfo
    $linkname   && $self->elink_linkname($linkname);    

    # for esearch, to hold original search query
    $term       && $self->esearch_query($term);
    # for esearch, holds total hits if present
    $total      && $self->esearch_total($total);

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

sub esearch_total {
    my $self = shift;
    return $self->{'_esearch_total'} = shift if @_;
    return $self->{'_esearch_total'};
}

sub esearch_query {
    my $self = shift;
    return $self->{'_esearch_query'} = shift if @_;
    return $self->{'_esearch_query'};
}

sub elink_dbfrom {
    my $self = shift;
    return $self->{'_elink_dbfrom'} = shift if @_;
    return $self->{'_elink_dbfrom'};
}

sub elink_queryids {
    my $self = shift;
    return $self->{'_query_ids'} = shift if @_;
    return @{ $self->{'_query_ids'} } if wantarray;
    return $self->{'_query_ids'};
}

sub elink_linkname {
    my $self = shift;
    return $self->{'_elink_linkname'} = shift if @_;
    return $self->{'_elink_linkname'};
}

1;
__END__