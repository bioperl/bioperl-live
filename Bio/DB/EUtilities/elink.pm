# Let the code begin...

package Bio::DB::EUtilities::elink;
use strict;
use warnings;
use Bio::DB::EUtilities;
use Bio::DB::EUtilities::Cookie;
use URI::Escape qw(uri_unescape);
use vars qw(@ISA $EUTIL $VERSION %CMD);

@ISA = qw(Bio::DB::EUtilities);

BEGIN {
    #set as default
    $EUTIL = 'elink';
    $VERSION = '1';
    # cmd parameter options; these haven't been mapped yet
    %CMD = ('prlinks'   => 1,
            'llinks'    => 1,
            'llinkslib' => 1,
            'lcheck'    => 1,
            'ncheck'    => 1,
            'neighbor'  => 1,
            'neighbor_history'  => 1,
            'acheck'    => 1,
           );
}

sub _initialize {
    my ($self, @args ) = @_;
    $self->SUPER::_initialize(@args);
	my ($term, $field, $reldate, $mindate, $maxdate, $datetype, $retstart,
        $retmax, $report, $dbfrom, $cmd, $holding, $version, $retmode, $linkname) = 
	  $self->_rearrange([qw(TERM FIELD RELDATE MINDATE MAXDATE DATETYPE
        RETSTART RETMAX REPORT DBFROM CMD HOLDING VERSION LINKNAME)], @args);
    # set by default
    $self->eutil($EUTIL);
    # defaults which can be overridden
    # Note : retmode should be 'xml' for all elink queries except when cmd=prlinks
    $datetype ||= 'mdat';
    $self->datetype($datetype);
    $version ||= $VERSION; # DTD to use, should leave alone
    $self->version($version);
    # normal settings
    $term       && $self->term($term);
    $field      && $self->field($field);
    $reldate    && $self->reldate($reldate);
    $mindate    && $self->mindate($mindate);
    $maxdate    && $self->maxdate($maxdate);
    $retstart   && $self->retstart($retstart);
    $retmax     && $self->retmax($retmax);
    $report     && $self->report($report);
    $dbfrom     && $self->dbfrom($dbfrom);
    $cmd        && $self->cmd($cmd);
    $holding    && $self->holding($holding);
    $linkname   && $self->linkname($linkname);
}

=head2 parse_response

 Title   : parse_response
 Usage   : $db->_parse_response($content)
 Function: parse out response for cookie
 Returns : response
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
    my $webenv;
    my @querykey;
    # go through to make sure this catches errors
    if (my ($warning) = $content =~ m!<ErrorList>(.+)</ErrorList>!s) {
        warn "Warning(s) from GenBank: $warning\n";
    }
    if (my ($error) = $content =~ /<OutputMessage>([^<]+)/) {
        $self->throw("Error from Genbank: $error");
    }
    if ($self->cmd && $self->cmd eq 'neighbor_history') {
        # use cookie objects here; need to switch to XML parsing
        # grab descriptions and other info
        ($webenv)    = ($content =~ m!<WebEnv>(\S+)</WebEnv>!);
        # for db=all, get multiple query keys!
        @querykey = ($content =~ m!<QueryKey>(\d+)!g );
        # if, for some reason, the query key is missing from the XML output,
        # defaults to 1
        push @querykey, 1 if !@querykey;
        for my $qk (@querykey) {
            my $cookie =
            Bio::DB::EUtilities::Cookie->new(
                                             -webenv    => $webenv,
                                             -querykey  => $qk,
                                             -eutil     => 'elink',
                                             -description   => 'test'
                                            );
        $self->add_cookie($cookie);
        }
    }
}

1;
__END__