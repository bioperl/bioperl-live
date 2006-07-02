# $Id$
# POD to come...

# Let the code begin...

package Bio::DB::EUtilities::elink;

use strict;
use warnings;
use Bio::DB::EUtilities;
use Bio::DB::EUtilities::Cookie;
use URI::Escape qw(uri_unescape);
use XML::Simple;
use Data::Dumper;

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
    $self->_eutil($EUTIL);
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
	# to add: parsing for dbfrom/dbto ids, tagging cookies with databases
    my $self    = shift;
    my $response = shift if @_;
    if (!$response || !$response->isa("HTTP::Response")) {
        $self->throw("Need HTTP::Response object");
    }
    my $xs = XML::Simple->new();
    my $simple = $xs->XMLin($response->content,
                            forcearray => [qw(LinkSetDb Link)]);
    # check for errors
    if ($simple->{ERROR}) {
        $self->throw("NCBI elink nonrecoverable error: ".$simple->{ERROR});
    }
	$self->debug("Response dumper:\n".Dumper($simple));
    my ($webenv, $dbfrom);
    my $cmd = $self->cmd;
    # go through to make sure this catches errors
    $dbfrom = $simple->{LinkSet}->{DbFrom};

    if ($cmd && $cmd eq 'neighbor_history') {
        $webenv   = $simple->{LinkSet}->{WebEnv};
        $self->throw("No links!") unless $simple->{LinkSet}->{LinkSetDbHistory};
        for my $linkset  (@{ $simple->{LinkSet}->{LinkSetDbHistory} }) {
            my $query_key = $linkset->{QueryKey};
            my $db = $linkset->{DbTo};
            my $cookie =
            Bio::DB::EUtilities::Cookie->new(
                                             -webenv    => $webenv,
                                             -querykey  => $query_key,
                                             -eutil     => 'elink',
                                             -database  => $db,
                                             -dbfrom    => $dbfrom,
                                            );
            $self->add_cookie($cookie);
        }
        
        return;
    }
    else { # this sets internal hash of databases and ids
        for my $linkset  ( @{ $simple->{LinkSet}->{LinkSetDb}}) {
            my $id_ref = [];

            my $db = $linkset->{DbTo};
            for my $id (@{ $linkset->{Link} }) {
                push @{ $id_ref }, $id->{Id};
                $self->_add_relevancy_ids($id->{Id},$id->{Score})
                    if ($id->{Score});
            }
    
            $self->throw("Missing database and/or id; something wrong with XML")
                if (!$db && !$id_ref);
                
            $self->_add_db_ids($db, $id_ref);
        }
        # toss the dbfrom ids (you should have these already!)
        return;
    }
}

1;
__END__