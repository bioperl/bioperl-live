# $Id$
#
# BioPerl module for Bio::DB::EUtilities::esearch
#
# Cared for by Chris Fields
#
# Copyright Chris Fields
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code
# 
# Part of the EUtilities BioPerl package

=head1 NAME

Bio::DB::EUtilities::ESearch - Base interface class for handling web queries and data
retrieval from Entrez Utilities from NCBI.  You shouldn't use this class directly.

=head1 SYNOPSIS

=head1 DESCRIPTION

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the
evolution of this and other Bioperl modules. Send
your comments and suggestions preferably to one
of the Bioperl mailing lists. Your participation
is much appreciated.

  bioperl-l@bioperl.org                  - General discussion
  http://bioperl.org/wiki/Mailing_lists  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to
help us keep track the bugs and their resolution.
Bug reports can be submitted via the web.

  http://bugzilla.bioperl.org/

=head1 AUTHOR 

Email cjfields at uiuc dot edu

=head1 APPENDIX

The rest of the documentation details each of the
object methods. Internal methods are usually
preceded with a _

=cut

# Let the code begin...

package Bio::DB::EUtilities::esearch;
use strict;
use warnings;
use Bio::DB::EUtilities;
use Bio::DB::EUtilities::Cookie;

use vars qw(@ISA $EUTIL);

@ISA = qw(Bio::DB::EUtilities);

BEGIN {
    $EUTIL = 'esearch';
}

sub _initialize {
    my ($self, @args ) = @_;
    $self->SUPER::_initialize(@args);
	my ($term, $field, $reldate, $mindate, $maxdate, $datetype, $rettype, $retstart, 
        $retmax, $sort, $usehistory) = 
	  $self->_rearrange([qw(TERM FIELD RELDATE MINDATE MAXDATE DATETYPE RETTYPE
        RETSTART RETMAX SORT USEHISTORY)],
		@args);    
    # set by default
    $self->_eutil($EUTIL);
    $datetype ||= 'mdat';
    $self->datetype($datetype) if $datetype;
	$term			&& $self->term($term);
	$field			&& $self->field($field);
	$reldate		&& $self->reldate($reldate);
	$mindate		&& $self->mindate($mindate);
	$maxdate		&& $self->maxdate($maxdate);
    $retstart       && $self->retstart($retstart);
    $retmax         && $self->retmax($retmax);
    $rettype        && $self->rettype($rettype);
    $sort           && $self->sort_results($sort);
	$usehistory		&& $self->usehistory($usehistory);
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
    my ($webenv, $querykey, $count);
    if (my ($warning) = $content =~ m!<ERROR>(.+)</ERROR>!s) {
        $self->warn("Warning(s) from GenBank: $warning\n");
    }
    if (my ($error) = $content =~ /<OutputMessage>([^<]+)/) {
        $self->throw("Error from Genbank: $error");
    }
    if ($self->usehistory && $self->usehistory eq 'y') {
		($count) = $content =~  /<Count>(\d+)/;
		$self->count($count);
        ($webenv) = $content =~ m!<WebEnv>(\S+)</WebEnv>!;
        ($querykey) = $content =~ m!<QueryKey>(\d+)!g;
		my $cookie = Bio::DB::EUtilities::Cookie->new(
										 -webenv    => $webenv,
										 -querykey  => $querykey,
										 -eutil     => 'esearch',
										 -description   => $self->term,
										 -total		=> $count
										);
        $self->add_cookie($cookie);
	}
}

1;
__END__