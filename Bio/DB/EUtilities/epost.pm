# $Id$
#
# BioPerl module for Bio::DB::EUtilities::epost
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

Bio::DB::EUtilities::epost

=head1 SYNOPSIS

=head1 DESCRIPTION

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the
evolution of this and other Bioperl modules. Send
your comments and suggestions preferably to one
of the Bioperl mailing lists. Your participation
is much appreciated.

  bioperl-l@lists.open-bio.org               - General discussion
  http://www.bioperl.org/wiki/Mailing_lists  - About the mailing lists
  
=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to
help us keep track the bugs and their resolution.
Bug reports can be submitted via the web.

  http://bugzilla.open-bio.org/

=head1 AUTHOR 

Email cjfields at uiuc dot edu

=head1 APPENDIX

The rest of the documentation details each of the
object methods. Internal methods are usually
preceded with a _

=cut

# Let the code begin...

package Bio::DB::EUtilities::epost;
use strict;
use warnings;
use Bio::DB::EUtilities;
use Bio::DB::EUtilities::Cookie;
use XML::Simple;
#use Data::Dumper;

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
    $self->_eutil($EUTIL);
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
    my $xs = XML::Simple->new();
    my $simple = $xs->XMLin($response->content);
    #$self->debug("Response dumper:\n".Dumper($simple));
    # check for errors
    if ($simple->{ERROR}) {
        $self->throw("NCBI epost nonrecoverable error: ".$simple->{ERROR});
    }
    if ($simple->{InvalidIdList}) {
        $self->warn("NCBI epost error: Invalid ID List".$simple->{InvalidIdList});
    }
    my $db = $self->db;
    my $webenv    = $simple->{WebEnv};
    my $querykey  = $simple->{QueryKey};
    my $cookie = Bio::DB::EUtilities::Cookie->new(-webenv   => $webenv,
                                                  -querykey => $querykey,
                                                  -eutil    => 'epost',
                                                  -database => $db,
                                                  );
    $self->add_cookie($cookie);
    return $response;
}

1;
__END__