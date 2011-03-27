#
# BioPerl module for Bio::Tools::EUtilities::History
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
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

Bio::Tools::EUtilities::History - lightweight implementation of HistoryI
interface (not bound to filehandles, extraneous methods, etc).

=head1 SYNOPSIS

  #should work for any class which is-a HistoryI

  if ($obj->has_History) {
      # do something here
  }

  ($webenv, $querykey) = $obj->history;

  $obj->history($webenv, $querykey);

  $webenv = $obj->get_webenv;

  $query_key = $obj->get_query_key;

=head1 DESCRIPTION

This class extends methods for any EUtilDataI implementation allow instances to
dealwith NCBI history data (WebEnv and query_key).  These can be used as
parameters for further queries against data sets stored on the NCBI server, much
like NCBI's Entrez search history. These are important when one wants to run
complex queries using esearch, retrieve related data using elink, and retrieve
large datasets using epost/efetch.

This class is the simplest implementation and merely holds data for future
queries from any HistoryI. See also Bio::Tools::EUtilities::Query (esearch) and
Bio::Tools::EUtilities::LinkSet (elink), which also implement HistoryI.

=cut

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other Bioperl
modules. Send your comments and suggestions preferably to one of the Bioperl
mailing lists. Your participation is much appreciated.

  bioperl-l@lists.open-bio.org               - General discussion
  http://www.bioperl.org/wiki/Mailing_lists  - About the mailing lists

=head2 Support 

Please direct usage questions or support issues to the mailing list:

I<bioperl-l@bioperl.org>

rather than to the module maintainer directly. Many experienced and 
reponsive experts will be able look at the problem and quickly 
address it. Please include a thorough description of the problem 
with code and data examples if at all possible.

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track the bugs
and their resolution. Bug reports can be submitted via the web.

  https://redmine.open-bio.org/projects/bioperl/

=head1 AUTHOR 

Email cjfields at bioperl dot org

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal
methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::Tools::EUtilities::History;
use strict;
use warnings;

use base qw(Bio::Root::Root Bio::Tools::EUtilities::HistoryI);

sub new {
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(@args);
    my ($eutil) = $self->_rearrange([qw(eutil)],@args);
    $eutil || $self->throw('eutil not defined');
    $self->eutil($eutil);
    $self->datatype('history');
    return $self;
}

=head2 history

 Title    : history
 Usage    : my ($webenv, $qk) = $hist->history
 Function : Get/Set two-element list of webenv() and query_key()
 Returns  : array
 Args     : two-element list of webenv, query key

=cut

=head2 get_webenv

 Title    : get_webenv
 Usage    : my $webenv = $hist->get_webenv
 Function : returns web environment key needed to retrieve results from
            NCBI server
 Returns  : string (encoded key)
 Args     : none

=cut

=head2 get_query_key

 Title    : get_query_key
 Usage    : my $qk = $hist->get_query_key
 Function : returns query key (integer) for the history number for this session
 Returns  : integer
 Args     : none

=cut

sub _add_data {
    my ($self, $simple) = @_;
    if (!exists $simple->{WebEnv} || !exists $simple->{QueryKey}) {
        $self->debug("Data:",Dumper($simple));
        $self->throw("Missing webenv/query key in history output");
    }
    $self->{_webenv} = $simple->{WebEnv} && delete $simple->{WebEnv};
    $self->{_querykey} = $simple->{QueryKey} && delete $simple->{QueryKey};
}

=head2 to_string

 Title    : to_string
 Usage    : $foo->to_string()
 Function : converts current object to string
 Returns  : none
 Args     : (optional) simple data for text formatting
 Note     : Used generally for debugging and for the print_* methods

=cut

sub to_string {
    my $self = shift;
    my $string;
    my %map = (
        'get_webenv' => 'WebEnv',
        'get_query_key' => 'Key'
        );
    for my $m (qw(get_webenv get_query_key)) {
        $string .= sprintf("%-20s:%s\n", $map{$m}, $self->$m);
    }
    return $string;
}

1;
