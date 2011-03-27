#
# BioPerl module for Bio::Tools::EUtilities::HistoryI
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

Bio::Tools::EUtilities::HistoryI - simple extension of EUtilDataI interface
class for classes which hold NCBI server history data

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

The simplest implementation is Bio::Tools::EUtilities::History, which holds the
history data for epost.  See also Bio::Tools::EUtilities::Query (esearch) and
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

package Bio::Tools::EUtilities::HistoryI;
use strict;
use warnings;

use base qw(Bio::Tools::EUtilities::EUtilDataI);

=head2 history

 Title    : history
 Usage    : my ($webenv, $qk) = $hist->history
 Function : Get/Set two-element list of webenv() and query_key()
 Returns  : array
 Args     : two-element list of webenv, querykey

=cut

sub history {
    my $self = shift;
    $self->parse_data if ($self->can('parse_data') && !$self->data_parsed);
    if (@_) {
        my ($webenv, $querykey) = (shift, shift);
        $self->throw("Missing part of cookie!") if (!$webenv || !$querykey);
        ($self->{'_webenv'}, $self->{'_querykey'}) = ($webenv, $querykey);
    }
    return ($self->get_webenv, $self->get_query_key);
}

=head2 get_webenv

 Title    : get_webenv
 Usage    : my $webenv = $hist->get_webenv
 Function : returns web environment key needed to retrieve results from
            NCBI server
 Returns  : string (encoded key)
 Args     : none

=cut

sub get_webenv {
    my $self = shift;
    $self->parse_data if ($self->can('parse_data') && !$self->data_parsed);
    return $self->{'_webenv'};
}

=head2 get_query_key

 Title    : get_query_key
 Usage    : my $qk = $hist->get_query_key
 Function : returns query key (integer) for the history number for this session
 Returns  : integer
 Args     : none

=cut

sub get_query_key {
    my $self = shift;
    $self->parse_data if ($self->can('parse_data') && !$self->data_parsed);    
    return $self->{'_querykey'};
}

1;
__END__
