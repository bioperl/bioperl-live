# $Id$
#
# BioPerl module for Bio::Tools::EUtilities::History
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

Bio::Tools::EUtilities::History - simple HistoryI implementation for holding
NCBI search history.

=head1 SYNOPSIS

  #### should not create instance directly; Bio::Tools::EUtilities does this ####
  
  my $hist = Bio::Tools::EUtilities->new(-eutil => 'epost',
                                         -file => 'epost.xml');
  
  if ($hist->has_History) {
      # do something here
  }
  
  ($webenv, $querykey) = $hist->get_history;
  
  $webenv = $hist->get_webenv;

  $query_key = $hist->get_query_key;

=head1 DESCRIPTION

This class is the simplest Bio::Tools::EUtilities::HistoryI implementation and
is primarily used for epost. This class adds methods to objects for dealing with
NCBI history data (WebEnv and query_key). These can be used as parameters for
futher queries against data sets stored on the NCBI server, much like NCBI's
Entrez search history.

NOTE: I plan on testing how long this data is stored remotely. I believe NCBI
deletes remote data after 24 hours.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other Bioperl
modules. Send your comments and suggestions preferably to one of the Bioperl
mailing lists. Your participation is much appreciated.

  bioperl-l@lists.open-bio.org               - General discussion
  http://www.bioperl.org/wiki/Mailing_lists  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track the bugs
and their resolution. Bug reports can be submitted via the web.

  http://bugzilla.open-bio.org/

=head1 AUTHOR 

Email cjfields at uiuc dot edu

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal
methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::Tools::EUtilities::History;
use strict;
use warnings;
#use URI::Escape qw(uri_unescape); 

use base qw(Bio::Tools::EUtilities);

=head1 Bio::Tools::EUtilities::History relevant methods

These are inherited through the Bio::Tools::EUtilities interface
Bio::Tools::EUtilities::HistoryI

=head2 history

 Title    : history
 Usage    : my ($webenv, $qk) = $hist->history
 Function : Get/Set two-element list of webenv() and query_key()
 Returns  : array
 Args     : two-element list of webenv, querykey

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

=head2 has_History

 Title    : has_History
 Usage    : if ($hist->has_History) {...}
 Function : returns TRUE if full history (webenv, query_key) is present 
 Returns  : BOOLEAN, value eval'ing to TRUE or FALUE
 Args     : none

=cut

# private method

sub _add_data {
    my ($self, $simple) = @_;
    map { $self->{'_'.lc $_} = $simple->{$_} } keys %$simple;
}

1;
__END__
