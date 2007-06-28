# $Id$
#
# BioPerl module for Bio::DB::EUtilities::DocSum
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

Bio::DB::EUtilities::Summary

=head1 SYNOPSIS

  #### should not create instance directly; Bio::Tools::EUtilities does this ####
    
  my $esum = Bio::Tools::EUtilities->new(-eutil => 'esummary',
                                         -file => 'summary.xml');
  # can also use '-response' (for HTTP::Response objects) or '-fh' (for filehandles)
  
  while (my $docsum = $esum->next_DocSum) {
      my $id = $docsum->get_ids;  # EUtilDataI compliant method, returns docsum ID
      my @names = $docsum->get_item_names;
  }

=head1 DESCRIPTION

This class handles data output (XML) from esummary.  

esummary retrieves information in the form of document summaries (docsums) when
passed a list of primary IDs or if using a previous search history.

This module breaks down the returned data from esummary into individual document
summaries per ID (using a DocSum object). As the data in a docsum can be nested,
subclasses of DocSums (Item, ListItem, Structure) are also present. 

Further documentation for Link and Field subclass methods is included below.

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

package Bio::Tools::EUtilities::Summary;
use strict;
use warnings;

use base qw(Bio::Tools::EUtilities Bio::Tools::EUtilities::EUtilDataI);

use Bio::Tools::EUtilities::Summary::DocSum;

# private EUtilDataI method

sub _add_data {
    my ($self, $data) = @_;
    if (!exists $data->{DocSum}) {
        $self->warn('No returned docsums.');
        return;
    }
    
    my @docs;
    for my $docsum (@{ $data->{DocSum} }) {
        my $ds = Bio::Tools::EUtilities::Summary::DocSum->new(-datatype => 'docsum',
                                                              -verbose => $self->verbose);
        $ds->_add_data($docsum);
        push @{ $self->{'_docsums'} }, $ds;
    }
}

1;

__END__