#
# BioPerl module for Bio::Tools::EUtilities::DocSum
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

Bio::Tools::EUtilities::Summary - class for handling data output (XML) from
esummary.

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

=head2 Support 

Please direct usage questions or support issues to the mailing list:

I<bioperl-l@bioperl.org>

rather than to the module maintainer directly. Many experienced and 
reponsive experts will be able look at the problem and quickly 
address it. Please include a thorough description of the problem 
with code and data examples if at all possible.

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to
help us keep track the bugs and their resolution.
Bug reports can be submitted via the web.

  https://redmine.open-bio.org/projects/bioperl/

=head1 AUTHOR 

Email cjfields at bioperl dot org

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
    my %data = (
        'DB'    => [1, join(', ',$self->get_databases) || ''],
    );
    my $string = $self->SUPER::to_string."\n";
    for my $k (sort {$data{$a}->[0] <=> $data{$b}->[0]} keys %data) {
        $string .= sprintf("%-20s:%s\n\n",$k, $self->_text_wrap('',' 'x 20 .':', $data{$k}->[1]));
    }
    while (my $ds = $self->next_DocSum) {
        $string .= $ds->to_string."\n";
    }
    return $string;
}

1;

__END__
