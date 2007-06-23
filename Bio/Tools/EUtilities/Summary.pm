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

=head2 get_ids

 Title    : get_ids
 Usage    : @ids = $esum->get_ids
 Function : returns array or array ref of IDs
 Returns  : array or array ref of IDs (depending on wantarray)
 Args     : none

=cut

sub get_ids {
    my $self = shift;
    unless (exists $self->{'_id'}) {
        push @{$self->{'_id'}}, map {$_->get_id } $self->get_DocSums;
    }
    return wantarray ? @{$self->{'_id'}} : $self->{'_id'};
}

=head2 next_DocSum

 Title    : next_DocSum
 Usage    : while (my $ds = $esum->next_DocSum) {...}
 Function : iterate through DocSum instances
 Returns  : single Bio::Tools::EUtilities::Summary::DocSum
 Args     : none yet

=cut

# add an option (?) for lazy parsing via fh (which allows tempfile or piping)
# add data to Simple object in chunks, create DocSum, return

# maybe allow callback to only return interesting DocSums?

sub next_DocSum {
    my $self = shift;
    $self->parse_data unless $self->data_parsed;
    unless ($self->{"_docsums_it"}) {
        my @docsums = $self->get_DocSums;
        $self->{"_docsums_it"} = sub {return shift @docsums}
    }
    $self->{'_docsums_it'}->();
}

=head2 get_DocSums

 Title    : get_DocSums
 Usage    : my @docsums = $esum->get_DocSums
 Function : retrieve a list of DocSum instances
 Returns  : array of Bio::Tools::EUtilities::Summary::DocSum
 Args     : none

=cut

sub get_DocSums {
    my $self = shift;
    $self->parse_data unless $self->data_parsed;
    # always return a list for consistency...
    return ref $self->{'_docsums'} ? @{ $self->{'_docsums'} } : return ();
}

=head2 rewind

 Title    : rewind
 Usage    : $esum->rewind()
            $esum->rewind('recursive')
 Function : retrieve a list of DocSum instances
 Returns  : array of Bio::Tools::EUtilities::Summary::DocSum
 Args     : [optional]
           'recursive' - rewind all DocSum object layers
                         (Items, ListItems, StructureItems)

=cut

sub rewind {
    my ($self, $request) = @_;
    if ($request && $request eq 'recursive') {
        map {$_->rewind('recursive') } $self->get_DocSums;
    }
    delete $self->{"_docsums_it"};
}

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
    $self->{'_parsed'} = 1;
}

1;

__END__