# $Id$
#
# BioPerl module for Bio::DB::EUtilities::Summary::DocSum
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

Bio::DB::EUtilities::Summary::DocSum - data object for document summary data
from esummary

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

=head1 AUTHOR Chris Fields

Email cjfields at uiuc dot edu

=head1 APPENDIX

The rest of the documentation details each of the
object methods. Internal methods are usually
preceded with a _

=cut

# Let the code begin...

package Bio::Tools::EUtilities::Summary::DocSum;

use strict;
use warnings;
use base qw(Bio::Root::Root Bio::Tools::EUtilities::EUtilDataI);

use Bio::Tools::EUtilities::Summary::Item;

=head2 new

 Title    : new
 Usage    : 
 Function : 
 Returns  : 
 Args     : 

=cut

sub new {
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(@args);
    my ($type) = $self->_rearrange(['DATATYPE'],@args);
    $type ||= 'docsum';
    $self->eutil('esummary');
    $self->datatype($type);
    return $self;
}

=head2 get_ids

 Title    : get_ids
 Usage    : my ($id) = $item->get_ids;
 Function : returns array or array ref with id
 Returns  : array or array ref
 Args     : none
 Note     : the behavior of this method remains consistent with other
            implementations of get_ids(). To retrieve the single DocSum ID
            use get_id()

=cut

sub get_ids {
    my $self = shift;
    return wantarray ? $self->{'_id'} : [$self->{'_id'}];
}

=head2 get_id

 Title    : get_id
 Usage    : my ($id) = $item->get_id;
 Function : returns UID of record
 Returns  : integer
 Args     : none

=cut

sub get_id {
    my $self = shift;
    return $self->{'_id'};
}

=head2 next_Item

 Title    : next_Item
 Usage    : while (my $item = $docsum->next_Item) {...}
 Function : iterates through Items (nested layer of Item)
 Returns  : single Item
 Args     : [optional] single arg (string)
            'flatten' - iterates through a flattened list ala
                          get_all_DocSum_Items()

=cut

sub next_Item {
    my ($self, $request) = @_;
    unless ($self->{"_items_it"}) {
        #my @items = $self->get_Items;
        my @items = ($request && $request eq 'flatten') ?
                    $self->get_all_Items :
                    $self->get_Items ;
        $self->{"_items_it"} = sub {return shift @items}
    }
    $self->{'_items_it'}->();
}

=head2 get_Items

 Title    : get_Items
 Usage    : my @items = $docsum->get_Items
 Function : returns list of, well, Items
 Returns  : array of Items
 Args     : none

=cut

sub get_Items {
    my $self = shift;
    return ref $self->{'_items'} ? @{ $self->{'_items'} } : return ();
}

=head2 get_all_Items

 Title    : get_all_Items
 Usage    : my @items = $docsum->get_all_Items
 Function : returns flattened list of all Item objects (Items, ListItems,
            StructureItems)
 Returns  : array of Items
 Args     : none
 Note     : items are added top-down (similar order to using nested calls)
            in original list order.

             1         2        7        8
           Item  -   Item  -  Item  -  Item ...
                     |
                    | 3        6
                 ListItem - ListItem
                   |
                  | 4          5
               Structure - Structure

=cut

sub get_all_Items {
    my $self = shift;
    unless ($self->{'_ordered_items'}) {
        for my $item ($self->get_Items) {
            push @{$self->{'_ordered_items'}}, $item;
            for my $ls ($item->get_ListItems) {
                push @{$self->{'_ordered_items'}}, $ls;
                for my $st ($ls->get_StructureItems) {
                    push @{$self->{'_ordered_items'}}, $st;                
                } 
            }
        }
    }
    return @{$self->{'_ordered_items'}};
}

=head2 get_content_by_name

 Title    : get_content_by_Item_name
 Usage    : my $data = get_content_by_name('CreateDate')
 Function : Returns scalar content for named Item in DocSum (indicated by
            passed argument)
 Returns  : scalar value (string) if present
 Args     : string (Item name)
 Warns    : If Item with name is not found 

=cut

sub get_content_by_name {
    my ($self, $key) = @_;
    return unless $key;
    my ($it) = grep {$_->get_name eq $key} $self->get_all_Items;
    return $it->get_content;
}

=head2 get_type_by_name

 Title    : get_type_by_name
 Usage    : my $data = get_type_by_name('CreateDate')
 Function : Returns data type for named Item in DocSum (indicated by
            passed argument)
 Returns  : scalar value (string) if present
 Args     : string (Item name)
 Warns    : If Item with name is not found 

=cut

sub get_type_by_name {
    my ($self, $key) = @_;
    return unless $key;
    my ($it) = grep {$_->get_name eq $key} $self->get_all_Items;
    return $it->get_type;
}

=head2 rewind

 Title    : rewind
 Usage    : $docsum->rewind();
 Function : rewinds DocSum iterator
 Returns  : none
 Args     : [optional]
           'recursive' - rewind all DocSum object layers
                         (Items, ListItems, StructureItems)

=cut

sub rewind {
    my ($self, $request) = @_;
    if ($request && $request eq 'all') {
        map {$_->rewind('all') } $self->get_Items;
    }
    delete $self->{"_items_it"};
}

# private EUtilDataI method

sub _add_data {
    my ($self, $data) = @_;
    if ($data->{Item}) {
        $self->{'_id'} = $data->{Id} if exists $data->{Id};
        for my $sd (@{ $data->{Item} } ) {
            $sd->{Id} = $data->{Id} if exists $data->{Id};
            my $subdoc = 
                Bio::Tools::EUtilities::Summary::Item->new(-datatype => 'item',
                                                  -verbose => $self->verbose);
            $subdoc->_add_data($sd);
            push @{ $self->{'_items'} }, $subdoc;
        }
    }
    $self->{'_id'} = $data->{Id} if exists $data->{Id};
}

1;