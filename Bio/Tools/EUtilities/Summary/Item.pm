# $Id$
#
# BioPerl module for Bio::Tools::EUtilities::Summary::Item
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

Bio::Tools::EUtilities::Summary::Item - simple layered object for DocSum item data

=head1 SYNOPSIS

  # Items can be nested up to three levels at this time. These levels can be
  # accessed via Item, ListItem, or StructureItem methods:

  while (my $item = $docsum->next_Item) {
     print "Name: ",$item->get_name,"\n";
     print "Data: ",$item->get_content,"\n";
     print "Type: ",$item->get_type,"\n";
     while (my $ls = $item->next_ListItem) {
        # do same here
        while (my $struct = $ls->next_StructureItem) {
           # do more stuff here
        }
     }
  }


=head1 DESCRIPTION

DocSum data, as returned from esummary, normally is a simple list of
item-content-content_type groups. However, items can also contain nested data to
represent more complex data (such as structural data). Up to three nested layers
may appear in any document summary.

This class contains methods to access data that can appear in a docsum for any
individual item as well as describes methods to traverse the hierarchy of items
present in a document summary.

The unique name for items are accessed via get_name(), the content by
get_content() (if present), and the data type by get_type(). Items can have
ListItems (Item objects with a datatype() 'list'), which in turn can have
StructureItems (Item objects with a datatype of 'structure'). Items are
initially traversed via a DocSum object using next_Item() or obtained all at
once with get_Items(). Similarly, nested Items can be accessed by using
next_ListItem/get_ListItems and next_StructureItem/get_StructureItem.  A
flattened list of items can be accessed with get_all_Items().

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

Email cjfields at bioperl dot org

=head1 APPENDIX

The rest of the documentation details each of the
object methods. Internal methods are usually
preceded with a _

=cut

# Let the code begin...

package Bio::Tools::EUtilities::Summary::Item;
use strict;
use warnings;

use base qw(Bio::Root::Root Bio::Tools::EUtilities::EUtilDataI);

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
    $type ||= 'item';
    $self->datatype($type);
    $self->eutil('esummary');
    $self->rewind('recursive');
    return $self;
}

=head2 get_ids

 Title    : get_ids
 Usage    : my ($id) = $item->get_ids;
 Function : returns array or array ref with id
 Returns  : array or array ref
 Args     : none
 Note     : the behavior of this method is to remain consistent with other 
            implementations of get_ids(). To retrieve the single DocSum ID use
            get_id()

=cut

sub get_ids {
    my $self = shift;
    return ($self->{'_id'});
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

=head2 next_ListItem

 Title    : next_ListItem
 Usage    : while (my $ls = $item->next_ListItem) {...}
 Function : iterates through ListItems (nested layer of Item)
 Returns  : single ListItem
 Args     : none

=cut

sub next_ListItem {
    my $self = shift;
    unless ($self->{'_lists_it'}) {
        my @lists = $self->get_ListItems;
        # reset the structure iterator (required!)
        delete $self->{'_structures_it'} if $self->{'_structures_it'};
        $self->{'_lists_it'} = sub {return shift @lists}
    }        
    return $self->{'_lists_it'}->();
}

=head2 get_ListItems

 Title    : get_ListItems
 Usage    : my @ls = $item->get_ListItems
 Function : returns list of, well, List Items
 Returns  : array of List Items
 Args     : none

=cut

sub get_ListItems {
    my $self = shift;
    my @items = $self->get_type eq 'List' ? $self->get_subItems : ();
    return @items;
}

=head2 next_StructureItem

 Title    : next_StructureItem
 Usage    : while (my $struc = $ls->next_StructureItem) {...}
 Function : iterates through StructureItems (nested layer of ListItem)
 Returns  : single StructureItems
 Args     : none

=cut

sub next_StructureItem {
    my $self = shift;
    unless ($self->{'_structures_it'}) {
        my @structs = $self->get_StructureItems;
        $self->{'_structures_it'} = sub {return shift @structs}
    }
    return $self->{'_structures_it'}->();
}

=head2 get_StructureItems

 Title    : get_StructureItems
 Usage    : my @structs = $ls->get_StructureItems
 Function : returns list of Structure Items
 Returns  : array of StructureItems
 Args     : none

=cut

sub get_StructureItems {
    my $self = shift;
    my @items = $self->get_type eq 'Structure' ? $self->get_subItems : ();
    return @items;
}

=head2 next_subItem

 Title    : next_subItem
 Usage    : while (my $it = $ls->next_subItem) {...}
 Function : iterates through the next layer of Items
 Returns  : single Item
 Args     : none
 Notes    : unlike next_ListItem and next_Structureitem, this generically
            accesses any sub Items (useful for recursive calls, for example)

=cut

sub next_subItem {
    my $self = shift;
    unless ($self->{'_subitem_it'}) {
        my @structs = $self->get_subItems;
        $self->{'_subitem_it'} = sub {return shift @structs}
    }
    return $self->{'_subitem_it'}->();
}

=head2 get_subItems

 Title    : get_subItems
 Usage    : my @items = $ls->get_subItems
 Function : returns list of sub Items
 Returns  : array of Items
 Args     : none
 Notes    : unlike get_ListItems and get_StructureItems, this generically
            accesses any sub Items (useful for recursive calls, for example)

=cut

sub get_subItems {
    my $self = shift;
    ref $self->{'_items'}  ? return @{ $self->{'_items'} } : return ();
}

=head2 get_name

 Title    : get_name
 Usage    : my $nm = $ls->get_name
 Function : retrieves Item/ListItem/StructureItem name
 Returns  : string
 Args     : none

=cut

sub get_name {
    my $self = shift;
    return $self->{'_itemname'};    
}

=head2 get_type

 Title    : get_type
 Usage    : my $type = $ls->get_type
 Function : retrieves Item/ListItem/StructureItem type 
 Returns  : string
 Args     : none
 Note     : this is not the same as the datatype(), which describes the
            group this Item object belongs to

=cut

sub get_type {
    my $self = shift;
    return $self->{'_itemtype'};
}

=head2 get_content

 Title    : get_content
 Usage    : my $data = $ls->get_content
 Function : retrieves Item/ListItem/StructureItem content (if any)
 Returns  : string
 Args     : none

=cut

sub get_content {
    my $self = shift;
    return $self->{'_itemcontent'};    
}

=head2 rewind

 Title    : rewind
 Usage    : $item->rewind()
 Function : rewinds iterators
 Returns  : none
 Args     : [optional] No arg - only rewinds current layer
           'all' - rewind all DocSum object layers
                (Items, ListItems, StructureItems)

=cut

sub rewind {
    my ($self, $request) = @_;
    if ($request && $request eq 'all') {
        map {$_->rewind()} $self->get_ListItems;
    }
    delete $self->{"_lists_it"} if $self->{"_lists_it"};
    delete $self->{"_structures_it"} if $self->{"_structures_it"};
}

# private data method

sub _add_data {
    my ($self, $data) = @_;
    if ($data->{Item}) {
        my $objtype = lc $data->{Type}.'_item';
        $self->{'_id'} = $data->{Id} if exists $data->{Id};
        for my $sd (@{ $data->{Item} } ) {
            $sd->{Id} = $data->{Id} if exists $data->{Id};
            my $subdoc = Bio::Tools::EUtilities::Summary::Item->new(
                                -datatype => $objtype,
                                -verbose => $self->verbose);
            $subdoc->_add_data($sd);
            push @{ $self->{'_items'} }, $subdoc;
        }
    }
    for my $nm (qw(Type content Name)) {
        $self->{'_item'.lc $nm} = $data->{$nm} if defined $data->{$nm};
    }
    $self->{'_id'} = $data->{Id} if exists $data->{Id};    
}

=head2 to_string

 Title    : to_string
 Usage    : $foo->to_string()
 Function : converts current object to string
 Returns  : none
 Args     : (optional) simple data for text formatting.  This implementation
            passes an argument for layering Items/subItems
 Note     : Used generically for debugging and print_DocSums methods

=cut

# recursively called to grab subitems, then layer

sub to_string {
    my $self = shift;
    my $level = shift || 0;
    # this is the field length for the initial data (spaces are padded in front)
    my $pad = 20 - $level;
    my $content = $self->get_content || '';
    my $string .= sprintf("%-*s%-*s%s\n",
        $level, '',
        $pad, $self->get_name(),
        $self->_text_wrap(':',
             ' ' x ($pad).':',
             $content));
    for my $sub ($self->get_subItems) {
        $string .= $sub->to_string(4 + $level);
    }
    return $string;
}

1;

