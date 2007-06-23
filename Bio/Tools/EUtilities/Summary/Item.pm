# $Id$
#
# BioPerl module for Bio::DB::EUtilities::Summary::Item
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

Bio::DB::EUtilities::Summary::Item - simple layered object for DocSum item data

=head1 SYNOPSIS

  # Items can be nested up to three levels at this time (Item, ListItem,
  # StructureItem).
  
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
next_ListItem/get_ListItems and next_StructureItem/get_StructureItem.

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
            implementations of get_ids().  To retrieve the single ID use get_item_id()

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
 Function : returns list of, well, ListItems
 Returns  : array of ListItems
 Args     : none

=cut

sub get_ListItems {
    my $self = shift;
    ref $self->{'_lists'}  ? return @{ $self->{'_lists'} } : return ();
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
 Function : returns list of StructureItems
 Returns  : array of StructureItems
 Args     : none

=cut

sub get_StructureItems {
    my $self = shift;
    ref $self->{'_structures'}  ? return @{ $self->{'_structures'} } : return ();
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
            group this Item ojbect belongs to
 
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
 Args     : [optional]
           'recursive' - rewind all DocSum object layers
                         (Items, ListItems, StructureItems)

=cut

sub rewind {
    my ($self, $request) = @_;
    if ($request && $request eq 'recursive') {
        map {$_->rewind()} $self->get_ListItems;
    }
    delete $self->{"_lists_it"} if $self->{"_lists_it"};
    delete $self->{"_structures_it"} if $self->{"_structures_it"};
}

# private data method

sub _add_data {
    my ($self, $data) = @_;
    if ($data->{Item}) {
        my $objtype = lc $data->{Type};
        $self->{'_id'} = $data->{Id} if exists $data->{Id};
        for my $sd (@{ $data->{Item} } ) {
            $sd->{Id} = $data->{Id} if exists $data->{Id};
            my $subdoc = Bio::Tools::EUtilities::Summary::Item->new(
                                -datatype => $objtype,
                                -verbose => $self->verbose);
            $subdoc->_add_data($sd);
            push @{ $self->{'_'.lc $objtype.'s'} }, $subdoc;
        }
    }
    for my $nm (qw(Type content Name)) {
        $self->{'_item'.lc $nm} = $data->{$nm} if $data->{$nm};
    }
    $self->{'_id'} = $data->{Id} if exists $data->{Id};    
}

1;

