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

Bio::DB::EUtilities::DocSum 

=head1 SYNOPSIS

# use only in conjunction with Bio::DB::EUtilities::esummary

    my $esum = Bio::DB::EUtilities->new(-verbose => 1,
                                        -eutil  => 'esummary',
                                        -cookie => $esearch->next_cookie,
                                        -retmax => 20
                                         );
    
    $esum->get_response;
    
    # get docsum objects
    while (my $docsum  = $esum->next_docsum) {
        # do stuff here
    }

=head1 DESCRIPTION

This is a remedial object that acts as a container for DocSum data from
Bio::DB::EUtilities::esummary.  It is in the very early stages of
development, so don't too be offended if the API changes.  It is possible
the various EUtilities container objects will be reorganized to have a
more consistent API; however, note that, due to the differences in the
actual data this may be next to impossible.

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

package Bio::DB::EUtilities::DocSum;
use strict;
use warnings;
use Data::Dumper;

use base qw(Bio::Root::Root);

sub new {
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(@args);
    $self->{'_itemindex'} = 0;
    $self->{'_docdata'} = [];
    return $self;
}

# this should make a deep copy of the XML data for each docsum

sub _add_data {
    my ($self, $docsum) = @_;
    if (!$docsum || ref($docsum) !~ /array/i) {
        $self->throw("Bad ESummary DocSum");
    }
    #$self->debug("Docsum: ".Dumper($docsum));
    for my $item (@ {$docsum} ) {
        my ($name, $type, $content) = ($item->{Name},
                                       $item->{Type},
                                       $item->{content});
        # layered docsum data
        if ($type eq 'List' && exists $item->{Item}) {
            my $ds = $self->new(-verbose => $self->verbose);
            $ds->_add_data($item->{Item});
            $content = $ds;
        }
        push @{$self->{'_docdata'}}, {'Name' => $name,
                                      'Type' => $type,
                                      'Content' => $content || ''};
    }
    return;
}

=head2 esummary_id

 Title   : esummary_id
 Usage   : $id = $docsum->esummary_id();
 Function: get/set ID value for DocSum object
 Returns : UID for DocSum object
 Args    : OPTIONAL : UID to set docsum object 

=cut

sub esummary_id {
    my $self = shift;
    return $self->{'_esum_id'} = shift if @_;
    return $self->{'_esum_id'};
}

=head2 get_all_names

 Title   : get_all_names
 Usage   : @names = $docsum->get_all_names;
 Function: get array of DocSum item names
 Returns : array of names for the items in DocSum object
 Args    : none

=cut

sub get_all_names {
    my $self = shift;
    my @names = map {$_->{Name}} @{ $self->{'_docdata'} };
    return @names;
}

=head2 get_item_by_name

 Title   : get_item_by_name
 Usage   : %item = $docsum->get_item_by_name($name);
 Function: retrieve docsum item hash by item name
           (retrieved via get_all_names())
 Returns : hash containing all information for the DocSum item
 Args    : REQUIRED: name of item to be retrieved

=cut

# the item should have a unique name, so a grep should work fine
sub get_item_by_name {
    my ($self, $name) = @_;
    $self->throw('Must supply name for get_data_by_name') if !$name;
    my ($data) = grep {$_->{Name} eq $name} @{ $self->{'_docdata'} };
    return %{ $data } if $data;
    return;
}

=head2 get_Type_by_name

 Title   : get_Type_by_name
 Usage   : $type = $docsum->get_Type_by_name($name);
 Function: retrieve the type of information held by docsum item 
           (retrieved via get_all_names())
 Returns : string 
 Args    : REQUIRED: name of item to be retrieved

=cut

sub get_Type_by_name {
    my ($self, $name) = @_;
    $self->throw('Must supply name for get_Type_by_name') if !$name;
    my ($data) = grep {$_->{Name} eq $name} @{ $self->{'_docdata'} };
    return $data->{Type} if exists $data->{Type};
    return;
}

=head2 get_Content_by_name

 Title   : get_Content_by_name
 Usage   : $type = $docsum->get_Content_by_name($name);
 Function: retrieve the Content held by docsum item
           (retrieved via get_all_names())
           
           NOTE:  the content is normally a string; however, items
           that have the Type 'List' will have more complex content
           containing subitems, subnames, and subcontent;
           these are represented by sub-DocSum objects.
           
           Unless you know that a particular DocSum item always
           has a 'String' or 'Integer' Type, it might
           be safe to preceed this method with calls to
           get_Type_by_name() to verify their content type
           
 Returns : string or DocSum object reference
 Args    : REQUIRED: name of item to be retrieved

=cut

sub get_Content_by_name {
    my ($self, $name) = @_;
    $self->throw('Must supply name for get_Content_by_name') if !$name;
    my ($data) = grep {$_->{Name} eq $name} @{ $self->{'_docdata'} };
    return $data->{Content} if exists $data->{Content};
    return;
}

=head2 next_docsum_item

 Title   : next_docsum_item
 Usage   : while ($docsum->next_docsum_item) {;
 Function: set the index value for the next item in the DocSum list
 Returns : hash containing docsum item data (Name, Type, Content)
 Args    : none

=cut

sub next_docsum_item {
    my $self = shift;
    my $index = $self->_next_item_index;
    if (exists $self->{'_docdata'}->[$index]) {
        return %{ $self->{'_docdata'}->[$index] };
    } else {
        return;
    }
}

=head2 rewind_docsum_items

 Title   : rewind_docsum_items
 Usage   : $docsum->rewind_docsum_items();
 Function: rewind the item index to the beginning
           (iterated via next_docsum_item) 
 Returns : none
 Args    : none

=cut

sub rewind_docsum_items{
    my $self = shift;
    $self->{'_itemindex'} = 0;
    return;
}

sub _next_item_index{
    my $self = shift;
    return $self->{'_itemindex'}++;
}

1;

__END__