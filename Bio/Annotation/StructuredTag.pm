# $Id: StructuredTag.pm 11693 2007-09-17 20:54:04Z cjfields $
#
# BioPerl module for Bio::Annotation::StructuredTag
#
# Cared for Chris Fields
#
# You may distribute this module under the same terms as perl itself.
# Refer to the Perl Artistic License (see the license accompanying this
# software package, or see http://www.perl.com/language/misc/Artistic.html)
# for the terms under which you may use, modify, and redistribute this module.
#
# THIS PACKAGE IS PROVIDED "AS IS" AND WITHOUT ANY EXPRESS OR IMPLIED
# WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED WARRANTIES OF
# MERCHANTIBILITY AND FITNESS FOR A PARTICULAR PURPOSE.
#
# POD documentation - main docs before the code

=head1 NAME

Bio::Annotation::StructuredTag - AnnotationI with simple tag-value relationships

=head1 SYNOPSIS

   use Bio::Annotation::StructuredTag;
   use Bio::Annotation::Collection;

   my $col = Bio::Annotation::Collection->new();
   
   # data structure can be an array reference with a data structure
   # corresponding to that defined by Data::Stag:
   
   my $sv = Bio::Annotation::StructuredTag->new(-tagname => 'mytag1',
                                                -value => $data_structure);
   $col->add_Annotation($sv);
   
   # regular text passed is parsed based on the tagformat().
   my $sv2 = Bio::Annotation::StructuredTag->new(-tagname => 'mytag2',
                                                -tagformat => 'xml',
                                                -value => $xmltext);
   $col->add_Annotation($sv2);
   
   ###### MORE TO FOLLOW ######

=head1 DESCRIPTION

For now, this is an experimental AnnotationI stub to determine whether this can
serve as an adequate replacement for Bio::Annotation::StructuredValue with
BioSQL. YMMV.

This takes tagged data values and stores them in a simple structured tag-based
hierarchy (complements of Chris Mungall's Data::Stag module). Data can then be
represented as text using a variety of output formats (indention, itext, xml,
spxr).  See L<Data::Stag> for details.

Data passed in using value() or the '-value' parameter upon instantiation
can either be:

1) an array reference corresponding to the data structure for Data::Stag;

2) a text string in 'xml', 'itext', 'spxr', or 'indent' format. The default
format is 'xml'; this can be changed using tagformat() prior to using value() or
by passing in the proper format using '-tagformat' upon instantiation.

For now, beyond checking for an array reference no format guessing occurs (so,
for roundtrip tests ensure that the IO formats correspond). Therefore, we
recommend when using text input to set tagformat() to one of these formats prior
to data loading to ensure the proper Data::Stag parser is selected. After data
loading, the tagformat() can be changed to change the text string format
returned by value().

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists. Your participation is much appreciated.

  bioperl-l@bioperl.org                  - General discussion
  http://bioperl.org/wiki/Mailing_lists  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via
or the web:

  http://bugzilla.open-bio.org/

=head1 AUTHOR - Chris Fields

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal
methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::Annotation::StructuredTag;
use strict;

# Object preamble - inherits from Bio::Root::Root

use base qw(Bio::Root::Root Bio::AnnotationI);
use Data::Stag;

=head2 new

 Title   : new
 Usage   : my $sv = Bio::Annotation::StructuredTag->new();
 Function: Instantiate a new StructuredTag object
 Returns : Bio::Annotation::StructuredTag object
 Args    : -value => $value to initialize the object data field [optional]
           -tagname => $tag to initialize the tagname [optional]
           -tagformat => format for output [optional]
                      (types 'xml', 'itext', 'sxpr', 'indent', default = 'itext')

=cut

sub new{
   my ($class,@args) = @_;
   my $self = $class->SUPER::new(@args);   
   my ($node, $value,$tag, $format) = $self->_rearrange([qw(
                                       NODE
                                       VALUE
                                       TAGNAME
                                       TAGFORMAT)], @args);
   defined $tag    && $self->tagname($tag);
   $format ||= 'xml';
   $self->tagformat($format);
   defined $value  && $self->value($value);
   defined $node   && $self->node($node);
   return $self;
}

=head1 AnnotationI implementing functions

=cut

=head2 as_text

 Title   : as_text
 Usage   : my $text = $obj->as_text
 Function: return the string "Value: $v" where $v is the value
 Returns : string
 Args    : none

=cut

sub as_text{
   my ($self) = @_;
   return "Value: ".$self->value;
}

=head2 display_text

 Title   : display_text
 Usage   : my $str = $ann->display_text();
 Function: returns a string. Unlike as_text(), this method returns a string
           formatted as would be expected for the specific implementation.
           
           One can pass a callback as an argument which allows custom text
           generation; the callback is passed the current instance and any text
           returned
 Example :
 Returns : a string
 Args    : [optional] callback

=cut

{
  my $DEFAULT_CB = sub { $_[0]->value || ''};

  sub display_text {
    my ($self, $cb) = @_;
    $cb ||= $DEFAULT_CB;
    $self->throw("Callback must be a code reference") if ref $cb ne 'CODE';
    return $cb->($self);
  }

}

=head2 hash_tree

 Title   : hash_tree
 Usage   : my $hashtree = $value->hash_tree
 Function: For supporting the AnnotationI interface just returns the value
           as a hashref with the key 'value' pointing to the value
 Returns : hashrf
 Args    : none

=cut

sub hash_tree{
   my ($self) = @_;
   my $h = {};
   $h->{'value'} = $self->value;
}

=head2 tagname

 Title   : tagname
 Usage   : $obj->tagname($newval)
 Function: Get/set the tagname for this annotation value.

           Setting this is optional. If set, it obviates the need to provide
           a tag to AnnotationCollection when adding this object.
 Example :
 Returns : value of tagname (a scalar)
 Args    : new value (a scalar, optional)

=cut

################# NEED TO ADD TOP-LEVEL TAGNAME CHANGE #################

sub tagname{
   my ($self,$value) = @_;
   if( defined $value) {
   $self->{'tagname'} = $value;
   }
   return $self->{'tagname'};
}

=head1 Specific accessors for StructuredTag

=cut

=head2 value

 Title   : value
 Usage   : $obj->value($newval)
 Function: Get/set the value for this annotation.
 Returns : value of value
 Args    : newvalue (optional)

=cut

sub value{
   my ($self,$value) = @_;
   # set mode? This resets the entire tagged database
   my $format = $self->tagformat;
   if ($value) {
      if (ref $value eq 'ARRAY') {
         my $tagname = $self->tagname ? $self->tagname : 'unknown';
         eval { $self->{db} = Data::Stag->nodify([$tagname => $value]) };
      } else {
         # not trying to guess here for now; we go by the tagformat() setting
         my $h = Data::Stag->getformathandler($format);
         eval {$self->{db} = Data::Stag->from($format.'str',$value)};
      }
      $self->throw("Data::Stag error:\n$@") if $@;
   }
   # get mode?
   # How do we return a data structure?
   # for now, we use the output (if there is a Data::Stag node present)
   # may need to an an eval {} to catch Data::Stag output errors
   $self->{db} ? $self->{db}->$format : '';
}

=head2 tagformat

 Title   : tagformat
 Usage   : $obj->tagformat($newval)
 Function: Get/set the output tag format for this annotation.
 Returns : value of tagformat
 Args    : newvalue (optional) - format for the data passed into value
           must be of values 'xml', 'indent', 'sxpr', 'itext'

=cut

my %IS_VALID_FORMAT = map {$_ => 1} qw(xml indent sxpr itext);

sub tagformat{
   my ($self,$value) = @_;
   if( defined $value && $IS_VALID_FORMAT{$value}) {
      $self->{'tagformat'} = $value;
   }
   return $self->{'tagformat'};
}

=head2 node

 Title   : node
 Usage   : $obj->node()
 Function: Get/set the Data::Stag node used for this annotation
 Returns : Data::Stag node implementation instance
           (default is Data::Stag::StagImpl)
 Args    : (optional) Data::Stag node implementation instance
 
=cut

sub node{
   my ($self,$value) = @_;
   if( defined $value && ref $value && $value->isa('Data::Stag::StagI')) {
      $self->{'db'} = $value;
   }
   return $self->{'db'};
}

################# NEED TO REIMPLEMENT USING DATA::STAG #################

#=head2 get_values
#
# Title   : get_values
# Usage   :
# Function: Get the top-level array of values. Each of the elements will
#           recursively be a reference to an array or a scalar, depending
#           on the depth of this structured value annotation.
# Example :
# Returns : an array
# Args    : none
#
#
#=cut
#
#sub get_values{
#    my $self = shift;
#
#    return @{$self->{'values'}};
#}

#=head2 get_all_values
#
# Title   : get_all_values
# Usage   :
# Function: Flattens all values in this structured annotation and
#           returns them as an array.
# Example :
# Returns : the (flat) array of values
# Args    : none
#
#
#=cut
#
#sub get_all_values{
#    my ($self) = @_;
#    # we code lazy here and just take advantage of value()
#    my $txt = $self->value(-joins => ['@!@'], -brackets => ['','']);
#    return split(/\@!\@/, $txt);
#}
#
#=head2 add_value
#
# Title   : add_value
# Usage   :
# Function: Adds the given value to the structured annotation at the
#           given index.
#
#           The index is multi-dimensional, with the first dimension
#           applying to the first level, and so forth. If a particular
#           dimension or a particular index does not exist yet, it will
#           be created. If it does exist and adding the value would
#           mean replacing a scalar with an array reference, we throw
#           an exception to prevent unintended damage. An index of -1
#           at any dimension means append.
#
#           If an array of values is to be added, it will create an
#           additional dimension at the index specified, unless the
#           last index value is -1, in which case they will all be
#           appended to the last dimension.
#
# Example :
# Returns : none
# Args    : the index at which to add (a reference to an array)
#           the value(s) to add
#
#=cut
#
#sub add_value{
#    my ($self,$index,@values) = @_;
#
#    my $tree = $self->{'values'};
#    my $lastidx = pop(@$index);
#    foreach my $i (@$index) {
#   if($i < 0) {
#       my $subtree = [];
#       push(@$tree, $subtree);
#       $tree = $subtree;
#   } elsif((! $tree->[$i]) || (ref($tree->[$i]) eq "ARRAY")) {
#       $tree->[$i] = [] unless ref($tree->[$i]) eq "ARRAY";
#       $tree = $tree->[$i];
#   } else {
#       $self->throw("element $i is a scalar but not in last dimension");
#   }
#    }
#    if($lastidx < 0) {
#   push(@$tree, @values);
#    } elsif(@values < 2) {
#   $tree->[$lastidx] = shift(@values);
#    } else {
#   $tree->[$lastidx] = [@values];
#    }
#
#}

1;
