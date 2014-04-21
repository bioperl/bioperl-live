#
# BioPerl module for BaseSAXHandler
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Juguang Xiao, juguang@tll.org.sg
#
# Copyright Juguang Xiao 
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::OntologyIO::Handlers::BaseSAXHandler base class for SAX Handlers

=head1 SYNOPSIS

See description.

=head1 DESCRIPTION

This module is an abstract module, serving as the base of any SAX Handler 
implementation. It tries to offer the framework that SAX handlers generally 
need, such as tag_stack, char_store, etc.

In the implementation handler, you can take advantage of this based module by
the following suggestions.

1) In start_element,

 sub start_element {
     my $self=shift;
     my $tag=$_[0]->{Name};
     my %args=%{$_[0]->{Attributes}};
     # Your code here.

     # Before you conclude the method, write these 2 line.
     $self->_visited_count_inc($tag);
     $self->_push_tag($tag);
 }

2) In end_element,

 sub end_element {
     my $self=shift;
     my $tag=shift->{Name};
     # Your code here.

     # Before you conclude the method, write these 2 lines.
     $self->_visited_count_dec($tag);
     $self->_pop_tag;
 }

3) In characters, or any other methods where you may use the tag
stack or count

 sub characters {
     my $self=shift;
     my $text=shift->{Data};

     $self->_chars_hash->{$self->_top_tag} .= $text;

 }
 $count = $self->_visited_count('myTag');
 $tag = $self->_top_tag;


=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists.

Your participation is much appreciated.

  bioperl-l@bioperl.org                  - General discussion
  http://bioperl.org/wiki/Mailing_lists  - About the mailing lists

=head2 Support 

Please direct usage questions or support issues to the mailing list:

I<bioperl-l@bioperl.org>

rather than to the module maintainer directly. Many experienced and 
reponsive experts will be able look at the problem and quickly 
address it. Please include a thorough description of the problem 
with code and data examples if at all possible.

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via the
web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR

Juguang Xiao, juguang@tll.org.sg

=head2 APPENDIX

The rest of the documentation details each of the object methods.
Interal methods are usually preceded with a _

=cut

package Bio::OntologyIO::Handlers::BaseSAXHandler;
use strict;
use base qw(Bio::Root::Root);


sub new {
    my ($class, @args) = @_;
    my $self=$class->SUPER::new(@args);
    $self->_initialize(@args);
    return $self;
}

sub _initialize {
    my $self = shift;
    $self->{_tag_stack} = [];
    $self->{_visited_count} = {};
    $self->{_chars_hash} = {};
    $self->{_current_hash} = {};
}

=head2 _tag_stack

  Title   : _tag_stack
  Usage   : @tags = $self->_tag_stack;
  Function: Get an array of tags that have been accessed but not enclosed.
  Return  : 
  Args    :    

=cut

sub _tag_stack {
    return @{shift->{_tag_stack}};
}

=head2 _push_tag

=cut

sub _push_tag {
    my($self,$tag)=@_;
    push @{$self->{_tag_stack}}, $tag;
}

=head2 _pop_tag

=cut

sub _pop_tag {
    my $self=shift;
    return pop @{$self->{_tag_stack}};
}

=head2 _top_tag

  Title   : _top_tag
  Usage   : $top = $self->_top_tag;
  Function: get the top tag in the tag stack.
  Return  : a tag name
  Args    : [none]   

=cut

sub _top_tag {
    my $self = shift;
    my @stack=@{$self->{_tag_stack}};
    return $stack[-1];
# get the last element in an array while remaining it in. There are few  ways
# 1) $stack[-1]
# 2) $stack[$#stack]
# 3) $stack[@stack-1]
}


=head2 _chars_hash

  Title   : _chars_hash
  Usage   : $hash= $self->_chars_hash;
  Function: return the character cache for the specific tag
  Return  : a hash reference, which is intent for character storage for tags
  Args    : [none]

=cut

sub _chars_hash {
    return shift->{_chars_hash};
}

=head2 _current_hash

=cut

sub _current_hash {
    return  shift->{_current_hash};
}

=head2 _visited_count_inc

  Title   : _vistied_count_inc
  Usage   : $self->vistied_count_inc($tag); # the counter for the tag increase
  Function: the counter for the tag increase
  Return  : the current count after this increment
  Args    : the tag name [scalar]

=cut

sub _visited_count_inc {
    my ($self, $tag) = @_;
    my $visited_count=$self->{_visited_count};
    if(exists $visited_count->{$tag}){
        $visited_count->{$tag}++;
    }else{
        $visited_count->{$tag}=1;
    }
    return $visited_count->{$tag};
}

=head2 _visited_count_dec

  Title   : _visited_count_dec
  Usage   : $self->_visited_count_dec($tag);
  Function: the counter for the tag decreases by one
  Return  : the current count for the specific tag after the decrement
  Args    : the tag name [scalar]

=cut

sub _visited_count_dec {
    my ($self, $tag) = @_;
    my $visited_count=$self->{_visited_count};
    if(exists $visited_count->{$tag}){
        $visited_count->{$tag}--;
    }else{
        $self->throw("'$tag' has not been visited yet. How to decrease it?!");
    }
    return $visited_count->{$tag};
}

=head2 _visited_count

  Title   : _visited_count
  Usage   : $count = $self->_visited_count
  Function: return the counter for the tag
  Return  : the current counter for the specific tag
  Args    : the tag name [scalar]

=cut

sub _visited_count {
    my ($self, $tag) = @_;
    return $self->{_visited_count}->{$tag};
}

1;
