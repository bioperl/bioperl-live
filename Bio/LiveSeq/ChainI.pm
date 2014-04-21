#
# bioperl module for Bio::LiveSeq::ChainI
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Joseph Insana <insana@ebi.ac.uk> <jinsana@gmx.net>
#
# Copyright Joseph Insana
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=head1 NAME

Bio::LiveSeq::ChainI - Double linked chain data structure

=head1 SYNOPSIS

  #documentation needed

=head1 DESCRIPTION

This class generates and manipulates generic double linked list, chain,
that can be used to manage biological sequences.

The advantages over strings or plain arrays is the ease of tracking
changes (mutations) in the elements (sequence). The other side of the
coin is that these structures need consideraly more memory, but that
is cheap and constantly inceasing resource in computers.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists.  Your participation is much appreciated.

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

=head1 AUTHOR - Joseph A.L. Insana

Email:  Insana@ebi.ac.uk, jinsana@gmx.net

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::LiveSeq::ChainI;

use Carp qw(croak);
use strict; # this will be moved before when strict enforced in Chain.pm

use Bio::LiveSeq::Chain; # package where all the subroutines are defined


=head2 new

  Title   : new
  Usage   : $chain = Bio::LiveSeq::ChainI->new(-string => "thequickbrownfoxjumpsoverthelazydog",
					     -offset => 3 );
         OR $chain = Bio::LiveSeq::ChainI->new(-array => \@array,
					     -offset => 3 );
  Function: generates a new Bio::LiveSeq:ChainI
  Returns : a new Chain
  Args    : string
         OR arrayreference
        AND optional offset to create element labels

=cut

sub new {
  my ($thing, %args) = @_;
  my $class = ref($thing) || $thing;
  my $obj;

  if ($args{-string}) {
    $obj = $thing->string2chain($args{-string}, $args{-offset});
  } elsif ($args{-array}) {
    $obj = $thing->array2chain($args{-array}, $args{-offset});
  } else {
    croak "$class not initialized properly";
  }

  $obj = bless $obj, $class;
  return $obj;
}

# added as of 1.9
sub string2chain {
  shift @_; # so that it doesn't pass the object reference
  return Bio::LiveSeq::Chain::string2chain(@_);
}
sub array2chain {
  shift @_; # so that it doesn't pass the object reference
  return Bio::LiveSeq::Chain::array2chain(@_);
}
#
sub chain2string {
  croak "ambiguous method call. Explicit down_ or up_";
}
sub down_chain2string {
  return Bio::LiveSeq::Chain::down_chain2string(@_);
}
sub up_chain2string {
  return Bio::LiveSeq::Chain::up_chain2string(@_);
}
sub chain2string_verbose {
  croak "ambiguous method call. Explicit down_ or up_";
}
sub down_chain2string_verbose {
  return Bio::LiveSeq::Chain::down_chain2string_verbose(@_);
}
sub up_chain2string_verbose {
  return Bio::LiveSeq::Chain::up_chain2string_verbose(@_);
}
sub invert_chain {
  return Bio::LiveSeq::Chain::invert_chain(@_);
}
sub mutate_element {
  croak "Old method name, please update code to: set_value_at_label";
}

# new as of version 2.33 of Chain.pm
sub down_labels {
  return Bio::LiveSeq::Chain::down_labels(@_);
}
sub up_labels {
  return Bio::LiveSeq::Chain::up_labels(@_);
}

sub start {
  return Bio::LiveSeq::Chain::start(@_);
}
sub end {
  return Bio::LiveSeq::Chain::end(@_);
}
sub label_exists {
  return Bio::LiveSeq::Chain::label_exists(@_);
}

sub get_value_at_pos {
  croak "ambiguous method call. Explicit down_ or up_";
}
sub down_get_value_at_pos {
  return Bio::LiveSeq::Chain::down_get_value_at_pos(@_);
}
sub up_get_value_at_pos {
  return Bio::LiveSeq::Chain::up_get_value_at_pos(@_);
}
sub set_value_at_pos {
  croak "ambiguous method call. Explicit down_ or up_";
}
sub down_set_value_at_pos {
  return Bio::LiveSeq::Chain::down_set_value_at_pos(@_);
}
sub up_set_value_at_pos {
  return Bio::LiveSeq::Chain::up_set_value_at_pos(@_);
}
sub get_value_at_label {
  return Bio::LiveSeq::Chain::get_value_at_label(@_);
}
sub set_value_at_label {
  return Bio::LiveSeq::Chain::set_value_at_label(@_);
}
sub get_label_at_pos {
  croak "ambiguous method call. Explicit down_ or up_";
}
sub up_get_label_at_pos {
  return Bio::LiveSeq::Chain::up_get_label_at_pos(@_);
}
sub down_get_label_at_pos {
  return Bio::LiveSeq::Chain::down_get_label_at_pos(@_);
}
sub get_pos_of_label {
  croak "ambiguous method call. Explicit down_ or up_";
}
sub up_get_pos_of_label {
  return Bio::LiveSeq::Chain::up_get_pos_of_label(@_);
}
sub down_get_pos_of_label {
  return Bio::LiveSeq::Chain::down_get_pos_of_label(@_);
}
#

sub preinsert_string {
  return Bio::LiveSeq::Chain::praeinsert_string(@_);
}
sub preinsert_array {
  return Bio::LiveSeq::Chain::praeinsert_array(@_);
}
sub praeinsert_string {
  return Bio::LiveSeq::Chain::praeinsert_string(@_);
}
sub postinsert_string {
  return Bio::LiveSeq::Chain::postinsert_string(@_);
}
sub praeinsert_array {
  return Bio::LiveSeq::Chain::praeinsert_array(@_);
}
sub postinsert_array {
  return Bio::LiveSeq::Chain::postinsert_array(@_);
}
sub down_element{
  return Bio::LiveSeq::Chain::down_element(@_);
}
sub up_element {
  return Bio::LiveSeq::Chain::up_element(@_);
}
sub is_downstream {
  return Bio::LiveSeq::Chain::is_downstream(@_);
}
sub is_upstream {
  return Bio::LiveSeq::Chain::is_upstream(@_);
}
sub check_chain {
  return Bio::LiveSeq::Chain::check_chain(@_);
}
sub chain_length {
  return Bio::LiveSeq::Chain::chain_length(@_);
}
sub splice_chain {
  return Bio::LiveSeq::Chain::splice_chain(@_);
}
sub pos_of_element {
  croak "ambiguous and old method name. use: down_pos_of_label";
}
sub up_pos_of_element {
  croak "old method name. use: down_pos_of_label";
  return Bio::LiveSeq::Chain::up_pos_of_element(@_);
}
sub down_pos_of_element {
  croak "old method name. use: up_pos_of_label";
  return Bio::LiveSeq::Chain::down_pos_of_element(@_);
}
sub subchain_length {
  croak "ambiguous method call. Explicit down_ or up_";
}
sub down_subchain_length {
  return Bio::LiveSeq::Chain::down_subchain_length(@_);
}
sub up_subchain_length {
  return Bio::LiveSeq::Chain::up_subchain_length(@_);
}

# these have to be deleted and changed names to conform to terminology
sub elements {
  return Bio::LiveSeq::Chain::down_elements(@_);
}
sub up_elements {
  return Bio::LiveSeq::Chain::up_elements(@_);
}
sub down_elements {
  return Bio::LiveSeq::Chain::down_elements(@_);
}

1;
