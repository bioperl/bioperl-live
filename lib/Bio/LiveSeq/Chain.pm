#!/usr/bin/perl
#
# bioperl module for Bio::LiveSeq::Chain
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
#

=head1 NAME

Bio::LiveSeq::Chain - DoubleChain DataStructure for Perl

=head1 SYNOPSIS

  #documentation needed

=head1 DESCRIPTION

This is a general purpose module (that's why it's not in object-oriented
form) that introduces a novel datastructure in PERL. It implements
the "double linked chain". The elements of the chain can contain basically
everything. From chars to strings, from object references to arrays or hashes.
It is used in the LiveSequence project to create a dynamical DNA sequence,
easier to manipulate and change. It's use is mainly for sequence variation
analysis but it could be used - for example - in e-cell projects.
The Chain module in itself doesn't have any biological bias, so can be
used for any programming purpose.

Each element of the chain (with the exclusion of the first and the last of the
chain) is connected to other two elements (the PREVious and the NEXT one).
There is no absolute position (like in an array), hence if positions are
important, they need to be computed (methods are provided).
Otherwise it's easy to keep track of the elements with their "LABELs".
There is one LABEL (think of it as a pointer) to each ELEMENT. The labels
won't change after insertions or deletions of the chain. So it's
always possible to retrieve an element even if the chain has been
modified by successive insertions or deletions.
From this the high potential profit for bioinformatics: dealing with
sequences in a way that doesn't have to rely on positions, without
the need of constantly updating them if the sequence changes, even
dramatically.

=head1 AUTHOR - Joseph A.L. Insana

Email:  Insana@ebi.ac.uk, jinsana@gmx.net

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

# DoubleChain Data Structure for PERL
# by Joseph A.L. Insana - Deathson - Filius Mortis - Fal Mortais
# insana@ebi.ac.uk, jinsana@gmx.net

package Bio::LiveSeq::Chain;
# TODO_list:
# **** cleanup code
# **** performance concerns
# *??* create hash2dchain ???? (with hashkeys used for label)
# **????** how about using array of arrays instead than hash of arrays??
#
# further strict complaints:
# in verbose $string assignment around line 721 ???

# TERMINOLOGY update, naming convention:
# "chain" the datastructure
# "element" the individual units that compose a chain
# "label" the unique name of a single element
# "position" the position of an element into the chain according to a
#            particular coordinate system (e.g. counting from the start)
# "value" what is stored in a single element

use Carp qw(croak cluck carp);
use Bio::Root::Version;
use strict; 
use integer; # WARNING: this is to increase performance
             # a little bit of attention has to be given if float need to
             # be stored as elements of the array
             # the use of this "integer" affects all operations but not
             # assignments. So float CAN be assigned as elements of the chain
             # BUT, if you assign $z=-1.8;, $z will be equal to -1 because
             # "-" counts as a unary operation!

=head2 _updown_chain2string

  Title   : chain2string
  Usage   : $string = Bio::LiveSeq::Chain::chain2string("down",$chain,6,9)
  Function: reads the contents of the chain, outputting a string
  Returns : a string
  Examples:
          : down_chain2string($chain) -> all the chain from begin to end
          : down_chain2string($chain,6) -> from 6 to the end
          : down_chain2string($chain,6,4) -> from 6, going on 4 elements
          : down_chain2string($chain,6,"",10) -> from 6 to 10
          : up_chain2string($chain,10,"",6) -> from 10 to 6 upstream
  Defaults: start=first element; if len undef, goes to last
            if last undef, goes to end
            if last defined, it overrides len (undefining it)
  Error code: -1
  Args    : "up"||"down" as first argument to specify the reading direction 
            reference (to the chain)
            [first] [len] [last] optional integer arguments to specify how
            much and from (and to) where to read

=cut

# methods rewritten 2.61
sub up_chain2string {
  _updown_chain2string("up",@_);
}
sub down_chain2string {
  _updown_chain2string("down",@_);
}

sub _updown_chain2string {
  my ($direction,$chain,$first,$len,$last)=@_;
  unless($chain) { cluck "no chain input"; return (-1); }
  my $begin=$chain->{'begin'}; # the label of the BEGIN element
  my $end=$chain->{'end'}; # the label of the END element
  my $flow;

  if ($direction eq "up") {
    $flow=2; # used to determine the direction of chain navigation
    unless ($first) { $first=$end; } # if undef or 0, use $end
  } else { # defaults to "down"
    $flow=1; # used to determine the direction of chain navigation
    unless ($first) { $first=$begin; } # if undef or 0, use $begin
  }

  unless($chain->{$first}) {
    cluck "label for first not defined"; return (-1); }
  if ($last) { # if last is defined, it gets priority and len is not used
    unless($chain->{$last}) {
      cluck "label for last not defined"; return (-1); }
    if ($len) {
      warn "Warning chain2string: argument LAST:$last overriding LEN:$len!";
      undef $len;
    }
  } else {
    if ($direction eq "up") {
      $last=$begin; # if last not defined, go 'till begin (or upto len elements)
    } else {
      $last=$end; # if last not defined, go 'till end (or upto len elements)
    }
  }

  my ($string,@array);
  my $label=$first; my $i=1;
  my $afterlast=$chain->{$last}[$flow]; # if last=end, afterlast is undef
  unless (defined $afterlast) { $afterlast=0; } # keep strict happy

  # proceed for len elements or until last, whichever comes first
  # if $len undef goes till end
  while (($label) && ($label != $afterlast) && ($i <= ($len || $i + 1))) {
    @array=@{$chain->{$label}};
    $string .= $array[0];
    $label = $array[$flow];
    $i++;
  }
  return ($string); # if chain is interrupted $string won't be complete
}

=head2 _updown_labels

 Title   : labels
 Usage   : @labels = Bio::LiveSeq::Chain::_updown_labels("down",$chain,4,16)
 Function: returns all the labels in a chain or those between two
           specified ones (termed "first" and "last")
 Returns : a reference to an array containing the labels
 Args    : "up"||"down" as first argument to specify the reading direction 
           reference (to the chain)
           [first] [last] (integer for the starting and eneding labels)

=cut


# arguments: CHAIN_REF [FIRSTLABEL] [LASTLABEL]
# returns: reference to array containing the labels
sub down_labels {
  my ($chain,$first,$last)=@_;
  _updown_labels("down",$chain,$first,$last);
}
sub up_labels {
  my ($chain,$first,$last)=@_;
  _updown_labels("up",$chain,$first,$last);
}
# arguments: "up"||"down" CHAIN_REF [FIRSTLABEL] [LASTLABEL]
# returns: reference to array containing the labels
sub _updown_labels {
  my ($direction,$chain,$first,$last)=@_;
  unless($chain) { cluck "no chain input"; return (0); }
  my $begin=$chain->{'begin'}; # the label of the BEGIN element
  my $end=$chain->{'end'}; # the label of the END element
  my $flow;
  if ($direction eq "up") { $flow=2;
    unless ($first) { $first=$end; }
    unless ($last) { $last=$begin; }
  } else { $flow=1;
    unless ($last) { $last=$end; }
    unless ($first) { $first=$begin; }
  }
  unless($chain->{$first}) { warn "not existing label $first"; return (0); }
  unless($chain->{$last}) { warn "not existing label $last"; return (0); }

  my $label=$first; my @labels;
  my $afterlast=$chain->{$last}[$flow]; # if last=end, afterlast is undef
  unless (defined $afterlast) { $afterlast=0; } # keep strict happy

  while (($label)&&($label != $afterlast)) {
    push(@labels,$label);
    $label=$chain->{$label}[$flow];
  }
  return (\@labels); # if chain is interrupted @labels won't be complete
}


=head2 start

 Title   : start
 Usage   : $start = Bio::LiveSeq::Chain::start()
 Returns : the label marking the start of the chain
 Errorcode: -1
 Args    : none

=cut

sub start {
  my $chain=$_[0];
  unless($chain) { cluck "no chain input"; return (-1); }
  return ($chain->{'begin'});
}

=head2 end

 Title   : end
 Usage   : $end = Bio::LiveSeq::Chain::end()
 Returns : the label marking the end of the chain
 Errorcode: -1
 Args    : none

=cut

sub end {
  my $chain=$_[0];
  unless($chain) { cluck "no chain input"; return (-1); }
  return ($chain->{'end'});
}

=head2 label_exists

 Title   : label_exists
 Usage   : $check = Bio::LiveSeq::Chain::label_exists($chain,$label)
 Function: It checks if a label is defined, i.e. if an element is there or
           is not there anymore
 Returns : 1 if the label exists, 0 if it is not there, -1 error
 Errorcode: -1
 Args    : reference to the chain, integer

=cut

sub label_exists {
  my ($chain,$label)=@_;
  unless($chain) { cluck "no chain input"; return (-1); }
  if ($label && $chain->{$label}) { return (1); } else { return (0) };
}


=head2 down_get_pos_of_label

 Title   : down_get_pos_of_label
 Usage   : $position = Bio::LiveSeq::Chain::down_get_pos_of_label($chain,$label,$first)
 Function: returns the position of $label counting from $first, i.e. taking
           $first as 1 of coordinate system. If $first is not specified it will
           count from the start of the chain.
 Returns : 
 Errorcode: 0
 Args    : reference to the chain, integer (the label of interest)
           optional: integer (a different label that will be taken as the
           first one, i.e. the one to count from)
 Note:     It counts "downstream". To proceed backward use up_get_pos_of_label

=cut

sub down_get_pos_of_label {
  #down_chain2string($_[0],$_[2],undef,$_[1],"counting");
  my ($chain,$label,$first)=@_;
  _updown_count("down",$chain,$first,$label);
}
sub up_get_pos_of_label {
  #up_chain2string($_[0],$_[2],undef,$_[1],"counting");
  my ($chain,$label,$first)=@_;
  _updown_count("up",$chain,$first,$label);
}

=head2 down_subchain_length

 Title   : down_subchain_length
 Usage   : $length = Bio::LiveSeq::Chain::down_subchain_length($chain,$first,$last)
 Function: returns the length of the chain between the labels "first" and "last", included
 Returns : integer
 Errorcode: 0
 Args    : reference to the chain, integer, integer
 Note:     It counts "downstream". To proceed backward use up_subchain_length

=cut

# arguments: chain_ref [first] [last]
# returns the length of the chain between first and last (included)
sub down_subchain_length {
  #down_chain2string($_[0],$_[1],undef,$_[2],"counting");
  my ($chain,$first,$last)=@_;
  _updown_count("down",$chain,$first,$last);
}
sub up_subchain_length {
  #up_chain2string($_[0],$_[1],undef,$_[2],"counting");
  my ($chain,$first,$last)=@_;
  _updown_count("up",$chain,$first,$last);
}

# arguments: DIRECTION CHAIN_REF FIRSTLABEL LASTLABEL
# errorcode 0
sub _updown_count {
  my ($direction,$chain,$first,$last)=@_;
  unless($chain) { cluck "no chain input"; return (0); }
  my $begin=$chain->{'begin'}; # the label of the BEGIN element
  my $end=$chain->{'end'}; # the label of the END element
  my $flow;
  if ($direction eq "up") { $flow=2;
    unless ($first) { $first=$end; }
    unless ($last) { $last=$begin; }
  } else { $flow=1;
    unless ($last) { $last=$end; }
    unless ($first) { $first=$begin; }
  }
  unless($chain->{$first}) { warn "not existing label $first"; return (0); }
  unless($chain->{$last}) { warn "not existing label $last"; return (0); }

  my $label=$first; my $count;
  my $afterlast=$chain->{$last}[$flow]; # if last=end, afterlast is undef
  unless (defined $afterlast) { $afterlast=0; } # keep strict happy

  while (($label)&&($label != $afterlast)) {
    $count++;
    $label=$chain->{$label}[$flow];
  }
  return ($count); # if chain is interrupted, $i will be up to the breaking point
}

=head2 invert_chain

 Title   : invert_chain
 Usage   : $errorcode=Bio::LiveSeq::Chain::invert_chain($chain)
 Function: completely inverts the order of the chain elements; begin is swapped with end and all links updated (PREV&NEXT fields swapped)
 Returns : 1 if all OK, 0 if errors
 Errorcode: 0
 Args    : reference to the chain

=cut

sub invert_chain {
  my $chain=$_[0];
  unless($chain) { cluck "no chain input"; return (0); }
  my $begin=$chain->{'begin'}; # the name of the first element
  my $end=$chain->{'end'}; # the name of the last element
  my ($label,@array);
  $label=$begin; # starts from the beginning
  while ($label) { # proceed with linked elements, swapping PREV and NEXT
    @array=@{$chain->{$label}};
    ($chain->{$label}[1],$chain->{$label}[2])=($array[2],$array[1]); # swap
    $label = $array[1]; # go to the next one
  }
  # now swap begin and end fields
  ($chain->{'begin'},$chain->{'end'})=($end,$begin);
  return (1); # that's it
}

# warning that method has changed name
#sub mutate_element {
  #croak "Warning: old method name. Please update code to 'set_value_at_label'\n";
  # &set_value_at_label;
#}

=head2 down_get_value_at_pos

 Title   : down_get_value_at_pos
 Usage   : $value = Bio::LiveSeq::Chain::down_get_value_at_pos($chain,$position,$first)
 Function: used to access the value of the chain at a particular position instead than directly with a label pointer. It will count the position from the start of the chain or from the label $first, if $first is specified
 Returns : whatever is stored in the element of the chain
 Errorcode: 0
 Args    : reference to the chain, integer, [integer]
 Note:     It works "downstream". To proceed backward use up_get_value_at_pos

=cut

#sub get_value_at_pos {
  #croak "Please use instead: down_get_value_at_pos";
  ##&down_get_value_at_pos;
#}
sub down_get_value_at_pos {
  my ($chain,$position,$first)=@_;
  my $label=down_get_label_at_pos($chain,$position,$first);
  # check place of change
  if (($label eq -1)||($label eq 0)) { # complain if label doesn't exist
    warn "not existing element $label"; return (0); }
  return _get_value($chain,$label);
}
sub up_get_value_at_pos {
  my ($chain,$position,$first)=@_;
  my $label=up_get_label_at_pos($chain,$position,$first);
  # check place of change
  if (($label eq -1)||($label eq 0)) { # complain if label doesn't exist
    warn "not existing element $label"; return (0); }
  return _get_value($chain,$label);
}

=head2 down_set_value_at_pos

 Title   : down_set_value_at_pos
 Usage   : $errorcode = Bio::LiveSeq::Chain::down_set_value_at_pos($chain,$newvalue,$position,$first)
 Function: used to store a new value inside an element of the chain at a particular position instead than directly with a label pointer. It will count the position from the start of the chain or from the label $first, if $first is specified
 Returns : 1
 Errorcode: 0
 Args    : reference to the chain, newvalue, integer, [integer]
           (newvalue can be: integer, string, object reference, hash ref)
 Note:     It works "downstream". To proceed backward use up_set_value_at_pos
 Note2:    If the $newvalue is undef, it will delete the contents of the
           element but it won't remove the element from the chain.

=cut

#sub set_value_at_pos {
  #croak "Please use instead: down_set_value_at_pos";
  ##&down_set_value_at_pos;
#}
sub down_set_value_at_pos {
  my ($chain,$value,$position,$first)=@_;
  my $label=down_get_label_at_pos($chain,$position,$first);
  # check place of change
  if (($label eq -1)||($label eq 0)) { # complain if label doesn't exist
    warn "not existing element $label"; return (0); }
  _set_value($chain,$label,$value);
  return (1);
}
sub up_set_value_at_pos {
  my ($chain,$value,$position,$first)=@_;
  my $label=up_get_label_at_pos($chain,$position,$first);
  # check place of change
  if (($label eq -1)||($label eq 0)) { # complain if label doesn't exist
    warn "not existing element $label"; return (0); }
  _set_value($chain,$label,$value);
  return (1);
}


=head2 down_set_value_at_label

 Title   : down_set_value_at_label
 Usage   : $errorcode = Bio::LiveSeq::Chain::down_set_value_at_label($chain,$newvalue,$label)
 Function: used to store a new value inside an element of the chain defined by its label.
 Returns : 1
 Errorcode: 0
 Args    : reference to the chain, newvalue, integer
           (newvalue can be: integer, string, object reference, hash ref)
 Note:     It works "downstream". To proceed backward use up_set_value_at_label
 Note2:    If the $newvalue is undef, it will delete the contents of the
           element but it won't remove the element from the chain.

=cut

sub set_value_at_label {
  my ($chain,$value,$label)=@_;
  unless($chain) { cluck "no chain input"; return (0); }

  # check place of change
  unless($chain->{$label}) { # complain if label doesn't exist
    warn "not existing element $label"; return (0); }
  _set_value($chain,$label,$value);
  return (1);
}

=head2 down_get_value_at_label

 Title   : down_get_value_at_label
 Usage   : $value = Bio::LiveSeq::Chain::down_get_value_at_label($chain,$label)
 Function: used to access the value of the chain from one element defined by its label.
 Returns : whatever is stored in the element of the chain
 Errorcode: 0
 Args    : reference to the chain, integer
 Note:     It works "downstream". To proceed backward use up_get_value_at_label

=cut

sub get_value_at_label {
  my $chain=$_[0];
  unless($chain) { cluck "no chain input"; return (0); }
  my $label = $_[1]; # the name of the element

  # check place of change
  unless($chain->{$label}) { # complain if label doesn't exist
    warn "not existing label $label"; return (0); }
  return _get_value($chain,$label);
}

# arguments: CHAIN_REF LABEL VALUE
sub _set_value {
  my ($chain,$label,$value)=@_;
  $chain->{$label}[0]=$value;
}
# arguments: CHAIN_REF LABEL
sub _get_value {
  my ($chain,$label)=@_;
  return $chain->{$label}[0];
}

=head2 down_get_label_at_pos

 Title   : down_get_label_at_pos
 Usage   : $label = Bio::LiveSeq::Chain::down_get_label_at_pos($chain,$position,$first)
 Function: used to retrieve the label of an an element of the chain at a particular position. It will count the position from the start of the chain or from the label $first, if $first is specified
 Returns : integer
 Errorcode: 0
 Args    : reference to the chain, integer, [integer]
 Note:     It works "downstream". To proceed backward use up_get_label_at_pos

=cut

# arguments: CHAIN_REF POSITION [FIRST]
# returns: LABEL of element found counting from FIRST
sub down_get_label_at_pos {
  _updown_get_label_at_pos("down",@_);
}
sub up_get_label_at_pos {
  _updown_get_label_at_pos("up",@_);
}

# arguments: [DIRECTION] CHAIN_REF POSITION [FIRST]
# Default DIRECTION="down"
# if FIRST is undefined, FIRST=START (if DIRECTION=down) or FIRST=END (up)

sub _updown_get_label_at_pos {
  my ($direction,$chain,$position,$first)=@_;
  unless($chain) { cluck "no chain input"; return (0); }
  my $begin=$chain->{'begin'}; # the label of the BEGIN element
  my $end=$chain->{'end'}; # the label of the END element
  my $flow;
  if ($direction eq "up") { $flow=2; unless ($first) { $first=$end; }
  } else { $flow=1; unless ($first) { $first=$begin; } }
  unless($chain->{$first}) { warn "not existing label $first"; return (0); }

  my $label=$first;
  my $i=1;
  while ($i < $position) {
    $label=$chain->{$label}[$flow];
    $i++;
    unless ($label) { return (0); } # chain ended before position reached
  }
  return ($label);
}

# for english_concerned, latin_unconcerned people
sub preinsert_string { &praeinsert_string }
sub preinsert_array { &praeinsert_array }

# praeinsert_string CHAIN_REF STRING [POSITION]
# the chars of STRING are passed to praeinsert_array
# the chars are inserted in CHAIN, before POSITION
# if POSITION is undef, default is to prepend the string to the beginning
# i.e. POSITION is START of CHAIN
sub praeinsert_string {
  my @string=split(//,$_[1]);
  praeinsert_array($_[0],\@string,$_[2]);
}

# postinsert_string CHAIN_REF STRING [POSITION]
# the chars of STRING are passed to postinsert_array
# the chars are inserted in CHAIN, after POSITION
# if POSITION is undef, default is to append the string to the end
# i.e. POSITION is END of CHAIN
sub postinsert_string {
  my @string=split(//,$_[1]);
  postinsert_array($_[0],\@string,$_[2]);
}

# praeinsert_array CHAIN_REF ARRAY_REF [POSITION]
# the elements of ARRAY are inserted in CHAIN, before POSITION
# if POSITION is undef, default is to prepend the elements to the beginning
# i.e. POSITION is START of CHAIN
sub praeinsert_array {
  _praepostinsert_array($_[0],"prae",$_[1],$_[2]);
}

# postinsert_array CHAIN_REF ARRAY_REF [POSITION]
# the elements of ARRAY are inserted in CHAIN, after POSITION
# if POSITION is undef, default is to append the elements to the end
# i.e. POSITION is END of CHAIN
sub postinsert_array {
  _praepostinsert_array($_[0],"post",$_[1],$_[2]);
}


=head2 _praepostinsert_array

 Title   : _praepostinsert_array
 Usage   : ($insbegin,$insend) = Bio::LiveSeq::Chain::_praepostinsert_array($chainref,"post",$arrayref,$position)
 Function: the elements of the array specified by $arrayref are inserted (creating a new subchain) in the chain specified by $chainref, before or after (depending on the "prae"||"post" keyword passed as second argument) the specified position.
 Returns : two labels: the first and the last of the inserted subchain
 Defaults: if no position is specified, the new chain will be inserted after
 (post) the first element of the chain
 Errorcode: 0
 Args    : chainref, "prae"||"post", arrayref, integer (position)

=cut

# returns: 0 if errors, otherwise returns references of begin and end of
# the insertion
sub _praepostinsert_array {
  my $chain=$_[0];
  unless($chain) { cluck "no chain input"; return (0); }
  my $praepost=$_[1] || "post"; # defaults to post
  my ($prae,$post);
  my $position=$_[3];
  my $begin=$chain->{'begin'}; # the name of the first element of the chain
  my $end=$chain->{'end'}; # the name of the the last element of the chain
  # check if prae or post insertion and prepare accordingly
  if ($praepost eq "prae") {
    $prae=1;
    unless (($position eq 0)||($position)) { $position=$begin; } # if undef, use $begin
  } else {
    $post=1;
    unless (($position eq 0)||($position)) { $position=$end; } # if undef, use $end
  }
  # check place of insertion
  unless($chain->{$position}) { # complain if position doesn't exist
    warn ("Warning _praepostinsert_array: not existing element $position");
    return (0);
  }

  # check if there are elements to insert
  my $elements=$_[2]; # reference to the array containing the new elements
  my $elements_count=scalar(@{$elements});
  unless ($elements_count) {
    warn ("Warning _praepostinsert_array: no elements input"); return (0); }

  # create new chainelements with offset=firstfree(chain)
  my ($insertbegin,$insertend)=_create_chain_elements($chain,$elements);

  # DEBUGGING
  #print "Executing ${praepost}insertion of $elements_count elements ('@{$elements}') at position: $position\n";

  # attach the new chain to the old chain
  # 4 cases: prae@begin, prae@middle, post@middle, post@end
  # NOTE: in case of double joinings always join wisely so not to
  # delete the PREV/NEXT attribute before it is needed
  my $noerror=1;
  if ($prae) {
    if ($position==$begin) { # 1st case: prae@begin
      $noerror=_join_chain_elements($chain,$insertend,$begin);
      $chain->{'begin'}=$insertbegin;
    } else { # 2nd case: prae@middle
      $noerror=_join_chain_elements($chain,up_element($chain,$position),$insertbegin);
      $noerror=_join_chain_elements($chain,$insertend,$position);
    }
  } elsif ($post) {
    if ($position==$end) { # 4th case: post@end
      $noerror=_join_chain_elements($chain,$end,$insertbegin);
      $chain->{'end'}=$insertend;
    } else { # 3rd case: post@middle # note the order of joins (important)
      $noerror=_join_chain_elements($chain,$insertend,down_element($chain,$position));
      $noerror=_join_chain_elements($chain,$position,$insertbegin);
    }
  } else { # this should never happen
    die "_praepostinsert_array: Something went very wrong";
  }

  # check for errors and return begin,end of insertion
  if ($noerror) {
    return ($insertbegin,$insertend);
  } else { # something went wrong with the joinings
    warn "Warning _praepostinsert_array: Joining of insertion failed";
    return (0);
  }
}

# create new chain elements with offset=firstfree
# arguments: CHAIN_REF ARRAY_REF
# returns: pointers to BEGIN and END of new chained elements created
# returns 0 if error(s) encountered
sub _create_chain_elements {
  my $chain=$_[0];
  unless($chain) {
    warn ("Warning _create_chain_elements: no chain input"); return (0); }
  my $arrayref=$_[1];
  my $array_count=scalar(@{$arrayref});
  unless ($array_count) {
    warn ("Warning _create_chain_elements: no elements input"); return (0); }
  my $begin=$chain->{'firstfree'};
  my $i=$begin-1;
  my $element;
  foreach $element (@{$arrayref}) {
    $i++;
    $chain->{$i}=[$element,$i+1,$i-1];
  }
  my $end=$i;
  $chain->{'firstfree'}=$i+1; # what a new added element should be called
  $chain->{'size'} += $end-$begin+1; # increase size of chain
  # leave sticky edges (to be joined by whoever called this subroutine)
  $chain->{$begin}[2]=undef;
  $chain->{$end}[1]=undef;
  return ($begin,$end); # return pointers to first and last of the newelements
}

# argument: CHAIN_REF ELEMENT
# returns: name of DOWN/NEXT element (the downstream one)
# returns -1 if error encountered (e.g. chain or elements undefined)
# returns 0 if there's no DOWN element
sub down_element {
  _updown_element("down",@_);
}
# argument: CHAIN_REF ELEMENT
# returns: name of UP/PREV element (the upstream one)
# returns -1 if error encountered (e.g. chain or elements undefined)
# returns 0 if there's no UP element
sub up_element {
  _updown_element("up",@_);
}

# used by both is_up_element and down_element
sub _updown_element {
  my $direction=$_[0] || "down"; # defaults to downstream
  my $flow;
  if ($direction eq "up") {
    $flow=2; # used to determine the direction of chain navigation
  } else {
    $flow=1; # used to determine the direction of chain navigation
  }
  my $chain=$_[1];
  unless($chain) {
    warn ("Warning ${direction}_element: no chain input"); return (-1); }
  my $me = $_[2]; # the name of the element
  my $it = $chain->{$me}[$flow]; # the prev||next one, upstream||downstream
  if ($it) {
    return ($it); # return the name of prev||next element
  } else {
    return (0); # there is no prev||next element ($it is undef)
  }
}

# used by both is_downstream and is_upstream
sub _is_updownstream {
  my $direction=$_[0] || "down"; # defaults to downstream
  my $flow;
  if ($direction eq "up") {
    $flow=2; # used to determine the direction of chain navigation
  } else {
    $flow=1; # used to determine the direction of chain navigation
  }
  my $chain=$_[1];
  unless($chain) {
    warn ("Warning is_${direction}stream: no chain input"); return (-1); }
  my $first=$_[2]; # the name of the first element
  my $second=$_[3]; # the name of the first element
  if ($first==$second) {
    warn ("Warning is_${direction}stream: first==second!!"); return (0); }
  unless($chain->{$first}) {
    warn ("Warning is_${direction}stream: first element not defined"); return (-1); }
  unless($chain->{$second}) {
    warn ("Warning is_${direction}stream: second element not defined"); return (-1); }
  my ($label,@array);
  $label=$first;
  my $found=0;
  while (($label)&&(!($found))) { # searches till the end or till found
    if ($label==$second) {
      $found=1;
    }
    @array=@{$chain->{$label}};
    $label = $array[$flow]; # go to the prev||next one, upstream||downstream
  }
  return $found;
}

=head2 is_downstream

  Title   : is_downstream
  Usage   : Bio::LiveSeq::Chain::is_downstream($chainref,$firstlabel,$secondlabel)
  Function: checks if SECONDlabel follows FIRSTlabel
            It runs downstream the elements of the chain from FIRST searching
            for SECOND.
  Returns : 1 if SECOND is found /after/ FIRST; 0 otherwise (i.e. if it
            reaches the end of the chain without having found it)
  Errorcode -1
  Args    : two labels (integer)

=cut

sub is_downstream {
  _is_updownstream("down",@_);
}

=head2 is_upstream

  Title   : is_upstream
  Usage   : Bio::LiveSeq::Chain::is_upstream($chainref,$firstlabel,$secondlabel)
  Function: checks if SECONDlabel follows FIRSTlabel
            It runs upstream the elements of the chain from FIRST searching
            for SECOND.
  Returns : 1 if SECOND is found /after/ FIRST; 0 otherwise (i.e. if it
            reaches the end of the chain without having found it)
  Errorcode -1
  Args    : two labels (integer)

=cut

sub is_upstream {
  _is_updownstream("up",@_);
}

=head2 check_chain

 Title   : check_chain
 Usage   : @errorcodes = Bio::LiveSeq::Chain::check_chain()
 Function: a wraparound to a series of check for consistency of the chain
           It will check for boundaries, size, backlinking and forwardlinking
 Returns : array of 4 warn codes, each can be 1 (all ok) or 0 (something wrong)
 Errorcode: 0
 Args    : none
 Note    : this is slow and through. It is not really needed. It is mostly
           a code-developer tool.

=cut

sub check_chain {
  my $chain=$_[0];
  unless($chain) {
    warn ("Warning check_chain: no chain input"); return (-1); }
  my ($warnbound,$warnsize,$warnbacklink,$warnforlink);
  $warnbound=&_boundcheck; # passes on the arguments of the subroutine
  $warnsize=&_sizecheck;
  $warnbacklink=&_downlinkcheck;
  $warnforlink=&_uplinkcheck;
  return ($warnbound,$warnsize,$warnbacklink,$warnforlink);
}

# consistency check for forwardlinks walking upstream
# argument: a chain reference
# returns: 1 all OK 0 problems
sub _uplinkcheck {
  _updownlinkcheck("up",@_);
}

# consistency check for backlinks walking downstream
# argument: a chain reference
# returns: 1 all OK 0 problems
sub _downlinkcheck {
  _updownlinkcheck("down",@_);
}

# consistency check for links, common to _uplinkcheck and _downlinkcheck
# argument: "up"||"down", check_ref
# returns: 1 all OK 0 problems
sub _updownlinkcheck {
  my $direction=$_[0] || "down"; # defaults to downstream
  my ($flow,$wolf);
  my $chain=$_[1];
  unless($chain) {
    warn ("Warning _${direction}linkcheck: no chain input"); return (0); }
  my $begin=$chain->{'begin'}; # the name of the first element
  my $end=$chain->{'end'}; # the name of the last element
  my ($label,@array,$me,$it,$itpoints);
  if ($direction eq "up") {
    $flow=2; # used to determine the direction of chain navigation
    $wolf=1;
    $label=$end; # start from end
  } else {
    $flow=1; # used to determine the direction of chain navigation
    $wolf=2;
    $label=$begin; # start from beginning
  }
  my $warncode=1;

  while ($label) { # proceed with linked elements, checking neighbours
    $me=$label;
    @array=@{$chain->{$label}};
    $label = $array[$flow]; # go to the next one
    $it=$label;
    if ($it) { # no sense in checking if next one not defined (END element)
      @array=@{$chain->{$label}};
      $itpoints=$array[$wolf];
      unless ($me==$itpoints) {
	warn "Warning: ${direction}LinkCheck: LINK wrong in $it, that doesn't point back to me ($me). It points to $itpoints\n";
	$warncode=0;
      }
    }
  }
  return $warncode;
}

# consistency check for size of chain
# argument: a chain reference
# returns: 1 all OK 0 wrong size
sub _sizecheck {
  my $chain=$_[0];
  unless($chain) {
    warn ("Warning _sizecheck: no chain input"); return (0); }
  my $begin=$chain->{'begin'}; # the name of the first element
  my $warncode=1;
  my ($label,@array);
  my $size=$chain->{'size'};
  my $count=0;
  $label=$begin;
  while ($label) { # proceed with linked elements, counting
    @array=@{$chain->{$label}};
    $label = $array[1]; # go to the next one
    $count++;
  }
  if ($size != $count) {
    warn "Size check reports error: assumed size: $size, real size: $count ";
    $warncode=0;
  }
  return $warncode;
}


# consistency check for begin and end (boundaries)
# argument: a chain reference
# returns: 1 all OK 0 problems
sub _boundcheck {
  my $chain=$_[0];
  unless($chain) {
    warn ("Warning _boundcheck: no chain input"); return (0); }
  my $begin=$chain->{'begin'}; # the name of the first element
  my $end=$chain->{'end'}; # the name of the (supposedly) last element
  my $warncode=1;

  # check SYNC of beginning
  if (($begin)&&($chain->{$begin})) { # if the BEGIN points to existing element
    if ($chain->{$begin}[2]) { # if BEGIN element has PREV not undef
      warn "Warning: BEGIN element has PREV field defined \n";
      warn "\tWDEBUG begin: $begin\t";
      warn "\tWDEBUG begin's PREV: $chain->{$begin}[2] \n";
      $warncode=0;
    }
  } else {
    warn "Warning: BEGIN key of chain does not point to existing element!\n";
    warn "\tWDEBUG begin: $begin\n";
    $warncode=0;
  }
  # check SYNC of end
  if (($end)&&($chain->{$end})) { # if the END points to an existing element
    if ($chain->{$end}[1]) { # if END element has NEXT not undef
      warn "Warning: END element has NEXT field defined \n";
      warn "\tWDEBUG end: $end\t";
      warn "\tWDEBUG end's NEXT: $chain->{$end}[1] \n";
      $warncode=0;
    }
  } else {
    warn "Warning: END key of chain does not point to existing element!\n";
    warn "\tWDEBUG end: $end\n";
    $warncode=0;
  }
  return $warncode;
}

# arguments: chain_ref
# returns: the size of the chain (the number of elements)
# return code -1: unexistant chain, errors...
sub chain_length {
  my $chain=$_[0];
  unless($chain) {
    warn ("Warning chain_length: no chain input"); return (-1); }
  my $size=$chain->{'size'};
  if ($size) {
    return ($size);
  } else {
    return (-1);
  }
}

# arguments: chain ref, first element name, second element name
# returns: 1 or 0 (1 ok, 0 errors)
sub _join_chain_elements {
  my $chain=$_[0];
  unless($chain) {
    warn ("Warning _join_chain_elements: no chain input"); return (0); }
  my $leftelem=$_[1];
  my $rightelem=$_[2];
  unless(($leftelem)&&($rightelem)) {
    warn ("Warning _join_chain_elements: element arguments??"); return (0); }
  if (($chain->{$leftelem})&&($chain->{$rightelem})) { # if the elements exist
    $chain->{$leftelem}[1]=$rightelem;
    $chain->{$rightelem}[2]=$leftelem;
    return 1;
  } else {
    warn ("Warning _join_chain_elements: elements not defined");
    return 0;
  }
}

=head2 splice_chain

 Title   : splice_chain
 Usage   : @errorcodes = Bio::LiveSeq::Chain::splice_chain($chainref,$first,$length,$last)
 Function: removes the elements designated by FIRST and LENGTH from a chain.
           The chain shrinks accordingly. If LENGTH is omitted, removes
           everything from FIRST onward.
           If END is specified, LENGTH is ignored and instead the removal
           occurs from FIRST to LAST.
 Returns : the elements removed as a string
 Errorcode: -1
 Args    : chainref, integer, integer, integer

=cut

sub splice_chain {
  my $chain=$_[0];
  unless($chain) {
    warn ("Warning splice_chain: no chain input"); return (-1); }
  my $begin=$chain->{'begin'}; # the name of the first element
  my $end=$chain->{'end'}; # the name of the (supposedly) last element
  my $first=$_[1];
  unless (($first eq 0)||($first)) { $first=$begin; } # if undef, use $begin
  my $len=$_[2];
  my $last=$_[3];
  my (@array, $string);
  my ($beforecut,$aftercut);

  unless($chain->{$first}) {
    warn ("Warning splice_chain: first element not defined"); return (-1); }
  if ($last) { # if last is defined, it gets priority and len is not used
    unless($chain->{$last}) {
      warn ("Warning splice_chain: last element not defined"); return (-1); }
    if ($len) {
      warn ("Warning splice_chain: argument LAST:$last overriding LEN:$len!");
      undef $len;
    }
  } else {
    $last=$end; # if last not defined, go 'till end (or to len, whichever 1st)
  }

  $beforecut=$chain->{$first}[2]; # what's the element before 1st deleted?
  # if it is undef then it means we are splicing since the beginning

  my $i=1;
  my $label=$first;
  my $afterlast=$chain->{$last}[1]; # if $last=$end $afterlast should be undef
  unless (defined $afterlast) { $afterlast=0; } # keep strict happy

  # proceed for len elements or until the end, whichever comes first
  # if len undef goes till last
  while (($label)&&($label != $afterlast) && ($i <= ($len || $i + 1))) {
    @array=@{$chain->{$label}};
    $string .= $array[0];
    $aftercut = $array[1]; # what's the element next last deleted?
			   # also used as savevar to change label posdeletion
    delete $chain->{$label}; # this can be deleted now
    $label=$aftercut; # label is updated using the savevar
    $i++;
  }
  
  # Now fix the chain (sticky edges, fields)
  # 4 cases: cut in the middle, cut from beginning, cut till end, cut all
    #print "\n\tstickyDEBUG beforecut: $beforecut "; # DEBUG
    #print "\taftercut: $aftercut \n"; # DEBUG
  if ($beforecut) {
    if ($aftercut) { # 1st case, middle cut
      _join_chain_elements($chain,$beforecut,$aftercut);
    } else { # 3rd case, end cut
      $chain->{'end'}=$beforecut; # update the END field
      $chain->{$beforecut}[1]=undef; # since we cut till the end
    }
  } else {
    if ($aftercut) { # 2nd case, begin cut
      $chain->{'begin'}=$aftercut; # update the BEGIN field
      $chain->{$aftercut}[2]=undef; # since we cut from beginning
    } else { # 4th case, all has been cut
      $chain->{'begin'}=undef;
      $chain->{'end'}=undef;
    }
  }
  $chain->{'size'}=($chain->{'size'}) - $i + 1; # update the SIZE field

  return $string;
}


# arguments: CHAIN_REF POSITION [FIRST]
# returns: element counting POSITION from FIRST or from START if FIRST undef
# i.e. returns the element at POSITION counting from FIRST
#sub element_at_pos {
  #croak "Warning: old method name. Please update code to 'down_get_label_at_position'\n";
  ##&down_element_at_pos;
#}
#sub up_element_at_pos {
  ## old wraparound
  ##my @array=up_chain2string($_[0],$_[2],$_[1],undef,"elements");
  ##return $array[-1];
  #croak "old method name. Update code to: up_get_label_at_position";
  ##&up_get_label_at_pos;
#}
#sub down_element_at_pos {
  ## old wraparound
  ##my @array=down_chain2string($_[0],$_[2],$_[1],undef,"elements");
  ##return $array[-1];
  #croak "old method name. Update code to: down_get_label_at_position";
  ##&down_get_label_at_pos;
#}

# arguments: CHAIN_REF ELEMENT [FIRST]
# returns: the position of ELEMENT counting from FIRST or from START
#i         if FIRST is undef
# i.e. returns the Number of elements between FIRST and ELEMENT
# i.e. returns the position of element taking FIRST as 1 of coordinate system
#sub pos_of_element {
  #croak ("Warning: old and ambiguous method name. Please update code to 'down_get_pos_of_label'\n");
  ##&down_pos_of_element;
#}
#sub up_pos_of_element {
  #croak ("Warning: old method name. Please update code to 'up_get_pos_of_label'\n");
  ##up_chain2string($_[0],$_[2],undef,$_[1],"counting");
#}
#sub down_pos_of_element {
  #croak ("Warning: old method name. Please update code to 'down_get_pos_of_label'\n");
  ##down_chain2string($_[0],$_[2],undef,$_[1],"counting");
#}

# wraparounds to calculate length of subchain from first to last
# arguments: chain_ref [first] [last]
#sub subchain_length {
  #croak "Warning: old method name. Please update code to 'down_subchain_length'\n";
  ##&down_subchain_length;
#}

# wraparounds to have elements output
# same arguments as chain2string
# returns label|name of every element
#sub elements {
  #croak ("Warning: method no more supported. Please update code to 'down_labels' (NB: now it returns ref to array and doesn't allow length argument!)\n");
  ##&down_elements;
#}
#sub up_elements {
  #croak ("Warning: method no more supported. Please update code to 'up_labels' (NB: now it returns ref to array and doesn't allow length argument!)\n");
  ##up_chain2string($_[0],$_[1],$_[2],$_[3],"elements");
#}
#sub down_elements {
  #croak ("Warning: method no more supported. Please update code to 'down_labels' (NB: now it returns ref to array and doesn't allow length argument!)\n");
  ##down_chain2string($_[0],$_[1],$_[2],$_[3],"elements");
#}

# wraparounds to have verbose output
# same arguments as chain2string
# returns the chain in a very verbose way
sub chain2string_verbose {
  carp "Warning: method no more supported.\n";
  &old_down_chain2string_verbose;
}
sub up_chain2string_verbose {
  carp "Warning: method no more supported.\n";
  old_up_chain2string($_[0],$_[1],$_[2],$_[3],"verbose");
}
sub down_chain2string_verbose {
  carp "Warning: method no more supported.\n";
  old_down_chain2string($_[0],$_[1],$_[2],$_[3],"verbose");
}

#sub chain2string {
  #croak ("Warning: old method name. Please update code to 'down_chain2string'\n");
  ##&down_chain2string;
#}
sub old_up_chain2string {
  old_updown_chain2string("up",@_);
}
sub old_down_chain2string {
  old_updown_chain2string("down",@_);
}

# common to up_chain2string and down_chain2string
# arguments: "up"||"down" chain_ref [first] [len] [last] [option]
# [option] can be any of "verbose", "counting", "elements"
# error: return -1
# defaults: start = first element; if len undef, goes to last
#           if last undef, goes to end
#           if last def it overrides len (that gets undef)
# returns: a string
# example usage: down_chain2string($chain) -> all the chain from begin to end
# example usage: down_chain2string($chain,6) -> from 6 to the end
# example usage: down_chain2string($chain,6,4) -> from 6, going on 4 elements
# example usage: down_chain2string($chain,6,"",10) -> from 6 to 10
# example usage: up_chain2string($chain,10,"",6) -> from 10 to 6 upstream
sub old_updown_chain2string {
  my ($direction,$chain,$first,$len,$last,$option)=@_;
  unless($chain) {
    warn ("Warning chain2string: no chain input"); return (-1); }
  my $begin=$chain->{'begin'}; # the name of the BEGIN element
  my $end=$chain->{'end'}; # the name of the END element
  my $flow;
  if ($direction eq "up") {
    $flow=2; # used to determine the direction of chain navigation
    unless ($first) { $first=$end; } # if undef or 0, use $end
  } else { # defaults to "down"
    $flow=1; # used to determine the direction of chain navigation
    unless ($first) { $first=$begin; } # if undef or 0, use $begin
  }

  unless($chain->{$first}) {
    warn ("Warning chain2string: first element not defined"); return (-1); }
  if ($last) { # if last is defined, it gets priority and len is not used
    unless($chain->{$last}) {
      warn ("Warning chain2string: last element not defined"); return (-1); }
    if ($len) {
      warn ("Warning chain2string: argument LAST:$last overriding LEN:$len!");
      undef $len;
    }
  } else {
    if ($direction eq "up") {
      $last=$begin; # if last not defined, go 'till begin (or upto len elements)
    } else {
      $last=$end; # if last not defined, go 'till end (or upto len elements)
    }
  }
  my (@array, $string, $count);
  # call for verbosity (by way of chain2string_verbose);
  my $verbose=0; my $elements=0; my @elements; my $counting=0;
  if ($option) { # keep strict happy
    if ($option eq "verbose") { $verbose=1; }
    if ($option eq "elements") { $elements=1; }
    if ($option eq "counting") { $counting=1; }
  }

  if ($verbose) {
    print "BEGIN=$begin"; print " END=$end"; print " SIZE=$chain->{'size'}";
    print " FIRSTFREE=$chain->{'firstfree'} \n";
  }

  my $i=1;
  my $label=$first;
  my $afterlast=$chain->{$last}[$flow]; # if $last=$end $afterlast should be undef
  unless (defined $afterlast) { $afterlast=0; } # keep strict happy

  # proceed for len elements or until last, whichever comes first
  # if $len undef goes till end
  while (($label)&&($label != $afterlast) && ($i <= ($len || $i + 1))) {
    @array=@{$chain->{$label}};
    if ($verbose) {
      $string .= "$array[2]_${label}_$array[1]=$array[0] ";
      $count++;
    } elsif ($elements) {
      push (@elements,$label); # returning element names/references/identifiers
    } elsif ($counting) {
      $count++;
    } else {
      $string .= $array[0];   # returning element content
    }
    $label = $array[$flow]; # go to next||prev i.e. downstream||upstream
    $i++;
  }
#DEBUG#print "len: $len, first: $first, last: $last, afterlast=$afterlast \n";
  if ($verbose) { print "TOTALprinted: $count\n"; } 
  if ($counting) {
    return $count;
  } elsif ($elements) {
    return @elements;
  } else {
    return $string;
  }
}

# sub string2schain
# --------> deleted, no more supported <--------
# creation of a single linked list/chain from a string
# basically could be recreated by taking the *2chain methods and
# omitting to set the 3rd field (label 2) containing the back links


# creation of a double linked list/chain from a string
# returns reference to a hash containing the chain
# arguments: STRING [OFFSET]
# defaults: OFFSET defaults to 1 if undef
# the chain will contain as elements the single characters in the string
sub string2chain {
  my @string=split(//,$_[0]);
  array2chain(\@string,$_[1]);
}

=head2 array2chain

  Title   : array2chain
  Usage   : $chainref = Bio::LiveSeq::Chain::array2chain($arrayref,$offset)
  Function: creation of a double linked chain from an array
  Returns : reference to a hash containing the chain
  Defaults: OFFSET defaults to 1 if undef
  Error code: 0
  Args    : a reference to an array containing the elements to be chainlinked
            an optional integer > 0 (this will be the starting count for
            the chain labels instead than having them begin from "1")

=cut

sub array2chain {
  my $arrayref=$_[0];
  my $array_count=scalar(@{$arrayref});
  unless ($array_count) {
    warn ("Warning array2chain: no elements input"); return (0); }
  my $begin=$_[1];
  if (defined $begin) {    
    if ($begin < 1) {
      warn "Warning array2chain: Zero or Negative offsets not allowed"; return (0); }
  } else {
    $begin=1;
  }
  my ($element,%hash);
  $hash{'begin'}=$begin;
  my $i=$begin-1;
  foreach $element (@{$arrayref}) {
    $i++;
    # hash with keys begin..end pointing to the arrays
    $hash{$i}=[$element,$i+1,$i-1];
  }
  my $end=$i;
  $hash{'end'}=$end;
  $hash{firstfree}=$i+1; # what a new added element should be called
  $hash{size}=$end-$begin+1; # how many elements in the chain

  # eliminate pointers to unexisting elements
  $hash{$begin}[2]=undef;
  $hash{$end}[1]=undef;

  return (\%hash);
}

1; # returns 1
