# Bio::Tools::Alignment::Trim.pm
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Chad Matsalla
#
# Copyright Chad Matsalla
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code 

=head1 NAME 

Bio::Tools::Alignment::Trim - A kludge to do specialized trimming of
	sequence based on quality.

=head1 SYNOPSIS

  use Bio::Tools::Alignment::Trim;
  $o_trim = Bio::Tools::Alignment::Trim->new();
  $o_trim->set_reverse_designator("R");
  $o_trim->set_forward_designator("F");


=head1 DESCRIPTION

This is a specialized module designed by Chad for Chad to trim sequences
based on a highly specialized list of requirements. In other words, write
something that will trim sequences 'just like the people in the lab would
do manually'.

I settled on a sliding-window-average style of search which is ugly and
slow but does _exactly_ what I want it to do.

Mental note: rewrite this.

It is very important to keep in mind the context in which this module was
written: strictly to support the projects for which Consed.pm was
designed.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists.  Your participation is much appreciated.

  bioperl-l@bioperl.org                  - General discussion
  http://bioperl.org/wiki/Mailing_lists     - About the mailing
lists

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

=head1 AUTHOR - Chad Matsalla

Email bioinformatics-at-dieselwurks.com

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

package Bio::Tools::Alignment::Trim;

use strict;
use Dumpvalue;

use vars qw(%DEFAULTS);

use base qw(Bio::Root::Root);

BEGIN {
    %DEFAULTS = ( 'f_designator' => 'f',
		  'r_designator' => 'r',
                  'windowsize' => '10',
                  'phreds' => '20');
}

=head2 new()

 Title   : new()
 Usage   : $o_trim = Bio::Tools::Alignment::Trim->new();
 Function: Construct the Bio::Tools::Alignment::Trim object. No parameters
	   are required to create this object. It is strictly a bundle of
	   functions, as far as I am concerned.
 Returns : A reference to a Bio::Tools::Alignment::Trim object.
 Args    : (optional)
           -windowsize (default 10)
           -phreds (default 20)


=cut 

sub new {
    my ($class,@args) = @_;    
    my $self = $class->SUPER::new(@args);
    my($windowsize,$phreds) =
        $self->_rearrange([qw(
                    WINDOWSIZE
                    PHREDS
                              )],
                          @args);
    $self->{windowsize} = $windowsize || $DEFAULTS{'windowsize'};
    $self->{phreds} = $phreds || $DEFAULTS{'phreds'};
          # print("Constructor set phreds to ".$self->{phreds}."\n") if $self->verbose > 0;
    $self->set_designators($DEFAULTS{'f_designator'},
			   $DEFAULTS{'r_designator'});
    return $self;
}

=head2 set_designators($forward_designator,$reverse_designator)

 Title   : set_designators(<forward>,<reverse>)
 Usage   : $o_trim->set_designators("F","R")
 Function: Set the string by which the system determines whether a given
	sequence represents a forward or a reverse read.
 Returns : Nothing.
 Args    : two scalars: one representing the forward designator and one
	representing the reverse designator

=cut 

sub set_designators {
    my $self = shift;
    ($self->{'f_designator'},$self->{'r_designator'}) = @_;
}

=head2 set_forward_designator($designator)

 Title   : set_forward_designator($designator)
 Usage   : $o_trim->set_forward_designator("F")
 Function: Set the string by which the system determines if a given
	sequence is a forward read.
 Returns : Nothing.
 Args    : A string representing the forward designator of this project.

=cut 

sub set_forward_designator {
    my ($self,$desig) = @_;
    $self->{'f_designator'} = $desig;
}

=head2 set_reverse_designator($reverse_designator)

 Title   : set_reverse_designator($reverse_designator)
 Function: Set the string by which the system determines if a given
	sequence is a reverse read.
 Usage   : $o_trim->set_reverse_designator("R")
 Returns : Nothing.
 Args    : A string representing the forward designator of this project.

=cut 

sub set_reverse_designator {
    my ($self,$desig) = @_;
    $self->{'r_designator'} = $desig;
}

=head2 get_designators()

 Title   : get_designators()
 Usage   : $o_trim->get_designators()
 Returns : A string describing the current designators.
 Args    : None
 Notes   : Really for informational purposes only. Duh.

=cut 

sub get_designators {
    my $self = shift;
    return("forward: ".$self->{'f_designator'}." reverse: ".$self->{'r_designator'}); 
}	

=head2 trim_leading_polys()

 Title   : trim_leading_polys()
 Usage   : $o_trim->trim_leading_polys()
 Function: Not implemented. Does nothing.
 Returns : Nothing.
 Args    : None.
 Notes   : This function is not implemented. Part of something I wanted to
	do but never got around to doing.

=cut 

sub trim_leading_polys {
	my ($self, $sequence) = @_;
}

=head2 dump_hash()

 Title   : dump_hash()
 Usage   : $o_trim->dump_hash()
 Function: Unimplemented.
 Returns : Nothing.
 Args    : None.
 Notes   : Does nothing.

=cut 

sub dump_hash {
	my $self = shift;
	my %hash = %{$self->{'qualities'}};
} # end dump_hash

=head2 trim_singlet($sequence,$quality,$name,$class)

 Title   : trim_singlet($sequence,$quality,$name,$class)
 Usage   : ($r_trim_points,$trimmed_sequence) =
	@{$o_trim->trim_singlet($sequence,$quality,$name,$class)};
 Function: Trim a singlet based on its quality.
 Returns : a reference to an array containing the forward and reverse
	trim points and the trimmed sequence.
 Args    : $sequence : A sequence (SCALAR, please)
	   $quality : A _scalar_ of space-delimited quality values.
	   $name : the name of the sequence
	   $class : The class of the sequence. One of qw(singlet
		singleton doublet pair multiplet)
 Notes   : At the time this was written the bioperl objects SeqWithQuality
	and PrimaryQual did not exist. This is what is with the clumsy
	passing of references and so on. I will rewrite this next time I
	have to work with it. I also wasn't sure whether this function
	should return just the trim points or the points and the sequence.
	I decided that I always wanted both so that's how I implemented
	it.
     - Note that the size of the sliding windows is set during construction of
       the Bio::Tools::Alignment::Trim object.

=cut 

sub trim_singlet {
    my ($self,$sequence,$quality,$name,$class) = @_;
    # this split is done because I normally store quality values in a
    # space-delimited scalar rather then in an array.
    # I do this because serialization of the arrays is tough.
    my @qual = split(' ',$quality);
    my @points;
    my $sequence_length = length($sequence);
    my ($returnstring,$processed_sequence);
    # smooth out the qualities
    my $r_windows = &_sliding_window(\@qual,$self->{windowsize});
    # find out the leading and trailing trimpoints
    my $start_base = $self->_get_start($r_windows,$self->{windowsize},$self->{phreds});
    my (@new_points,$trimmed_sequence);
    # do you think that any sequence shorter then 100 should be
    # discarded? I don't think that this should be the decision of this
    # module.
    # removed, 020926
    $points[0] = $start_base;
    # whew! now for the end base
    # required parameters: reference_to_windows,windowsize,$phredvalue,start_base
    my $end_base = &_get_end($r_windows,$self->{windowsize},
			     $self->{phreds},$start_base);
    $points[1] = $end_base;
    # now do the actual trimming
    # CHAD : I don't think that it is a good idea to call chop_sequence here
    # because chop_sequence also removes X's and N's and things
    # and that is not always what is wanted
    return \@points;
}

=head2 trim_doublet($sequence,$quality,$name,$class)

 Title   : trim_doublet($sequence,$quality,$name,$class) 
 Usage   : ($r_trim_points,$trimmed_sequence) =
	    @{$o_trim->trim_singlet($sequence,$quality,$name,$class)};
 Function: Trim a singlet based on its quality.
 Returns : a reference to an array containing the forward and reverse
 Args    : $sequence : A sequence
	   $quality : A _scalar_ of space-delimited quality values.
	   $name : the name of the sequence
	   $class : The class of the sequence. One of qw(singlet
		singleton doublet pair multiplet)
 Notes   : At the time this was written the bioperl objects SeqWithQuality
	and PrimaryQual did not exist. This is what is with the clumsy
	passing of references and so on. I will rewrite this next time I
	have to work with it. I also wasn't sure whether this function
	should return just the trim points or the points and the sequence.
	I decided that I always wanted both so that's how I implemented
	it.

=cut 

#'
sub trim_doublet {
    my ($self,$sequence,$quality,$name,$class) = @_;
    my @qual = split(' ',$quality);
    my @points;
    my $sequence_length = length($sequence);
    my ($returnstring,$processed_sequence);
          # smooth out the qualities
    my $r_windows = &_sliding_window(\@qual,$self->{windowsize});
          # determine where the consensus sequence starts
    my $offset = 0;
    for (my $current = 0; $current<$sequence_length;$current++) {
          if ($qual[$current] != 0) {
               $offset = $current;
               last;
          }
    }
          # start_base required: r_quality,$windowsize,$phredvalue
    my $start_base = $self->_get_start($r_windows,$self->{windowsize},$self->{phreds},$offset);
    if ($start_base > ($sequence_length - 100)) {
          $points[0] = ("FAILED");
	     $points[1] = ("FAILED");
          return @points;
     }
    $points[0] = $start_base;
         #
         # whew! now for the end base
         # 
         # required parameters: reference_to_windows,windowsize,$phredvalue,start_base
         #								    |	
         # 010420 NOTE: We will no longer get the end base to avoid the Q/--\___/-- syndrome
    my $end_base = $sequence_length;
    my $start_of_trailing_zeros = &count_doublet_trailing_zeros(\@qual);
    $points[1] = $end_base;
          # CHAD : I don't think that it is a good idea to call chop_sequence here
          # because chop_sequence also removes X's and N's and things
          # and that is not always what is wanted
     return @points;
}				# end trim_doublet

=head2 chop_sequence($name,$class,$sequence,@points)

 Title   : chop_sequence($name,$class,$sequence,@points)
 Usage   : ($start_point,$end_point,$chopped_sequence) = 
	$o_trim->chop_sequence($name,$class,$sequence,@points);
 Function: Chop a sequence based on its name, class, and sequence.
 Returns : an array containing three scalars:
	1- the start trim point
	2- the end trim point
	3- the chopped sequence
 Args    :
	   $name : the name of the sequence
	   $class : The class of the sequence. One of qw(singlet
		singleton doublet pair multiplet)
	   $sequence : A sequence
	   @points : An array containing two elements- the first contains
		the start trim point and the second conatines the end trim
		point.

=cut

sub chop_sequence {
    my ($self,$name,$class,$sequence,@points) = @_;
     print("Coming into chop_sequence, \@points are @points\n");
    my $fdesig = $self->{'f_designator'};
    my $rdesig = $self->{'r_designator'};
    if (!$points[0] && !$points[1]) {
	$sequence = "junk";
	return $sequence;
    }
    if ($class eq "singlet" && $name =~ /$fdesig$/) {
	$sequence = substr($sequence,$points[0],$points[1]-$points[0]);
    }
    elsif ($class eq "singlet" && $name =~ /$rdesig$/) {
	$sequence = substr($sequence,$points[0],$points[1]-$points[0]);
    }		
    elsif ($class eq "singleton" && $name =~ /$fdesig$/) {
	$sequence = substr($sequence,$points[0],$points[1]-$points[0]);
    }
    elsif ($class eq "singleton" && $name =~ /$rdesig$/) {
	$sequence = substr($sequence,$points[0],$points[1]-$points[0]);
    }
    elsif ($class eq "doublet") {
	$sequence = substr($sequence,$points[0],$points[1]-$points[0]);
    }
    # this is a _terrible_ to do this! i couldn't seem to find a better way
    # i thought something like s/(^.*[Xx]{5,})//g; might work, but no go
    # no time to find a fix!
    my $length_before_trimming = length($sequence);
    my $subs_Xs = $sequence =~ s/^.*[Xx]{5,}//g;
    if ($subs_Xs) {
	my $length_after_trimming = length($sequence);
	my $number_Xs_trimmed = $length_before_trimming - $length_after_trimming;
	$points[0] += $number_Xs_trimmed;
    }
    $length_before_trimming = length($sequence);
    my $subs_Ns = $sequence =~ s/[Nn]{1,}$//g;
    if ($subs_Ns) {
	my $length_after_trimming = length($sequence);
	my $number_Ns_trimmed = $length_before_trimming - $length_after_trimming;
	$points[1] -= $number_Ns_trimmed;
	$points[1] -= 1;
    }
     push @points,$sequence;
     print("chop_sequence \@points are @points\n");
    return @points;
}

=head2 _get_start($r_quals,$windowsize,$phreds,$offset)

 Title   : _get_start($r_quals,$windowsize,$phreds,$offset)
 Usage   : $start_base = $self->_get_start($r_windows,5,20);
 Function: Provide the start trim point for this sequence.
 Returns : a scalar representing the start of the sequence
 Args    : 
	$r_quals : A reference to an array containing quality values. In
		context, this array of values has been smoothed by then
		sliding window-look ahead algorithm.
	$windowsize : The size of the window used when the sliding window
		look-ahead average was calculated.
	$phreds : <fill in what this does here>
	$offset : <fill in what this does here>

=cut 

sub _get_start {
    my ($self,$r_quals,$windowsize,$phreds,$offset) = @_;
     print("Using $phreds phreds\n")  if $self->verbose > 0;
          # this is to help determine whether the sequence is good at all
    my @quals = @$r_quals;
    my ($count,$count2,$qualsum);
    if ($offset) { $count = $offset; } else { $count = 0; }
          # search along the length of the sequence
    for (; ($count+$windowsize) <= scalar(@quals); $count++) {
               # sum all of the quality values in this window.
          my $cumulative=0;
          for($count2 = $count; $count2 < $count+$windowsize; $count2++) {
               if (!$quals[$count2]) {
                         # print("Quals don't exist here!\n");
               }
               else {
                    $qualsum += $quals[$count2]; 
                         # print("Incremented qualsum to ($qualsum)\n");
               }
               $cumulative++;
          }
               # print("The sum of this window (starting at $count) is $qualsum. I counted $cumulative bases.\n");
               # if the total of windowsize * phreds is 
          if ($qualsum && $qualsum >= $windowsize*$phreds) { return $count; }
	     $qualsum = 0;
    }
    # if ($count > scalar(@quals)-$windowsize) { return; }
    return $count;
}	

=head2 _get_end($r_qual,$windowsize,$phreds,$count)

 Title   : _get_end($r_qual,$windowsize,$phreds,$count)
 Usage   : my $end_base = &_get_end($r_windows,20,20,$start_base);
 Function: Get the end trim point for this sequence.
 Returns : A scalar representing the end trim point for this sequence.
 Args    : 
	$r_qual : A reference to an array containing quality values. In
		context, this array of values has been smoothed by then
		sliding window-look ahead algorithm.
	$windowsize : The size of the window used when the sliding window
		look-ahead average was calculated.
	$phreds : <fill in what this does here>
	$count : Start looking for the end of the sequence here.

=cut 

sub _get_end {
    my ($r_qual,$windowsize,$phreds,$count) = @_;
    my @quals = @$r_qual;
    my $total_bases = scalar(@quals);
    my ($count2,$qualsum,$end_of_quals,$bases_counted);
    if (!$count) { $count=0; }
  BASE: for (; $count < $total_bases; $count++) {
      $bases_counted = 0;
      $qualsum = 0;
    POSITION: for($count2 = $count; $count2 < $total_bases; $count2++) {
	$bases_counted++;

	if ($count2 == $total_bases-1) {
	    $qualsum += $quals[$count2];
	    $bases_counted++;
	    last BASE;
	}
	elsif ($bases_counted == $windowsize) {
	    $qualsum += $quals[$count2];
	    if ($qualsum < $bases_counted*$phreds) {
		return $count+$bases_counted+$windowsize;
	    }
	    next BASE;
	}
	else {
	    $qualsum += $quals[$count2];
	}
    }
      if ($qualsum < $bases_counted*$phreds) {
	  return $count+$bases_counted+$windowsize;
      }
      else { }
      $qualsum = 0;
  }				# end for
    if ($end_of_quals) {
	my $bases_for_average = $total_bases-$count2;
	return $count2;
    }
    else { }
    if ($qualsum) { } # print ("$qualsum\n");
    return $total_bases;
} # end get_end

=head2 count_doublet_trailing_zeros($r_qual)

 Title   : count_doublet_trailing_zeros($r_qual)
 Usage   : my $start_of_trailing_zeros = &count_doublet_trailing_zeros(\@qual);
 Function: Find out when the trailing zero qualities start.
 Returns : A scalar representing where the zeros start.
 Args    : A reference to an array of quality values.
 Notes   : Again, this should be rewritten to use PrimaryQual objects.
	A more detailed explanation of why phrap puts these zeros here should
	be written and placed here. Please email and hassle the author.


=cut 

sub count_doublet_trailing_zeros {
	my ($r_qual) = shift;
	my $number_of_trailing_zeros = 0;
	my @qualities = @$r_qual;
	for (my $current=scalar(@qualities);$current>0;$current--) {
		if ($qualities[$current] && $qualities[$current] != 0) {
			$number_of_trailing_zeros = scalar(@qualities)-$current;
			return $current+1;
		}
	}
	return scalar(@qualities);
} # end count_doublet_trailing_zeros

=head2 _sliding_window($r_quals,$windowsize)

 Title   : _sliding_window($r_quals,$windowsize)
 Usage   : my $r_windows = &_sliding_window(\@qual,$windowsize);
 Function: Create a sliding window, look-forward-average on an array
	of quality values. Used to smooth out differences in qualities.
 Returns : A reference to an array containing the smoothed values.
 Args    : $r_quals: A reference to an array containing quality values.
	   $windowsize : The size of the sliding window.
 Notes   : This was written before PrimaryQual objects existed. They
	   should use that object but I haven't rewritten this yet.

=cut 

#'
sub _sliding_window {
    my ($r_quals,$windowsize) = @_;
    my (@window,@quals,$qualsum,$count,$count2,$average,@averages,$bases_counted);
    @quals = @$r_quals;    
    my $size_of_quality = scalar(@quals);
          # do this loop for all of the qualities
     for ($count=0; $count <= $size_of_quality; $count++) {
          $bases_counted = 0;
          BASE: for($count2 = $count; $count2 < $size_of_quality; $count2++) {
               $bases_counted++;
                    # if the search hits the end of the averages, stop
                    # this is for the case near the end where bases remaining < windowsize
               if ($count2 == $size_of_quality) {
                    $qualsum += $quals[$count2];
                    last BASE;
               }				
                    # if the search hits the size of the window
               elsif ($bases_counted == $windowsize) {
                    $qualsum += $quals[$count2];
                    last BASE;
               }
                    # otherwise add the quality value
               unless (!$quals[$count2]) {
                    $qualsum += $quals[$count2];
               }
          }
          unless (!$qualsum || !$windowsize) {
              $average = $qualsum / $bases_counted;
               if (!$average) { $average = "0"; }
     	     push @averages,$average;
          }
	     $qualsum = 0;
     }
          # 02101 Yes, I repaired the mismatching numbers between averages and windows.
          # print("There are ".scalar(@$r_quals)." quality values. They are @$r_quals\n");
          # print("There are ".scalar(@averages)." average values. They are @averages\n");
    return \@averages;
     
}

=head2 _print_formatted_qualities

 Title   : _print_formatted_qualities(\@quals)
 Usage   : &_print_formatted_qualities(\@quals);
 Returns : Nothing. Prints.
 Args    : A reference to an array containing quality values.
 Notes   : An internal procedure used in debugging. Prints out an array nicely.

=cut 

sub _print_formatted_qualities {
    my $rquals = shift;
    my @qual = @$rquals;
    for (my $count=0; $count<scalar(@qual) ; $count++) {
	if (($count%10)==0) { print("\n$count\t"); }
	if ($qual[$count]) { print ("$qual[$count]\t");}
	else { print("0\t"); }
    }
    print("\n");
}

=head2 _get_end_old($r_qual,$windowsize,$phreds,$count)

 Title   : _get_end_old($r_qual,$windowsize,$phreds,$count)
 Usage   : Deprecated. Don't use this!
 Returns : Deprecated. Don't use this!
 Args    : Deprecated. Don't use this!

=cut 

#'
sub _get_end_old {
    my ($r_qual,$windowsize,$phreds,$count) = @_;
    warn("Do Not Use this function (_get_end_old)");
    my $target = $windowsize*$phreds;
    my @quals = @$r_qual;
    my $total_bases = scalar(@quals);
    my ($count2,$qualsum,$end_of_quals);
    if (!$count) { $count=0; }
  BASE: for (; $count < $total_bases; $count++) {
      for($count2 = $count; $count2 < $count+$windowsize; $count2++) {
	  if ($count2 == scalar(@quals)-1) {
	      $qualsum += $quals[$count2];
	      $end_of_quals = 1;
	      last BASE;

	  }
	  $qualsum += $quals[$count2];
      }
      if ($qualsum < $windowsize*$phreds) {
	  return $count+$windowsize;
      }
      $qualsum = 0;
  }				# end for
}				# end get_end_old


# Autoload methods go after =cut, and are processed by the autosplit program.

1;
__END__
