# Bio::Tools::Alignment::Trim.pm
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
  $o_trim = new Bio::Tools::Alignment::Trim;
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

  bioperl-l@bioperl.org                 - General discussion
  http://bio.perl.org/MailList.html     - About the mailing
lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via
email or the web:

  bioperl-bugs@bio.perl.org
  http://bio.perl.org/bioperl-bugs/

=head1 AUTHOR - Chad Matsalla

Email bioinformatics@dieselwurks.com

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

package Bio::Tools::Alignment::Trim;

use strict;
use vars qw($VERSION @ISA);

$VERSION = '0.01';

# Preloaded methods go here.

my $f_designator = "f";
my $r_designator = "r";
#my $acefile;

=head2 new()

 Title   : new()
 Usage   : $o_trim = Bio::Tools::Alignment::Trim->new();
 Function: Construct the Bio::Tools::Alignment::Trim object. No parameters
	are required to create this object. It is strictly a bundle of
	functions, as far as I am concerned.
 Returns : A reference to a Bio::Tools::Alignment::Trim object.
 Args    : (none)

=cut 

sub new {
    my $class = shift;
    my $self = {};
    bless ($self,$class);
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
	($self->{f_designator},$self->{r_designator}) = @_;
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
	$self->{f_designator} = $desig;
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
	$self->{r_designator} = $desig;
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
	return("forward: ".$self->{f_designator}." reverse: ".$self->{r_designator}); 
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
	my %hash = %{$self->{qualities}};
} # end dump_hash

=head2 trim_singleton()

 Title   : trim_singleton()
 Usage   : $o_trim->trim_singleton()
 Function: Not implemented.
 Returns : Nothing.
 Args    : None.
 Notes   : Does nothing. Unimplemented.

=cut 

sub trim_singleton {
}

=head2 trim_singlet($sequence,$quality,$name,$class)

 Title   : trim_singlet($sequence,$quality,$name,$class)
 Usage   : ($r_trim_points,$trimmed_sequence) =
	@o_trim->trim_singlet($sequence,$quality,$name,$class);
 Function: Trim a singlet based on its quality.
 Returns : a reference to an array containing the forward and reverse
	trim points and the trimmed sequence.
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

sub trim_singlet {
    my ($self,$sequence,$quality,$name,$class) = @_;
    my @qual = split(' ',$quality);
    # print("\@qual: @qual\n");
    my @points;
    my $sequence_length = length($sequence);
    my ($returnstring,$processed_sequence);
    # print("Trim: \$sequence: $sequence\nTrim: \$name: $name\nTrim: \$f_designator: $self->{f_designator}\n");
    # print("Trim: \$r_designator: $self->{r_designator}\n");My own designators are".$self->{f_designator}."\n");
    # $points[0] = "$name: singlet that ends in the forward designator";
    # find out the leading and trailing trimpoints
    # for now, the rule for trailing points will be a run of $windowsize each less then 10phreds
    # print("Creating sliding window averages...\n");
    my $windowsize = 10;
    my $r_windows = &_sliding_window(\@qual,$windowsize);
    # print("Those sliding windows look like:\n");
    # &_print_formatted_qualities($r_windows);
    # my $r_windows = \@qual;
    my $windowtrail = 10;
    my $phreds = 20;
    # start_base required: r_quality,$windowsize,$phredvalue
    my $start_base = &_get_start($r_windows,5,20);
    if ($start_base > ($sequence_length - 100)) {
	# print("find_start: FAILED\n");
	$points[0] = ("FAILED");
	$points[1] = ("FAILED");
	return @points;
    }
    $points[0] = $start_base;
    #
    # whew! now for the end base
    # 
    # required parameters: reference_to_windows,windowsize,$phredvalue,start_base
    my $end_base = &_get_end($r_windows,20,20,$start_base);
    $points[1] = $end_base;
    # now do the actual trimming
    my @new_points = &chop_sequence($self,$name,$class,$sequence,@points);
    # substr EXPR,OFFSET,LEN,REPLACEMENT
    my $trimmed_sequence = pop(@new_points);
    return \@new_points,$trimmed_sequence;
}

=head2 trim_doublet($sequence,$quality,$name,$class)

 Title   : trim_doublet($sequence,$quality,$name,$class) 
 Usage   : ($r_trim_points,$trimmed_sequence) =
	@o_trim->trim_singlet($sequence,$quality,$name,$class);
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

sub trim_doublet {
    #
    my ($self,$sequence,$quality,$name,$class) = @_;
    my @qual = split(' ',$quality);
    # some debugging lines
    # print("CSM::Consed::trim_doublet: there are ".scalar(@qual)." quality values here.\n");
    # print("CSM::Consed::trim_doublet: working on $name.<br>\n");
    my @points;
    my $sequence_length = length($sequence);
    my ($returnstring,$processed_sequence);
    # print("Trim: \$sequence: $sequence\nTrim: \$name: $name\nTrim: \$f_designator: $self->{f_designator}\n");
    # print("Trim: \$r_designator: $self->{r_designator}\n");My own designators are".$self->{f_designator}."\n");
    # $points[0] = "$name: singlet that ends in the forward designator";
    # find out the leading and trailing trimpoints
    # for now, the rule for trailing points will be a run of $windowsize each less then 10phreds
    # print("Creating sliding window averages...\n");
    my $windowsize = 10;
    my $r_windows = &_sliding_window(\@qual,$windowsize);
    # print("Those sliding windows look like:\n");
    # &_print_formatted_qualities($r_windows);
    # my $r_windows = \@qual;
    my $windowtrail = 10;
    my $phreds = 20;
    # determine where the consensus sequence starts
    my $offset = 0;
    for (my $current = 0; $current<$sequence_length;$current++) {
	# print("\$current is $current\n");
	if ($qual[$current] != 0) {
	    $offset = $current;
	    last;
	}
    }
    # print("The offset for this doublet is $offset\n");
    # start_base required: r_quality,$windowsize,$phredvalue
    my $start_base = &_get_start($r_windows,5,20,$offset);
    if ($start_base > ($sequence_length - 100)) {
	# print("find_start: FAILED\n");
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
    # my $end_base = &_get_end($r_windows,20,20,$start_base);
    my $end_base = $sequence_length;
    my $start_of_trailing_zeros = &count_doublet_trailing_zeros(\@qual);
    # print("The trailing zero point at the end of this sequence ($name) is $start_of_trailing_zeros\n");
    # if ($start_of_trailing_zeros < $end_base) { $points[1] = $start_of_trailing_zeros; }
    $points[1] = $end_base;
    # now do the actual trimming
    my @new_points = &chop_sequence($self,$name,$class,$sequence,@points);
    # substr EXPR,OFFSET,LEN,REPLACEMENT
    my $trimmed_sequence = pop(@new_points);
    return @new_points,$trimmed_sequence;
}				# end trim_doublet

=head2 chop_sequence($name,$class,$sequence,@points)

 Title   : chop_sequence($name,$class,$sequence,@points)
 Usage   : ($start_point,$end_point,$chopped_sequence) = 
	$o_trim->chop_sequence($name,$class,$sequence,@points);
 Function: Chop a sequence based on its name, class, and sequence.
 Returns : an array containing three elements:
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
    my $fdesig = $self->{f_designator};
    my $rdesig = $self->{r_designator};
    if (!$points[0] && !$points[1]) {
	$sequence = "junk";
	return $sequence;
    }
    # print("\tchop_sequence: before trimming at $points[0] and $points[1], the sequence length was ".length($sequence)."\n");
    if ($class eq "singlet" && $name =~ /$fdesig$/) {
	# print("\tchop_sequence: chopping singlet sequence for $name with F designator\n");
	# EXPR,OFFSET,LEN
	$sequence = substr($sequence,$points[0],$points[1]-$points[0]);
    }
    elsif ($class eq "singlet" && $name =~ /$rdesig$/) {
	# print("\tchop_sequence: chopping singlet sequence for $name with R designator\n");
	$sequence = substr($sequence,$points[0],$points[1]-$points[0]);
    }		
    elsif ($class eq "singleton" && $name =~ /$fdesig$/) {
	# print("\tchop_sequence: chopping singleton sequence for $name with F designator\n");
	$sequence = substr($sequence,$points[0],$points[1]-$points[0]);
    }
    elsif ($class eq "singleton" && $name =~ /$rdesig$/) {
	# print("\tchop__sequence: chopping singleton sequence for $name with R designator\n");
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
	# print("\tTrimmed $number_Ns_trimmed N's from the end of the sequence\n");
	$points[1] -= $number_Ns_trimmed;
	$points[1] -= 1;
    }
    # print("\tchop_sequence: I chopped $subs_Xs for X\'s and $subs_Ns for N\'s\n");
    # print("\tchop_sequence: after trimming, the sequence length was ".length($sequence)."\n");
    push @points,$sequence;
    return @points;
}				# end chop_sequence

=head2 _get_start($r_quals,$windowsize,$phreds,$offset)

 Title   : _get_start($r_quals,$windowsize,$phreds,$offset)
 Usage   : $start_base = &_get_start($r_windows,5,20);
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
    my ($r_quals,$windowsize,$phreds,$offset) = @_;
    # this is to help determine whether the sequence is good at all
    # print("r_quals,windowsize,phreds = #$r_quals#,#$windowsize#,#$phreds#\n");
    my @quals = @$r_quals;
    my ($count,$count2,$qualsum);
    if ($offset) { $count = $offset; } else { $count = 0; }
    for (; ($count+$windowsize) <= scalar(@quals); $count++) {
	# print("On window starting at base $count\n");
	for($count2 = $count; $count2 < $count+$windowsize-1; $count2++) {
				# print("Adding base $count2 ($quals[$count2]). Running total:");
	    unless (!$quals[$count2]) {
		$qualsum += $quals[$count2];
	    }
				# print("qualsum: $qualsum\n");

	}
	if ($qualsum && $qualsum >= $windowsize*$phreds) {
				# print("\tset_start: $qualsum is greater then windowsize $windowsize times $phreds phreds in window starting at $count\n");
	    return $count;
	}
	$qualsum = 0;
    }
    # if ($count > scalar(@quals)-$windowsize) { return; }
    return $count;
}				# end get_start

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
		# print("rqual,windowsize,phreds: #$r_qual#,#$windowsize#,#$phreds#\n");
		# print("CSM::Consed::Trim::get_end: Starting the search for the end at $count\n");
	my @quals = @$r_qual;
	my $total_bases = scalar(@quals);
		# print("CSM::Consed::Trim::get_end: \@quals are @quals\n");
		# print("CSM::Consed::Trim::get_end: \$total_bases: $total_bases\n");
	my ($count2,$qualsum,$end_of_quals,$bases_counted);
	if (!$count) { $count=0; }
	BASE: for (; $count < $total_bases; $count++) {
			# print("Looking for quality window starting at base $count\n");
		$bases_counted = 0;
		$qualsum = 0;
		POSITION: for($count2 = $count; $count2 < $total_bases; $count2++) {
			$bases_counted++;
				# print("CSM::Consed::Trim::get_end: \$bases_counted is $bases_counted\n");
			if ($count2 == $total_bases-1) {
				$qualsum += $quals[$count2];
					# print("CSM::Consed::Trim::get_end: Hit the end of the quals. Current start base: $count \$bases_counted: $bases_counted\n");
					# print("CSM::Consed::Trim::get_end: Hit the end of the quals.\n");
				$bases_counted++;
				last BASE;
			}
			elsif ($bases_counted == $windowsize) {
						# print("The number of bases counted equals the windowsize. Moving to the next window\n");
						# print("adding quality $quals[$count2] at position $count2 to the total for the window starting at $count2\n");
				$qualsum += $quals[$count2];
				if ($qualsum < $bases_counted*$phreds) {
					print("CSM::Consed::Trim::get_end: Found a quality problem at $count:\n\$qualsum: $qualsum\n\$bases_counted: $bases_counted\n");
					return $count+$bases_counted+$windowsize;
				}
				next BASE;
			}
			else {
				$qualsum += $quals[$count2];
			}
		}
		print("Looking for the end in a window starting at $count: \$qualsum there was  $qualsum\n");
		if ($qualsum < $bases_counted*$phreds) {
			print("CSM::Consed::Trim::get_end: Found a quality problem at $count:\n\$qualsum: $qualsum\n\$bases_counted: $bases_counted\n");
			return $count+$bases_counted+$windowsize;
		}
		else {
			print("CSM::Consed::Trim::get_end: looks ok here at bast $count: \$bases_counted times \$phreds is ".$bases_counted*$phreds." which is greater then $qualsum\n"); 
		}
		$qualsum = 0;
	} # end for
	if ($end_of_quals) {
		my $bases_for_average = $total_bases-$count2;
		print("CSM::Consed::Trim::get_end: Ran out of quals! count2 being $count2 and \$qualsum being $qualsum and \$bases_for_average is $bases_for_average\n");
		return $count2;
	}
	else {
		print("Nope. Not the end of the quals yet\n");
	}
	print("I didn't find an end anywhere. Search ended at base number $count where windowsize*phreds was ".$windowsize*$phreds." but qualsum was only ");
		if ($qualsum) { print("$qualsum\n"); }
		# else { print("0\n"); }
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
				# print("CSM::Consed::Trim::count_doublet_trailing_zeros:  The current quality $current is not zero. It is $qualities[$current]\n");
				# print("returning $number_of_trailing_zeros\n");
				# return $number_of_trailing_zeros;
				# print("CSM::Consed::Trim::count_doublet_trailing_zeros: returning $current\n");
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
	   should use that object but I haven\'t rewritten this yet.

=cut 

sub _sliding_window {
    my ($r_quals,$windowsize) = @_;
    my (@window,@quals,$qualsum,$count,$count2,$average,@averages,$bases_counted);
    @quals = @$r_quals;
    # print("@quals\n");
    my $size_of_quality = scalar(@quals);
    # print("CSM::Consed::Trim::sliding_window: there are $size_of_quality quality values to build windows for\n");
    for ($count=0; $count <= $size_of_quality; $count++) {
				# print("SlidingWindow: Averaging bases from $count\n");
				# original: 010105
				# for($count2 = $count; $count2 < $count+$windowsize; $count2++) {
	$bases_counted = 0;
      BASE: for($count2 = $count; $count2 < $size_of_quality; $count2++) {
	  $bases_counted++;
				# print("Adding base at $count2 ($quals[$count2]) \$qualsum is $qualsum\n");
				# if the search hits the end of the averages, stop
				# this is for the case near the end where bases remaining < windowsize
	  if ($count2 == $size_of_quality) {
				# print("Hit the end of the quality values.\n");
	      $qualsum += $quals[$count2];
	      last BASE;
	  }				
				# if the search hits the size of the window
				# 010116 this is wrong! elsif ($count2 == $windowsize) {
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
	}
	if (!$average) { $average = "0"; }
	# print("At position $count, \$qualsum is $qualsum and \$average is $average\n");
	# print("Pushing average $average onto the averages array\n");
	# print("Counted $bases_counted to get an average of $average.\n");
	push @averages,$average;
	$qualsum = 0;
	# print("Done with count2: final value: $count2\n");
    }
    # print("CSM::Consed::Trim::sliding_window: there are ".scalar(@averages)." averages here\n");
    return \@averages;
}				# end sliding_window

=head2 _print_formatted_qualities(\@quals)

 Title   : _print_formatted_qualities(\@quals)
 Usage   : &_print_formatted_qualities(\@quals);
 Returns : Nothing. Prints.
 Args    : A reference to an array containing quality values.
 Notes   : An internal procedure used in debugging. Prints out an array nicely.

=cut 

sub _print_formatted_qualities {
        my $rquals = shift;
		# print("print_formatted_qualities: \$rquals is $rquals\n");
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
 Usage   : Deprecated. Don\'t use this!
 Returns : Deprecated. Don\'t use this!
 Args    : Deprecated. Don\'t use this!

=cut 

sub _get_end_old {
    my ($r_qual,$windowsize,$phreds,$count) = @_;
    # print("rqual,windowsize,phreds: #$r_qual#,#$windowsize#,#$phreds#\n");
    print("Starting the serch for the end at $count\n");
    my $target = $windowsize*$phreds;
    my @quals = @$r_qual;
    my $total_bases = scalar(@quals);
    my ($count2,$qualsum,$end_of_quals);
    if (!$count) { $count=0; }
  BASE: for (; $count < $total_bases; $count++) {
      print("Looking for quality window starting at base $count\n");
      for($count2 = $count; $count2 < $count+$windowsize; $count2++) {
	  if ($count2 == scalar(@quals)-1) {
	      print("Hit the end of the quals. Summing...\n");
	      $qualsum += $quals[$count2];
	      $end_of_quals = 1;
	      last BASE;

	  }
				# print("adding position $quals[$count2] $count2 to the total for the window starting at $count2\n");
	  $qualsum += $quals[$count2];
      }
      # print("Looking for the end in a window starting at $count: \$qualsum there was  $qualsum\n");
      if ($qualsum < $windowsize*$phreds) {
	  print("\tget_end: $qualsum is less then a window of $windowsize times $phreds phreds in window starting at $count\n");
	  return $count+$windowsize;
      }
      $qualsum = 0;
  }				# end for


}				# end get_end_old


# Autoload methods go after =cut, and are processed by the autosplit program.

1;
__END__
