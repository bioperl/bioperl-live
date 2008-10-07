# $Id$
#
# BioPerl module for Bio::LocatableSeq
#
# Cared for by Ewan Birney <birney@ebi.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::LocatableSeq - A Bio::PrimarySeq object with start/end points on it
that can be projected into a MSA or have coordinates relative to
another seq.

=head1 SYNOPSIS


    use Bio::LocatableSeq;
    my $seq = Bio::LocatableSeq->new(-seq => "CAGT-GGT",
				    -id  => "seq1",
				    -start => 1,
				    -end   => 7);


=head1 DESCRIPTION

    # a normal sequence object
    $locseq->seq();
    $locseq->id();

    # has start,end points
    $locseq->start();
    $locseq->end();

    # inherits off RangeI, so range operations possible

=head1 FEEDBACK


=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists.  Your participation is much appreciated.

  bioperl-l@bioperl.org                  - General discussion
  http://bioperl.org/wiki/Mailing_lists  - About the mailing lists

The locatable sequence object was developed mainly because the
SimpleAlign object requires this functionality, and in the rewrite
of the Sequence object we had to decide what to do with this.

It is, to be honest, not well integrated with the rest of bioperl, for
example, the trunc() function does not return a LocatableSeq object,
as some might have thought. There are all sorts of nasty gotcha's
about interactions between coordinate systems when these sort of
objects are used.


=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via the
web:

  http://bugzilla.open-bio.org/


=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

#'
# Let the code begin...

package Bio::LocatableSeq;
use strict;

use Bio::Location::Simple;
use Bio::Location::Fuzzy;
use vars qw($GAP_SYMBOLS $OTHER_SYMBOLS $MATCHPATTERN);

# should we change these to non-globals?  (I can see this
# causing problems down the road...) - cjfields

$GAP_SYMBOLS = '\-\.=~';

$OTHER_SYMBOLS = '\*\?';

$MATCHPATTERN = '0-9A-Za-z'.$GAP_SYMBOLS.$OTHER_SYMBOLS;

use base qw(Bio::PrimarySeq Bio::RangeI);

sub new {
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(@args);

    my ($start,$end,$strand, $tuple) =
	$self->_rearrange( [qw(START END STRAND TUPLE)],
			   @args);

    defined $start && $self->start($start);
    defined $end   && $self->end($end);
    defined $strand && $self->strand($strand);
    defined $tuple && $self->tuple($tuple);

    return $self; # success - we hope!
}

=head2 start

 Title   : start
 Usage   : $obj->start($newval)
 Function: Get/set the 1-based start position of this sequence in the original
           sequence. '0' means before the original sequence starts.
 Returns : value of start
 Args    : newvalue (optional)

=cut

sub start{
	my $self = shift;
	if( @_ ) {
		my $value = shift;
		$self->{'start'} = $value;
	}
	return $self->{'start'} if defined $self->{'start'};
	return 1                if $self->seq;
	return;
}

=head2 end

 Title   : end
 Usage   : $obj->end($newval)
 Function: Get/set the 1-based end position of this sequence in the original
           sequence. '0' means before the original sequence starts.
 Returns : value of end
 Args    : newvalue (optional)

=cut

sub end {
	my $self = shift;
	if( @_ ) {
		my $value = shift;
		my $string = $self->seq;
        # start of 0 usually means the sequence is blank (all gaps)
		if ($self->seq && $self->start != 0 ) {
			my $len = $self->_ungapped_len;
			my $id = $self->id;
			# changed 9/14/08
			if ($len != $value) {
				$self->warn("In sequence $id residue count gives end value ".
                "$len.  \nOverriding value [$value] with value $len for ".
                "Bio::LocatableSeq::end().\n$string");
				$value = $len;
			}
		}

		$self->{'end'} = $value;
    }

	return defined $self->{'end'} ? $self->{'end'} : $self->_ungapped_len;
}

sub _ungapped_len {
    my $self = shift;
    my $string = $self->seq || '';
    $string =~ s/[$GAP_SYMBOLS]+//g;
    $self->seq ? (return $self->start + CORE::length($string) - 1 ) : undef;
}

=head2 strand

 Title   : strand
 Usage   : $obj->strand($newval)
 Function: return or set the strandedness
 Returns : the value of the strandedness (-1, 0 or 1)
 Args    : the value of the strandedness (-1, 0 or 1)

=cut

sub strand{
   my $self = shift;
   if( @_ ) {
		my $value = shift;
		$self->{'strand'} = $value;
    }
    return $self->{'strand'};
}

=head2 get_nse

 Title   : get_nse
 Usage   :
 Function: read-only name of form id/start-end
 Example :
 Returns :
 Args    :

=cut

sub get_nse{
   my ($self,$char1,$char2) = @_;

   $char1 ||= "/";
   $char2 ||= "-";

   $self->throw("Attribute id not set") unless defined($self->id());
   $self->throw("Attribute start not set") unless defined($self->start());
   $self->throw("Attribute end not set") unless defined($self->end());

   return $self->id() . $char1 . $self->start . $char2 . $self->end ;

}


=head2 no_gap

 Title   : no_gaps
 Usage   :$self->no_gaps('.')
 Function:

           Gets number of gaps in the sequence. The count excludes
           leading or trailing gap characters.

           Valid bioperl sequence characters are [A-Za-z\-\.\*]. Of
           these, '.' and '-' are counted as gap characters unless an
           optional argument specifies one of them.

 Returns : number of internal gaps in the sequnce.
 Args    : a gap character (optional)

=cut

sub no_gaps {
    my ($self,$char) = @_;
    my ($seq, $count) = (undef, 0);

    # default gap characters
    $char ||= $GAP_SYMBOLS;

    $self->warn("I hope you know what you are doing setting gap to [$char]")
	unless $char =~ /[-.]/;

    $seq = $self->seq;
    return 0 unless $seq; # empty sequence does not have gaps

    $seq =~ s/^([$char]+)//;
    $seq =~ s/([$char]+)$//;
    $count++ while $seq =~ /[$char]+/g;

    return $count;

}


=head2 column_from_residue_number

 Title   : column_from_residue_number
 Usage   : $col = $seq->column_from_residue_number($resnumber)
 Function:

           This function gives the position in the alignment
           (i.e. column number) of the given residue number in the
           sequence. For example, for the sequence

  	     Seq1/91-97 AC..DEF.GH

           column_from_residue_number(94) returns 5.

           An exception is thrown if the residue number would lie
           outside the length of the aligment
           (e.g. column_from_residue_number( "Seq2", 22 )

 Returns : A column number for the position of the
           given residue in the given sequence (1 = first column)
 Args    : A residue number in the whole sequence (not just that
           segment of it in the alignment)

=cut

sub column_from_residue_number {
    my ($self, $resnumber) = @_;

    $self->throw("Residue number has to be a positive integer, not [$resnumber]")
	unless $resnumber =~ /^\d+$/ and $resnumber > 0;

    if ($resnumber >= $self->start() and $resnumber <= $self->end()) {
	my @residues = split //, $self->seq;
	my $count = $self->start();
	my $i;
	my ($start,$end,$inc,$test);
        my $strand = $self->strand || 0;
	# the following bit of "magic" allows the main loop logic to be the
	# same regardless of the strand of the sequence
	($start,$end,$inc,$test)= ($strand == -1)?
            (scalar(@residues-1),0,-1,sub{$i >= $end}) :
                (0,scalar(@residues-1),1,sub{$i <= $end});

	for ($i=$start; $test->(); $i+= $inc) {
	    if ($residues[$i] ne '.' and $residues[$i] ne '-') {
		$count == $resnumber and last;
		$count++;
	    }
	}
	# $i now holds the index of the column.
        # The actual column number is this index + 1

	return $i+1;
    }

    $self->throw("Could not find residue number $resnumber");

}

=head2 location_from_column

 Title   : location_from_column
 Usage   : $loc = $ali->location_from_column($column_number)
 Function:

           This function gives the residue number for a given position
           in the alignment (i.e. column number) of the given. Gaps
           complicate this process and force the output to be a
           L<Bio::Range> where values can be undefined. For example,
           for the sequence:

  	     Seq/91-97 .AC..DEF.G.

           location_from_column( 3 ) position 93
           location_from_column( 2 ) position 92^93
           location_from_column(10 ) position 97^98
           location_from_column( 1 ) position undef

           An exact position returns a Bio::Location::Simple object
           where where location_type() returns 'EXACT', if a position
           is between bases location_type() returns 'IN-BETWEEN'.
           Column before the first residue returns undef. Note that if
           the position is after the last residue in the alignment,
           that there is no guarantee that the original sequence has
           residues after that position.

           An exception is thrown if the column number is not within
           the sequence.

 Returns : Bio::Location::Simple or undef
 Args    : A column number
 Throws  : If column is not within the sequence

See L<Bio::Location::Simple> for more.

=cut

sub location_from_column {
    my ($self, $column) = @_;

    $self->throw("Column number has to be a positive integer, not [$column]")
	unless $column =~ /^\d+$/ and $column > 0;
    $self->throw("Column number [$column] is larger than".
		 " sequence length [". $self->length. "]")
	unless $column <= $self->length;

    my ($loc);
    my $s = $self->subseq(1,$column);
    $s =~ s/[^a-zA-Z\*]//g;

    my $pos = CORE::length $s;

    my $start = $self->start || 0 ;
    my $strand = $self->strand() || 1;
    my $relative_pos = ($strand == -1)
        ? ($self->end - $pos + 1)
	: ($pos + $start - 1);
    if ($self->subseq($column, $column) =~ /[a-zA-Z\*]/ ) {
	$loc = Bio::Location::Simple->new
	    (-start  => $relative_pos,
	     -end    => $relative_pos,
	     -strand => 1,
	     );
    } elsif ($pos == 0 and $self->start == 1) {
    } else {
      my ($start,$end) = ($relative_pos, $relative_pos + $strand);
      if ($strand == -1) {
	($start,$end) = ($end,$start);
      }
	$loc = Bio::Location::Simple->new
	    (-start         => $start,
	     -end           => $end,
	     -strand        => 1,
	     -location_type => 'IN-BETWEEN'
	     );
    }
    return $loc;
}

=head2 revcom

 Title   : revcom
 Usage   : $rev = $seq->revcom()
 Function: Produces a new Bio::LocatableSeq object which
           has the reversed complement of the sequence. For protein
           sequences this throws an exception of "Sequence is a
           protein. Cannot revcom"

 Returns : A new Bio::LocatableSeq object
 Args    : none

=cut

sub revcom {
    my ($self) = @_;

    my $new = $self->SUPER::revcom;
    $new->strand($self->strand * -1) if $self->strand;
    $new->start($self->start) if $self->start;
    $new->end($self->end) if $self->end;
    return $new;
}


=head2 trunc

 Title   : trunc
 Usage   : $subseq = $myseq->trunc(10,100);
 Function: Provides a truncation of a sequence,

 Example :
 Returns : a fresh Bio::PrimarySeqI implementing object
 Args    : Two integers denoting first and last columns of the
           sequence to be included into sub-sequence.


=cut

sub trunc {
    my ($self, $start, $end) = @_;
    my $new = $self->SUPER::trunc($start, $end);
    $new->strand($self->strand);

    # end will be automatically calculated
    $start = $end if $self->strand == -1;

    $start = $self->location_from_column($start);
    $start ? ($start = $start->end) : ($start = 1);
    $new->start($start) if $start;

    return $new;
}

=head2 validate_seq

 Title   : validate_seq
 Usage   : if(! $seq->validate_seq($seq_str) ) {
                print "sequence $seq_str is not valid for an object of
                alphabet ",$seq->alphabet, "\n";
	   }
 Function: Validates a given sequence string. A validating sequence string
           must be accepted by seq(). A string that does not validate will
           lead to an exception if passed to seq().

           The implementation provided here does not take alphabet() into
           account. Allowed are all letters (A-Z), numbers [0-9] 
           and '-','.','*','?','=',and '~'.

 Example :
 Returns : 1 if the supplied sequence string is valid for the object, and
           0 otherwise.
 Args    : The sequence string to be validated.


=cut

sub validate_seq {
	my ($self,$seqstr) = @_;
	if( ! defined $seqstr ){ $seqstr = $self->seq(); }
	return 0 unless( defined $seqstr);
	
	if((CORE::length($seqstr) > 0) &&
	   ($seqstr !~ /^([$MATCHPATTERN]+)$/)) {
	    $self->warn("seq doesn't validate, mismatch is " .
			join(",",($seqstr =~ /([^$MATCHPATTERN]+)/g)));
		return 0;
	}
	return 1;
}

1;
