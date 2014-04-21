#---------------------------------------------------------------------------
# PACKAGE    : Bio::Tools::SeqWords
# PURPOSE    : To count n-mers in any sequence of characters
# AUTHOR     : Derek Gatherer (d.gatherer@vir.gla.ac.uk)
# SOURCE     : 
# CREATED    : 21st March 2000
# MODIFIED   : 11th November 2003 (DG - new method, count_overlap_words)
# LICENCE    : You may distribute this module under the same terms 
#	          : as the rest of BioPerl.
#---------------------------------------------------------------------------

=head1 NAME

Bio::Tools::SeqWords - Object holding n-mer statistics for a sequence

=head1 SYNOPSIS

  # Create the SeqWords object, e.g.:

  my $inputstream = Bio::SeqIO->new(-file => "seqfile", 
	                                 -format => 'Fasta');
  my $seqobj = $inputstream->next_seq();
  my $seq_word = Bio::Tools::SeqWords->new(-seq => $seqobj);

  # Or:
  my $seqobj = Bio::PrimarySeq->new(-seq => "agggtttccc",
                                    -alphabet => 'dna',
                                    -id => 'test');
  my $seq_word  =  Bio::Tools::SeqWords->new(-seq => $seqobj);

  # obtain a hash of word counts, eg:
  my $hash_ref = $seq_stats->count_words($word_length);

  # display hash table, eg:
  my %hash = %$hash_ref;
  foreach my $key(sort keys %hash)
  {
    print "\n$key\t$hash{$key}";
  }

  # Or:

  my $hash_ref =
     Bio::Tools::SeqWords->count_words($seqobj,$word_length);

=head1 DESCRIPTION

L<Bio::Tools::SeqWords> is a featherweight object for the calculation
of n-mer word occurrences in a single sequence.  It is envisaged that
the object will be useful for construction of scripts which use n-mer
word tables as the raw material for statistical calculations; for
instance, hexamer frequency for the calculation of coding protential,
or the calculation of periodicity in repetitive DNA.  Triplet
frequency is already handled by L<Bio::Tools::SeqStats> (author: Peter
Schattner).

There are a few possible applications for protein, e.g. hypothesised
amino acid 7-mers in heat shock proteins, or proteins with multiple
simple motifs.  Sometimes these protein periodicities are best seen
when the amino acid alphabet is truncated, e.g. Shulman alphabet.
Since there are quite a few of these shortened alphabets, this module
does not specify any particular alphabet.

See Synopsis above for object creation code.

=head2 Rationale

Take a sequence object and create an object for the purposes of
holding n-mer word statistics about that sequence. The sequence can be
nucleic acid or protein.

In count_words() the words are counted in a non-overlapping manner,
ie. in the style of a codon table, but with any word length.

In count_overlap_words() the words are counted in an overlapping
manner.

For counts on opposite strand (DNA/RNA), a reverse complement method
should be performed, and then the count repeated.

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

=head1 AUTHOR

Derek Gatherer, in the loosest sense of the word 'author'.  The
general shape of the module is lifted directly from the SeqStat module
of Peter Schattner. The central subroutine to count the words is
adapted from original code provided by Dave Shivak, in response to a
query on the bioperl mailing list.  At least 2 other people provided
alternative means (equally good but not used in the end) of performing
the same calculation.  Thanks to all for your assistance.

=head1 CONTRIBUTORS

Jason Stajich, jason-at-bioperl.org

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

package Bio::Tools::SeqWords;
use strict;


use base qw(Bio::Root::Root);

sub new {
    my($class,@args) = @_;
    # our new standard way of instantiation
    my $self = $class->SUPER::new(@args);

    my ($seqobj) = $self->_rearrange([qw(SEQ)],@args);
    if((! defined($seqobj)) && @args && ref($args[0])) {
	# parameter not passed as named parameter?
	$seqobj = $args[0];
    }
    
    if(! $seqobj->isa("Bio::PrimarySeqI")) { 
	$self->throw(ref($self) . " works only on PrimarySeqI objects\n");
    }
	
    $self->{'_seqref'} = $seqobj;
    return $self; 
}


=head2 count_words

 Title   : count_words
 Usage   : $word_count = $seq_stats->count_words($word_length)
                or 
           $word_count = $seq_stats->Bio::Tools::SeqWords->($seqobj,$word_length);
 Function: Counts non-overlapping words within a string, any alphabet is 
           used
 Example : a sequence ACCGTCCGT, counted at word length 4, will give the hash
           {ACCG => 1, TCCG => 1}
 Returns : Reference to a hash in which keys are words (any length) of the
           alphabet used and values are number of occurrences of the word 
           in the sequence.
 Args    : Word length as scalar and, reference to sequence object if
           required

           Throws an exception word length is not a positive integer
           or if word length is longer than the sequence.

=cut

sub count_words
{
    my ($self,$seqobj,$word_length) = @_;

    # check how we were called, and if necessary rearrange arguments
    if(ref($seqobj)) {
	# call as SeqWords->count_words($seq, $wordlen)
	if(! $seqobj->isa("Bio::PrimarySeqI")) { 
	    $self->throw("SeqWords works only on PrimarySeqI objects\n");
	}
    } else {
	# call as $obj->count_words($wordlen)
	$word_length = $seqobj;
	$seqobj = undef;
    }

    if(! defined($seqobj)){
	      $seqobj =  $self->{'_seqref'};
    }
    
    if($word_length eq "" || $word_length =~ /[a-z]/i){
	      $self->throw("SeqWords cannot accept non-numeric characters".
		     " or a null value in the \$word_length variable\n");
    }elsif ($word_length <1 || ($word_length - int($word_length)) >0){
	      $self->throw("SeqWords requires the word length to be a ".
		     "positive integer\n");
    }

    my $seqstring = uc $seqobj->seq();

    if($word_length > length($seqstring)){
	      $self->throw("die in _count, \$word_length is bigger ".
		    "than sequence length\n");
    }

    my $type = "non-overlap";
    my $words = _count($seqobj, $word_length, $type);
    return $words;   # ref. to a hash
}

=head2 count_overlap_words

 Title   : count_overlap_words
 Usage   : $word_count = $word_obj->count_overlap_words($word_length);
 Function: Counts overlapping words within a string, any alphabet is used
 Example : A sequence ACCAACCA, counted at word length 4, will give the hash
	        {ACCA=>2, CCAA=>1, CAAC=>1, AACC=>1}
 Returns : Reference to a hash in which keys are words (any length) of the 
           alphabet used and values are number of occurrences of the word in 
           the sequence.
 Args    : Word length as scalar

           Throws an exception if word length is not a positive integer
           or if word length is longer than the sequence.

=cut

sub count_overlap_words
{
    my ($self,$seqobj,$word_length) = @_;
 # check how we were called, and if necessary rearrange arguments
    if(ref($seqobj)){
	# call as SeqWords->count_words($seq, $wordlen)
	      if(! $seqobj->isa("Bio::PrimarySeqI")){
	          $self->throw("SeqWords works only on PrimarySeqI objects\n");
	      }
    }else{
	# call as $obj->count_words($wordlen)
	      $word_length = $seqobj;
	      $seqobj = undef;
    }

    if(! defined($seqobj)) {
	$seqobj =  $self->{'_seqref'};
    }
    my $seqstring = uc $seqobj->seq();

    if($word_length > length($seqstring)){
	      $self->throw("die in _count, \$word_length is bigger ".
		    "than sequence length\n");
    }
    
    my $type = "overlap";
    my $words = _count($seqobj, $word_length, $type);
    return $words;   # ref. to a hash
}

# the actual counting routine
# used by both count_words and count_overlap_words
sub _count {
    my ($seqobj, $word_length, $type) = @_;
    my %codon = ();

    # now the real business
    # JS - remove DNA assumption

    my $seqstring = uc $seqobj->seq();
    if($type eq "non-overlap")
    {
	while($seqstring =~ /((\w){$word_length})/gim){
	    $codon{uc($1)}++;
	}
    } elsif($type eq "overlap"){
	my $seqlen = $seqobj->length();	# measure length
	for (my $frame = 1; $frame <= $word_length; $frame++) { 
            # run through frames
	    my $seqstring = uc($seqobj->subseq($frame,$seqlen));
            # take the relevant substring
	    while($seqstring =~ /((\w){$word_length})/gim){
		$codon{uc($1)}++; # keep adding to hash
	    }
	}
    } else {
	Bio::Root::Root->throw("\nSomething badly wrong here. \$type: $type can only be overlap or non-overlap");
      }
    return \%codon;
}

1;
