# $Id$

#---------------------------------------------------------------------------
# PACKAGE    : SeqWords.pm
# PURPOSE    : To count n-mers in any sequence of characters
# AUTHOR     : Derek Gatherer (D.Gatherer@organon.nhe.akzonobel.nl)
# SOURCE     : 
# CREATED    : 21st March 2000
# MODIFIED   : 
# DISCLAIMER : I am employed in the pharmaceutical industry but my 
#	     : employers do not endorse or sponsor this module
#	     : in any way whatsoever.  The above email address is
#	     : given purely for the purpose of easy communication
#            : with the author, and does not imply any connection
#	     : between my employers and anything written below.
# LICENCE    : You may distribute this module under the same terms 
#	     : as the rest of BioPerl.
#---------------------------------------------------------------------------

=head1 NAME

Bio::Tools::SeqWords - Object holding n-mer statistics for one sequence

=head1 SYNOPSIS

Take a sequence object from eg, an inputstream, and creates an object
for the purposes of holding n-mer word statistics about that sequence.
The sequence can be nucleic acid or protein, but the module is
probably most relevant for DNA.  The words are counted in a
non-overlapping manner, ie. in the style of a codon table, but with
any word length.  For overlapping word counts, a sequence can be
'shifted' to remove the first character and then the count repeated.
For counts on opposite strand (DNA/RNA), a reverse complement method
should be performed, and then the count repeated.

Creating the SeqWords object, eg:

  my $inputstream = Bio::SeqIO->new( -file => "seqfile", 
	                             -format => 'Fasta');
  my $seqobj = $inputstream->next_seq();
  my $seq_word = Bio::Tools::SeqWords->new(-seq => $seqobj);

or:

  my $seqobj = Bio::PrimarySeq->new(-seq=>'[cut and paste a sequence here]', 
                                    -alphabet => 'dna', 
                                    -id => 'test');
  my $seq_word  =  Bio::Tools::SeqWords->new(-seq => $seqobj);

obtain a hash of word counts, eg:

  my $hash_ref = $seq_stats->count_words($word_length);

display hash table, eg:

  my %hash = %$hash_ref;
  foreach my $key(sort keys %hash)
  {
    print "\n$key\t$hash{$key}";
  }

or	

  my $hash_ref = Bio::SeqWords->count_words($seqobj,$word_length);


=head1 DESCRIPTION

Bio:SeqWords is a featherweight object for the calculation of n-mer
word occurrences in a single sequence.  It is envisaged that the
object will be useful for construction of scripts which use n-mer word
tables as the raw material for statistical calculations; for instance,
hexamer frequency for the calculation of coding protential, or the
calculation of periodicity in repetitive DNA.  Triplet frequency is
already handled by Bio::SeqStats.pm (author: Peter Schattner).  There
are a few possible applications for protein, eg: hypothesised amino
acid 7-mers in heat shock proteins, or proteins with multiple simple
motifs.  Sometimes these protein periodicities are best seen when the
amino acid alphabet is truncated, eg Shulman alphabet.  Since there
are quite a few of these shortened alphabets, this module does not
specify any particular alphabet.

See Synopsis above for object creation code.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this
and other Bioperl modules. Send your comments and suggestions preferably
to one of the Bioperl mailing lists.
Your participation is much appreciated.

  bioperl-l@bioperl.org                 - General discussion
  http://bio.perl.org/MailList.html     - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.
Bug reports can be submitted via email or the web:

  bioperl-bugs@bio.perl.org
  http://bugzilla.bioperl.org/

=head1 AUTHOR

Derek Gatherer, in the loosest sense of the word 'author'.  The
general shape of the module is lifted directly from Peter Schattner's
SeqStats.pm module.  The central subroutine to count the words is
adapted from original code provided by Dave Shivak, in response to a
query on the bioperl mailing list.  At least 2 other people provided
alternative means (equally good but not used in the end) of performing
the same calculation.  Thanks to all for your assistance.

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

package Bio::Tools::SeqWords;
use vars qw(@ISA);
use strict;

use Bio::Root::Root;

@ISA = qw(Bio::Root::Root);


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
 Usage   : $word_count = $seq_stats->count_words($word_length); 
 or 	 : $word_count = $seq_stats->Bio::SeqWords->($seqobj,$word_length);
 Function: Counts non-overlapping words within a string
	 : any alphabet is used
 Example : a sequence ACCGTCCGT, counted at word length 4,
	 : will give the hash
	 : ACCG 1, TCCG 1
 Returns : Reference to a hash in which keys are words (any length) of the
         : alphabet used and values are number of occurrences of the word 
         : in the sequence.
 Args    : Word length as scalar and, reference to sequence object if
         : required

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

	if($word_length eq "" || $word_length =~ /[a-z]/i)
	{
	    $self->throw("SeqWords cannot accept non-numeric characters".
                         " or a null value in the \$word_length variable\n");
	}
	elsif ($word_length <1 || ($word_length - int($word_length)) >0)
	{
	    $self->throw("SeqWords requires the word length to be a ".
                         "positive integer\n");
    	}
	
	if(! defined($seqobj)) {
	    $seqobj =  $self->{'_seqref'};
	}
	my $seqstring = uc $seqobj->seq();

	if($word_length > length($seqstring))
	{
	    $self->throw("die in count words, \$word_length is bigger ".
                         "than sequence length\n");
	}
	
	my %codon = ();

# now the real business

	while($seqstring =~ /(([ACGT]){$word_length})/gim)
	{
		$codon{uc($1)}++;
	}
	return \%codon;

# and that's it
}

1;
