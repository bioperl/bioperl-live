#$Id$
#-----------------------------------------------------------------------------
# PACKAGE    : OddCodes.pm
# PURPOSE    : To write amino acid sequences in alternative alphabets
# AUTHOR     : Derek Gatherer (D.Gatherer@organon.nhe.akzonobel.nl)
# SOURCE     :
# CREATED    : 8th July 2000
# MODIFIED   :
# DISCLAIMER : I am employed in the pharmaceutical industry but my
#            : employers do not endorse or sponsor this module
#            : in any way whatsoever.  The above email address is
#            : given purely for the purpose of easy communication
#            : with the author, and does not imply any connection
#	     : between my employers and anything written below.
# LICENCE    : You may distribute this module under the same terms
#	     : as the rest of BioPerl.
#----------------------------------------------------------------------------

=head1 NAME

Bio::Tools::OddCodes - Object holding alternative alphabet coding for
one protein sequence

=head1 SYNOPSIS

  # Take a sequence object from eg, an inputstream, and creates an
  # object for the purposes of rewriting that sequence in another
  # alphabet.  These are abbreviated amino acid sequence alphabets,
  # designed to simplify the statistical aspects of analysing protein
  # sequences, by reducing the combinatorial explosion of the
  # 20-letter alphabet.  These abbreviated alphabets range in size
  # from 2 to 8.

  # Creating the OddCodes object, eg:

	my $inputstream = Bio::SeqIO->new( '-file' => "seqfile",
                                           '-format' => 'Fasta');
	my $seqobj = $inputstream->next_seq();
	my $oddcode_obj = Bio::Tools::Oddcodes->new(-seq => $seqobj);

  # or:

	my $seqobj = Bio::PrimarySeq->new
              (-seq=>'[cut and paste a sequence here]',
               -alphabet => 'protein',
               -id => 'test');
	my $oddcode_obj  =  Bio::Tools::OddCodes->new(-seq => $seqobj);

  # do the alternative coding, returning the answer as a reference to
  # a string

	my $output = $oddcode_obj->structural();
	my $output = $oddcode_obj->chemical();
	my $output = $oddcode_obj->functional();
	my $output = $oddcode_obj->charge();
	my $output = $oddcode_obj->hydrophobic();
	my $output = $oddcode_obj->Dayhoff();
	my $output = $oddcode_obj->Sneath();
	my $output = $oddcode_obj->Stanfel();


  # display sequence in new form, eg:

	my $new_coding = $$output;
	print "\n$new_coding";

=head1 DESCRIPTION

Bio::Tools::Oddcodes is a welterweight object for rewriting a protein
sequence in an alternative alphabet.  Eight of these are provided, ranging
from the the 2-letter hydrophobic alphabet, to the 8-letter chemical
alphabet.  These are useful for the statistical analysis of protein
sequences since they can partially avoid the combinatorial explosion
produced by the full 20-letter alphabet (eg. 400 dimers, 8000 trimers
etc.)

The objects will print out a warning if the input sequence is not a
protein. If you know what you are doing, you can silence the warning
by setting verbose() to a negative value.

See SYNOPSIS above for object creation code.

=head1 REFERENCES

Stanfel LE (1996) A new approach to clustering the amino acids.  J. theor.
Biol. 183, 195-205.

Karlin S, Ost F and Blaisdell BE (1989)  Patterns in DNA and amino acid
sequences and their statistical significance.  Chapter 6 of: Mathematical
Methods for DNA Sequences.  Waterman MS (ed.)  CRC Press, Boca Raton , FL.

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

Derek Gatherer

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

package Bio::Tools::OddCodes;
use strict;


use base qw(Bio::Root::Root);

sub new
{
    my($class,@args) = @_;

    my $self = $class->SUPER::new(@args);

    my ($seqobj) = $self->_rearrange([qw(SEQ)],@args);
    if((! defined($seqobj)) && @args && ref($args[0])) {
	# parameter not passed as named parameter?
	$seqobj = $args[0];
    }
    unless  ($seqobj->isa("Bio::PrimarySeqI"))
    {
        $self->throw("Bio::Tools::OddCodes only works on PrimarySeqI objects");
    }

    $self->{'_seqref'} = $seqobj;

    return $self;
}

=head2 structural

 Title   : structural
 Usage   : $output = $oddcode_obj->structural();
 Function: turns amino acid sequence into 3-letter structural alphabet
	 : A (ambivalent), E (external), I (internal)
 Example : a sequence ACDEFGH will become AAEEIAE
 Returns : Reference to the new sequence string
 Args    : none

=cut

sub structural()
{
	my $self = $_[0];
	my $seqstring = &_pullseq($self);	# see _pullseq() below

# now the real business

	$seqstring =~ tr/[ACGPSTWY]/1/;
	$seqstring =~ tr/[RNDQEHK]/2/;
	$seqstring =~ tr/[ILMFV]/3/;
	$seqstring =~ tr/1/A/;
	$seqstring =~ tr/2/E/;
	$seqstring =~ tr/3/I/;

	return \$seqstring;

# and that's that one
}

=head2 functional

 Title   : functional
 Usage   : $output = $oddcode_obj->functional();
 Function: turns amino acid sequence into 4-letter functional alphabet
	 : A (acidic), C (basic), H (hydrophobic), P (polar)
 Example : a sequence ACDEFGH will become HPAAHHC
 Returns : Reference to the new sequence string
 Args    : none

=cut

sub functional()
{
	my $self = $_[0];
	my $seqstring = &_pullseq($self);

# now the real business

	$seqstring =~ tr/[DE]/1/;
	$seqstring =~ tr/[HKR]/2/;
	$seqstring =~ tr/[AFILMPVW]/3/;
	$seqstring =~ tr/[CGNQSTY]/4/;
	$seqstring =~ tr/1/A/;
	$seqstring =~ tr/2/C/;
	$seqstring =~ tr/3/H/;
	$seqstring =~ tr/4/P/;

	return \$seqstring;

# and that's that one
}

=head2 hydrophobic

 Title   : hydrophobic
 Usage   : $output = $oddcode_obj->hydrophobic();
 Function: turns amino acid sequence into 2-letter hydrophobicity alphabet
	 : O (hydrophobic), I (hydrophilic)
 Example : a sequence ACDEFGH will become OIIIOII
 Returns : Reference to the new sequence string
 Args    : none

=cut

sub hydrophobic()
{
	my $self = $_[0];
	my $seqstring = &_pullseq($self);

# now the real business

	$seqstring =~ tr/[AFILMPVW]/1/;
	$seqstring =~ tr/[CDEGHKNQRSTY]/2/;
	$seqstring =~ tr/1/I/;
	$seqstring =~ tr/2/O/;

	return \$seqstring;

# and that's that one
}

=head2 Dayhoff

 Title   : Dayhoff
 Usage   : $output = $oddcode_obj->Dayhoff();
 Function: turns amino acid sequence into 6-letter Dayhoff alphabet
 Example : a sequence ACDEFGH will become CADDGCE
         : A (=C),   C (=AGPST), D (=DENQ),
         : E (=HKR), F (=ILMV),  G (=FWY)
 Returns : Reference to the new sequence string
 Args    : none

=cut

sub Dayhoff()
{
	my $self = $_[0];
	my $seqstring = &_pullseq($self);

# now the real business

	$seqstring =~ tr/[C]/1/;
	$seqstring =~ tr/[AGPST]/2/;
	$seqstring =~ tr/[DENQ]/3/;
	$seqstring =~ tr/[HKR]/4/;
	$seqstring =~ tr/[ILMV]/5/;
	$seqstring =~ tr/[FWY]/6/;
	$seqstring =~ tr/1/A/;
	$seqstring =~ tr/2/C/;
	$seqstring =~ tr/3/D/;
	$seqstring =~ tr/4/E/;
	$seqstring =~ tr/5/F/;
	$seqstring =~ tr/6/G/;

	return \$seqstring;

# and that's that one
}

=head2 Sneath

 Title   : Sneath
 Usage   : $output = $oddcode_obj->Sneath();
 Function: turns amino acid sequence into 7-letter Sneath alphabet
 Example : a sequence ACDEFGH will become CEFFHCF
         : A (=ILV), C (=AGP), D (=MNQ), E (=CST),
         : F (=DE),  G (=KR),  H (=FHWY)
 Returns : Reference to the new sequence string
 Args    : none

=cut

sub Sneath()
{
	my $self = $_[0];
	my $seqstring = &_pullseq($self);

# now the real business

	$seqstring =~ tr/[ILV]/1/;
	$seqstring =~ tr/[AGP]/2/;
	$seqstring =~ tr/[MNQ]/3/;
	$seqstring =~ tr/[CST]/4/;
	$seqstring =~ tr/[DE]/5/;
	$seqstring =~ tr/[KR]/6/;
	$seqstring =~ tr/[FHWY]/7/;
	$seqstring =~ tr/1/A/;
	$seqstring =~ tr/2/C/;
	$seqstring =~ tr/3/D/;
	$seqstring =~ tr/4/E/;
	$seqstring =~ tr/5/F/;
	$seqstring =~ tr/6/G/;
	$seqstring =~ tr/7/H/;

	return \$seqstring;

# and that's that one
}

=head2 Stanfel

 Title   : Stanfel
 Usage   : $output = $oddcode_obj->Stanfel();
 Function: turns amino acid sequence into 4-letter Stanfel alphabet
 Example : a sequence ACDEFGH will become AACCDAE
         : A (=ACGILMPSTV), C (=DENQ), D (=FWY), E (=HKR)
 Returns : Reference to the new sequence string
 Args    : none

=cut

sub Stanfel()
{
	my $self = $_[0];
	my $seqstring = &_pullseq($self);

# now the real business

	$seqstring =~ tr/[ACGILMPSTV]/1/;
	$seqstring =~ tr/[DENQ]/2/;
	$seqstring =~ tr/[FWY]/3/;
	$seqstring =~ tr/[HKR]/4/;
	$seqstring =~ tr/1/A/;
	$seqstring =~ tr/2/C/;
	$seqstring =~ tr/3/D/;
	$seqstring =~ tr/4/E/;

	return \$seqstring;

# and that's that one
}

=head2 chemical

 Title   : chemical
 Usage   : $output = $oddcode_obj->chemical();
 Function: turns amino acid sequence into 8-letter chemical alphabet
	 : A (acidic), L (aliphatic), M (amide), R (aromatic)
	 : C (basic),  H (hydroxyl),  I (imino), S (sulphur)
 Example : a sequence ACDEFGH will become LSAARAC
 Returns : Reference to the new sequence string
 Args    : none

=cut

sub chemical()
{
	my $self = $_[0];
	my $seqstring = &_pullseq($self);

# now the real business

	$seqstring =~ tr/[DE]/1/;
	$seqstring =~ tr/[AGILV]/2/;
	$seqstring =~ tr/[NQ]/3/;
	$seqstring =~ tr/[FWY]/4/;
	$seqstring =~ tr/[RHK]/5/;
	$seqstring =~ tr/[ST]/6/;
	$seqstring =~ tr/P/7/;
	$seqstring =~ tr/[CM]/8/;
	$seqstring =~ tr/1/A/;
	$seqstring =~ tr/2/L/;
	$seqstring =~ tr/3/M/;
	$seqstring =~ tr/4/R/;
	$seqstring =~ tr/5/C/;
	$seqstring =~ tr/6/H/;
	$seqstring =~ tr/7/I/;
	$seqstring =~ tr/8/S/;

	return \$seqstring;

# and that's that one
}

=head2 charge

 Title   : charge
 Usage   : $output = $oddcode_obj->charge();
 Function: turns amino acid sequence into 3-letter charge alphabet
 Example : a sequence ACDEFGH will become NNAANNC
         : A (negative; NOT anode), C (positive; NOT cathode), N (neutral)
 Returns : Reference to the new sequence string
 Args    : none

=cut

sub charge()
{
	my $self = $_[0];
	my $seqstring = &_pullseq($self);

# now the real business

	$seqstring =~ tr/[DE]/1/;
	$seqstring =~ tr/[HKR]/2/;
	$seqstring =~ tr/[ACFGILMNPQSTVWY]/3/;
	$seqstring =~ tr/1/A/;
	$seqstring =~ tr/2/C/;
	$seqstring =~ tr/3/N/;

	return \$seqstring;

# and that's that one
}

# _pullseq is called within each of the subroutines
# it just checks a few things and returns the sequence

sub _pullseq
{
	my $self = $_[0];

	my $seqobj =  $self->{'_seqref'};

	unless  ($seqobj->isa("Bio::PrimarySeqI"))
	{
		$self->throw("die, OddCodes works only on PrimarySeqI objects\n");
    	}
        $self->warn("\tAll OddCode alphabets need a protein sequence,\n".
                    "\tbut BioPerl thinks this is not: [". $seqobj->id. "]")
            unless $seqobj->alphabet eq 'protein' or $self->verbose < 0;;

	my $seqstring = uc $seqobj->seq();

	if(length($seqstring)<1)
	{
		$self->throw("$seqstring: die, sequence has zero length\n");
	}
	return $seqstring;
}

1;
