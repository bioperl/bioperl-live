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

Take a sequence object from eg, an inputstream, and creates an object
for the purposes of rewriting that sequence in another alphabet.
These are abbreviated amino acid sequence alphabets, designed to
simplify the statistical aspects of analysing protein sequences, by
reducing the combinatorial explosion of the 20-letter alphabet.  These
abbreviated alphabets range in size from 2 to 8.

Creating the OddCodes object, eg:

	my $inputstream = Bio::SeqIO->new( '-file' => "seqfile", 
                                           '-format' => 'Fasta');
	my $seqobj = $inputstream->next_seq();
	my $oddcode_obj = Bio::Tools::Oddcodes->new(-seq => $seqobj);

or:

	my $seqobj = Bio::PrimarySeq->new
              (-seq=>'[cut and paste a sequence here]', 
               -alphabet = 'protein', 
               -id = 'test');
	my $oddcode_obj  =  Bio::Tools::OddCodes->new(-seq => $seqobj);

do the alternative coding, returning the answer as a reference to a string

	my $output = $oddcode_obj->structural();
	my $output = $oddcode_obj->chemical();
	my $output = $oddcode_obj->functional();
	my $output = $oddcode_obj->charge();
	my $output = $oddcode_obj->hydrophobic();
	my $output = $oddcode_obj->Dayhoff();
	my $output = $oddcode_obj->Sneath();
	my $output = $oddcode_obj->Stanfel();


display sequence in new form, eg:

	my $new_coding = $$output;
	print "\n$new_coding";

=head1 DESCRIPTION

Bio::Tools::Oddcodes is a welterweight object for rewriting a protein
sequence in an alternative alphabet.  8 of these are provided, ranging
from the the 2-letter hydrophobic alphabet, to the 8-letter chemical
alphabet.  These are useful for the statistical analysis of protein
sequences since they can partially avoid the combinatorial explosion
produced by the full 20-letter alphabet (eg. 400 dimers, 8000 trimers
etc.)

See Synopsis above for object creation code.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this
and other Bioperl modules. Send your comments and suggestions preferably
to one of the Bioperl mailing lists.
Your participation is much appreciated.

  bioperl-l@bioperl.org                - General discussion
  http://www.bioperl.org/MailList.html - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via email
or the web:

  bioperl-bugs@bioperl.org
  http://www.bioperl.org/bioperl-bugs/

=head1 AUTHOR

Derek Gatherer

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

#'

package Bio::Tools::OddCodes;
use vars qw(@ISA);
use strict;

use Bio::Root::Root;

@ISA = qw(Bio::Root::Root);


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
	die("die in _init, OddCodes works only on PrimarySeqI
objects\n");
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

=head2 chemical()

 Title   : chemical
 Usage   : $output = $oddcode_obj->chemical(); 
 Function: turns amino acid sequence into 8-letter chemical alphabet
	 : A (acidic), L (aliphatic), M (amide), R (aromatic)
	 : C (basic), H (hydroxyl), I (imino), S (sulphur)
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
		die("die, OddCodes works only on PrimarySeqI objects\n");
    	}

	my $seqstring = uc $seqobj->seq();

	if(length($seqstring)<1)
	{
		die("$seqstring: die, sequence has zero length\n");
	}
	return $seqstring;
}

1;
