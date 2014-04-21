#
# BioPerl module for Bio::Variation::AAReverseMutate
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Heikki Lehvaslaiho <heikki-at-bioperl-dot-org>
#
# Copyright Heikki Lehvaslaiho
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Variation::AAReverseMutate - point mutation and codon
  information from single amino acid changes

=head1 SYNOPSIS

    $aamut = Bio::Variation::AAReverseMutate->new
                 (-aa_ori => 'F', 
		  -aa_mut => 'S',
		  -codon_ori => 'ttc', # optional
		  -codon_table => '3'  # defaults to 1
		  );

   @points = $aamut->each_Variant;

   if (scalar @points > 0 ) {
       foreach $rnachange ( @points ) {
   	   # $rnachange is a Bio::Variation::RNAChange object
   	   print " ", $rnachange->allele_ori->seq, ">", 
	   $rnachange->allele_mut->seq, " in ",
   	   $rnachange->codon_ori, ">", $rnachange->codon_mut,
   	   " at position ", $rnachange->codon_pos, "\n";
       }
   } else  {
       print "No point mutations possible\n",
   }

=head1 DESCRIPTION

Bio::Variation::AAReverseMutate objects take in reference and mutated
amino acid information and deduces potential point mutations at RNA
level leading to this change. The choice can be further limited by
letting the object know what is the the codon in the reference
sequence. The results are returned as L<Bio::Variation::RNAChange>
objects.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to the 
Bioperl mailing lists  Your participation is much appreciated.

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

=head1 AUTHOR - Heikki Lehvaslaiho

Email:  heikki-at-bioperl-dot-org

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...

package Bio::Variation::AAReverseMutate;

use strict;

# Object preamble - inheritance
use Bio::Tools::CodonTable;
use Bio::Variation::RNAChange;
use Bio::Variation::Allele;

use base qw(Bio::Root::Root);

sub new {
    my($class,@args) = @_;
    my $self = $class->SUPER::new(@args);

    my ($aa_ori, $aa_mut, $codon_ori, $codon_table) =
	$self->_rearrange([qw(AA_ORI
			      AA_MUT
			      CODON
			      CODON_TABLE
			      )],@args);

  $aa_ori && $self->aa_ori($aa_ori);
  $aa_mut  && $self->aa_mut($aa_mut);
  $codon_ori  && $self->codon_ori($codon_ori);
  $codon_table && $self->codon_table($codon_table);

  return $self; # success - we hope!

}


=head2 aa_ori

 Title   : aa_ori
 Usage   : $obj->aa_ori();
 Function: 

            Sets and returns original aa sequence.  If value is not
            set, returns false. 

            Amino acid sequences are stored in upper case characters,
            others in lower case.

 Example : 
 Returns : string
 Args    : single character amino acid code

=cut

sub aa_ori {
    my ($self,$value) = @_;
    if( defined $value) {
	if (  uc($value) !~ /^[ARNDCQEGHILKMFPSTWYVBZX*]$/ ) {
	    $self->throw("'$value' is not a valid one letter amino acid symbol\n");
	} else {
	    $self->{'aa_ori'} = uc $value;
	}
    }
    return $self->{'aa_ori'};    
}


=head2 aa_mut

 Title   : aa_mut
 Usage   : $obj->aa_mut();
 Function: 

            Sets and returns the mutated allele sequence.  If value is not
            set, returns false. 

 Example : 
 Returns : string
 Args    : single character amino acid code

=cut


sub aa_mut {
    my ($self,$value) = @_;
    if( defined $value) {
	if (  uc($value) !~ /^[ARNDCQEGHILKMFPSTWYVBZX*]$/ ) {
	    $self->throw("'$value' is not a valid one letter amino acid symbol\n");
	} else {
	    $self->{'aa_mut'} = uc $value;
	}
    }
    return $self->{'aa_mut'};
}


=head2 codon_ori

 Title   : codon_ori
 Usage   : $obj->codon_ori();
 Function: 

            Sets and returns codon_ori triplet.  If value is not set,
            returns false.  The string has to be three characters
            long. The chracter content is not checked.

 Example : 
 Returns : string
 Args    : string

=cut

sub codon_ori {
    my ($self,$value) = @_;
    if( defined $value) {
	if (length $value != 3 or lc $value =~ /[^atgc]/) {
	    $self->warn("Codon string \"$value\" is not valid unique codon");
	}
	$self->{'codon_ori'} = lc $value;
    }
    return $self->{'codon_ori'};
}

=head2 codon_table

 Title   : codon_table
 Usage   : $obj->codon_table();
 Function: 

            Sets and returns the codon table id of the RNA
            If value is not set, returns 1, 'universal' code, as the default.

 Example : 
 Returns : integer
 Args    : none if get, the new value if set

=cut


sub codon_table {
    my ($self,$value) = @_;
    if( defined $value) {
	if (  not $value =~ /^\d+$/ ) {
	    $self->throw("'$value' is not a valid codon table ID\n".
			"Has to be a positive integer. Defaulting to 1\n");
	} else {
	    $self->{'codon_table'} = $value;
	}
    }
    if( ! exists $self->{'codon_table'} ) {
	return 1;
    } else {
	return $self->{'codon_table'};
    }
}


=head2 each_Variant

 Title   : each_Variant
 Usage   : $obj->each_Variant();
 Function: 

            Returns a list of Variants.

 Example : 
 Returns : list of Variants
 Args    : none

=cut

sub each_Variant{
   my ($self,@args) = @_;
   
   $self->throw("aa_ori is not defined\n") if not defined $self->aa_ori;
   $self->throw("aa_mut is not defined\n") if not defined $self->aa_mut;

   my (@points, $codon_pos, $allele_ori, $allele_mut);
   my $ct  = Bio::Tools::CodonTable->new( '-id' => $self->codon_table );
   foreach my $codon_ori ($ct->revtranslate($self->aa_ori)) {
       next if $self->codon_ori and $self->codon_ori ne $codon_ori;
       foreach my $codon_mut ($ct->revtranslate($self->aa_mut)) {
	   my $k = 0;
	   my $length = 0; 
	   $codon_pos = $allele_ori = $allele_mut = undef;
	   while ($k<3) {
	       my $nt_ori = substr ($codon_ori, $k, 1);
	       my $nt_mut = substr ($codon_mut, $k, 1);
	       if ($nt_ori ne $nt_mut) {
		   $length++;
		   $codon_pos = $k+1;
		   $allele_ori = $nt_ori;
		   $allele_mut = $nt_mut;
	       }
	       $k++;
	   }
	   if ($length == 1) {
	       my $rna = Bio::Variation::RNAChange->new 
		   ('-length'        => '1',     
		    '-codon_ori'     => $codon_ori,
		    '-codon_mut'     => $codon_mut,
		    '-codon_pos'     => $codon_pos,
		    '-isMutation'    => 1
		    );
	       my $all_ori = Bio::Variation::Allele->new('-seq'=>$allele_ori);
	       $rna->allele_ori($all_ori);
	       my $all_mut = Bio::Variation::Allele->new('-seq'=>$allele_mut);
	       $rna->allele_mut($all_mut);
	       push @points, $rna;
	   }
       }
   }
   return @points;
}

1;
