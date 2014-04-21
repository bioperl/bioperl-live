# bioperl module for Bio::Variation::SeqDiff
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Heikki Lehvaslaiho <heikki-at-bioperl-dot-org>
#
# Copyright Heikki Lehvaslaiho
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

# cds_end definition?

=head1 NAME

Bio::Variation::SeqDiff - Container class for mutation/variant descriptions

=head1 SYNOPSIS

  $seqDiff = Bio::Variation::SeqDiff->new (
                                           -id => $M20132,
					   -alphabet => 'rna',
                                           -gene_symbol => 'AR'
                                           -chromosome => 'X',
                                           -numbering => 'coding'
                                           );
  # get a DNAMutation object somehow
  $seqDiff->add_Variant($dnamut);
  print  $seqDiff->sys_name(), "\n"; 

=head1 DESCRIPTION

SeqDiff stores Bio::Variation::VariantI object references and
descriptive information common to all changes in a sequence. Mutations
are understood to be any kind of sequence markers and are expected to
occur in the same chromosome. See L<Bio::Variation::VariantI> for details.

The methods of SeqDiff are geared towards describing mutations in
human genes using gene-based coordinate system where 'A' of the
initiator codon has number 1 and the one before it -1. This is
according to conventions of human genetics.

There will be class Bio::Variation::Genotype to describe markers in
different chromosomes and diploid genototypes.

Classes implementing Bio::Variation::VariantI interface are 
Bio::Variation::DNAMutation, Bio::Variation::RNAChange, and
Bio::Variation::AAChange. See L<Bio::Variation::VariantI>,
L<Bio::Variation::DNAMutation>, L<Bio::Variation::RNAChange>, and
L<Bio::Variation::AAChange> for more information.

Variant objects can be added using two ways: an array passed to the
constructor or as individual Variant objects with add_Variant
method.


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

=head1 CONTRIBUTORS

Eckhard Lehmann, ecky@e-lehmann.de

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::Variation::SeqDiff;

use strict;
use Bio::Tools::CodonTable;
use Bio::PrimarySeq;

use base qw(Bio::Root::Root);


=head2 new

  Title   : new
  Usage   : $seqDiff = Bio::Variation::SeqDiff->new;
  Function: generates a new Bio::Variation::SeqDiff
  Returns : reference to a new object of class SeqDiff
  Args    : 

=cut

sub new {
    my($class,@args) = @_;
    my $self = $class->SUPER::new(@args);

    my($id, $sysname, $trivname, $chr, $gene_symbol, 
       $desc, $alphabet, $numbering, $offset, $rna_offset, $rna_id, $cds_end,
       $dna_ori, $dna_mut, $rna_ori, $rna_mut, $aa_ori, $aa_mut
       #@variants, @genes
       ) =
	   $self->_rearrange([qw(ID
				 SYSNAME
				 TRIVNAME
				 CHR
				 GENE_SYMBOL
				 DESC
				 ALPHABET
				 NUMBERING
				 OFFSET
				 RNA_OFFSET
				 RNA_ID
				 CDS_END
				 DNA_ORI
				 DNA_MUT
				 RNA_ORI
				 AA_ORI
				 AA_MUT
				 )],
			    @args);
    
    #my $make = $self->SUPER::_initialize(@args);
    
    $id        && $self->id($id);           
    $sysname   && $self->sysname($sysname); 
    $trivname  && $self->trivname($trivname);
    $chr       && $self->chromosome($chr);  
    $gene_symbol && $self->gene_symbol($chr);
    $desc      && $self->description($desc);
    $alphabet   && $self->alphabet($alphabet);
    $numbering && $self->numbering($numbering);
    $offset    && $self->offset($offset);   
    $rna_offset && $self->rna_offset($rna_offset);   
    $rna_id    && $self->rna_id($rna_id);   
    $cds_end   && $self->cds_end($cds_end);   

    $dna_ori   && $self->dna_ori($dna_ori); 
    $dna_mut   && $self->dna_mut($dna_mut); 
    $rna_ori   && $self->rna_ori($rna_ori); 
    $rna_mut   && $self->rna_mut($rna_mut); 
    $aa_ori    && $self->aa_ori ($aa_ori);  
    $aa_mut    && $self->aa_mut ($aa_mut);  

    $self->{ 'variants' } = [];
    #@variants && push(@{$self->{'variants'}},@variants);

    $self->{ 'genes' } = [];
    #@genes && push(@{$self->{'genes'}},@genes);

    return $self; # success - we hope!
}


=head2 id

 Title   : id
 Usage   : $obj->id(H0001); $id = $obj->id();
 Function: 

           Sets or returns the id of the seqDiff.
           Should be used to give the collection of variants a UID
           without semantic associations.

 Example : 
 Returns : value of id, a scalar
 Args    : newvalue (optional)

=cut


sub id {
  my ($self,$value) = @_;
  if (defined $value) {
    $self->{'id'} = $value;
  }
  else {
      return $self->{'id'};
  }
}


=head2 sysname

 Title   : sysname
 Usage   : $obj->sysname('5C>G'); $sysname = $obj->sysname();
 Function: 

           Sets or returns the systematic name of the seqDiff.  The
           name should follow the HUGO Mutation Database Initiative
           approved nomenclature. If called without first setting the
           value, will generate it from L<Bio::Variation::DNAMutation>
           objects attached.

 Example : 
 Returns : value of sysname, a scalar
 Args    : newvalue (optional)

=cut


sub sysname {
    my ($self,$value) = @_;
    if (defined $value) {
	$self->{'sysname'} = $value;
    }
    elsif (not defined $self->{'sysname'}) {

	my $sysname = ''; 
	my $c = 0;
	foreach my $mut ($self->each_Variant) {
	    if( $mut->isa('Bio::Variation::DNAMutation') ) {
		$c++;
		if ($c == 1 ) {
		    $sysname = $mut->sysname ;
		}
		else {
		    $sysname .= ";". $mut->sysname;
		}
	    }
	}
	$sysname  = "[". $sysname. "]" if $c > 1;
	$self->{'sysname'} = $sysname;
    }
    return $self->{'sysname'};
}


=head2 trivname

 Title   : trivname
 Usage   : $obj->trivname('[A2G;T56G]'); $trivname = $obj->trivname();
 Function: 

           Sets or returns the trivial name of the seqDiff.
           The name should follow the HUGO Mutation Database Initiative
           approved nomenclature. If called without first setting the
           value, will generate it from L<Bio::Variation::AAChange>
           objects attached.

 Example : 
 Returns : value of trivname, a scalar
 Args    : newvalue (optional)

=cut


sub trivname {
    my ($self,$value) = @_;
    if (defined $value) {
	$self->{'trivname'} = $value;
    }
    elsif (not defined $self->{'trivname'}) {
	
	my $trivname = ''; 
	my $c = 0;
	foreach my $mut ($self->each_Variant) {
	    if( $mut->isa('Bio::Variation::AAChange') ) {
		$c++;
		if ($c == 1 ) {
		    $trivname = $mut->trivname ;
		}
		else {
		    $trivname .= ";". $mut->trivname;
		}
	    }
	}
	$trivname  = "[". $trivname. "]" if $c > 1;
	$self->{'trivname'} = $trivname;
    }

  else {
      return $self->{'trivname'};
  }
}


=head2 chromosome

 Title   : chromosome
 Usage   : $obj->chromosome('X'); $chromosome = $obj->chromosome();
 Function: 

           Sets or returns the chromosome ("linkage group") of the seqDiff.

 Example : 
 Returns : value of chromosome, a scalar
 Args    : newvalue (optional)

=cut


sub chromosome {
  my ($self,$value) = @_;
  if (defined $value) {
    $self->{'chromosome'} = $value;
  }
  else {
      return $self->{'chromosome'};
  }
}


=head2 gene_symbol

 Title   : gene_symbol
 Usage   : $obj->gene_symbol('FOS'); $gene_symbol = $obj->gene_symbol;
 Function: 

           Sets or returns the gene symbol for the studied CDS.

 Example : 
 Returns : value of gene_symbol, a scalar
 Args    : newvalue (optional)

=cut


sub gene_symbol {
  my ($self,$value) = @_;
  if (defined $value) {
    $self->{'gene_symbol'} = $value;
  }
  else {
      return $self->{'gene_symbol'};
  }
}



=head2 description

 Title   : description
 Usage   : $obj->description('short description'); $descr = $obj->description();
 Function: 

           Sets or returns the short description of the seqDiff.

 Example : 
 Returns : value of description, a scalar
 Args    : newvalue (optional)

=cut


sub description {
  my ($self,$value) = @_;
  if (defined $value) {
    $self->{'description'} = $value;
  }
  else {
      return $self->{'description'};
  }
}


=head2 alphabet

 Title   : alphabet
 Usage   : if( $obj->alphabet eq 'dna' ) { /Do Something/ }
 Function: Returns the type of primary reference sequence being one of 
           'dna', 'rna' or 'protein'. This is case sensitive.

 Returns : a string either 'dna','rna','protein'. 
 Args    : none


=cut

sub alphabet {
   my ($self,$value) = @_;
   my %type = (dna => 1,
	       rna => 1,
	       protein => 1);
   if( defined $value ) {
       if ($type{$value}) {
	   $self->{'alphabet'} = $value;
       } else {
	   $self->throw("$value is not valid alphabet value!");
       }
   }
   return $self->{'alphabet'};
}


=head2 numbering

 Title   : numbering
 Usage   : $obj->numbering('coding'); $numbering = $obj->numbering();
 Function: 

           Sets or returns the string giving the numbering schema used
           to describe the variants.

 Example : 
 Returns : value of numbering, a scalar
 Args    : newvalue (optional)

=cut



sub numbering {
  my ($self,$value) = @_;
  if (defined $value) {
    $self->{'numbering'} = $value;
  }
  else {
      return $self->{'numbering'};
  }
}


=head2 offset

 Title   : offset
 Usage   : $obj->offset(124); $offset = $obj->offset();
 Function: 

           Sets or returns the offset from the beginning of the DNA sequence 
           to the coordinate start used to describe variants. Typically
           the beginning of the coding region of the gene. 
           The cds_start should be 1 + offset.

 Example : 
 Returns : value of offset, a scalar
 Args    : newvalue (optional)

=cut



sub offset {
  my ($self,$value) = @_;
  if (defined $value) {
    $self->{'offset'} = $value;
  }
  elsif (not defined $self->{'offset'} ) {
      return $self->{'offset'} = 0;
  }
  else {
      return $self->{'offset'};
  }
}


=head2 cds_start

 Title   : cds_start
 Usage   : $obj->cds_start(123); $cds_start = $obj->cds_start();
 Function: 

           Sets or returns the cds_start from the beginning of the DNA
           sequence to the coordinate start used to describe
           variants. Typically the beginning of the coding region of
           the gene. Needs to be and is implemented as 1 + offset.

 Example : 
 Returns : value of cds_start, a scalar
 Args    : newvalue (optional)

=cut



sub cds_start {
  my ($self,$value) = @_;
  if (defined $value) {
    $self->{'offset'} = $value - 1;
  }
  else {
      return $self->{'offset'} + 1;
  }
}


=head2 cds_end

 Title   : cds_end
 Usage   : $obj->cds_end(321); $cds_end = $obj->cds_end();
 Function: 

           Sets or returns the position of the last nucleotitide of the
           termination codon. The coordinate system starts from cds_start.

 Example : 
 Returns : value of cds_end, a scalar
 Args    : newvalue (optional)

=cut



sub cds_end {
  my ($self,$value) = @_;
  if (defined $value) {
    $self->{'cds_end'} = $value;
  }
  else {
      return $self->{'cds_end'};
      #$self->{'cds_end'} = CORE::length($self->SeqDiff->rna_ori)/3;
  }
}


=head2 rna_offset

 Title   : rna_offset
 Usage   : $obj->rna_offset(124); $rna_offset = $obj->rna_offset();
 Function: 

           Sets or returns the rna_offset from the beginning of the RNA sequence 
           to the coordinate start used to describe variants. Typically
           the beginning of the coding region of the gene. 

 Example : 
 Returns : value of rna_offset, a scalar
 Args    : newvalue (optional)

=cut



sub rna_offset {
  my ($self,$value) = @_;
  if (defined $value) {
    $self->{'rna_offset'} = $value;
  }
  elsif (not defined $self->{'rna_offset'} ) {
      return $self->{'rna_offset'} = 0;
  }
  else {
      return $self->{'rna_offset'};
  }
}


=head2 rna_id

 Title   : rna_id
 Usage   : $obj->rna_id('transcript#3'); $rna_id = $obj->rna_id();
 Function: 

	    Sets or returns the ID for original RNA sequence of the seqDiff.

 Example : 
 Returns : value of rna_id, a scalar
 Args    : newvalue (optional)

=cut


sub rna_id {
  my ($self,$value) = @_;
  if (defined $value) {
    $self->{'rna_id'} = $value;
  }
  else {
      return $self->{'rna_id'};
  }
}



=head2 add_Variant

 Title   : add_Variant
 Usage   : $obj->add_Variant($variant)
 Function: 

           Pushes one Bio::Variation::Variant into the list of variants.
           At the same time, creates a link from the Variant to SeqDiff
           using its SeqDiff method.

 Example : 
 Returns : 1 when succeeds, 0 for failure.
 Args    : Variant object

=cut

sub add_Variant {
  my ($self,$value) = @_;
  if (defined $value) {
      if( ! $value->isa('Bio::Variation::VariantI') ) {
	  $self->throw("Is not a VariantI complying  object but a [$self]");
	  return 0;
      }
      else {
	  push(@{$self->{'variants'}},$value);
	  $value->SeqDiff($self);
	  return 1;
      }
  }
  else {
      return 0;
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
   
   return @{$self->{'variants'}}; 
}



=head2 add_Gene

 Title   : add_Gene
 Usage   : $obj->add_Gene($gene)
 Function: 

           Pushes one L<Bio::LiveSeq::Gene> into the list of genes.

 Example : 
 Returns : 1 when succeeds, 0 for failure.
 Args    : Bio::LiveSeq::Gene object

See L<Bio::LiveSeq::Gene> for more information.

=cut


sub add_Gene {
  my ($self,$value) = @_;
  if (defined $value) {
      if( ! $value->isa('Bio::LiveSeq::Gene') ) {
	  $value->throw("Is not a Bio::LiveSeq::Gene object but a  [$value]");
	  return 0;
      }
      else {
	  push(@{$self->{'genes'}},$value);
	  return 1;
      }
  }
  else {
      return 0;
  }
}


=head2 each_Gene

 Title   : each_Gene
 Usage   : $obj->each_Gene();
 Function: 

            Returns a list of L<Bio::LiveSeq::Gene>s.

 Example : 
 Returns : list of Genes
 Args    : none

=cut

sub each_Gene{
   my ($self,@args) = @_;

   return @{$self->{'genes'}}; 
}


=head2 dna_ori

 Title   : dna_ori
 Usage   : $obj->dna_ori('atgctgctgctgct'); $dna_ori = $obj->dna_ori();
 Function: 

	    Sets or returns the original DNA sequence string of the seqDiff.

 Example : 
 Returns : value of dna_ori, a scalar
 Args    : newvalue (optional)

=cut


sub dna_ori {
  my ($self,$value) = @_;
  if (defined $value) {
      $self->{'dna_ori'} = $value;
  }
  else {
      return $self->{'dna_ori'};
  }
}


=head2 dna_mut

 Title   : dna_mut
 Usage   : $obj->dna_mut('atgctggtgctgct'); $dna_mut = $obj->dna_mut();
 Function: 

	    Sets or returns the mutated DNA sequence of the seqDiff.
            If sequence has not been set generates it from the
            original sequence and DNA mutations.

 Example : 
 Returns : value of dna_mut, a scalar
 Args    : newvalue (optional)

=cut


sub dna_mut {
  my ($self,$value) = @_;
  if (defined $value) {
      $self->{'dna_mut'} = $value;
  }
  else {
      $self->_set_dnamut() unless $self->{'dna_mut'};
      return $self->{'dna_mut'};
  }
}

sub _set_dnamut {
    my $self = shift;

    return unless $self->{'dna_ori'}  && $self->each_Variant;

    $self->{'dna_mut'} = $self->{'dna_ori'};
    foreach ($self->each_Variant) {
	next unless $_->isa('Bio::Variation::DNAMutation');
	next unless $_->isMutation;

	my ($s, $la, $le);
	#lies the mutation less than 25 bases after the start of sequence?
	if ($_->start < 25) {
	    $s = 0; $la = $_->start - 1;
	} else {
	    $s = $_->start - 25; $la = 25;
	}

	#is the mutation an insertion?
	$_->end($_->start) unless $_->allele_ori->seq;

	#does the mutation end greater than 25 bases before the end of
	#sequence?
	if (($_->end + 25) > length($self->{'dna_mut'})) {
	    $le = length($self->{'dna_mut'}) - $_->end;
	} else {
	    $le = 25;
	}

	$_->dnStreamSeq(substr($self->{'dna_mut'}, $s, $la));
	$_->upStreamSeq(substr($self->{'dna_mut'}, $_->end, $le));

	my $s_ori = $_->dnStreamSeq . $_->allele_ori->seq . $_->upStreamSeq;
	my $s_mut = $_->dnStreamSeq . $_->allele_mut->seq . $_->upStreamSeq;

	(my $str = $self->{'dna_mut'}) =~ s/$s_ori/$s_mut/;
	$self->{'dna_mut'} = $str;
    }
}


=head2 rna_ori

 Title   : rna_ori
 Usage   : $obj->rna_ori('atgctgctgctgct'); $rna_ori = $obj->rna_ori();
 Function: 

	    Sets or returns the original RNA sequence of the seqDiff.

 Example : 
 Returns : value of rna_ori, a scalar
 Args    : newvalue (optional)

=cut


sub rna_ori {
  my ($self,$value) = @_;
  if (defined $value) {
    $self->{'rna_ori'} = $value;
  }
  else {
      return $self->{'rna_ori'};
  }
}


=head2 rna_mut

 Title   : rna_mut
 Usage   : $obj->rna_mut('atgctggtgctgct'); $rna_mut = $obj->rna_mut();
 Function: 

	    Sets or returns the mutated RNA sequence of the seqDiff.

 Example : 
 Returns : value of rna_mut, a scalar
 Args    : newvalue (optional)

=cut


sub rna_mut {
  my ($self,$value) = @_;
  if (defined $value) {
    $self->{'rna_mut'} = $value;
  }
  else {
      return $self->{'rna_mut'};
  }
}


=head2 aa_ori

 Title   : aa_ori
 Usage   : $obj->aa_ori('MAGVLL*'); $aa_ori = $obj->aa_ori();
 Function: 

	    Sets or returns the original protein sequence of the seqDiff.

 Example : 
 Returns : value of aa_ori, a scalar
 Args    : newvalue (optional)

=cut


sub aa_ori {
  my ($self,$value) = @_;
  if (defined $value) {
    $self->{'aa_ori'} = $value;
  }
  else {
      return $self->{'aa_ori'};
  }
}


=head2 aa_mut

 Title   : aa_mut
 Usage   : $obj->aa_mut('MA*'); $aa_mut = $obj->aa_mut();
 Function: 

	    Sets or returns the mutated protein sequence of the seqDiff.

 Example : 
 Returns : value of aa_mut, a scalar
 Args    : newvalue (optional)

=cut


sub aa_mut {
  my ($self,$value) = @_;
  if (defined $value) {
    $self->{'aa_mut'} = $value;
  }
  else {
      return $self->{'aa_mut'};
  }
}


=head2 seqobj

 Title   : seqobj
 Usage   : $dnaobj = $obj->seqobj('dna_mut');
 Function: 

	    Returns the any original or mutated sequences as a
	    Bio::PrimarySeq object.

 Example : 
 Returns : Bio::PrimarySeq object for the requested sequence
 Args    : string, method name for the sequence requested

See L<Bio::PrimarySeq> for more information.

=cut

sub seqobj {
  my ($self,$value) = @_;
  my $out;
  my %valid_obj = 
      map {$_, 1} qw(dna_ori rna_ori aa_ori dna_mut rna_mut aa_mut);
  $valid_obj{$value} ||
      $self->throw("Sequence type '$value' is not a valid type (".
                  join(',', map "'$_'", sort keys %valid_obj) .") lowercase");
  my ($alphabet) = $value =~ /([^_]+)/;
  my $id =  $self->id;
  $id =  $self->rna_id if $self->rna_id;
  $alphabet = 'protein' if $alphabet eq 'aa';
  $out = Bio::PrimarySeq->new
      ( '-seq' => $self->{$value},
	'-display_id'  => $id,
	'-accession_number' => $self->id,
	'-alphabet' => $alphabet
	) if   $self->{$value} ;
  return $out;
}

=head2 alignment

 Title   : alignment
 Usage   : $obj->alignment
 Function: 

           Returns a pretty RNA/AA sequence alignment from linked
           objects.  Under construction: Only simple coding region
           point mutations work.

 Example : 
 Returns : 
 Args    : none

=cut


sub alignment {
    my $self = shift;
    my (@entry, $text);

    my $maxflanklen = 12;

    foreach my $mut ($self->each_Variant) {
	if( $mut->isa('Bio::Variation::RNAChange') ) {

	    my $upflank = $mut->upStreamSeq;
	    my $dnflank = $mut->dnStreamSeq;
	    my $cposd = $mut->codon_pos;
	    my $rori = $mut->allele_ori->seq;
	    my $rmut =  $mut->allele_mut->seq;
	    my $rseqoriu = '';
	    my $rseqmutu = '';
	    my $rseqorid = '';
	    my $rseqmutd = '';
	    my $aaseqmutu = '';
	    my (@rseqori, @rseqmut );

	    #  point
	    if ($mut->DNAMutation->label =~ /point/) {
		if ($cposd == 1 ) {
		    my $nt2d = substr($dnflank, 0, 2);
		    push @rseqori, $rori. $nt2d;
		    push @rseqmut, uc ($rmut). $nt2d;
		    $dnflank = substr($dnflank, 2);
		}
		elsif ($cposd == 2) { 
		    my $ntu = chop $upflank;
		    my $ntd = substr($dnflank, 0, 1);
		    push @rseqori, $ntu. $rori. $ntd;
		    push @rseqmut,  $ntu. uc ($rmut). $ntd;
		    $dnflank =  substr($dnflank, 1);
		}
		elsif ($cposd == 3) {
		    my $ntu1 = chop $upflank;
		    my $ntu2 = chop $upflank;
		    push (@rseqori, $ntu2. $ntu1. $rori);
		    push (@rseqmut, $ntu2. $ntu1. uc $rmut);
		}		
	    }
	    #deletion
	    elsif ($mut->DNAMutation->label =~ /deletion/) {
		if ($cposd == 2 ) {
		    $rseqorid = chop $upflank;
		    $rseqmutd = $rseqorid;
		}
		for (my $i=1; $i<=$mut->length; $i++) {
		    my $ntd .= substr($mut->allele_ori, $i-1, 1);
		    $rseqorid .= $ntd;
		    if  (length($rseqorid) == 3 ) {
			push (@rseqori, $rseqorid);
			push (@rseqmut, "   ");
			$rseqorid = '';
		    }		    
		}

		if ($rseqorid) {
		    $rseqorid .= substr($dnflank, 0, 3-$rseqorid);
		    push (@rseqori, $rseqorid);
		    push (@rseqmut, "   ");
		    $dnflank = substr($dnflank,3-$rseqorid);
		} 
	    }
	    $upflank = reverse $upflank;
	    # loop throught the flanks
	    for (my $i=1; $i<=length($dnflank); $i++) {
		
		last if  $i > $maxflanklen;

		my $ntd .= substr($dnflank, $i-1, 1);
		my $ntu .= substr($upflank, $i-1, 1);

		$rseqmutd .= $ntd;
		$rseqorid .= $ntd;
		$rseqmutu = $ntu. $rseqmutu;
		$rseqoriu = $ntu. $rseqoriu;
		
		if  (length($rseqorid) == 3  and length($rseqorid) == 3) {
		    push (@rseqori, $rseqorid);
		    push (@rseqmut, $rseqmutd);
		    $rseqorid =  $rseqmutd ='';
		}
		if  (length($rseqoriu) == 3  and length($rseqoriu) == 3) {
		    unshift (@rseqori, $rseqoriu);
		    unshift (@rseqmut, $rseqmutu);
		    $rseqoriu =  $rseqmutu ='';
		}

		#print "|i=$i,  $cposd, $rseqmutd, $rseqorid\n";
		#print "|i=$i,  $cposu, $rseqmutu, $rseqoriu\n\n";

	    }

	    push (@rseqori, $rseqorid);
	    unshift (@rseqori, $rseqoriu);
	    push (@rseqmut, $rseqmutd);
	    unshift (@rseqmut, $rseqmutu);
	    
	    return unless $mut->AAChange;
	    #translate
	    my $tr = Bio::Tools::CodonTable->new('-id' => $mut->codon_table);
	    my $apos =  $mut->AAChange->start;
	    my $aposmax = CORE::length($self->aa_ori); #terminator codon no 
	    my $rseqori;
	    my $rseqmut;
	    my $aaseqori;
	    my $aaseqmut = "";
	    for (my $i = 0; $i <= $#rseqori; $i++) {
		 my $a = '';

		 $a =  $tr->translate($rseqori[$i]) if length($rseqori[$i]) == 3;
		 
		 if (length($a) != 1 or 
		     $apos - ( $maxflanklen/2 -1) + $i < 1 or 
		     $apos - ( $maxflanklen/2 -1) + $i > $aposmax ) {
		     $aaseqori .= "    ";
		 } else {
		     $aaseqori .= " ". $a. "  ";
		 }
		 my $b = '';
		 if (length($rseqmut[$i]) == 3) {
		     if ($rseqmut[$i] eq '   ') {
			 $b = "_";
		     } else {
			 $b = $tr->translate($rseqmut[$i]);
		     }
		 }
		 if (( $b ne $a and
		       length($b) == 1 and 
		       $apos - ( $maxflanklen/2 -1) + $i >= 1 ) or
		     ( $apos - ( $maxflanklen/2 -1) + $i >= $aposmax and 
		       $mut->label =~ 'termination')
		     ) {
		     $aaseqmut .= " ". $b. "  ";
		 } else {
		     $aaseqmut .= "    ";
		 }
		 
		 if ($i == 0 and length($rseqori[$i]) != 3) {
		     my $l = 3 - length($rseqori[$i]);
		     $rseqori[$i] = (" " x $l). $rseqori[$i];
		     $rseqmut[$i] = (" " x $l). $rseqmut[$i];
		 }
		 $rseqori .= $rseqori[$i]. " " if $rseqori[$i] ne '';
		 $rseqmut .= $rseqmut[$i]. " " if $rseqmut[$i] ne '';
	     }
	    
	    # collect the results
	    push (@entry, 
		  "\n"
		  );   	    
	    $text = "           ". $aaseqmut; 
	    push (@entry, 
		  $text
		  );   	    
	    $text = "Variant  : ". $rseqmut;
	    push (@entry, 
		  $text
		  );   	    
	    $text = "Reference: ". $rseqori;
	    push (@entry, 
		  $text
		  );   	    
	    $text = "           ". $aaseqori;
	    push (@entry, 
		  $text
		  );   
	    push (@entry, 
		  "\n"
		  );   
	}

    }

    my $res;
    foreach my $line (@entry) {
       $res .=  "$line\n";
    }
    return $res;
}

1;
