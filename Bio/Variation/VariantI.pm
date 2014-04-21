#
# BioPerl module for Bio::Variation::VariantI
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

Bio::Variation::VariantI - Sequence Change SeqFeature abstract class

=head1 SYNOPSIS

  #get Bio::Variant::VariantI somehow
  print $var->restriction_changes, "\n";
  foreach $allele ($var->each_Allele) {
      #work on Bio::Variation::Allele objects
  }

=head1 DESCRIPTION

This superclass defines common methods to basic sequence changes.  The
instantiable classes Bio::Variation::DNAMutation,
Bio::Variation::RNAChange and Bio::Variation::AAChange use them.
See L<Bio::Variation::DNAMutation>, L<Bio::Variation::RNAChange>,
and L<Bio::Variation::AAChange> for more information.

These classes store information, heavy computation to detemine allele
sequences is done elsewhere.

The database cross-references are implemented as
Bio::Annotation::DBLink objects. The methods to access them are
defined in Bio::DBLinkContainerI. See L<Bio::Annotation::DBLink>
and L<Bio::DBLinkContainerI> for details.

Bio::Variation::VariantI redifines and extends
Bio::SeqFeature::Generic for sequence variations. This class
describes specific sequence change events. These events are always
from a specific reference sequence to something different. See
L<Bio::SeqFeature::Generic> for more information.

IMPORTANT: The notion of reference sequence permeates all
Bio::Variation classes. This is especially important to remember when
dealing with Alleles. In a polymorphic site, there can be a large
number of alleles. One of then has to be selected to be the reference
allele (allele_ori). ALL the rest has to be passed to the Variant
using the method add_Allele, including the mutated allele in a
canonical mutation. The IO modules and generated attributes depend on
it. They ignore the allele linked to using allele_mut and circulate
each Allele returned by each_Allele into allele_mut and calculate
the changes between that and allele_ori.


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


package Bio::Variation::VariantI;
use strict;
# Object preamble - inheritance

use base qw(Bio::Root::Root Bio::SeqFeature::Generic Bio::DBLinkContainerI);

=head2 id

 Title   : id
 Usage   : $obj->id
 Function:

           Read only method. Returns the id of the variation object.
           The id is the id of the first DBLink object attached to this object.

 Example :
 Returns : scalar
 Args    : none

=cut

sub id {
   my ($self) = @_;
   my @ids = $self->each_DBLink;
   my $id = $ids[0] if scalar @ids > 0;
   return $id->database. "::". $id->primary_id if $id;
}


=head2 add_Allele

 Title   : add_Allele
 Usage   : $self->add_Allele($allele)
 Function: 

	    Adds one Bio::Variation::Allele into the list of alleles.
            Note that the method forces the convention that nucleotide
            sequence is in lower case and amino acds are in upper
            case.

 Example : 
 Returns : 1 when succeeds, 0 for failure.
 Args    : Allele object

=cut


sub add_Allele {
  my ($self,$value) = @_;
  if (defined $value) {
      if( ! $value->isa('Bio::Variation::Allele') ) {
	  my $com = ref $value;
	  $self->throw("Is not a Allele object but a  [$com]");
	  return 0;
      } else {
	  if ( $self->isa('Bio::Variation::AAChange') ) {
	      $value->seq( uc $value->seq) if $value->seq;
	  } else {
	      $value->seq( lc $value->seq) if $value->seq;
	  } 
	  push(@{$self->{'alleles'}},$value); 
	  $self->allele_mut($value); #????
	  return 1;
      }
  } else {
      return 0;
  }
}


=head2 each_Allele

 Title   : alleles
 Usage   : $obj->each_Allele();
 Function: 

	     Returns a list of Bio::Variation::Allele objects

 Example : 
 Returns : list of Alleles
 Args    : none

=cut

sub each_Allele{
   my ($self,@args) = @_;
   return @{$self->{'alleles'}};
}


=head2 isMutation

 Title   : isMutation
 Usage   : print join('/', $obj->each_Allele) if not $obj->isMutation;
 Function:

           Returns or sets the boolean value indicating that the
           variant descibed is a canonical mutation with two alleles
           assinged to be the original (wild type) allele and mutated
           allele, respectively. If this value is not set, it is
           assumed that the Variant descibes polymorphisms.

 Returns : a boolean

=cut

sub isMutation {
    my ($self,$value) = @_;
    if (defined $value) {
        if ($value ) {
            $self->{'isMutation'} = 1;
        } else {
            $self->{'isMutation'} = 0;
        }
    }
    return $self->{'isMutation'};
} 


=head2 allele_ori

 Title   : allele_ori
 Usage   : $obj->allele_ori();
 Function: 

            Links to and returns the Bio::Variation::Allele object.
            If value is not set, returns false. All other Alleles are
            compared to this.

            Amino acid sequences are stored in upper case characters,
            others in lower case.

 Example : 
 Returns : string
 Args    : string

See L<Bio::Variation::Allele> for more.

=cut

sub allele_ori {
   my ($self,$value) = @_;
   if( defined $value) {
       if ( ! ref $value || ! $value->isa('Bio::Variation::Allele')) {
	   $self->throw("Value is not Bio::Variation::Allele but [$value]");
       } else {
	   if ( $self->isa('Bio::Variation::AAChange') ) {
	       $value->seq( uc $value->seq) if $value->seq;
	   } else {
	       $value->seq( lc $value->seq) if $value->seq;
	   } 
	   $self->{'allele_ori'} = $value;
       }
   }
   return $self->{'allele_ori'};
}


=head2 allele_mut

 Title   : allele_mut
 Usage   : $obj->allele_mut();
 Function: 

             Links to and returns the Bio::Variation::Allele
             object.  Sets and returns the mutated allele sequence.
             If value is not set, returns false.

             Amino acid sequences are stored in upper case characters,
             others in lower case.

 Example : 
 Returns : string
 Args    : string

See L<Bio::Variation::Allele> for more.

=cut


sub allele_mut {
   my ($self,$value) = @_;
   if( defined $value) {
       if ( ! ref $value || ! $value->isa('Bio::Variation::Allele')) {
	   $self->throw("Value is not Bio::Variation::Allele but [$value]");
       } else {
	   if ( $self->isa('Bio::Variation::AAChange') ) {
	       $value->seq( uc $value->seq) if $value->seq;
	   } else {
	       $value->seq( lc $value->seq) if $value->seq;
	   } 
	   $self->{'allele_mut'} = $value;
       }
   }
   return $self->{'allele_mut'};
}

=head2 length

 Title   : length
 Usage   : $obj->length();
 Function: 

            Sets and returns the length of the affected original
            allele sequence.  If value is not set, returns false == 0.

            Value 0 means that the variant position is before the
            start=end sequence position. (Value 1 would denote a point
            mutation). This follows the convension to report an
            insertion (2insT) in equivalent way to a corresponding
            deletion (2delT) (Think about indel polymorpism ATC <=> AC
            where the origianal state is not known ).

 Example : 
 Returns : string
 Args    : string

=cut


sub length {
   my ($self,$value) = @_;
   if ( defined $value) {
       $self->{'length'} = $value;
  }
   if ( ! exists $self->{'length'} ) {
       return 0;
   } 
   return $self->{'length'};
}

=head2 upStreamSeq

 Title   : upStreamSeq
 Usage   : $obj->upStreamSeq();
 Function: 

            Sets and returns upstream flanking sequence string.  If
            value is not set, returns false. The sequence should be
            >=25 characters long, if possible.

 Example : 
 Returns : string or false
 Args    : string

=cut


sub upStreamSeq {
    my ($self,$value) = @_;
    if( defined $value) {
	$self->{'upstreamseq'} = $value;
    }
   return $self->{'upstreamseq'};
}


=head2 dnStreamSeq

 Title   : dnStreamSeq
 Usage   : $obj->dnStreamSeq();
 Function: 

            Sets and returns dnstream flanking sequence string.  If
            value is not set, returns false. The sequence should be
            >=25 characters long, if possible.

 Example : 
 Returns : string or false
 Args    : string

=cut


sub dnStreamSeq {
    my ($self,$value) = @_;
    if( defined $value) {
	$self->{'dnstreamseq'} = $value;
    }
    return $self->{'dnstreamseq'};
    
}


=head2 label

 Title   : label
 Usage   : $obj->label();
 Function: 

            Sets and returns mutation event label(s).  If value is not
            set, or no argument is given returns false.  Each
            instantiable class needs to implement this method. Valid
            values are listed in 'Mutation event controlled vocabulary' in
            http://www.ebi.ac.uk/mutations/recommendations/mutevent.html.

 Example : 
 Returns : string
 Args    : string

=cut


sub label {
    my ($self,$value) = @_;
    $self->throw_not_implemented();
}



=head2 status

 Title   : status
 Usage   : $obj->status()
 Function: 

           Returns the status of the sequence change object.
           Valid values are: 'suspected' and 'proven'

 Example : $obj->status('proven');
 Returns : scalar
 Args    : valid string (optional, for setting)


=cut


sub status {
   my ($self,$value) = @_;
   my %status = (suspected => 1,
		 proven => 1
		 );

   if( defined $value) {
       $value = lc $value;
       if ($status{$value}) {
	   $self->{'status'} = $value;
       } 
       else {
	   $self->throw("$value is not valid status value!");
       }
    }
   if( ! exists $self->{'status'} ) {
       return "$self";
   }
   return $self->{'status'};
}


=head2 proof

 Title   : proof
 Usage   : $obj->proof()
 Function: 

           Returns the proof of the sequence change object.
           Valid values are: 'computed' and 'experimental'.

 Example : $obj->proof('computed');
 Returns : scalar
 Args    : valid string (optional, for setting)


=cut


sub proof {
    my ($self,$value) = @_;
    my %proof = (computed => 1,
		 experimental => 1
		 );

    if( defined $value) {
	$value = lc $value;
	if ($proof{$value}) {
	    $self->{'proof'} = $value;
	} else {
	    $self->throw("$value is not valid proof value!");
	}
    }
    return $self->{'proof'};
}


=head2 region

 Title   : region
 Usage   : $obj->region();
 Function: 

            Sets and returns the name of the sequence region type or
            protein domain at this location.  If value is not set,
            returns false.

 Example : 
 Returns : string
 Args    : string

=cut


sub region {
    my ($self,$value) = @_;
    if( defined $value) {
	$self->{'region'} = $value;
    }
    return $self->{'region'};
}


=head2 region_value

 Title   : region_value
 Usage   : $obj->region_value();
 Function: 

            Sets and returns the name of the sequence region_value or
            protein domain at this location.  If value is not set,
            returns false.

 Example : 
 Returns : string
 Args    : string

=cut


sub region_value {
    my ($self,$value) = @_;
    if( defined $value) {
	$self->{'region_value'} = $value;
    }
    return $self->{'region_value'};
}

=head2 region_dist

 Title   : region_dist
 Usage   : $obj->region_dist();
 Function: 

            Sets and returns the distance tot the closest region
            (i.e. intro/exon or domain) boundary. If distance is not
            set, returns false.

 Example : 
 Returns : integer
 Args    : integer

=cut


sub region_dist {
    my ($self,$value) = @_;
    if( defined $value) {
       if (  not $value =~ /^[+-]?\d+$/ ) {
	   $self->throw("[$value] for region_dist has to be an integer\n");
        } else {
	    $self->{'region_dist'} = $value;
        }
    }
    return $self->{'region_dist'};
}


=head2 numbering

 Title   : numbering
 Usage   : $obj->numbering()
 Function: 

           Returns the numbering chema used locating sequnce features.
           Valid values are: 'entry' and 'coding'

 Example : $obj->numbering('coding');
 Returns : scalar
 Args    : valid string (optional, for setting)


=cut


sub numbering {
   my ($self,$value) = @_;
   my %numbering = (entry => 1,
		    coding => 1
		    );

   if( defined $value) {
       $value = lc $value;
       if ($numbering{$value}) {
	   $self->{'numbering'} = $value;
       } 
       else {
	   $self->throw("'$value' is not a valid for numbering!");
       }
    }
   if( ! exists $self->{'numbering'} ) {
       return "$self";
   }
   return $self->{'numbering'};
}

=head2 mut_number

 Title   : mut_number
 Usage   : $num = $obj->mut_number;
         : $num = $obj->mut_number($number);
 Function: 

           Returns or sets the number identifying the order in which the
           mutation has been issued. Numbers shouldstart from 1.
           If the number has never been set, the method will return ''

           If you want the output from IO modules look nice and, for
           multivariant/allele variations, make sense you better set
           this attribute.

 Returns : an integer

=cut


sub mut_number {
    my ($self,$value) = @_;
    if (defined $value) {
	$self->{'mut_number'} = $value;
    }
    unless (exists $self->{'mut_number'}) {
	return ('');
    } else {
	return $self->{'mut_number'};
    }
}       


=head2 SeqDiff

 Title   : SeqDiff
 Usage   : $mutobj = $obj->SeqDiff;
         : $mutobj = $obj->SeqDiff($objref);
 Function: 

           Returns or sets the link-reference to the umbrella
           Bio::Variation::SeqDiff object.  If there is no link,
           it will return undef

           Note: Adding a variant into a SeqDiff object will
           automatically set this value.

 Returns : an obj_ref or undef

See L<Bio::Variation::SeqDiff> for more information.

=cut

sub SeqDiff {
    my ($self,$value) = @_;
    if (defined $value) {
	if( ! $value->isa('Bio::Variation::SeqDiff') ) {
	    $self->throw("Is not a Bio::Variation::SeqDiff object but a [$value]");
	    return;
	}
	else {
	    $self->{'seqDiff'} = $value;
	}
    }
    unless (exists $self->{'seqDiff'}) {
	return;
    } else {
	return $self->{'seqDiff'};
    }
}

=head2 add_DBLink

 Title   : add_DBLink
 Usage   : $self->add_DBLink($ref)
 Function: adds a link object
 Example :
 Returns : 
 Args    :


=cut


sub add_DBLink{
   my ($self,$com) = @_;
   if( $com && ! $com->isa('Bio::Annotation::DBLink') ) {
       $self->throw("Is not a link object but a  [$com]");
   }
   $com && push(@{$self->{'link'}},$com);
}

=head2 each_DBLink

 Title   : each_DBLink
 Usage   : foreach $ref ( $self->each_DBlink() )
 Function: gets an array of DBlink of objects
 Example :
 Returns : 
 Args    :


=cut

sub each_DBLink{
   my ($self) = @_;
   
   return @{$self->{'link'}}; 
}

=head2 restriction_changes

 Title   : restriction_changes
 Usage   : $obj->restriction_changes();
 Function: 

            Returns a string containing a list of restriction
            enzyme changes of form +EcoRI, separated by
            commas. Strings need to be valid restriction enzyme names
            as stored in REBASE. allele_ori and allele_mut need to be assigned.

 Example : 
 Returns : string
 Args    : string

=cut

sub restriction_changes { 
    my ($self) = @_;

    if (not $self->{'re_changes'}) { 
	my %re = &_enzymes;
	
	# complain if used on AA data
	if ($self->isa('Bio::Variation::AAChange')) {
	    $self->throw('Restriction enzymes do not bite polypeptides!');
	}
	
	#sanity checks
	$self->warn('Upstream sequence is empty!')
	    if $self->upStreamSeq eq '';
	$self->warn('Downstream sequence is empty!')
	    if $self->dnStreamSeq eq '';
#	 $self->warn('Original allele sequence is empty!')
#	     if $self->allele_ori eq '';
#	 $self->warn('Mutated allele sequence is empty!')
#	     if $self->allele_mut eq '';
	
	#reuse the non empty DNA level list at RNA level if the flanks are identical
	#Hint: Check DNAMutation object first
	if ($self->isa('Bio::Variation::RNAChange') and  $self->DNAMutation and
	    $self->upStreamSeq eq $self->DNAMutation->upStreamSeq  and 
	    $self->dnStreamSeq eq $self->DNAMutation->dnStreamSeq and
	    $self->DNAMutation->restriction_changes ne '' ) {
	    $self->{'re_changes'} = $self->DNAMutation->restriction_changes;
	} else {
	    
	    #maximum length of a type II restriction site in the current REBASE
	    my ($le_dn) = 15; 
	    my ($le_up) = $le_dn;

	    #reduce the flank lengths if the desired length is not available
	    $le_dn = CORE::length ($self->dnStreamSeq) if $le_dn > CORE::length ($self->dnStreamSeq);
	    $le_up = CORE::length ($self->upStreamSeq) if $le_up > CORE::length ($self->upStreamSeq);

	    #Build sequence strings to compare
	    my ($oriseq, $mutseq);    
	    $oriseq  = $mutseq = substr($self->upStreamSeq, -$le_up, $le_up);
	    $oriseq .= $self->allele_ori->seq if $self->allele_ori->seq;
	    $mutseq .= $self->allele_mut->seq if $self->allele_mut->seq;
	    $oriseq .= substr($self->dnStreamSeq, 0, $le_dn);
	    $mutseq .= substr($self->dnStreamSeq, 0, $le_dn);
	    
	    # ... and their reverse complements
	    my $oriseq_rev = _revcompl ($oriseq);
	    my $mutseq_rev = _revcompl ($mutseq);

	    # collect results into a string
	    my $rec = '';
	    foreach my $enz (sort keys (%re)) {
		my $site = $re{$enz};
		my @ori = ($oriseq=~ /$site/g);
		my @mut = ($mutseq=~ /$site/g);
		my @ori_r = ($oriseq_rev =~ /$site/g);
		my @mut_r = ($mutseq_rev =~ /$site/g);
		
		$rec .= '+'. $enz. ", " 
		    if (scalar @ori < scalar @mut) or (scalar @ori_r < scalar @mut_r);
		$rec .= '-'. $enz. ", " 		    
		    if (scalar @ori > scalar @mut) or (scalar @ori_r > scalar @mut_r);
		
	    }
	    $rec = substr($rec, 0, CORE::length($rec) - 2) if $rec ne '';
	    $self->{'re_changes'} =  $rec;
	}
    }
    return $self->{'re_changes'}
}


sub _revcompl { 
    # side effect: lower case letters
    my ($seq) = shift;

    $seq = lc $seq;
    $seq =~ tr/acgtrymkswhbvdnx/tgcayrkmswdvbhnx/;
    return CORE::reverse $seq;
}


sub _enzymes {
 #REBASE version 005   type2.005
 my %enzymes =  (
         'AarI' => 'cacctgc',
         'AatII' => 'gacgtc',
         'AccI' => 'gt[ac][gt]ac',
         'AceIII' => 'cagctc',
         'AciI' => 'ccgc',
         'AclI' => 'aacgtt',
         'AcyI' => 'g[ag]cg[ct]c',
         'AflII' => 'cttaag',
         'AflIII' => 'ac[ag][ct]gt',
         'AgeI' => 'accggt',
         'AhaIII' => 'tttaaa',
         'AloI' => 'gaac[acgt][acgt][acgt][acgt][acgt][acgt]tcc',
         'AluI' => 'agct',
         'AlwNI' => 'cag[acgt][acgt][acgt]ctg',
         'ApaBI' => 'gca[acgt][acgt][acgt][acgt][acgt]tgc',
         'ApaI' => 'gggccc',
         'ApaLI' => 'gtgcac',
         'ApoI' => '[ag]aatt[ct]',
         'AscI' => 'ggcgcgcc',
         'AsuI' => 'gg[acgt]cc',
         'AsuII' => 'ttcgaa',
         'AvaI' => 'c[ct]cg[ag]g',
         'AvaII' => 'gg[at]cc',
         'AvaIII' => 'atgcat',
         'AvrII' => 'cctagg',
         'BaeI' => 'ac[acgt][acgt][acgt][acgt]gta[ct]c',
         'BalI' => 'tggcca',
         'BamHI' => 'ggatcc',
         'BbvCI' => 'cctcagc',
         'BbvI' => 'gcagc',
         'BbvII' => 'gaagac',
         'BccI' => 'ccatc',
         'Bce83I' => 'cttgag',
         'BcefI' => 'acggc',
         'BcgI' => 'cga[acgt][acgt][acgt][acgt][acgt][acgt]tgc',
         'BciVI' => 'gtatcc',
         'BclI' => 'tgatca',
         'BetI' => '[at]ccgg[at]',
         'BfiI' => 'actggg',
         'BglI' => 'gcc[acgt][acgt][acgt][acgt][acgt]ggc',
         'BglII' => 'agatct',
         'BinI' => 'ggatc',
         'BmgI' => 'g[gt]gccc',
         'BplI' => 'gag[acgt][acgt][acgt][acgt][acgt]ctc',
         'Bpu10I' => 'cct[acgt]agc',
         'BsaAI' => '[ct]acgt[ag]',
         'BsaBI' => 'gat[acgt][acgt][acgt][acgt]atc',
         'BsaXI' => 'ac[acgt][acgt][acgt][acgt][acgt]ctcc',
         'BsbI' => 'caacac',
         'BscGI' => 'cccgt',
         'BseMII' => 'ctcag',
         'BsePI' => 'gcgcgc',
         'BseRI' => 'gaggag',
         'BseSI' => 'g[gt]gc[ac]c',
         'BsgI' => 'gtgcag',
         'BsiI' => 'cacgag',
         'BsiYI' => 'cc[acgt][acgt][acgt][acgt][acgt][acgt][acgt]gg',
         'BsmAI' => 'gtctc',
         'BsmI' => 'gaatgc',
         'Bsp1407I' => 'tgtaca',
         'Bsp24I' => 'gac[acgt][acgt][acgt][acgt][acgt][acgt]tgg',
         'BspGI' => 'ctggac',
         'BspHI' => 'tcatga',
         'BspLU11I' => 'acatgt',
         'BspMI' => 'acctgc',
         'BspMII' => 'tccgga',
         'BsrBI' => 'ccgctc',
         'BsrDI' => 'gcaatg',
         'BsrI' => 'actgg',
         'BstEII' => 'ggt[acgt]acc',
         'BstXI' => 'cca[acgt][acgt][acgt][acgt][acgt][acgt]tgg',
         'BtrI' => 'cacgtc',
         'BtsI' => 'gcagtg',
         'Cac8I' => 'gc[acgt][acgt]gc',
         'CauII' => 'cc[cg]gg',
         'Cfr10I' => '[ag]ccgg[ct]',
         'CfrI' => '[ct]ggcc[ag]',
         'CjeI' => 'cca[acgt][acgt][acgt][acgt][acgt][acgt]gt',
         'CjePI' => 'cca[acgt][acgt][acgt][acgt][acgt][acgt][acgt]tc',
         'ClaI' => 'atcgat',
         'CviJI' => '[ag]gc[ct]',
         'CviRI' => 'tgca',
         'DdeI' => 'ct[acgt]ag',
         'DpnI' => 'gatc',
         'DraII' => '[ag]gg[acgt]cc[ct]',
         'DraIII' => 'cac[acgt][acgt][acgt]gtg',
         'DrdI' => 'gac[acgt][acgt][acgt][acgt][acgt][acgt]gtc',
         'DrdII' => 'gaacca',
         'DsaI' => 'cc[ag][ct]gg',
         'Eam1105I' => 'gac[acgt][acgt][acgt][acgt][acgt]gtc',
         'EciI' => 'ggcgga',
         'Eco31I' => 'ggtctc',
         'Eco47III' => 'agcgct',
         'Eco57I' => 'ctgaag',
         'EcoNI' => 'cct[acgt][acgt][acgt][acgt][acgt]agg',
         'EcoRI' => 'gaattc',
         'EcoRII' => 'cc[at]gg',
         'EcoRV' => 'gatatc',
         'Esp3I' => 'cgtctc',
         'EspI' => 'gct[acgt]agc',
         'FauI' => 'cccgc',
         'FinI' => 'gggac',
         'Fnu4HI' => 'gc[acgt]gc',
         'FnuDII' => 'cgcg',
         'FokI' => 'ggatg',
         'FseI' => 'ggccggcc',
         'GdiII' => 'cggcc[ag]',
         'GsuI' => 'ctggag',
         'HaeI' => '[at]ggcc[at]',
         'HaeII' => '[ag]gcgc[ct]',
         'HaeIII' => 'ggcc',
         'HaeIV' => 'ga[ct][acgt][acgt][acgt][acgt][acgt][ag]tc',
         'HgaI' => 'gacgc',
         'HgiAI' => 'g[at]gc[at]c',
         'HgiCI' => 'gg[ct][ag]cc',
         'HgiEII' => 'acc[acgt][acgt][acgt][acgt][acgt][acgt]ggt',
         'HgiJII' => 'g[ag]gc[ct]c',
         'HhaI' => 'gcgc',
         'Hin4I' => 'ga[cgt][acgt][acgt][acgt][acgt][acgt][acg]tc',
         'HindII' => 'gt[ct][ag]ac',
         'HindIII' => 'aagctt',
         'HinfI' => 'ga[acgt]tc',
         'HpaI' => 'gttaac',
         'HpaII' => 'ccgg',
         'HphI' => 'ggtga',
         'Hpy178III' => 'tc[acgt][acgt]ga',
         'Hpy188I' => 'tc[acgt]ga',
         'Hpy99I' => 'cg[at]cg',
         'KpnI' => 'ggtacc',
         'Ksp632I' => 'ctcttc',
         'MaeI' => 'ctag',
         'MaeII' => 'acgt',
         'MaeIII' => 'gt[acgt]ac',
         'MboI' => 'gatc',
         'MboII' => 'gaaga',
         'McrI' => 'cg[ag][ct]cg',
         'MfeI' => 'caattg',
         'MjaIV' => 'gt[acgt][acgt]ac',
         'MluI' => 'acgcgt',
         'MmeI' => 'tcc[ag]ac',
         'MnlI' => 'cctc',
         'MseI' => 'ttaa',
         'MslI' => 'ca[ct][acgt][acgt][acgt][acgt][ag]tg',
         'MstI' => 'tgcgca',
         'MwoI' => 'gc[acgt][acgt][acgt][acgt][acgt][acgt][acgt]gc',
         'NaeI' => 'gccggc',
         'NarI' => 'ggcgcc',
         'NcoI' => 'ccatgg',
         'NdeI' => 'catatg',
         'NheI' => 'gctagc',
         'NlaIII' => 'catg',
         'NlaIV' => 'gg[acgt][acgt]cc',
         'NotI' => 'gcggccgc',
         'NruI' => 'tcgcga',
         'NspBII' => 'c[ac]gc[gt]g',
         'NspI' => '[ag]catg[ct]',
         'PacI' => 'ttaattaa',
         'Pfl1108I' => 'tcgtag',
         'PflMI' => 'cca[acgt][acgt][acgt][acgt][acgt]tgg',
         'PleI' => 'gagtc',
         'PmaCI' => 'cacgtg',
         'PmeI' => 'gtttaaac',
         'PpiI' => 'gaac[acgt][acgt][acgt][acgt][acgt]ctc',
         'PpuMI' => '[ag]gg[at]cc[ct]',
         'PshAI' => 'gac[acgt][acgt][acgt][acgt]gtc',
         'PsiI' => 'ttataa',
         'PstI' => 'ctgcag',
         'PvuI' => 'cgatcg',
         'PvuII' => 'cagctg',
         'RleAI' => 'cccaca',
         'RsaI' => 'gtac',
         'RsrII' => 'cgg[at]ccg',
         'SacI' => 'gagctc',
         'SacII' => 'ccgcgg',
         'SalI' => 'gtcgac',
         'SanDI' => 'ggg[at]ccc',
         'SapI' => 'gctcttc',
         'SauI' => 'cct[acgt]agg',
         'ScaI' => 'agtact',
         'ScrFI' => 'cc[acgt]gg',
         'SduI' => 'g[agt]gc[act]c',
         'SecI' => 'cc[acgt][acgt]gg',
         'SexAI' => 'acc[at]ggt',
         'SfaNI' => 'gcatc',
         'SfeI' => 'ct[ag][ct]ag',
         'SfiI' => 'ggcc[acgt][acgt][acgt][acgt][acgt]ggcc',
         'SgfI' => 'gcgatcgc',
         'SgrAI' => 'c[ag]ccgg[ct]g',
         'SimI' => 'gggtc',
         'SmaI' => 'cccggg',
         'SmlI' => 'ct[ct][ag]ag',
         'SnaBI' => 'tacgta',
         'SnaI' => 'gtatac',
         'SpeI' => 'actagt',
         'SphI' => 'gcatgc',
         'SplI' => 'cgtacg',
         'SrfI' => 'gcccgggc',
         'Sse232I' => 'cgccggcg',
         'Sse8387I' => 'cctgcagg',
         'Sse8647I' => 'agg[at]cct',
         'SspI' => 'aatatt',
         'Sth132I' => 'cccg',
         'StuI' => 'aggcct',
         'StyI' => 'cc[at][at]gg',
         'SwaI' => 'atttaaat',
         'TaqI' => 'tcga',
         'TaqII' => 'gaccga',
         'TatI' => '[at]gtac[at]',
         'TauI' => 'gc[cg]gc',
         'TfiI' => 'ga[at]tc',
         'TseI' => 'gc[at]gc',
         'Tsp45I' => 'gt[cg]ac',
         'Tsp4CI' => 'ac[acgt]gt',
         'TspEI' => 'aatt',
         'TspRI' => 'ca[cg]tg[acgt][acgt]',
         'Tth111I' => 'gac[acgt][acgt][acgt]gtc',
         'Tth111II' => 'caa[ag]ca',
         'UbaGI' => 'cac[acgt][acgt][acgt][acgt]gtg',
         'UbaPI' => 'cgaacg',
         'VspI' => 'attaat',
         'XbaI' => 'tctaga',
         'XcmI' => 'cca[acgt][acgt][acgt][acgt][acgt][acgt][acgt][acgt][acgt]tgg',
         'XhoI' => 'ctcgag',
         'XhoII' => '[ag]gatc[ct]',
         'XmaIII' => 'cggccg',
         'XmnI' => 'gaa[acgt][acgt][acgt][acgt]ttc'
        );

    return %enzymes;
}

1;
