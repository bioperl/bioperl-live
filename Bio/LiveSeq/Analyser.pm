#
# Bioperl module for Bio::LiveSeq::Analyser
#
# Cared for by Heikki Lehvaslaiho <heikki@ebi.ac.uk>
#
# Copyright Heikki Lehvaslaiho
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=head1 NAME


  Bio::LiveSeq::Analyser - Package analysing mutated LiveSequences

=head1 SYNOPSIS


=head1 DESCRIPTION

  This package takes results in form of L<Bio::Variation::Haplotype>,
  typically created by L<Bio::LiveSeq::Mutator>, analyses the
  mutations further and adds the results into the Haplotype object.

=head1 FEEDBACK

=head2 Mailing Lists

  User feedback is an integral part of the evolution of this and other
  Bioperl modules. Send your comments and suggestions preferably to one
  of the Bioperl mailing lists.  Your participation is much appreciated.

  vsns-bcd-perl@lists.uni-bielefeld.de          - General discussion
  vsns-bcd-perl-guts@lists.uni-bielefeld.de     - Technically-oriented discussion
  http://bio.perl.org/MailList.html             - About the mailing lists

=head2 Reporting Bugs

  report bugs to the Bioperl bug tracking system to help us keep track
  the bugs and their resolution.  Bug reports can be submitted via
  email or the web:

  bioperl-bugs@bio.perl.org
  http://bio.perl.org/bioperl-bugs/

=head1 AUTHOR - Heikki Lehvaslaiho 

  Email:  heikki@ebi.ac.uk

  Address: 

     EMBL Outstation, European Bioinformatics Institute
     Wellcome Trust Genome Campus, Hinxton
     Cambs. CB10 1SD, United Kingdom 

=head1 APPENDIX

  The rest of the documentation details each of the object
  methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::LiveSeq::Analyser;
$VERSION=2.31;
# Version history:
#Sat May 13 23:33:42 BST 2000 v 1.0 begun

use strict;
use Carp qw(cluck croak carp);
use vars qw($VERSION @ISA);
use Bio::Variation::Haplotype;
use Bio::Variation::DNAMutation;
use Bio::Variation::RNAChange;
use Bio::Variation::AAChange;


=head2 analyse

 Title   : analyse
 Usage   : $results = Bio::LiveSeq::Analyser::analyse(-object => $results,
                                                -mutmatrix => \@mutmatrix)
 Function: 

 Args    : a reference to a LiveSeq object
           a reference to a mutmatrix array of arrays
 Returns : a reference to a hash
 Errorcode: 0

=cut



sub analyse {
    my $obj = shift;

#    foreach my $gene ($obj->each_Gene) {
#	 my $cds_start = $gene->upbound;
#	 my $cds_end   = $gene->downbound;
#	 
    
#    print  "--", substr($obj->aa_ori, -10, 10), "\n";
#    print  "--", $obj->aa_ori, "\n";

    foreach my $mut ($obj->each_Variant) {
	
	    
	if(  $mut->isa('Bio::Variation::DNAMutation') ) {
	    #label
	    $mut->label(_DNA_mutevent($mut->allele_ori,$mut->allele_mut) );
	    #print $mut->label, "\n";
	    
	}
	
	elsif(  $mut->isa('Bio::Variation::RNAChange') ) {
	    
	    #region
	    if ($mut->end < 0 ){
		$mut->region_type('5\'UTR');
	    }
	    elsif ($mut->end >= 1 and $mut->start <= ($obj->cds_end - $obj->offset)) {
		$mut->region_type('coding');
	    }
	    else {
		$mut->region_type('3\'UTR');
	    }
	    # translation table defaults to '1'
	    $mut->translation_table(1) if not $mut->translation_table;

	    #label
	    if (not $mut->region_type('coding')) {
		$mut->label('unknown') ;
		#label
		#$mut->label( _RNA_mutevent($mut->allele_ori,$mut->allele_mut) );
	    }
	    else {
		
	    }
	}
	
	if(  $mut->isa('Bio::Variation::AAChange') ) {
	    #label
	    my ($aalabel, $rnalabel);
	    ($aalabel, $rnalabel) = 
		_AA_mutevent($mut->allele_ori,$mut->allele_mut, $mut->start, $obj->aa_ori);
	    $mut->label($aalabel);
	    $mut->RNAChange->label($rnalabel);
	}	
    }
    return $obj;
}
#}

sub print_alignment {
    my $obj = shift;
    #my (@entry, $text);
    
    foreach my $mut ($obj->each_Variant) {
	if(  $mut->isa('Bio::Variation::AAChange') and $mut->lenght <= 1 ) {
	    
	}    
#	     push (@entry, 
#		   "\n"
#		   );   
#	     
#	     $text = "           ". $mut->allele_mut;
#	     push (@entry, 
#		   $text
#		   );   
#	     
#	     $text = "Variant  : ". $seqmut;
##    $text = "Sequence mut: ". substr($cdna,$ntnumber-1-10, 10).
#	     $ntmut. substr($cdna,$ntnumber, 10);
#	     push (@entry, 
#		   $text
#		   );   
#	     
#	     $text = "Reference: ". $seqori;
##    $text = "Sequence ori: ". substr($cdna,$ntnumber-1-10, 21);
#	     push (@entry, 
#		   $text
#		   );   
#	     
#	     $text = "            ". $aaseqori;
#	     push (@entry, 
#		   $text
#		   );   
#	 }
#
#
#
#	 }
#	 

##-------------------print seq alignment-----------------------
#    if ($query->param('mutation') eq 'point') {
#	 $seq = substr($cdna,$ntnumber-1-10, 21);
#	 for ($i=1; $i<=length($seq); $i++) {
#	     $seqori .=  ' ', $seqmut .= ' ', $aaseqmut .= ' '
#		 if $i%3 == 3-$query->param('codonposition') and 
#		     $atgnumber <= $ntnumber-11+$i and 
#		     $termnumber > $ntnumber-11+$i-2;
#	     $seqori .= substr($seq, $i-1,1);
#	     if ($i==11) {
#		 $seqmut .= '<font color="Red"><b>'. $ntmut. 
#		     '</b><font color="Black">';
#		 $aaseqmut .= ' ' if $query->param('codonposition')==1;
#		 if ($query->param('codonposition')==3) {
#		     chop $aaseqmut;
#		     substr($aaseqori, 0, 1) = '';
#		     $aaseqori = ' '. $aaseqori 
#			 if $query->param('aanumber') <=3;
#		 }
#		 $aaseqmut .= $query->param('aamut');
#	     }
#	     else {
#		 $seqmut .= substr($seq, $i-1,1);
#	     }
#	     if ($i%3 == (4-$query->param('codonposition'))%3) {
#		 if ($atgnumber <= $ntnumber-11+$i and 
#		     $termnumber > $ntnumber-11+$i-2) {
#		     $aaseqori .= 
#			 substr($ttables[$ttabid-1],$codons->{substr($seq,$i-2,3)},1);
#		 }
#	     }
#	     $aaseqori .= ' ';
#	     $aaseqmut .= ' ';
#	 }
#	 chop $aaseqori, chop $aaseqori if $query->param('codonposition')==1;
##    }
##    if ($query->param('mutation') eq 'deletion') {
##	 if (length($ntori) < 10){
##	     $seq = substr($cdna,$ntnumber-1-5, 11 + length($ntori)-1);
##	     for ($i=1; $i<=length($seq); $i++) {
##		 
##	     }
##	 }
##    }
#	 push (@entry, 
#	       "\n"
#	       );   
#	 
#	 $text = "           ". $aaseqmut; 
#	 push (@entry, 
#	       $text
#	       );   
#	 
#	 $text = "Variant  : ". $seqmut;
##    $text = "Sequence mut: ". substr($cdna,$ntnumber-1-10, 10).
#	 $ntmut. substr($cdna,$ntnumber, 10);
#	 push (@entry, 
#	       $text
#	       );   
#	 
#	 $text = "Reference: ". $seqori;
##    $text = "Sequence ori: ". substr($cdna,$ntnumber-1-10, 21);
#	 push (@entry, 
#	       $text
#	       );   
#	 
#	 $text = "            ". $aaseqori;
#	 push (@entry, 
#	       $text
#	       );   
#    }
#

}

}


sub print {
    my ($obj, @args) = @_;

    # HTML | ASCII
    my ($mode) = @args;


    print "\nNo of Variants: ", scalar ($obj->each_Variant),  "\n\n";

	     
    my %tag = 
	(
	 'ID'               => 'ID           ',
	 'Description'      => 'Description  ',
	 'FeatureKey'       => 'Feature      ',
	 'FeatureQual'      => 'Feature        ',
	 'ErrorComment'     => 'Comment      ',
	 'Comment'          => 'Comment      -!-',
	 'Commentline'      => 'Comment         ',
	 );
    
    my @entry =();
    
    my ($text, $tmp, $tmp2);
    my ($count) = 0;
    foreach my $mut ($obj->each_Variant) {
	
	 if(  $mut->isa('Bio::Variation::DNAMutation') ) {
	     
      	     
	     $text = $tag{ID};
	     
	     if ($mode eq 'HTML' ) {
		 $text .= '<A HREF=http://srs.ebi.ac.uk/srs6bin/cgi-bin/wgetz?-e+'.
		     '[embl-acc:'. $obj->id. ']>'. $obj->id. '</A>';
	     }
	     else {
		 $text .= $obj->id;
	     }
	     push (@entry, $text);
	     push (@entry,  
		   $tag{FeatureKey}. 'DNA'. "; ". $mut->mut_number
		   );

	     $text=$tag{FeatureQual}. '/label: '. $mut->label;
	     push (@entry, $text);

	     if ($mut->proof) {
		 $text = $tag{FeatureQual}. '/proof: '. $mut->proof;
		 push (@entry, $text) ;
	     }
     

	     $text = $tag{FeatureQual}. '/location: '. 
		 #$mut->id. '; '. $mut->start; 
		  $mut->start; 

	     if ($mut->end - $mut->start ) {
		 $text .= '..'. $mut->end;

		 $tmp2 = $mut->end - $obj->offset;
		 if ($tmp <= 0) {
		     $tmp -= 1;
		 }
		 
	     }
	     #$offset = $atgnumber-1;
	     $tmp = $mut->start + $obj->offset;
	     if ($tmp <= 0) {
		 $tmp -= 1;
	     }
	     
	     $text.= ' ('. $obj->id. "::". $tmp;
	     if ( ($mut->end - $mut->start) > 0 ) {
		 $text.= "..". $tmp2;
	     }
	     $text .= ')';
	     push (@entry, $text);
	     
	     push (@entry,  
		   $tag{FeatureQual}. '/change: '. $mut->allele_ori. 
		   '>'. $mut->allele_mut
		   );


	     if ($mut->region_type ) {
	         push (@entry,  
		       $tag{FeatureQual}. '/region: '. 
		       $mut->region_type. '; '. $mut->region_value
		       );
	     }

	     if ($mut->re_changes ne '') {
#		  my @enz = $mut->param('re');
#		  my ($rtmp, @tmp);
#		  foreach $rtmp (@enz) {
#		      $rtmp =~ s/([+-])(.*)/$2\]\>$1$2/;
#		      $rtmp = '<A HREF=http://srs6.ebi.ac.uk/srs6bin/cgi-bin/wgetz?-e+[rebase-id:'. 
#			  $rtmp. '</A>';
#		      push (@tmp, $rtmp);
#		  }
#		  $text = join ('; ', sort(@tmp));
#		  push (@entry,
#			$tag{FeatureQual}. '/re_site: '. $text
#			);
	     }
	     if ($mut->CpG) {
		 push (@entry,  
		       $tag{FeatureQual}. "/CpG"
		       );
	     }
    
	     

	     
	 }

	 if(  $mut->isa('Bio::Variation::RNAChange') ) {


	     push (@entry,  
		   $tag{FeatureKey}. 'RNA'. "; ". $mut->mut_number
		   );

	     $text=$tag{FeatureQual}. '/label: '. $mut->label;
	     push (@entry, $text);
	     
	     if ($mut->proof) {
		 $text = $tag{FeatureQual}. '/proof: '. $mut->proof;
		 push (@entry, $text) ;
	     }

	
	     $text = $tag{FeatureQual}. '/location: '. 
		 #$mut->id. '; '. $mut->start; 
		  $mut->start; 

	     if ($mut->end - $mut->start ) {
		 $text .= '..'. $mut->end;

		 $tmp2 = $mut->end - $obj->offset;
		 if ($tmp <= 0) {
		     $tmp -= 1;
		 }
		 
	     }
	     #$offset = $atgnumber-1;
	     $tmp = $mut->start + $obj->offset;
	     if ($tmp <= 0) {
		 $tmp -= 1;
	     }
	     
	     $text.= ' ('. $obj->id. '::'. $tmp;
	     if ( ($mut->end - $mut->start) > 0 ) {
		 $text.= "..". $tmp2;
	     }
	     $text .= ')';
	     push (@entry, $text);
	     
	     push (@entry,  
		   $tag{FeatureQual}. '/change: '. $mut->allele_ori. 
		   '>'. $mut->allele_mut
		   );

         
	     if ($mut->region_type eq 'coding') {
		 $text = $tag{FeatureQual}. '/codon: '. $mut->codon_ori. '>';
		 if (substr($mut->DNAMutation->label,0,5) eq 'point') {
		     $text .= $mut->codon_mut;		     
		 }
		 else {
		     $text .= '-';
		 }
		 $text .= "; ". $mut->codon_pos;
		 push (@entry, $text);
	     }

	     if ($mut->re_changes ne '') {
#		  my @enz = $mut->param('re');
#		  my ($rtmp, @tmp);
#		  foreach $rtmp (@enz) {
#		      $rtmp =~ s/([+-])(.*)/$2\]\>$1$2/;
#		      $rtmp = '<A HREF=http://srs6.ebi.ac.uk/srs6bin/cgi-bin/wgetz?-e+[rebase-id:'. 
#			  $rtmp. '</A>';
#		      push (@tmp, $rtmp);
#		  }
#		  $text = join ('; ', sort(@tmp));
#		  push (@entry,
#			$tag{FeatureQual}. '/re_site: '. $text
#			);
	     }

	     $text =  $tag{FeatureQual}. '/transl_table: ';
	     if ($mode eq 'HTML' ) {
		 $text .=  "<A HREF=http://www.ebi.ac.uk/cgi-bin/mutations/trtables.cgi?".
		     "id=". $mut->translation_table. "&Action=Show>".
		       $mut->translation_table.  '</A>';
	     }
	     else {
		 $text .= $mut->translation_table;
	     }
	     push (@entry, $text);

	

	     #	  
	     #	  if ($mut->param('molecule') eq "DNA") {
	     #	      if (abs $mut->param('splicedist') < 10 ) {
	     #		  push (@entry,  
	     #			$tag{FeatureQual}. '/splicedist: '. $mut->param('splicedist')
	     #			);
	     #	      }
	     #	  }
	     #
	     #
	     #

	     if ($mut->region_type ) {
	         push (@entry,  
		       $tag{FeatureQual}. '/region: '. $mut->region_type
		       );
	     }

	 }

	 if(  $mut->isa('Bio::Variation::AAChange') ) {

	     push (@entry,  
		   $tag{FeatureKey}. 'AA'. "; ". $mut->mut_number
		   );

	     $text=$tag{FeatureQual}. '/label: '. $mut->label;
	     push (@entry, $text) ;

	     if ($mut->proof) {
		 $text = $tag{FeatureQual}. '/proof: '. $mut->proof;
		 push (@entry, $text) ;
	     }

	     $text = $tag{FeatureQual}. '/location: '. 
		 #$mut->id. '; '. $mut->start; 
		  $mut->start; 

	     if ($mut->end - $mut->start ) {
		 $text .= '..'. $mut->end;

		 $tmp2 = $mut->end - $obj->offset;
		 if ($tmp <= 0) {
		     $tmp -= 1;
		 }
		 
	     }
	     #$offset = $atgnumber-1;
	     $tmp = $mut->start + $obj->offset;
	     if ($tmp <= 0) {
		 $tmp -= 1;
	     }
	     push (@entry, $text) ;
#	      $text.= ' ('. $obj->id. '::'. $tmp;
#	      if ( ($mut->end - $mut->start) > 0 ) {
#		  $text.= "..". $tmp2;
#	      }
#	      $text .= ')';
#	      push (@entry, $text);
	     

	     push (@entry,  
		   $tag{FeatureQual}. '/change: '. $mut->allele_ori. 
		   '>'. $mut->allele_mut
		   );


	     if ($mut->region_type ) {
	         push (@entry,  
		       $tag{FeatureQual}. '/region: '. $mut->region_type
		       );
	     }

	 #    if ( $mut->param('coding') and 
	 #	   !($mut->param('mutation') eq 'insertion' and 
	 #	     $mut->param('ntnumber') == $atgnumber+1)) {
	 #	  push (@entry,  
	 #		$tag{FeatureKey}. 'AA'
	 #		);
	 #	  $text = join(", ", $mut->param('aatype'));
	 #	  push (@entry,
	 #		$tag{FeatureQual}. '/label: '. $text
	 #		);
	 #	  
	 #	  $text = $tag{FeatureQual}. '/location: ';
	 #	  $text .= $mut->param('aanumber');
	 #	  if (length($mut->param('aaori')) != 1) {
	 #	      $text .= '..'. ($mut->param('aanumber')+length($mut->param('aaori'))-1);
	 #	  }
	 #	  push (@entry,  
	 #		$text
	 #		);
	 #	  
	 #	  if ($mut->param('indeltype') eq 'inframe' and 
	 #	      $mut->param('codonposition') == 1 ) {
	 #	      if ($mut->param('mutation') eq 'insertion') {
	 #		  push (@entry,  
	 #			$tag{FeatureQual}. '/change: '. '+'. $mut->param('aamut')
	 #			);
	 #	      }
	 #	      elsif ($mut->param('mutation') eq 'deletion') {
	 #		  push (@entry,  
	 #			$tag{FeatureQual}. '/change: '. '-'. $mut->param('aaori')
	 #			);
	 #	      }
	 #	  }
	 #	  else {
	 #	      push (@entry,  
	 #		    $tag{FeatureQual}. '/change: '. $mut->param('aaori'). 
	 #		    '>'. $mut->param('aamut')
	 #		    );
	 #	  }
	 #    }
	 #}


	 }

 

     }
    push (@entry, 
	  "//"
	  );  

    foreach my $line (@entry) {
	print $line, "\n";
    }
}



sub _DNA_mutevent {
    my ($n, $m) = @_;
    my ($type);

    if (length($n) == length($n) && length($n) == 1) {
	$type = 'point';
	$type .= ", ". _pointType($n, $m);
    }
    elsif (length($n) == 0 ) {
	$type = 'insertion';

    }
    elsif (length($m) == 0 ) {
	$type = 'deletion';
    }
    else {
	$type = 'complex';
    }
}

sub _pointType {
    my ($n, $m) = @_;
    my ($type);
    my %transition = ('a' => 'g',
		   'g' => 'a',
		   'c' => 't',
		   't' => 'c');
    $n = lc $n;
    $m = lc $m;
    if ($n eq $m) {
	$type = 'no change';
    }
    elsif ($transition{$n} eq $m ) {
	$type = 'transition';
    }
    else {
	$type = 'transversion';
    }
}


sub CpG {
    my ($ntori, $upflank, $dnflank) = @_;
    #my ($ntori, $ntmut, $upflank, $dnflank) = @_;
    # valid only for point mutations
    # CpG methylation-mediated deamination:
    #   CG -> TG | CG -> CA substitutions
    # implement only weather CpG dinucleotide was hit
    
    if ( ( ($ntori eq 'c') && (substr($upflank,0, 1) eq 'g') ) ||
	 ( ($ntori eq 'g') && (substr($dnflank, -1, 1) eq 'c') ) ) {
	return 1;
    }
    else {
	return 0;
      };
}

sub _RNA_mutevent {
    my ($n, $m, $start, $codn, $codm) = @_;
    my ($type);

    

#    if ($start >= 1 and $start <= 3 ) {
#	 $type='initiation codon';
#    }
#    if ($codn eq $codm ) {
#	 $type = 'silent';
#    }
#
#
#
#    if (length($n) == length($n) && length($n) == 1) {
#	 $type = 'point';
#	 $type .= ", ". _pointType($n, $m);
#    }
#    elsif (length($n) == 0 ) {
#	 $type = 'insertion';
#
#    }
#    elsif (length($m) == 0 ) {
#	 $type = 'deletion';
#    }
#    else {
#	 $type = 'complex';
#    }

    return $type;
}



sub _AA_mutevent {
    my ($n, $m, $start, $len) = @_;
    my ($rnatype, $aatype);
    
    
    if ($start == 1 ) {
	$rnatype='initiation codon';
	
	if (substr($n, 0, 1) ne substr($m, 0, 1)) {
	    $aatype = 'no translation';
	}
	elsif ($n eq $m ) {
	    $aatype = 'silent';
	}
	# more ...
    }
    elsif ($start == $len+1 ) {
	$rnatype='terminator codon';

	if (substr($n, 0, 1) ne substr($m, 0, 1)) {
	    $aatype = 'post-elongation';
	}
	elsif ($n eq $m ) {
	    $aatype = 'silent';
	}
    }
    elsif ($n eq $m) {
	$rnatype = 'silent';
	$aatype = 'silent';
    }
    elsif ($m eq '*') {
	$rnatype = 'nonsense';
	$aatype = 'truncation';
    }
    elsif  ($n ne $m) {
       	$rnatype = 'missense';
	$aatype = 'substitution';
    }
    

#    else { #------------------ indels -----------------------
#	 
#	 if ( $query->param('coding') and 
#	     $query->param('ntnumber') != $atgnumber) {
#	     
#	     $query->append(-name=>'rnatype',-value=>$indeltype);
#	     if ($query->param('indeltype') eq 'inframe'){ 
#		 if ($query->param('mutation') eq 'insertion') {
#		     $query->append(-name=>'aatype',-value=>'insertion');
#		 }
#		 elsif ($query->param('mutation') eq 'deletion') {
#		     $query->append(-name=>'aatype',-value=>'deletion');
#		 }
#	     }
#	     elsif ($indeltype eq 'frameshift') {
#		 if (length($query->param('aamut')) != 1) {
#		     $query->append(-name=>'aatype',-value=>
#				    'out-of-frame extension');
#		 }
#		 $query->append(-name=>'aatype',-value=>'truncation');
#	     }
#	 }
#	 else {
#	     $query->append(-name=>'rnatype',-value=>'unknown');
#	 }
#    }
#    if ($query->param('mutation') eq 'deletion') {
#	 if ($query->param('coding')) {
#	     $query->param('rnatype', $indeltype);
#	 }
#	 else {
#	     $query->param('rnatype', 'unknown');
#	 }
#    }
#    
##########################

    return ( $aatype, $rnatype) ;
       
}



sub trivialname {

#    $aamutterm = substr ($aamut, -1, 1);
#    
#    if  (($query->param('rnatype') eq 'initiation codon'
#	   and $aaori ne $aamut)) {
#	 $aamut = 'X';
#    } 
#    elsif ($query->param('mutation') eq 'point') {
#	 $aamut = 'X' if $aamut eq '*';
#    }
#    elsif ($query->param('mutation') eq 'deletion') {
#	 $aamutsymbol = 'del';
#	 
#	 if ($aamutterm eq '*') {
#	     $aatermnumber = $aanumber + length($aamut)-1;
#	     $aamut = 'X'. $aatermnumber;
#	 }
#	 
#	 $query->append(-name=>'aamutterm',-value=>$aamut); 
#	 
#	 if ($query->param('indeltype') eq 'inframe'){
#	     $aamut = '-'. length($query->param('ntori'))/3 ;
#	 }		
#    }
#    elsif ($query->param('mutation') eq 'insertion') {
#	 $aamutsymbol = 'ins';
#    
#	 if (($aamutterm eq '*') && (length($aamut)-1 != 0)) {
#	     $aatermnumber = $aanumber + length($aamut)-1;
#	     $aamut =  $aatermnumber. 'X';
#	 }
#	 $query->append(-name=>'aamutterm',-value=>$aamut); 
#	 
#	 
#	 if ($query->param('indeltype') eq 'inframe'){
#	     $aamut = '+'. length($query->param('ntmut'))/3 ;
#	 }		
#    }
#    else {
#	 $aamutsymbol = '';
#    }
#    
#    $aaori = substr ($aaori, 0, 1);
#    $pin = $aaori.$aanumber.$aamutsymbol.$aamut;
#    $query->append(-name=>'pin',-value=>$pin); 

}




1;








