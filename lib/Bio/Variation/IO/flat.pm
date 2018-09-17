# BioPerl module for Bio::Variation::IO::flat
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

=head1 NAME

Bio::Variation::IO::flat - flat file sequence variation input/output stream

=head1 SYNOPSIS

Do not use this module directly. Use it via the Bio::Variation::IO class.

=head1 DESCRIPTION

This object can transform Bio::Variation::SeqDiff objects to and from
flat file databases.

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

report bugs to the Bioperl bug tracking system to help us keep track
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

package Bio::Variation::IO::flat;

use strict;

use Text::Wrap;
use Bio::Variation::SeqDiff;
use Bio::Variation::DNAMutation;
use Bio::Variation::RNAChange;
use Bio::Variation::AAChange;
use Bio::Variation::Allele;


use base qw(Bio::Variation::IO);

sub new {
    my($class, @args) = @_;
    my $self = bless {}, $class;
    $self->_initialize(@args);
    return $self;
}

sub _initialize {
  my($self,@args) = @_;
  return unless $self->SUPER::_initialize(@args);
}

=head2 next


 Title   : next
 Usage   : $haplo = $stream->next()
 Function: returns the next seqDiff in the stream
 Returns : Bio::Variation::SeqDiff object
 Args    : NONE

=cut

sub next {
    my( $self ) = @_;
    local $/ = '//';
    return unless my $entry = $self->_readline;
    
    return if $entry =~ /^\s+$/;

    $entry =~ /\s*ID\s+\S+/ || $self->throw("We do need an ID!");

    my ($id, $offset, $alphabet) = $entry =~ /\s*ID +([^:]+)..(\d+)[^\)]*.\[?([cg])?/
	or $self->throw("Can't parse ID line");
#    $self->throw("$1|$2|$3");
    my $h =Bio::Variation::SeqDiff->new(-id         => $id,
					-offset     => $offset,
					  );
    if ($alphabet) { 
	if ($alphabet eq 'g') {
	    $alphabet = 'dna';
	} 
	elsif ($alphabet eq 'c') {
	    $alphabet = 'rna';
	}
	$h->alphabet($alphabet);
    }
    #
    # DNA 
    #
    my @dna = split ( / DNA;/, $entry );
    shift @dna;
    my $prevdnaobj;
    foreach my $dna (@dna) {
	$dna =~ s/Feature[ \t]+//g;
	($dna) = split "RNA; ", $dna; 
	#$self->warn("|$dna|") ;
	#exit;
	my ($mut_number, $proof, $location, $upflank, $change, $dnflank) = 
	    $dna =~ m|\W+([\d\.]+).+/proof: (\w+).+/location: ([^ \n]+).+/upflank: ([ \n\w]+).+/change: ([^ /]+).+/dnflank: ([ \n\w]+)|s;
	$change =~ s/[ \n]//g;	
	my ($ori, $mut) =  split /[>\|]/, $change;
	my ($variation_number, $change_number) = split /\./, $mut_number;
	#$self->warn("|$mut_number|>|$variation_number|$change_number|");
	my $dnamut;
	if ($change_number and $change_number > 1 ) {
	    my $a3 = Bio::Variation::Allele->new;
	    $a3->seq($mut) if $mut;
	    #$dnamut->add_Allele($a3);
	    $prevdnaobj->add_Allele($a3);
	} else {
	    $upflank =~ s/[ \n]//g;
	    $dnflank =~ s/[ \n]//g;
	    my ($region, $junk, $region_value, $junk2, $region_dist) =  
		$dna =~ m|.+/region: ([\w\']+)(; )?(\w+)?( ?\(\+?)?(-?\d+)?|s;
	    #my $s = join ("|", $mut_number, $proof, $location, $upflank, 
	    #	     $change, $dnflank, $region, $region_value, $region_dist, $1,$2,$3,$4,$5);
	    #$self->warn($s);
	    #exit;
	    my ($start, $sep, $end) = $location =~ /(-?\d+)(.)?\D?(-?\d+)?/;
	    $end = $start if not defined $end ;
	    my ($len) = $end - $start +1; 
	    $len = 0, $start = $end if defined $sep and $sep eq '^';
	    my $ismut = 0;
	    $ismut = 1 if $change =~ m/>/; 
	    
	    $dnamut = Bio::Variation::DNAMutation->new 
		('-start'         => $start,
		 '-end'           => $end,
		 '-length'        => $len,
		 '-upStreamSeq'   => $upflank,
		 '-dnStreamSeq'   => $dnflank,
		 '-proof'         => $proof,
		 '-mut_number'    => $mut_number
		 );
	    $prevdnaobj = $dnamut;
	    my $a1 = Bio::Variation::Allele->new;
	    $a1->seq($ori) if $ori;
	    $dnamut->allele_ori($a1);
	    my $a2 = Bio::Variation::Allele->new;
	    $a2->seq($mut) if $mut;
	    $dnamut->add_Allele($a2);
	    if ($ismut) {
		$dnamut->isMutation(1);
		$dnamut->allele_mut($a2);
	    }
	    $dnamut->region($region) if defined $region;
	    $dnamut->region_value($region_value) if defined $region_value;
	    $dnamut->region_dist($region_dist) if defined $region_dist;

	    $h->add_Variant($dnamut);
	    $dnamut->SeqDiff($h);
	}
    }

    #
    # RNA 
    #
    my @rna = split ( / RNA;/, $entry );
    shift @rna;
    my $prevrnaobj;
    foreach my $rna (@rna) {
	$rna = substr ($rna, 0, index($rna, 'Feature      AA'));
	$rna =~ s/Feature[ \t]+//g;
	($rna) = split "DNA; ", $rna; 
	#$self->warn("|$rna|") ;
	my ($mut_number, $proof, $location, $upflank, $change, $dnflank) = 
	    $rna =~ m|\W+([\d\.]+).+/proof: (\w+).+/location: ([^ \n]+).+/upflank: (\w+).+/change: ([^/]+).+/dnflank: (\w+)|s ;#'
	my ($region, $junk, $region_value, $junk2, $region_dist) =  
	    $rna =~ m|.+/region: ([\w\']+)(; )?(\w+)?( ?\(\+?)?(-?\d+)?|s;
	#my $s = join ("|", $mut_number, $proof, $location, $upflank, 
	#	      $change, $dnflank, $region, $region_value, $region_dist, $1,$2,$3,$4,$5);
	#$self->warn($s);
	#exit;
	$change =~ s/[ \n]//g;	
	my ($ori, $mut) =  split /[>\|]/, $change;
	my $rnamut;
	my ($variation_number, $change_number) = split /\./, $mut_number;
	if ($change_number and $change_number > 1 ) {
	    my $a3 = Bio::Variation::Allele->new;
	    $a3->seq($mut) if $mut;
	    #$rnamut->add_Allele($a3);
	    $prevrnaobj->add_Allele($a3);
	} else {
	    my ($start, $sep, $end) = $location =~ /(-?\d+)(.)?\D?(-?\d+)?/;
	    $end = $start if not defined $end ;
	    my ($len) = $end - $start + 1; 
	    $len = 0, $start = $end if defined $sep and $sep eq '^'; 
	    my $ismut;
	    $ismut = 1 if $change =~ m/>/; 
	    my ($codon_table) = $rna =~ m|.+/codon_table: (\d+)|s;
	    my ($codon_pos) = $rna =~ m|.+/codon:[^;]+; ([123])|s;

	    $rnamut = Bio::Variation::RNAChange->new 
		('-start'         => $start,
		 '-end'           => $end,
		 '-length'        => $len,
		 '-upStreamSeq'   => $upflank,
		 '-dnStreamSeq'   => $dnflank,
		 '-proof'         => $proof,
		 '-mut_number'    => $mut_number
		 
		 );
	    $prevrnaobj = $rnamut;
	    my $a1 = Bio::Variation::Allele->new;
	    $a1->seq($ori) if $ori;
	    $rnamut->allele_ori($a1);
	    my $a2 = Bio::Variation::Allele->new;
	    $a2->seq($mut) if $mut;
	    $rnamut->add_Allele($a2);
	    if ($ismut) {
		$rnamut->isMutation(1);
		$rnamut->allele_mut($a2);
	    }
	    $rnamut->region($region) if defined $region;
	    $rnamut->region_value($region_value) if defined $region_value;
	    $rnamut->region_dist($region_dist) if defined $region_dist;

	    $rnamut->codon_table($codon_table) if $codon_table;
	    $rnamut->codon_pos($codon_pos) if $codon_pos;
	    $h->add_Variant($rnamut);
	    foreach my $mut ($h->each_Variant) {
		if ($mut->isa('Bio::Variation::DNAMutation') ) {
		    if ($mut->mut_number == $rnamut->mut_number) {
			$rnamut->DNAMutation($mut);
			$mut->RNAChange($rnamut);
		    }
		}
	    }
	}
    }    
    #
    # AA 
    #
    my @aa = split ( / AA;/, $entry );
    shift @aa;
    my $prevaaobj;
    foreach my $aa (@aa) {
	$aa = substr ($aa, 0, index($aa, 'Feature      AA'));
	$aa =~ s/Feature[ \t]+//g;
	($aa) = split "DNA; ", $aa; 
	#$self->warn("|$aa|") ;
	my ($mut_number, $proof, $location, $change) = 
	    $aa =~ m|\W+([\d\.]+).+/proof: (\w+).+/location: ([^ \n]+)./change: ([^/;]+)|s;
	$change =~ s/[ \n]//g;	
	#my $s = join ("|", $mut_number, $proof, $location, $change);
	#$self->warn($s);
	#exit;
	$change =~ s/[ \n]//g;
	$change =~ s/DNA$//;
	my ($ori, $mut) =  split /[>\|]/, $change;
	#print "------$location----$ori-$mut-------------\n";
	my ($variation_number, $change_number) = split /\./, $mut_number;
	my $aamut;
	if ($change_number and $change_number > 1 ) {
	    my $a3 = Bio::Variation::Allele->new;
	    $a3->seq($mut) if $mut;
	    $prevaaobj->add_Allele($a3);
	} else {
	    my ($start, $sep, $end) = $location =~ /(-?\d+)(.)?\D?(-?\d+)?/;
	    $end = $start if not defined $end ;
	    my ($len) = $end - $start + 1; 
	    $len = 0, $start = $end if defined $sep and $sep eq '^'; 
	    my $ismut;
	    $ismut = 1 if $change =~ m/>/; 
	    my ($region) =  $aa =~ m|.+/region: (\w+)|s ;	
	    $aamut = Bio::Variation::AAChange->new 
		('-start'         => $start,
		 '-end'           => $end,
		 '-length'        => $len,
		 '-proof'         => $proof,
		 '-mut_number'    => $mut_number	     
		 );
	    $prevaaobj = $aamut;
	    my $a1 = Bio::Variation::Allele->new;
	    $a1->seq($ori) if $ori;
	    $aamut->allele_ori($a1);
	    my $a2 = Bio::Variation::Allele->new;
	    $a2->seq($mut) if $mut;
	    $aamut->add_Allele($a2);
	    if ($ismut) {
		$aamut->isMutation(1);
		$aamut->allele_mut($a2);
	    }
	    $region && $aamut->region($region);
	    $h->add_Variant($aamut); 
	    foreach my $mut ($h->each_Variant) {
		if ($mut->isa('Bio::Variation::RNAChange') ) {
		    if ($mut->mut_number == $aamut->mut_number) {
			$aamut->RNAChange($mut);
			$mut->AAChange($aamut);
		    }
		}
	    }

	}
    }
    return $h;
}

=head2 write

 Title   : write
 Usage   : $stream->write(@seqDiffs)
 Function: writes the $seqDiff object into the stream
 Returns : 1 for success and 0 for error
 Args    : Bio::Variation::SeqDiff object


=cut

sub write {
    my ($self,@h) = @_;

    #$columns = 75; #default for Text::Wrap
    my %tag = 
	(
	 'ID'               => 'ID           ',
	 'Description'      => 'Description  ',
	 'FeatureKey'       => 'Feature      ',
	 'FeatureQual'      => "Feature        ",
	 'FeatureWrap'      => "Feature         ",
	 'ErrorComment'     => 'Comment      '
	 #'Comment'          => 'Comment      -!-',
	 #'CommentLine'      => 'Comment         ',
	 );
    
    if( !defined $h[0] ) {
        $self->throw("Attempting to write with no information!");
    }

    foreach my $h (@h) {
	
	my @entry =();
    
	my ($text, $tmp, $tmp2, $sep);
	my ($count) = 0;

	
	$text = $tag{ID};
	
	$text .= $h->id;
	$text .= ":(". $h->offset;
	$text .= "+1" if $h->sysname =~ /-/;
	$text .= ")".  $h->sysname;
	$text .= "; ".  $h->trivname if $h->trivname;
	push (@entry, $text);

	#Variants need to be ordered accoding to mutation_number attribute
	#put them into a hash of arrays holding the Variant objects 
	#This is necessary for cases like several distict mutations present 
	# in the same sequence.
	my @allvariants = $h->each_Variant;
	my %variants = ();
	foreach my $mut ($h->each_Variant) {
	    push @{$variants{$mut->mut_number} }, $mut; 
	}
	#my ($variation_number, $change_number) = split /\./, $mut_number;
	foreach my $var (sort keys %variants) {
	    #print $var, ": ", join (" ", @{$variants{$var}}), "\n";	
	    
	    foreach my $mut (@{$variants{$var}}) {
		#
		# DNA
		#
		if ( $mut->isa('Bio::Variation::DNAMutation') ) {
		    #collect all non-reference alleles
		     $self->throw("allele_ori needs to be defined in [$mut]") 
			 if not $mut->allele_ori;
		     if ($mut->isMutation) {
			 $sep = '>';
		     } else {
			 $sep = '|';
		     }
		     my @alleles = $mut->each_Allele;
		     #push @alleles, $mut->allele_mut if $mut->allele_mut;
		     my $count = 0; # two alleles
		     foreach my $allele (@alleles) {
			 $count++;
			 my ($variation_number, $change_number) = split /\./, $mut->mut_number;
			 if ($change_number and $change_number != $count){
			     $mut->mut_number("$change_number.$count");
			 }
			 $mut->allele_mut($allele);
			 push (@entry,  
			       $tag{FeatureKey}. 'DNA'. "; ". $mut->mut_number
			       );
			 #label
			 $text=$tag{FeatureQual}. '/label: '. $mut->label;
			 push (@entry, $text);
			 
			 #proof
			 if ($mut->proof) {
			     $text = $tag{FeatureQual}. '/proof: '. $mut->proof;
			     push (@entry, $text) ;
			 }
			 #location
			 $text = $tag{FeatureQual}. '/location: '; 
			     #$mut->id. '; '. $mut->start; 
			 if ($mut->length > 1 ) {#	    if ($mut->end - $mut->start ) {
			     my $l = $mut->start + $mut->length -1;
			     $text .= $mut->start. '..'.  $l;
			 }
			 elsif ($mut->length == 0) {
			     my $tmp_start = $mut->start - 1;
			     $tmp_start-- if $tmp_start == 0;
			     $text .= $tmp_start. '^'. $mut->end;
			 } else {
			     $text .= $mut->start;
			 }

			 if ($h->alphabet && $h->alphabet eq 'dna') {
			     $tmp = $mut->start + $h->offset;
			     $tmp-- if $tmp <= 0;
			     $mut->start < 1 && $tmp++; 
			     #$text.= ' ('. $h->id. '::'. $tmp;
			     $tmp2 = $mut->end + $h->offset;
			     if ( $mut->length > 1 ) {
				 $mut->end < 1 && $tmp2++; 
				 $text.= ' ('. $h->id. '::'. $tmp. "..". $tmp2;
			     }
			     elsif ($mut->length == 0) {
				 $tmp--;
				 $tmp-- if $tmp == 0;
				 $text .= ' ('. $h->id. '::'. $tmp. '^'. $tmp2;
			     } else {
				 $text.= ' ('. $h->id. '::'. $tmp;
			     }
			     $text .= ')';
			 }
			 push (@entry, $text);
			 #sequence
			 push (@entry,  
			       $tag{FeatureQual}. '/upflank: '. $mut->upStreamSeq
			       );
			 $text = '';
			 $text = $mut->allele_ori->seq if $mut->allele_ori->seq;
			 $text .= $sep;
			 $text .= $mut->allele_mut->seq if $mut->allele_mut->seq;
			 push (@entry,  
			       wrap($tag{FeatureQual}. '/change: ', $tag{FeatureWrap}, 
				    $text)
			       );
			 
			 push (@entry,  
			       $tag{FeatureQual}. '/dnflank: '. $mut->dnStreamSeq
			       );
			 #restriction enzyme
			 if ($mut->restriction_changes ne '') {
			     $text = $mut->restriction_changes;
			     $text = wrap($tag{FeatureQual}. '/re_site: ', $tag{FeatureWrap}, $text); 
			     push (@entry,
				   $text
				   );
			 }
			 #region
			 if ($mut->region ) {
			     $text = $tag{FeatureQual}. '/region: '. $mut->region;
			     $text .= ';' if $mut->region_value or $mut->region_dist; 
			     $text .= ' '. $mut->region_value if $mut->region_value;
			     if ($mut->region_dist ) {
				 $tmp = '';
				 $tmp = '+' if $mut->region_dist > 1;
				 $text .= " (". $tmp. $mut->region_dist. ')';
			     }
			     push (@entry, $text);
			 }
			 #CpG
			 if ($mut->CpG) {
			     push (@entry,  
				   $tag{FeatureQual}. "/CpG"
				   );
			 }
		     }
		 }
		 #
		 # RNA
		 #	    
		 elsif ($mut->isa('Bio::Variation::RNAChange') ) {
		     #collect all non-reference alleles
		     $self->throw("allele_ori needs to be defined in [$mut]") 
			 if not $mut->allele_ori;
		     my @alleles = $mut->each_Allele;
		     #push @alleles, $mut->allele_mut if $mut->allele_mut;
		     if ($mut->isMutation) {
			 $sep = '>';
		     } else {
			 $sep = '|';
		     }

		     my $count = 0; # two alleles
		     foreach my $allele (@alleles) {
			 $count++;
			 my ($variation_number, $change_number) = split /\./, $mut->mut_number;
			 if ($change_number and $change_number != $count){
			     $mut->mut_number("$change_number.$count");
			 }
			 $mut->allele_mut($allele);
			 push (@entry,  
			       $tag{FeatureKey}. 'RNA'. "; ". $mut->mut_number
			       );
			 #label
			 $text=$tag{FeatureQual}. '/label: '. $mut->label;
			 push (@entry, $text);
			 #proof
			 if ($mut->proof) {
			     $text = $tag{FeatureQual}. '/proof: '. $mut->proof;
			     push (@entry, $text) ;
			 }
			 #location
			 $text = $tag{FeatureQual}. '/location: ' ; 
			 if ($mut->length > 1 ) {
			     $text .= $mut->start. '..'. $mut->end;
			     $tmp2 = $mut->end + $h->offset;
			 }
			 elsif ($mut->length == 0) {
			     my $tmp_start = $mut->start;
			     $tmp_start--;
			     $tmp_start-- if $tmp_start == 0;
			     $text .= $tmp_start. '^'. $mut->end;
			 } else {
			     $text .= $mut->start;
			 }

			 if ($h->alphabet && $h->alphabet eq 'rna') {
			     $tmp = $mut->start + $h->offset;
			     $tmp-- if $tmp <= 0;
			     #$mut->start < 1 && $tmp++;			     
			     #$text.= ' ('. $h->id. '::'. $tmp;
			     $tmp2 = $mut->end + $h->offset;
			     #$mut->end < 1 && $tmp2++; 
			     if ( $mut->length > 1 ) {
				 $text.= ' ('. $h->id. '::'. $tmp. "..". $tmp2;
			     }
			     elsif ($mut->length == 0) {
				 $tmp--;
				 $text .= ' ('. $h->id. '::'. $tmp. '^'. $tmp2;
			     } else {
				 $text.= ' ('. $h->id. '::'. $tmp;
			     }

			     $text .= ')';
			 }
			 push (@entry, $text);

			 #sequence
			 push (@entry,  
			       $tag{FeatureQual}. '/upflank: '. $mut->upStreamSeq
			       );
			 $text = '';
			 $text = $mut->allele_ori->seq if $mut->allele_ori->seq;
			 $text .= $sep;
			 $text .= $mut->allele_mut->seq if $mut->allele_mut->seq;
			 push (@entry,  
			       wrap($tag{FeatureQual}. '/change: ', $tag{FeatureWrap}, 
				    $text)
			       );
			 push (@entry,  
			       $tag{FeatureQual}. '/dnflank: '. $mut->dnStreamSeq
			       );
			 #restriction
			 if ($mut->restriction_changes ne '') {
			     $text = $mut->restriction_changes;
			     $text = wrap($tag{FeatureQual}. '/re_site: ', $tag{FeatureWrap}, $text); 
			     push (@entry,
				   $text
				   );
			 }
			 #coding
			 if ($mut->region eq 'coding') {
			     #codon table
			     $text =  $tag{FeatureQual}. '/codon_table: ';
			     $text .= $mut->codon_table;
			     push (@entry, $text);
			     #codon

			     $text = $tag{FeatureQual}. '/codon: '. $mut->codon_ori. $sep;
			     if ($mut->DNAMutation->label =~ /.*point/) {
				 $text .= $mut->codon_mut;		     
			     }
			     else {
				 $text .= '-';
			     }
			     $text .= "; ". $mut->codon_pos;
			     push (@entry, $text);
			 }
			 #region
			 if ($mut->region ) {
			     $text = $tag{FeatureQual}. '/region: '. $mut->region;
			     $text .= ';' if $mut->region_value or $mut->region_dist; 
			     $text .= ' '. $mut->region_value if $mut->region_value;
			     if ($mut->region_dist ) {
				 $tmp = '';
				 $tmp = '+' if $mut->region_dist > 1;
				 $text .= " (". $tmp. $mut->region_dist. ')';
			     }
			     push (@entry, $text);
			 }
		     }
		 }
		 #
		 # AA
		 #	    
		 elsif ($mut->isa('Bio::Variation::AAChange')) {
		     #collect all non-reference alleles
		     $self->throw("allele_ori needs to be defined in [$mut]") 
			 if not $mut->allele_ori;
		     if ($mut->isMutation) {
			 $sep = '>';
		     } else {
			 $sep = '|';
		     }
		     my @alleles = $mut->each_Allele;
		     #push @alleles, $mut->allele_mut if $mut->allele_mut;
		     my $count = 0; # two alleles		     
		     foreach my $allele (@alleles) {
			 $count++;
			 my ($variation_number, $change_number) = split /\./, $mut->mut_number;
			 if ($change_number and $change_number != $count){
			     $mut->mut_number("$change_number.$count");
			 }
			 $mut->allele_mut($allele);
			 push (@entry,  
			       $tag{FeatureKey}. 'AA'. "; ". $mut->mut_number
			       );
			 #label
			 $text=$tag{FeatureQual}. '/label: '. $mut->label;
			 push (@entry, $text) ;
			 #proof
			 if ($mut->proof) {
			     $text = $tag{FeatureQual}. '/proof: '. $mut->proof;
			     push (@entry, $text) ;
			 }
			 #location
			 $text = $tag{FeatureQual}. '/location: '. 
			     #$mut->id. '; '. $mut->start; 
			     $mut->start; 
			 if ($mut->length > 1 ) {
			     $tmp = $mut->start + $mut->length -1;
			     $text .= '..'. $tmp;
			 }
			 push (@entry, $text);
			 #sequence
			 $text = '';
			 $text = $mut->allele_ori->seq if $mut->allele_ori->seq;
			 $text .= $sep;
			 $text .= $mut->allele_mut->seq if $mut->allele_mut->seq;
			 push (@entry,  
			       wrap($tag{FeatureQual}. '/change: ', $tag{FeatureWrap}, 
				    $text)
			       );
			 #region
			 if ($mut->region ) {
			     $text = $tag{FeatureQual}. '/region: '. $mut->region;
			     $text .= ';' if $mut->region_value or $mut->region_dist; 
			     $text .= ' '. $mut->region_value if $mut->region_value;
			     if ($mut->region_dist ) {
				 $tmp = '';
				 $tmp = '+' if $mut->region_dist > 1;
				 $text .= " (". $tmp. $mut->region_dist. ')';
			     }
			     push (@entry, $text);
			 }
		     }
		  }
	     }
	}
	push (@entry, 
	      "//"
	      );  
	my $str = join ("\n", @entry). "\n";
	$str =~ s/\t/        /g;
	$self->_print($str);
    }
    return 1;
}

1;
