# $Id$
# BioPerl module for Bio::Variation::IO::xml
#
# Cared for by Heikki Lehvaslaiho <Heikki@ebi.ac.uk>
#
# Copyright Heikki Lehvaslaiho
#
# You may distribute this module under the same terms as perl itself
#

# POD documentation - main docs before the code

=head1 NAME

Bio::Variation::IO::xml - XML sequence variation input/output stream

=head1 SYNOPSIS

Do not use this module directly.  Use it via the Bio::Variation::IO class.

=head1 DESCRIPTION

This object can transform Bio::Variation::SeqDiff objects to and from XML
file databases.

The XML format, although consistent, is still evolving. The DTD is not
yet available. There is also a formatting problem: The module is not
able to write root tags around entries.

=head1 REQUIREMENTS

To use this code you need the CPAN module XML::Node to read XML and  
modules XML::Writer and IO::String to write XML out.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to the 
Bioperl mailing lists  Your participation is much appreciated.

  bioperl-l@bioperl.org                         - General discussion
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

package Bio::Variation::IO::xml;
my $VERSION=1.0;
use vars qw(@ISA  $h $id $moltype $offset $dna $start $end $len $ismut $number
	    $allele_ori $allele_mut $upFlank $dnFlank $proof $region $region_value 
	    $rna $codon_ori $codon_mut $codon_pos $codon_table $aa $upflank $dnflank
	    $prevdnaobj $prevrnaobj $prevaaobj);
use strict;
# Object preamble - inherits from Bio::Root::Object

#use XML::Parser;        
use XML::Node 0.10; 
use XML::Writer 0.4;
use IO::String;
use Bio::Variation::IO;
use Bio::Variation::SeqDiff;
use Bio::Variation::DNAMutation;
use Bio::Variation::RNAChange;
use Bio::Variation::AAChange;
use Bio::Variation::Allele;

# new() is inherited from Bio::Root::Object
@ISA = qw(Bio::Variation::IO);

# _initialize is where the heavy stuff will happen when new is called

sub new {
    my ($class,@args) = @_;
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


sub _seqDiff {
    $h->id($id);
    $h->moltype($moltype);
    $h->offset($offset);

    $id = $moltype = $offset = '';
}


sub _DNA {
    my ($variation_number, $change_number) = split /\./, $number;
    #$self->warn("|$mut_number|>|$variation_number|$change_number|");
    if ($change_number and $change_number > 1 ) {
	my $a3 = Bio::Variation::Allele->new;
	$a3->seq($allele_mut) if $allele_mut;
	$prevdnaobj->add_Allele($a3);
    } else {
	$dna = Bio::Variation::DNAMutation->new ('-start'         => $start,      
						 '-end'           => $end,        
						 '-lenght'        => $len,        
						 '-mut_number'    => $number,
						 '-upStreamSeq'   => $upFlank,    
						 '-dnStreamSeq'   => $dnFlank,    
						 '-proof'         => $proof,      
						 '-region'        => $region,     
						 '-region_value'  => $region_value
						 );
	
	my $a1 = Bio::Variation::Allele->new;
	$a1->seq($allele_ori) if $allele_ori;
	$dna->allele_ori($a1);
	my $a2 = Bio::Variation::Allele->new;
	$a2->seq($allele_mut) if $allele_mut;
	if ($ismut) {
	    $dna->isMutation(1);
	}
	$dna->allele_mut($a2);	
	$dna->add_Allele($a2);
	$dna->length($len);
	$h->add_Variant($dna);
	$prevdnaobj = $dna;
    }
    $upFlank = $dnFlank = '';
    $start = $end = $len = $ismut = $number = $allele_ori = $allele_mut = 
        $proof = $region = $region_value = '';
}

sub _RNA {
    my ($variation_number, $change_number) = split /\./, $number;
    #$self->warn("|$mut_number|>|$variation_number|$change_number|");
    if ($change_number and $change_number > 1 ) {
	my $a3 = Bio::Variation::Allele->new;
	$a3->seq($allele_mut) if $allele_mut;
	$prevrnaobj->add_Allele($a3);
    } else {
	$rna = Bio::Variation::RNAChange->new ('-start'         => $start,      
					       '-end'           => $end,        
					       '-lenght'        => $len,     
					       '-mut_number'    => $number,   
					       '-upStreamSeq'   => $upFlank,    
					       '-dnStreamSeq'   => $dnFlank,    
					       '-proof'         => $proof,      
					       '-region'        => $region
					       );

	my $a1 = Bio::Variation::Allele->new;
	$a1->seq($allele_ori) if $allele_ori;
	$rna->allele_ori($a1);
	my $a2 = Bio::Variation::Allele->new;
	$a2->seq($allele_mut) if $allele_mut;
	if ($ismut) {
	    $rna->isMutation(1);
	}
	$rna->allele_mut($a2);
	$rna->add_Allele($a2);
	$codon_table = 1 || $rna->codon_table($codon_table);
	$codon_pos  &&  $rna->codon_pos($codon_pos);
	$rna->length($len);
	$dna->RNAChange($rna);
	$rna->DNAMutation($dna);
	$h->add_Variant($rna);
	$prevrnaobj = $rna;
    }	
    $codon_table = $codon_ori = $codon_mut = $codon_pos ='';
    $upFlank = $dnFlank = '';
    $start = $end = $len = $ismut = $number = $allele_ori = $allele_mut = 
	$proof = $region  = '';
}


sub _AA {
    my ($variation_number, $change_number) = split /\./, $number;
    #$self->warn("|$mut_number|>|$variation_number|$change_number|");
    if ($change_number and $change_number > 1 ) {
	my $a3 = Bio::Variation::Allele->new;
	$a3->seq($allele_mut) if $allele_mut;
	$prevaaobj->add_Allele($a3);
    } else {    
	$aa = Bio::Variation::AAChange->new ('-start'         => $start,      
					     '-end'           => $end,        
					     '-lenght'        => $len,    
					     '-mut_number'    => $number,    
					     '-proof'         => $proof,      
					     '-region'        => $region
					     );
	my $a1 = Bio::Variation::Allele->new;
	$a1->seq($allele_ori) if $allele_ori;
	$aa->allele_ori($a1);
	my $a2 = Bio::Variation::Allele->new;
	$a2->seq($allele_mut) if $allele_mut;
	if ($ismut) {
	    $aa->isMutation(1);
	}
	$aa->allele_mut($a2);
	$aa->add_Allele($a2);
	$aa->length($len);
	$aa->region($region) if $region;

	$rna->AAChange($aa);
	$aa->RNAChange($rna);
	$h->add_Variant($aa);
	$prevaaobj = $aa;
    }
    $start = $end = $len = $ismut = $number = $allele_ori = $allele_mut = $upflank = 
	$dnflank = $proof = $region  = '';
}

sub next {
    my( $self ) = @_;

    #$self->_readline unless $self->_readline =~ /<seqDiff/;
    local $/ = '</seqDiff>';

    return unless my $entry = $self->_readline;
    #print  STDERR "|$entry|";
    return unless $entry =~ /^\W*<seqDiff/;

    $id = $offset = '';
    $start = $end = $len = $number = $allele_ori = $allele_mut = $upflank = 
	$dnflank = $proof = $region  = '';
    
    $h = Bio::Variation::SeqDiff->new;

    # create new parser object
    my $p = XML::Node->new();

    # tell object which elements and which attibutes to keep track on and what to do
    $p->register(">seqDiff:id","attr" => \$id);
    $p->register(">seqDiff:moltype","attr" => \$moltype);
    $p->register(">seqDiff:offset","attr" => \$offset);
    $p->register(">seqDiff","end" => \&_seqDiff);

    $p->register(">seqDiff>DNA:number","attr" => \$number);
    $p->register(">seqDiff>DNA:start","attr" => \$start);
    $p->register(">seqDiff>DNA:end","attr" => \$end);
    $p->register(">seqDiff>DNA:length","attr" => \$len);
    $p->register(">seqDiff>DNA:isMutation","attr" => \$ismut);

    $p->register(">seqDiff>RNA:number","attr" => \$number);
    $p->register(">seqDiff>RNA:start","attr" => \$start);
    $p->register(">seqDiff>RNA:end","attr" => \$end);
    $p->register(">seqDiff>RNA:length","attr" => \$len);
    $p->register(">seqDiff>RNA:isMutation","attr" => \$ismut);

    $p->register(">seqDiff>AA:number","attr" => \$number);
    $p->register(">seqDiff>AA:start","attr" => \$start);
    $p->register(">seqDiff>AA:end","attr" => \$end);
    $p->register(">seqDiff>AA:length","attr" => \$len);
    $p->register(">seqDiff>AA:isMutation","attr" => \$ismut);

    $p->register("proof","char" => \$proof);
    $p->register("upFlank","char" => \$upFlank);
    $p->register("dnFlank","char" => \$dnFlank);
    $p->register("allele_ori","char" => \$allele_ori);
    $p->register("allele_mut","char" => \$allele_mut);
    $p->register("region","char" => \$region);
    $p->register(">seqDiff>DNA>region:value","attr" => \$region_value);

    $p->register("codon_table","char" => \$codon_table);
    $p->register(">seqDiff>RNA>codon:codon_ori","attr" => \$codon_ori);
    $p->register(">seqDiff>RNA>codon:codon_mut","attr" => \$codon_mut);
    $p->register(">seqDiff>RNA>codon:codon_pos","attr" => \$codon_pos);

    $p->register("DNA","end" => \&_DNA);
    $p->register("RNA","end" => \&_RNA);
    $p->register("AA","end" => \&_AA);

    #parse the entry string
    $p->parse($entry); 
                                   
    return $h;
}

=head2 write

 Title   : write
 Usage   : $stream->write(@haplos)
 Function: writes the $seqDiff objects into the stream
 Returns : 1 for success and 0 for error
 Args    : Bio::Variation::SeqDiff object

=cut

sub write {
    my ($self,@h) = @_;

    if( !defined $h[0] ) {
        $self->throw("Attempting to write with no information!");
    }
    my $str;
    my $output = IO::String->new($str);     
    my $w = new XML::Writer(OUTPUT => $output, DATA_MODE => 1, DATA_INDENT => 4 ); 

    foreach my $h (@h) {
	#
	# seqDiff
	#
	$h->moltype || $self->throw("Moltype of the reference sequence is not set!");
	my $hasAA = 0;
	foreach my $mut ($h->each_Variant) {	    
	    $hasAA = 1 if  $mut->isa('Bio::Variation::AAChange'); 
	}
	if ($hasAA) {
	    $w->startTag("seqDiff",
			 "id" => $h->id,
			 "moltype" => $h->moltype,
			 "offset" => $h->offset,
			 "sysname" => $h->sysname,
			 "trivname" => $h->trivname
			 );
	} else {
	    $w->startTag("seqDiff",
			 "id" => $h->id,
			 "moltype" => $h->moltype,
			 "offset" => $h->offset,
			 "sysname" => $h->sysname
			 );
	}
	my @allvariants = $h->each_Variant;
	#print "allvars:", scalar @allvariants, "\n";
	my %variants = ();
	foreach my $mut ($h->each_Variant) {
	    #print STDERR  $mut->mut_number, "\t", $mut, "\t",
	    #$mut->proof, "\t", scalar $mut->each_Allele,  "\n";	    
	    push @{$variants{$mut->mut_number} }, $mut; 
	}
	foreach my $var (sort keys %variants) {
	    foreach my $mut (@{$variants{$var}}) {
		#
		# DNA
		#
		if( $mut->isa('Bio::Variation::DNAMutation') ) {
		    $mut->isMutation(0) if not $mut->isMutation;
		    my @alleles = $mut->each_Allele;
		    my $count = 0; 
		    foreach my $allele (@alleles) {
			$count++;
			my ($variation_number, $change_number) = split /\./, $mut->mut_number;
			if ($change_number and $change_number != $count){
			    $mut->mut_number("$change_number.$count");
			}
			$mut->allele_mut($allele);
			$w->startTag("DNA",
				     "number" => $mut->mut_number,
				     "start"  => $mut->start,
				     "end"    => $mut->end,
				     "length" => $mut->length,
				     "isMutation" => $mut->isMutation
				     #param('splicedist') < 10 $mut->param('splicedist')
				     );
			if ($mut->label) {
			    $w->startTag("label");
			    $w->characters($mut->label );
			    $w->endTag;
			}	
			if ($mut->proof) {
			    $w->startTag("proof");
			    $w->characters($mut->proof );
			    $w->endTag;
			}	
			if ($mut->upStreamSeq) {
			    $w->startTag("upFlank");
			    $w->characters($mut->upStreamSeq );
			    $w->endTag;
			}
			#if ( $mut->isMutation) {
			#if ($mut->allele_ori) {
			$w->startTag("allele_ori");
			$w->characters($mut->allele_ori->seq) if $mut->allele_ori->seq ;
			$w->endTag;
			#}	
			#if ($mut->allele_mut) {
			$w->startTag("allele_mut");
			$w->characters($mut->allele_mut->seq) if $mut->allele_mut->seq;
			$w->endTag;
			#}	
			#}
			if ($mut->dnStreamSeq) {
			    $w->startTag("dnFlank");
			    $w->characters($mut->dnStreamSeq );
			    $w->endTag;
			}
			if ($mut->region) {
			    if($mut->region_value) {
				$w->startTag("region",
					     "value" => $mut->region_value
					     );
			    } else {
				$w->startTag("region");
			    }
			    $w->characters($mut->region );
			    $w->endTag;
			}
			if ($mut->restriction_changes) {
			    $w->startTag("restriction_changes");
			    $w->characters($mut->restriction_changes);
			    $w->endTag;
			}	    
			$w->endTag; #DNA
		    }
		}
		#
		# RNA
		#
		elsif(  $mut->isa('Bio::Variation::RNAChange') ) {
		    $mut->isMutation(0) if not $mut->isMutation;
		    my @alleles = $mut->each_Allele;
		    my $count = 0; 
		    foreach my $allele (@alleles) {
			$count++;
			my ($variation_number, $change_number) = split /\./, $mut->mut_number;
			if ($change_number and $change_number != $count){
			    $mut->mut_number("$change_number.$count");
			}
			$mut->allele_mut($allele);
			$w->startTag("RNA",
				     "number" => $mut->mut_number,
				     "start"  => $mut->start,
				     "end"    => $mut->end,
				     "length" => $mut->length,
				     "isMutation" => $mut->isMutation
				     );

			if ($mut->label) {
			    $w->startTag("label");
			    $w->characters($mut->label );
			    $w->endTag;
			}	
			if ($mut->proof) {
			    $w->startTag("proof");
			    $w->characters($mut->proof );
			    $w->endTag;
			}	
			if ($mut->upStreamSeq) {
			    $w->startTag("upFlank");
			    $w->characters($mut->upStreamSeq );
			    $w->endTag;
			}	
			#if ( $mut->isMutation) {
			if ($mut->allele_ori) {
			    $w->startTag("allele_ori");
			    $w->characters($mut->allele_ori->seq) if $mut->allele_ori->seq ;
			    $w->endTag;
			}	
			if ($mut->allele_mut) {
			    $w->startTag("allele_mut");
			    $w->characters($mut->allele_mut->seq) if $mut->allele_mut->seq ;
			    $w->endTag;
			}	
			#}
			if ($mut->dnStreamSeq) {
			    $w->startTag("dnFlank");
			    $w->characters($mut->dnStreamSeq );
			    $w->endTag;
			}
			if ($mut->region eq 'coding') {
			    if (! $mut->codon_mut) {
				$w->startTag("codon",
					     "codon_ori" => $mut->codon_ori,
					     "codon_pos" => $mut->codon_pos
					     );
			    } else {
				$w->startTag("codon",
					     "codon_ori" => $mut->codon_ori,
					     "codon_mut" => $mut->codon_mut,
					     "codon_pos" => $mut->codon_pos
					     );
			    }
			    $w->endTag;
			}
			if ($mut->codon_table != 1) {
			    $w->startTag("codon_table");
			    $w->characters($mut->codon_table);
			    $w->endTag;
			}	
			
			if ($mut->restriction_changes) {
			    $w->startTag("restriction_changes");
			    $w->characters($mut->restriction_changes);
			    $w->endTag;
			}	    
			if ($mut->region) {
			    $w->startTag("region");
			    $w->characters($mut->region );
			    $w->endTag;
			}
			$w->endTag; #RNA
		    }
		}
		#
		# AA
		#
		elsif(  $mut->isa('Bio::Variation::AAChange') ) {
		    $mut->isMutation(0) if not $mut->isMutation;		
		    my @alleles = $mut->each_Allele;
		    my $count = 0; 
		    foreach my $allele (@alleles) {
			$count++;
			my ($variation_number, $change_number) = split /\./, $mut->mut_number;
			if ($change_number and $change_number != $count){
			    $mut->mut_number("$change_number.$count");
			}
			$mut->allele_mut($allele);
			$w->startTag("AA",
				     "number" => $mut->mut_number,
				     "start"  => $mut->start,
				     "end"    => $mut->end,
				     "length" => $mut->length,
				     "isMutation" => $mut->isMutation
				     );

			if ($mut->label) {
			    $w->startTag("label");
			    $w->characters($mut->label );
			    $w->endTag;
			}	
			if ($mut->proof) {
			    $w->startTag("proof");
			    $w->characters($mut->proof );
			    $w->endTag;
			}	
			#if ( $mut->isMutation) {
			if ($mut->allele_ori) {
			    $w->startTag("allele_ori");
			    $w->characters($mut->allele_ori->seq) if $mut->allele_ori->seq;
			    $w->endTag;
			}	
			if ($mut->allele_mut) {
			    $w->startTag("allele_mut");
			    $w->characters($mut->allele_mut->seq) if $mut->allele_mut->seq;
			    $w->endTag;
			}	
			#}
			if ($mut->region) {
			    $w->startTag("region");
			    $w->characters( $mut->region );
			    $w->endTag;
			}
			$w->endTag; #AA
		    }
		}
	    }
	}
    }
    $w->endTag;


    $w->end;
    $self->_print($str);
    $output->close;
    return 1;
}

1;
