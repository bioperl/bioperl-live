# BioPerl module for Bio::Variation::IO::xml
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

Bio::Variation::IO::xml - XML sequence variation input/output stream

=head1 SYNOPSIS

Do not use this module directly.  Use it via the Bio::Variation::IO class.

=head1 DESCRIPTION

This object can transform L<Bio::Variation::SeqDiff> objects to and from XML
file databases.

The XML format, although consistent, is still evolving. The current
DTD for it is at L<http://www.ebi.ac.uk/mutations/DTDE/seqDiff.dtd>.

=head1 REQUIREMENTS

To use this code you need the module L<XML::Twig> which creates an
interface to L<XML::Parser> to read XML and modules L<XML::Writer> and
L<IO::String> to write XML out.

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

package Bio::Variation::IO::xml;

use vars qw($seqdiff $var $prevdnaobj $prevrnaobj $prevaaobj);
use strict;

use XML::Twig;
use XML::Writer 0.4;
use IO::String;
use Bio::Variation::SeqDiff;
use Bio::Variation::DNAMutation;
use Bio::Variation::RNAChange;
use Bio::Variation::AAChange;
use Bio::Variation::Allele;

use base qw(Bio::Variation::IO);

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
    my ($t, $term)= @_;
    $seqdiff->id( $term->att('id') );
    $seqdiff->alphabet( $term->att('moltype') );
    $seqdiff->offset( $term->att('offset') );

    foreach my $child ($term->children) {
	 _variant($t, $child);
    }
}

sub _variant {
    my ($t, $term)= @_;
    my $var;
    my $att = $term->atts();
    my ($variation_number, $change_number) = split /\./,  $att->{number}; 
    
    # if more than two alleles
    if ($variation_number and $change_number and $change_number > 1 ) {
	my $a3 = Bio::Variation::Allele->new;
	$a3->seq( $term->first_child_text('allele_mut') ) 
	    if $term->first_child_text('allele_mut');
	if ($term->gi eq 'DNA') {
	    $prevdnaobj->add_Allele($a3);
	}
	elsif ($term->gi eq 'RNA') {
	    $prevrnaobj->add_Allele($a3);
	} else { # AA
	    $prevaaobj->add_Allele($a3);
	}
    } else { # create new variants
	if ($term->gi eq 'DNA') {
	    $var = Bio::Variation::DNAMutation->new();
	}
	elsif ($term->gi eq 'RNA') {
	    $var = Bio::Variation::RNAChange->new();
	} else { # AA
	    $var = Bio::Variation::AAChange->new();
	}

	# these are always present
	$var->start( $att->{start} );
	$var->end( $att->{end});
	$var->length($att->{len});
	$var->mut_number( $att->{number});
	$var->upStreamSeq($term->first_child_text('upFlank'));
	$var->dnStreamSeq($term->first_child_text('dnFlank'));
	$var->proof($term->first_child_text('proof'));

	# region
	my $region = $term->first_child('region');
	if ($region) {
	    $var->region($region->text);
	    my $region_atts = $region->atts;
	    $var->region_value( $region_atts->{value} ) 
		if $region_atts->{value};
	    $var->region_dist( $region_atts->{dist} ) 
		if $region_atts->{dist};
	}
	
	# alleles
	my $a1 = Bio::Variation::Allele->new;
	$a1->seq($term->first_child_text('allele_ori') ) 
	    if $term->first_child_text('allele_ori');
	$var->allele_ori($a1);
	my $a2 = Bio::Variation::Allele->new;
	$a2->seq($term->first_child_text('allele_mut') ) 
	    if $term->first_child_text('allele_mut');
	$var->isMutation(1) if $term->att('isMutation');
	$var->allele_mut($a2);	
	$var->add_Allele($a2);
	$var->length( $term->att('length') );
	$seqdiff->add_Variant($var);

	# variant specific code
	if ($term->gi eq 'DNA') {
	    $prevdnaobj = $var;
	}
	elsif ($term->gi eq 'RNA') {
	    my $codon = $term->first_child('codon');
	    if ($codon) {
		my $codon_atts = $codon->atts;
		$var->codon_table( $codon->att('codon_table') ) 
		    if $codon_atts->{codon_table} and $codon_atts->{codon_table} != 1;
		$var->codon_pos( $codon->att('codon_pos') ) 
		    if $codon_atts->{codon_pos};
	    }
	    $prevdnaobj->RNAChange($var);
	    $var->DNAMutation($prevdnaobj);
	    $prevrnaobj = $var;
	} else {
	    $prevrnaobj->AAChange($var);
	    $var->RNAChange($prevrnaobj);
	    $prevaaobj = $var;
	}
    }
}

sub next {
    my( $self ) = @_;

    local $/ = "</seqDiff>\n";
    return unless my $entry = $self->_readline;
#    print  STDERR "|$entry|";
    return unless $entry =~ /^\W*<seqDiff/;

    $seqdiff = Bio::Variation::SeqDiff->new;

    # create new parser object
    my $twig_handlers = {'seqDiff' =>  \&_seqDiff };
    my $t = XML::Twig->new ( TwigHandlers => $twig_handlers,
			    KeepEncoding => 1 );
    $t->parse($entry);

    return $seqdiff;
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
    my $w = XML::Writer->new(OUTPUT => $output, DATA_MODE => 1, DATA_INDENT => 4 );

    foreach my $h (@h) {
	#
	# seqDiff
	#
	$h->alphabet || $self->throw("Moltype of the reference sequence is not set!");
	my $hasAA = 0;
	foreach my $mut ($h->each_Variant) {	
	    $hasAA = 1 if  $mut->isa('Bio::Variation::AAChange');
	}
	if ($hasAA) {
	    $w->startTag("seqDiff",
			 "id" => $h->id,
			 "moltype" => $h->alphabet,
			 "offset" => $h->offset,
			 "sysname" => $h->sysname,
			 "trivname" => $h->trivname
			 );
	} else {
	    $w->startTag("seqDiff",
			 "id" => $h->id,
			 "moltype" => $h->alphabet,
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
				     );
			if ($mut->label) {
			    foreach my $label (split ', ', $mut->label) {
				$w->startTag("label");
				$w->characters($label);
				$w->endTag;
			    }
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
			if ($mut->restriction_changes) {
			    $w->startTag("restriction_changes");
			    $w->characters($mut->restriction_changes);
			    $w->endTag;
			}	
			if ($mut->region) {
			    if($mut->region_value and $mut->region_dist) {
				$w->startTag("region",
					     "value" => $mut->region_value,
					     "dist" => $mut->region_dist
					     );
			    }
			    elsif($mut->region_value) {
				$w->startTag("region",
					     "value" => $mut->region_value
					     );
			    }
			    elsif($mut->region_dist) {
				$w->startTag("region",
					     "dist" => $mut->region_dist
					     );
			    } else {
				$w->startTag("region");
			    }
			    $w->characters($mut->region );
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
			    foreach my $label (split ', ', $mut->label) {
				$w->startTag("label");
				$w->characters($label );
				$w->endTag;
			    }
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
			    if($mut->region_value and $mut->region_dist) {
				$w->startTag("region",
					     "value" => $mut->region_value,
					     "dist" => $mut->region_dist
					     );
			    }
			    elsif($mut->region_value) {
				$w->startTag("region",
					     "value" => $mut->region_value
					     );
			    }
			    elsif($mut->region_dist) {
				$w->startTag("region",
					     "dist" => $mut->region_dist
					     );
			    } else {
				$w->startTag("region");
			    }
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
			    foreach my $label (split ', ', $mut->label) {
				$w->startTag("label");
				$w->characters($label );
				$w->endTag;
			    }
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
			    if($mut->region_value and $mut->region_dist) {
				$w->startTag("region",
					     "value" => $mut->region_value,
					     "dist" => $mut->region_dist
					     );
			    }
			    elsif($mut->region_value) {
				$w->startTag("region",
					     "value" => $mut->region_value
					     );
			    }
			    elsif($mut->region_dist) {
				$w->startTag("region",
					     "dist" => $mut->region_dist
					     );
			    } else {
				$w->startTag("region");
			    }
			    $w->characters($mut->region );
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
    $output = undef;
    return 1;
}

1;
