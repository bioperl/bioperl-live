#
# BioPerl module for Bio::Tools::Spidey::Results
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Ryan Golhar <golharam@umdnj.edu>
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Tools::Spidey::Results - Results of a Spidey run

=head1 SYNOPSIS

   use Bio::Tools::Spidey::Results;
	my $spidey = Bio::Tools::Spidey::Results->new(-file => 'result.spidey' );

	# or

	my $spidey = Bio::Tools::Spidey::Results->new( -fh   => \*INPUT );

	# get the exons before doing anything else
	my $exonset = $spidey->next_exonset();

	# parse the results
	my @exons = $exonset->sub_SeqFeature();
	print "Total no of Exons: ", scalar(@exons), "\n";

	print "Genomic sequence length: ", $spidey->genomic_dna_length(), "\n";

	# $exonset is-a Bio::SeqFeature::Generic with Bio::Tools::Spidey::Exons
	# as sub features
	print "Delimited on sequence ", $exonset->seq_id(), " from ", 
		$exonset->start(), " to ", $exonset->end(), "\n";

	foreach my $exon ( $exonset->sub_SeqFeature() ) {
		# $exon is-a Bio::SeqFeature::FeaturePair
		print "Exon from ", $exon->start, " to ", $exon->end, 
			" on strand ", $exon->strand(), "\n";
		# you can get out what it matched using the est_hit attribute
		my $homol = $exon->est_hit();
		print "Matched to sequence ", $homol->seq_id, 
			" at ", $homol->start," to ", $homol->end, "\n";
	}

	# essential if you gave a filename at initialization (otherwise 
  	# the file stays open)
	$spidey->close();

=head1 DESCRIPTION

The spidey module provides a parser and results object for spidey 
output. The spidey results are specialised types of SeqFeatures, 
meaning you can add them to AnnSeq objects fine, and manipulate them 
in the "normal" seqfeature manner.

The spidey Exon objects are Bio::SeqFeature::FeaturePair inherited 
objects. The $esthit = $exon-E<gt>est_hit() is the alignment as a 
feature on the matching object (normally, a cDNA), in which the 
start/end points are where the hit lies.

To make this module work sensibly you need to run

     spidey -i genomic.fasta -m cDNA.fasta

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

=head1 AUTHOR - Ryan Golhar

Email golharam@umdnj.edu

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Tools::Spidey::Results;
use strict;

use File::Basename;
use Bio::Root::Root;
use Bio::Tools::Spidey::Exon;

use base qw(Bio::Tools::AnalysisResult);

sub _initialize_state {
    my($self,@args) = @_;

    # call the inherited method first
    my $make = $self->SUPER::_initialize_state(@args);

#    my ($est_is_first) = $self->_rearrange([qw(ESTFIRST)], @args);

#    delete($self->{'_est_is_first'});
#    $self->{'_est_is_first'} = $est_is_first if(defined($est_is_first));
    $self->analysis_method("Spidey");
}

=head2 analysis_method

 Usage     : $spidey->analysis_method();
 Purpose   : Inherited method. Overridden to ensure that the name matches
             /Spidey/i.
 Returns   : String
 Argument  : n/a

=cut

#-------------
sub analysis_method { 
#-------------
    my ($self, $method) = @_;  
    if($method && ($method !~ /Spidey/i)) {
	$self->throw("method $method not supported in " . ref($self));
    }
    return $self->SUPER::analysis_method($method);
}

=head2 parse_next_alignment

 Title   : parse_next_alignment
 Usage   : @exons = $spidey_result->parse_next_alignment;
           foreach $exon (@exons) {
               # do something
           }
 Function: Parses the next alignment of the Spidey result file and returns the
           found exons as an array of Bio::Tools::Spidey::Exon objects. Call
           this method repeatedly until an empty array is returned to get the
           results for all alignments.
 Example :
 Returns : An array of Bio::Tools::Spidey::Exon objects
 Args    :


=cut

sub parse_next_alignment {
    my ($self) = @_;
    # for strand 1 = plus, -1 = minus
    my ($started,$version,$strand, $exoncount) = (0,0,0,-1);
    my (%seq1props,%seq2props,@exons);
    
    # we refer to the properties of each seq by reference

    while(defined($_ = $self->_readline())) {
	chomp;
	#
	# bascially, parse a Spidey result...
	#
	# matches: --SPIDEY version 1.40--
	if( /^--SPIDEY\s+version\s+(\d+\.\d+)--/) {
	    if($started) {
		$self->_pushback($_);
		return \@exons;
	    }
	    $version = $1;
	    if ($version != 1.40) {
		$self->throw("Spidey parser only designed to work with Spidey v1.40\n");
	    }
	    $started = 1;
	} elsif (/^Genomic:\s+(\S+)\s.*,\s+(\d+)\sbp$/ ) {
	    # matches: Genomic: lcl|some_name other information, 1234 bp
	    # $seq1props{'filename'} = $1;
	    $seq1props{'seqname'} = $1;
	    $seq1props{'length'} = $2;	   
	    $self->genomic_dna_length($seq1props{'length'});

	} elsif( /^mRNA:\s+(\S+)\s.*,(?:\s+mRNA\s+sequence,)?\s(\d+)\sbp$/ ) {
	    # matches: mRNA:
	    # $seq2props{'filename'} = $1;
	    $seq2props{'seqname'} = $1;
	    $seq2props{'length'} = $2;

	} elsif( /^Strand:/ ) {
	    if (/plus/) {
		$strand = 1;
	    } else {
		$strand = -1;
	    }
	} elsif( /^Number of exons: (\d+)/ ) {
	    $exoncount = $1;

	    my ($genomic_start, $genomic_stop, $cdna_start, $cdna_stop,
		$id, $mismatches, $gaps, $splice_donor, 
		$splice_acceptor, $uncertain);

	    # the next $exoncount lines contains information
	    # about the matches of each exon.  we should parse
	    # this information here

	    for (my $ec = 1; $ec <= $exoncount; $ec++) {
		if (defined($_ = $self->_readline())) {
		    chomp;

		    if (/^Exon\s$ec[\(\)-]*:\s(\d+)-(\d+)\s\(gen\)\s+(\d+)-(\d+)\s\(mRNA\)\s+id\s([\d\.inf-]+)%\s+mismatches\s(\d+)\s+gaps\s(\d+)\s+splice\ssite\s\(d\s+a\):\s(\d+)\s+(\d+)\s*(\w*)/) {
			$genomic_start = $1;
			$genomic_stop = $2;
			$cdna_start = $3;
			$cdna_stop = $4;
			$id = $5;
			$mismatches = $6;
			$gaps = $7;
			$splice_donor = $8;
			$splice_acceptor = $9;
			$uncertain = $10;
		    } else {
			$self->throw( "Failed to match anything:\n$_\n");
		    }

		    my $exon = Bio::Tools::Spidey::Exon->new
			(-start  => $genomic_start,
			 -end    => $genomic_stop,
			 -strand => $strand);
		    $exon->seq_id($seq1props{'seqname'});

		    # feature1 is supposed to be initialized to a Similarity object, but we provide a safety net
		    if ($exon->feature1->can('seqlength')) {
			$exon->feature1->seqlength($seq1props{'length'});
		    } else {
			$exon->feature1->add_tag_value('seqlength', $seq1props{'length'});
		    }

		    # create and initialize the feature wrapping the 'hit' (the cDNA)
		    my $fea2 = Bio::SeqFeature::Similarity->new
			(-start   => $cdna_start,
			 -end     => $cdna_stop,
			 -strand  => $strand,
			 -seq_id  => $seq2props{'seqname'},
			 -primary => "aligning_cDNA");
		    $fea2->seqlength($seq2props{'length'});
		    # store
		    $exon->est_hit($fea2);	   

		    # general properties
		    $exon->source_tag($self->analysis_method());
		    $exon->percentage_id($5);
		    $exon->mismatches($6);
		    $exon->gaps($7);
		    $exon->donor($8);
		    $exon->acceptor($9);

		    # push onto array
		    push(@exons, $exon);
		} else {
		    $self->throw("Unexpected end of file reached\n");
		}
	    }
	} elsif( /^Number of splice sites:\s+(\d+)/ ) {
	    $self->splicesites($1);	
	} elsif( /^mRNA coverage:\s+(\d+)%/ ) { 
	    $self->est_coverage($1);
	} elsif(/^overall percent identity:\s+([\d\.]+)%/ ) {
	    $self->overall_percentage_id($1);
	} elsif(/^Missing mRNA ends:\s+(\w+)/ ) {
	    $self->missing_mrna_ends($1);
	} elsif( /^Exon (\d+): (\d+)-(\d+) \(gen\)\s+(\d+)-(\d+) \(mRNA\)/ ) {
	    my ($exon_num, $gen_start, $gen_stop, $cdna_start, $cdna_stop);						
	    $exon_num = $1;
	    $gen_start = $2;
	    $gen_stop = $3;
	    $cdna_start = $4;
	    $cdna_stop = $5;			
	} elsif( /No alignment found/ ) {
	    return [];
	} else {
	    #$self->debug("unmatched $_\n");
	}
    }
    # Typical format:
    # 	Exon 1: 36375798-36375691 (gen)  1-108 (mRNA)
    #
    #
    #	CCTCTTTTTCTTTGCAGGGTATATACCCAGTTACTTAGACAAGGATGAGCTATGTGTAGT
    #        	   |  ||||||||||||||||||||||||||||||||||||||||||||||
    #	          ATGTCAGGGTATATACCCAGTTACTTAGACAAGGATGAGCTATGTGTAGT
    #	           M  S  G  Y  I  P  S  Y  L  D  K  D  E  L  C  V  V 
    #
    #
    #	ATGTGGGGACAAAGCCACCGGATATCATTATCGCTGCATCACTTGTGAAGGTTGCAAGGT
    #	||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    #	ATGTGGGGACAAAGCCACCGGATATCATTATCGCTGCATCACTTGTGAAGGTTGCAAG
    #	  C  G  D  K  A  T  G  Y  H  Y  R  C  I  T  C  E  G  C  K 
    #
    #
    #	AAATGGCA
    #
    @exons ? return \@exons : return ;
}

=head2 next_exonset

  Title   : next_exonset
  Usage   : $exonset = $spidey_result->parse_next_exonset;
         print "Exons start at ", $exonset->start(), 
        "and end at ", $exonset->end(), "\n";
         for $exon ($exonset->sub_SeqFeature()) {
	    # do something
         }
  Function: Parses the next alignment of the Spidey result file and returns the
       set of exons as a container of features. The container is itself
       a Bio::SeqFeature::Generic object, with the Bio::Tools::Spidey::Exon
       objects as sub features. Start, end, and strand of the container
       will represent the total region covered by the exons of this set.

      See the documentation of parse_next_alignment() for further
      reference about parsing and how the information is stored.
 Example : 
 Returns : An Bio::SeqFeature::Generic object holding Bio::Tools::Spidey::Exon
          objects as sub features.
 Args    :

=cut

sub next_exonset {
    my $self = shift;
    my $exonset;

    # get the next array of exons
    my $exons = $self->parse_next_alignment();
    if( ! defined $exons ) {
        $self->warn("No exons returned");
        return;
    } 
    if( @$exons == 0 ) {
	return Bio::SeqFeature::Generic->new();
    }
    # create the container of exons as a feature object itself, with the
    # data of the first exon for initialization
    $exonset = Bio::SeqFeature::Generic->new('-start' => $exons->[0]->start(),
					     '-end' => $exons->[-1]->end(),
					     '-strand' => $exons->[0]->strand(),
					     '-primary' => "ExonSet");
    $exonset->source_tag($exons->[0]->source_tag());
    $exonset->seq_id($exons->[0]->seq_id());
    # now add all exons as sub features, with enabling EXPANsion of the region
    # covered in total
    foreach my $exon (@$exons) {
	$exonset->add_sub_SeqFeature($exon, 'EXPAND');
    }
    return $exonset;
}

=head2 next_feature

  Title   : next_feature
  Usage   : while($exonset = $spidey->next_feature()) {
	    # do something
	   }
  Function: Does the same as L<next_exonset()>. See there for documentation of
      the functionality. Call this method repeatedly until FALSE is
      returned.

      The returned object is actually a SeqFeatureI implementing object.
      This method is required for classes implementing the
      SeqAnalysisParserI interface, and is merely an alias for 
      next_exonset() at present.

  Example :
  Returns : A Bio::SeqFeature::Generic object.
  Args    :

=cut

sub next_feature {
    my ($self,@args) = @_;
    # even though next_exonset doesn't expect any args (and this method
    # does neither), we pass on args in order to be prepared if this changes
    # ever
    return $self->next_exonset(@args);
}

=head2 genomic_dna_length

    Title   : genomic_dna_length
    Usage   : $spidey->genomic_dna_length();
    Function: Returns the length of the genomic DNA used in this Spidey result
    Example :
    Returns : An integer value.
    Args    :

=cut

sub genomic_dna_length {
    my ($self, @args) = @_;
    my $val;

    if(@args) {
	$val = shift(@args);
	$self->{'genomic_dna_length'} = $val;
    } else {
	$val = $self->{'genomic_dna_length'};
    }
    return $val;
}

=head2 splicesites

    Title   : splicesites
    Usage   : $spidey->splicesites();
    Function: Returns the number of splice sites found in this Spidey result
    Example :
    Returns : An integer value.
    Args    :

=cut

sub splicesites {
    my ($self, @args) = @_;
    my $val;

    if(@args) {
	$val = shift(@args);
	$self->{'splicesites'} = $val;
    } else {
	$val = $self->{'splicesites'};
    }
    return $val;
}

=head2 est_coverage

    Title   : est_coverage
    Usage   : $spidey->est_coverage();
    Function: Returns the percent of est coverage in this Spidey result
    Example :
    Returns : An integer value.
    Args    :

=cut

sub est_coverage {
     my ($self, @args) = @_;
     my $val;
     
     if(@args) {
	 $val = shift(@args);
	 $self->{'est_coverage'} = $val;
     } else {
	 $val = $self->{'est_coverage'};
     }
     return $val;
 }

=head2 overall_percentage_id

    Title   : overall_percentage_id
    Usage   : $spidey->overall_percentage_id();
    Function: Returns the overall percent id in this Spidey result
    Example :
    Returns : An float value.
    Args    :

=cut

sub overall_percentage_id {
    my ($self, @args) = @_;
    my $val;

    if(@args) {
	$val = shift(@args);
	$self->{'overall_percentage_id'} = $val;
    } else {
	$val = $self->{'overall_percentage_id'};
    }
    return $val;
}

=head2 missing_mrna_ends

    Title   : missing_mrna_ends
    Usage   : $spidey->missing_mrna_ends();
    Function: Returns left/right/neither from Spidey
    Example :
    Returns : A string value.
    Args    :

=cut

sub missing_mrna_ends {
    my ($self, @args) = @_;
    my $val;

    if(@args) {
	$val = shift(@args);
	$self->{'missing_mrna_ends'} = $val;
    } else {
	$val = $self->{'missing_mrna_ends'};
    }
    return $val;
}

1;
