#
# BioPerl module for Bio::Tools::Est2Genome
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Jason Stajich <jason-at-bioperl.org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Tools::Est2Genome - Parse est2genome output, makes simple Bio::SeqFeature::Generic objects

=head1 SYNOPSIS

  use Bio::Tools::Est2Genome;

  my $featureiter = Bio::Tools::Est2Genome->new(-file => 'output.est2genome');

  # This is going to be fixed to use the SeqAnalysisI next_feature
  # Method eventually when we have the objects to put the data in
  # properly
  while( my $f = $featureiter->parse_next_gene ) {
   # process Bio::SeqFeature::Generic objects here
  }

=head1 DESCRIPTION

This module is a parser for C<est2genome> [EMBOSS] alignments of est/cdna
sequence to genomic DNA.  This is generally accepted as the best
program for predicting splice sites based on est/dnas (as far as I know).

This module currently does not try pull out the ungapped alignments
(Segment) but may in the future.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

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
of the bugs and their resolution. Bug reports can be submitted the
web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Jason Stajich

Email jason-at-bioperl.org

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Tools::Est2Genome;
use strict;

# Object preamble - inherits from Bio::Root::Root

use Bio::Root::Root;
use Bio::SeqFeature::Gene::Exon;
use Bio::SeqFeature::Gene::Transcript;
use Bio::SeqFeature::Gene::Intron;
use Bio::SeqFeature::Gene::GeneStructure;
use Bio::SeqFeature::SimilarityPair;

use base qw(Bio::Tools::AnalysisResult);

=head2 new

 Title   : new
 Usage   : my $obj = Bio::Tools::Est2Genome->new();
 Function: Builds a new Bio::Tools::Est2Genome object
 Returns : an instance of Bio::Tools::Est2Genome
 Args    : -file => 'output.est2genome' or
           -fh   => \*EST2GENOMEOUTPUT
           -genomefirst => 1  # genome was the first input (not standard)

=cut

sub _initialize_state {
    my($self,@args) = @_;

    # call the inherited method first
    my $make = $self->SUPER::_initialize_state(@args);

    my ($genome_is_first) = $self->_rearrange([qw(GENOMEFIRST)], @args);

    delete($self->{'_genome_is_first'});
    $self->{'_genome_is_first'} = $genome_is_first if(defined($genome_is_first));
    $self->analysis_method("est2genome");
}

=head2 analysis_method

 Usage     : $sim4->analysis_method();
 Purpose   : Inherited method. Overridden to ensure that the name matches
             /est2genome/i.
 Returns   : String
 Argument  : n/a

=cut

#-------------
sub analysis_method {
#-------------
    my ($self, $method) = @_;
    if($method && ($method !~ /est2genome/i)) {
	$self->throw("method $method not supported in " . ref($self));
    }
    return $self->SUPER::analysis_method($method);
}

=head2 parse_next_gene

 Title   : parse_next_gene
 Usage   : @gene = $est2genome_result->parse_next_gene;
           foreach $exon (@exons) {
               # do something
           }

 Function: Parses the next alignments of the est2genome result file and
           returns the found exons as an array of
           Bio::SeqFeature::SimilarityPair objects. Call
           this method repeatedly until an empty array is returned to get the
           results for all alignments.

           The $exon->seq_id() attribute will be set to the identifier of the
           respective sequence for both sequences.
           The length is accessible via the seqlength()
           attribute of $exon->query() and
           $exon->est_hit().
 Returns : An array (or array reference) of Bio::SeqFeature::SimilarityPair and Bio::SeqFeature::Generic objects
           or Bio::SeqFeature::Gene::GeneStructure
 Args    : flag(1/0) indicating to return Bio::SeqFeature::Gene::GeneStructure or Bio::SeqFeature::SimilarityPair
           defaults to 0

=cut

sub parse_next_gene {
   my ($self,$return_gene) = @_;
   return $self->_parse_gene_struct if $return_gene;
   my $seensegment = 0;
   my @features;
   my ($qstrand,$hstrand) = (1,1);
   my $lasthseqname;
   while( defined($_ = $self->_readline) ) {
       if( /Note Best alignment is between (reversed|forward) est and (reversed|forward) genome, (but|and) splice\s+sites imply\s+(forward gene|REVERSED GENE)/) {
	   if( $seensegment ) {
	       $self->_pushback($_);
	       return wantarray ? @features : \@features;
	   }
	   $hstrand = -1 if $1 eq 'reversed';
	   $qstrand = -1 if $4 eq 'REVERSED GENE';
	   #$self->debug( "1=$1, 2=$2, 4=$4\n");
       }
       elsif( /^Exon/ ) {
	   my ($name,$score,$perc_ident,$qstart,$qend,$qseqname,
	       $hstart,$hend, $hseqname) = split;
	   $lasthseqname = $hseqname;
	   my $query = Bio::SeqFeature::Similarity->new(-primary => $name,
						       -source  => $self->analysis_method,
						       -seq_id => $qseqname, # FIXME WHEN WE REDO THE GENERIC NAME CHANGE
						       -start   => $qstart,
						       -end     => $qend,
						       -strand  => $qstrand,
						       -score   => $score,
						       -tag => {
#							   'Location' => "$hstart..$hend",
							   'Sequence' => "$hseqname",
							   'identity' => $perc_ident,
							   }
						       );
	   my $hit = Bio::SeqFeature::Similarity->new(-primary => 'exon_hit',
						     -source  => $self->analysis_method,
						     -seq_id => $hseqname,
						     -start   => $hstart,
						     -end     => $hend,
						     -strand  => $hstrand,
						     -score   => $score,
						     -tag => {
#							 'Location' => "$qstart..$qend",
							 'Sequence' => "$qseqname",
  							 'identity' => $perc_ident,
						     }
						     );
	   push @features, Bio::SeqFeature::SimilarityPair->new
	       (-query => $query,
		-hit   => $hit,
		-source => $self->analysis_method);
       } elsif( /^([\-\+\?])(Intron)/) {
	   my ($name,$score,$perc_ident,$qstart,$qend,$qseqname) = split;
	   push @features, Bio::SeqFeature::Generic->new(-primary => $2,
							-source => $self->analysis_method,
							-start => $qstart,
							-end   => $qend,
							-strand => $qstrand,
							-score  => $score,
							-seq_id => $qseqname,
							-tag => {
							 'identity' => $perc_ident,
							    'Sequence' => $lasthseqname});
       } elsif( /^Span/ ) {
       } elsif( /^Segment/ ) {
	   $seensegment = 1;
       } elsif( /^\s+$/ ) { # do nothing
       } else {
	   $self->warn( "unknown line $_\n");
       }
   }
   return unless( @features );
   return wantarray ? @features : \@features;
}

sub _parse_gene_struct {
   my ($self) = @_;
   my $seensegment = 0;
   my @features;
   my ($qstrand,$hstrand) = (1,1);
   my $lasthseqname;
   my $gene = Bio::SeqFeature::Gene::GeneStructure->new(-source => $self->analysis_method);
   my $transcript = Bio::SeqFeature::Gene::Transcript->new(-source => $self->analysis_method);
   my @suppf;
   my @exon;
   while( defined($_ = $self->_readline) ) {
       if( /Note Best alignment is between (reversed|forward) est and (reversed|forward) genome, (but|and) splice\s+sites imply\s+(forward gene|REVERSED GENE)/) {
	      if( $seensegment ) {
	       $self->_pushback($_);
	       return $gene;
	      }
    	   $hstrand = -1 if $1 eq 'reversed';
    	   $qstrand = -1 if $4 eq 'REVERSED GENE';
       }
       elsif( /^Exon/ ) {
    	   my ($name,$score,$perc_ident,$qstart,$qend,$qseqname,$hstart,$hend, $hseqname) = split;
    	   $lasthseqname = $hseqname;
    	   my $exon = Bio::SeqFeature::Gene::Exon->new(-primary => $name,
						       -source  => $self->analysis_method,
						       -seq_id => $qseqname, # FIXME WHEN WE REDO THE GENERIC NAME CHANGE
						       -start   => $qstart,
						       -end     => $qend,
						       -strand  => $qstrand,
						       -score   => $score,
						       -tag => {
                            #'Location' => "$hstart..$hend",
  							 'identity' => $perc_ident,
           							   'Sequence' => "$hseqname",
							              }
						       );
          $transcript->seq_id($qseqname) unless $transcript->seq_id;
          $exon->add_tag_value('phase',0);
          push @exon, $exon;
              
       } elsif( /^([\-\+\?])(Intron)/) {
         next; #intron auto matically built from exons..hope thats ok..
       } elsif( /^Span/ ) {
       } elsif( /^Segment/ ) {
    	    my ($name,$score,$perc_ident,$qstart,$qend,$qseqname,$hstart,$hend, $hseqname) = split;
	         my $query = Bio::SeqFeature::Similarity->new(-primary => $name,
						       -source  => $self->analysis_method,
						       -seq_id => $qseqname, # FIXME WHEN WE REDO THE GENERIC NAME CHANGE
						       -start   => $qstart,
						       -end     => $qend,
						       -strand  => $qstrand,
						       -score   => $score,
						       -tag => {
#							   'Location' => "$hstart..$hend",
							   'Sequence' => "$hseqname",
  							 'identity' => $perc_ident,
							   }
						     );
      	   my $hit = Bio::SeqFeature::Similarity->new(-primary => 'exon_hit',
                                           						     -source  => $self->analysis_method,
                                          						     -seq_id => $hseqname,
                                          						     -start   => $hstart,
                                          						     -end     => $hend,
                                          						     -strand  => $hstrand,
                                          						     -score   => $score,
                                          						     -tag => {
                                                            #	'Location' => "$qstart..$qend",
                                               							 'Sequence' => "$qseqname",
  							 'identity' => $perc_ident,
						                                                }
						     );
        	   my $support =  Bio::SeqFeature::SimilarityPair->new(-query => $query,
                                                              		-hit   => $hit,
                                                              		-source => $self->analysis_method);
             push @suppf, $support;
       } elsif( /^\s+$/ ) { # do nothing
       } else {
      	   $self->warn( "unknown line $_\n");
       }
   }
   return unless $#exon >=0;
   foreach my $e(@exon){
    my @add;
    foreach my $sf(@suppf){
      if($sf->overlaps($e)){
          push @add,$sf;
      }
    }
    $e->add_tag_value('supporting_feature',@add);
    $transcript->add_exon($e);
  }
  
   $gene->add_transcript($transcript);
   $gene->seq_id($transcript->seq_id);
   return $gene;
}

=head2 next_feature

 Title   : next_feature
 Usage   : $seqfeature = $obj->next_feature();
 Function: Returns the next feature available in the analysis result, or
           undef if there are no more features.
 Example :
 Returns : A Bio::SeqFeatureI implementing object, or undef if there are no
           more features.
 Args    : none

=cut

sub next_feature {
    my ($self) = shift;
    $self->throw("We haven't really done this right, yet, use parse_next_gene");
}


1;
