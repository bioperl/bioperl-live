#
# BioPerl module for Bio::Tools::Sim4::Results
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Ewan Birney <birney-at-sanger.ac.uk>
#          and Hilmar Lapp <hlapp-at-gmx.net>
#
# Copyright Ewan Birney and Hilmar Lapp
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Tools::Sim4::Results - Results of one Sim4 run

=head1 SYNOPSIS

   # to preset the order of EST and genomic file as given on the sim4 
   # command line:
   my $sim4 = Bio::Tools::Sim4::Results->new(-file => 'result.sim4',
                                             -estfirst => 1);
   # to let the order be determined automatically (by length comparison):
   $sim4 = Bio::Tools::Sim4::Results->new( -file => 'sim4.results' );
   # filehandle:
   $sim4 = Bio::Tools::Sim4::Results->new( -fh   => \*INPUT );

   # parse the results
   while(my $exonset = $sim4->next_exonset()) {
       # $exonset is-a Bio::SeqFeature::Generic with Bio::Tools::Sim4::Exons
       # as sub features
       print "Delimited on sequence ", $exonset->seq_id(), 
             "from ", $exonset->start(), " to ", $exonset->end(), "\n";
       foreach my $exon ( $exonset->sub_SeqFeature() ) {
	  # $exon is-a Bio::SeqFeature::FeaturePair
	  print "Exon from ", $exon->start, " to ", $exon->end, 
                " on strand ", $exon->strand(), "\n";
          # you can get out what it matched using the est_hit attribute
          my $homol = $exon->est_hit();
          print "Matched to sequence ", $homol->seq_id, 
                " at ", $homol->start," to ", $homol->end, "\n";
      }
   }

   # essential if you gave a filename at initialization (otherwise the file
   # stays open)
   $sim4->close();

=head1 DESCRIPTION

The sim4 module provides a parser and results object for sim4 output. The
sim4 results are specialised types of SeqFeatures, meaning you can add them
to AnnSeq objects fine, and manipulate them in the "normal" seqfeature manner.

The sim4 Exon objects are Bio::SeqFeature::FeaturePair inherited objects. The 
$esthit = $exon-E<gt>est_hit() is the alignment as a feature on the matching 
object (normally, an EST), in which the start/end points are where the hit
lies.

To make this module work sensibly you need to run

     sim4 genomic.fasta est.database.fasta
or
     sim4 est.fasta genomic.database.fasta

To get the sequence identifiers recorded for the first sequence, too, use
A=4 as output option for sim4.

One fiddle here is that there are only two real possibilities to the matching
criteria: either one sequence needs reversing or not. Because of this, it
is impossible to tell whether the match is in the forward or reverse strand
of the genomic DNA. We solve this here by assuming that the genomic DNA is
always forward. As a consequence, the strand attribute of the matching EST is
unknown, and the strand attribute of the genomic DNA (i.e., the Exon object)
will reflect the direction of the hit.

See the documentation of parse_next_alignment() for abilities of the parser
to deal with the different output format options of sim4.

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

=head1 AUTHOR - Ewan Birney, Hilmar Lapp

Ewan Birney E<lt>birney-at-sanger.ac.ukE<gt>
Hilmar Lapp E<lt>hlapp-at-gmx.netE<gt> or E<lt>hilmar.lapp-at-pharma.novartis.comE<gt>.

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Tools::Sim4::Results;
use strict;


use File::Basename;
use Bio::Root::Root;
use Bio::Tools::Sim4::Exon;

use base qw(Bio::Tools::AnalysisResult);


sub _initialize_state {
    my($self,@args) = @_;

    # call the inherited method first
    my $make = $self->SUPER::_initialize_state(@args);

    my ($est_is_first) = $self->_rearrange([qw(ESTFIRST)], @args);

    delete($self->{'_est_is_first'});
    $self->{'_est_is_first'} = $est_is_first if(defined($est_is_first));
    $self->analysis_method("Sim4");
}

=head2 analysis_method

 Usage     : $sim4->analysis_method();
 Purpose   : Inherited method. Overridden to ensure that the name matches
             /sim4/i.
 Returns   : String
 Argument  : n/a

=cut

#-------------
sub analysis_method { 
#-------------
    my ($self, $method) = @_;  
    if($method && ($method !~ /sim4/i)) {
	$self->throw("method $method not supported in " . ref($self));
    }
    return $self->SUPER::analysis_method($method);
}

=head2 parse_next_alignment

 Title   : parse_next_alignment
 Usage   : @exons = $sim4_result->parse_next_alignment;
           foreach $exon (@exons) {
               # do something
           }
 Function: Parses the next alignment of the Sim4 result file and returns the
           found exons as an array of Bio::Tools::Sim4::Exon objects. Call
           this method repeatedly until an empty array is returned to get the
           results for all alignments.

           The $exon->seq_id() attribute will be set to the identifier of the
           respective sequence for both sequences if A=4 was used in the sim4
           run, and otherwise for the second sequence only. If the output does
           not contain the identifier, the filename stripped of path and 
           extension is used instead. In addition, the full filename 
           will be recorded for both features ($exon inherits off 
           Bio::SeqFeature::SimilarityPair) as tag 'filename'. The length
           is accessible via the seqlength() attribute of $exon->query() and
           $exon->est_hit().

           Note that this method is capable of dealing with outputs generated
           with format 0,1,3, and 4 (via the A=n option to sim4). It
           automatically determines which of the two sequences has been 
           reversed, and adjusts the coordinates for that sequence. It will
           also detect whether the EST sequence(s) were given as first or as
           second file to sim4, unless this has been specified at creation
           time of the object.

 Example :
 Returns : An array of Bio::Tools::Sim4::Exon objects
 Args    :


=cut

sub parse_next_alignment {
   my ($self) = @_;
   my @exons = ();
   my %seq1props = ();
   my %seq2props = ();
   # we refer to the properties of each seq by reference
   my ($estseq, $genomseq, $to_reverse);
   my $started = 0;
   my $hit_direction = 1;
   my $output_fmt = 3; # same as 0 and 1 (we cannot deal with A=2 produced
                       # output yet)
   
   while(defined($_ = $self->_readline())) {
       #chomp();
       #
       # bascially, each sim4 'hit' starts with seq1...
       #
       /^seq1/ && do {
	   if($started) {
	       $self->_pushback($_);
	       last;
	   }
	   $started = 1;

	   # filename and length of seq 1
	   /^seq1\s+=\s+(\S+)\,\s+(\d+)/ ||
	       $self->throw("Sim4 parsing error on seq1 [$_] line. Sorry!");
	   $seq1props{'filename'} = $1;
	   $seq1props{'length'} = $2;
	   next;
       };
       /^seq2/ && do {
	   # the second hit has also the database name in the >name syntax 
	   # (in brackets).
	   /^seq2\s+=\s+(\S+)\s+\(>?(\S+)\s*\)\,\s+(\d+)/||
	       $self->throw("Sim4 parsing error on seq2 [$_] line. Sorry!");
	   $seq2props{'filename'} = $1;
	   $seq2props{'seqname'} = $2;
	   $seq2props{'length'} = $3;
	   next;
       };
       if(/^>(\S+)\s*(.*)$/) {
	   # output option was A=4, which not only gives the complete
	   # description lines, but also causes the longer sequence to be
	   # reversed if the second file contained one (genomic) sequence
	   $seq1props{'seqname'} = $1;
	   $seq1props{'description'} = $2 if $2;
	   $output_fmt = 4;
	   # we handle seq1 and seq2 both here
	   if(defined($_ = $self->_readline()) && (/^>(\S+)\s*(.*)$/)) {
	       $seq2props{'seqname'} = $1; # redundant, since already set above
	       $seq2props{'description'} = $2 if $2;
	   }
	   next;
       }
       /^\(complement\)/ && do {
	   $hit_direction = -1;
	   next;
       };
       # this matches
       # start-end (start-end) pctid%
       if(/(\d+)-(\d+)\s+\((\d+)-(\d+)\)\s+(\d+)%/) {
 	   $seq1props{'start'} = $1;
 	   $seq1props{'end'} = $2;
 	   $seq2props{'start'} = $3;
 	   $seq2props{'end'} = $4;
	   my $pctid   = $5;
	   
	   if(! defined($estseq)) {
	       # for the first time here: need to set the references referring
	       # to seq1 and seq2 
	       if(! exists($self->{'_est_is_first'})) {
		   # detect which one is the EST by looking at the lengths,
		   # and assume that this holds throughout the entire result
		   # file (i.e., when this method is called for the next
		   # alignment, this will not be checked again)
		   if($seq1props{'length'} > $seq2props{'length'}) {
		       $self->{'_est_is_first'} = 0;
		   } else {
		       $self->{'_est_is_first'} = 1;
		   }
	       }
	       if($self->{'_est_is_first'}) {
		   $estseq = \%seq1props;
		   $genomseq = \%seq2props;
		   # if the EST is given first, A=4 selects the genomic
		   # seq for being reversed (reversing the EST is default)
		   $to_reverse = ($output_fmt == 4) ? $genomseq : $estseq;
	       } else {
		   $estseq = \%seq2props;
		   $genomseq = \%seq1props;
		   # if the EST is the second, A=4 does not change the
		   # seq being reversed (always the EST is reversed)
		   $to_reverse = $estseq;
	       }
	   }
	   if($hit_direction == -1) {
	       # we have to reverse the coordinates of one of both seqs
	       my $tmp = $to_reverse->{'start'};
	       $to_reverse->{'start'} =
		   $to_reverse->{'length'} - $to_reverse->{'end'} + 1;
	       $to_reverse->{'end'} = $to_reverse->{'length'} - $tmp + 1;
	   }
	   # create and initialize the exon object
	   my $exon = Bio::Tools::Sim4::Exon->new(
					    '-start' => $genomseq->{'start'},
					    '-end'   => $genomseq->{'end'},
					    '-strand' => $hit_direction);
	   if(exists($genomseq->{'seqname'})) {
	       $exon->seq_id($genomseq->{'seqname'});
	   } else {
	       # take filename stripped of path as fall back
	       my ($basename) = &File::Basename::fileparse($genomseq->{'filename'}, '\..*');
	       $exon->seq_id($basename);
	   }
	   $exon->feature1()->add_tag_value('filename',
					    $genomseq->{'filename'});
	   # feature1 is supposed to be initialized to a Similarity object,
           # but we provide a safety net
	   if($exon->feature1()->can('seqlength')) {
	       $exon->feature1()->seqlength($genomseq->{'length'});
	   } else {
	       $exon->feature1()->add_tag_value('SeqLength',
						$genomseq->{'length'});
	   }
	   # create and initialize the feature wrapping the 'hit' (the EST)
	   my $fea2 = Bio::SeqFeature::Similarity->new(
                                            '-start' => $estseq->{'start'},
					    '-end'   => $estseq->{'end'},
					    '-strand' => 0,
					    '-primary' => "aligning_EST");
	   if(exists($estseq->{'seqname'})) {
	       $fea2->seq_id($estseq->{'seqname'});
	   } else {
	       # take filename stripped of path as fall back
	       my ($basename) =
		   &File::Basename::fileparse($estseq->{'filename'}, '\..*');
	       $fea2->seq_id($basename);
	   }
	   $fea2->add_tag_value('filename', $estseq->{'filename'});
	   $fea2->seqlength($estseq->{'length'});
	   # store
	   $exon->est_hit($fea2);	   
	   # general properties
	   $exon->source_tag($self->analysis_method());
	   $exon->percentage_id($pctid);
	   $exon->score($exon->percentage_id());
	   # push onto array
	   push(@exons, $exon);
	   next; # back to while loop
       }
   }
   return @exons;
}

=head2 next_exonset

 Title   : next_exonset
 Usage   : $exonset = $sim4_result->parse_next_exonset;
           print "Exons start at ", $exonset->start(), 
                 "and end at ", $exonset->end(), "\n";
           foreach $exon ($exonset->sub_SeqFeature()) {
               # do something
           }
 Function: Parses the next alignment of the Sim4 result file and returns the
           set of exons as a container of features. The container is itself
           a Bio::SeqFeature::Generic object, with the Bio::Tools::Sim4::Exon
           objects as sub features. Start, end, and strand of the container
           will represent the total region covered by the exons of this set.

           See the documentation of parse_next_alignment() for further
           reference about parsing and how the information is stored.

 Example : 
 Returns : An Bio::SeqFeature::Generic object holding Bio::Tools::Sim4::Exon
           objects as sub features.
 Args    :

=cut

sub next_exonset {
    my $self = shift;
    my $exonset;

    # get the next array of exons
    my @exons = $self->parse_next_alignment();
    unless( @exons ) {
	return if eof($self->_fh);
	return $self->next_exonset;
    } 
    # create the container of exons as a feature object itself, with the
    # data of the first exon for initialization
    $exonset = Bio::SeqFeature::Generic->new('-start' => $exons[0]->start(),
					     '-end' => $exons[0]->end(),
					     '-strand' => $exons[0]->strand(),
					     '-primary' => "ExonSet");
    $exonset->source_tag($exons[0]->source_tag());
    $exonset->seq_id($exons[0]->seq_id());
    # now add all exons as sub features, with enabling EXPANsion of the region
    # covered in total
    foreach my $exon (@exons) {
	$exonset->add_sub_SeqFeature($exon, 'EXPAND');
    }
    return $exonset;
}

=head2 next_feature

 Title   : next_feature
 Usage   : while($exonset = $sim4->next_feature()) {
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

1;
