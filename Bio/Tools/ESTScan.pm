
#
# BioPerl module for Bio::Tools::ESTScan
#
# Cared for by Hilmar Lapp <hlapp@gmx.net>
#
# Copyright Hilmar Lapp
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Tools::ESTScan - Results of one ESTScan run

=head1 SYNOPSIS

   $estscan = Bio::Tools::ESTScan->new(-file => 'result.estscan');
   # filehandle:
   $estscan = Bio::Tools::ESTScan->new( -fh  => \*INPUT );

   # parse the results
   while($orf = $estscan->next_prediction()) {
       # $orf is an instance of Bio::Tools::Prediction::Exon and thus
       # implements SeqFeatureI
       
   }

   # essential if you gave a filename at initialization (otherwise the file
   # will stay open)
   $estscan->close();

=head1 DESCRIPTION

The ESTScan module provides a parser for ESTScan coding region prediction
output.

This module inherits off L<Bio::Tools::AnalysisResult> and therefore
implements the L<Bio::SeqAnalysisParserI> interface. 
See L<Bio::SeqAnalysisParserI>.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists.  Your participation is much appreciated.

  bioperl-l@bioperl.org          - General discussion
  http://bio.perl.org/MailList.html             - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via email
or the web:

  bioperl-bugs@bio.perl.org
  http://bio.perl.org/bioperl-bugs/

=head1 AUTHOR - Hilmar Lapp

Email hlapp@gmx.net (or hilmar.lapp@pharma.novartis.com)

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::Tools::ESTScan;
use vars qw(@ISA);
use strict;
use Symbol;

# Object preamble - inherits from Bio::Root::Object

use Bio::Root::Object;
use Bio::Tools::AnalysisResult;
use Bio::Tools::Prediction::Exon;

@ISA = qw(Bio::Tools::AnalysisResult);
# new() is inherited from Bio::Root::Object

# _initialize is where the heavy stuff will happen when new is called

# called by the inherited _initialize.
sub _initialize_state {
    my ($self,@args) = @_;

    # first call the inherited method!
    my $make = $self->SUPER::_initialize_state(@args);

    if(! $self->analysis_method()) {
	$self->analysis_method('ESTScan');
    }
}

=head2 analysis_method

 Usage     : $estscan->analysis_method();
 Purpose   : Inherited method. Overridden to ensure that the name matches
             /estscan/i.
 Returns   : String
 Argument  : n/a

=cut

#-------------
sub analysis_method { 
#-------------
    my ($self, $method) = @_;  
    if($method && ($method !~ /estscan/i)) {
	$self->throw("method $method not supported in " . ref($self));
    }
    return $self->SUPER::analysis_method($method);
}

=head2 next_feature

 Title   : next_feature
 Usage   : while($orf = $estscan->next_feature()) {
                  # do something
           }
 Function: Returns the next gene structure prediction of the ESTScan result
           file. Call this method repeatedly until FALSE is returned.

           The returned object is actually a SeqFeatureI implementing object.
           This method is required for classes implementing the
           SeqAnalysisParserI interface, and is merely an alias for 
           next_prediction() at present.

 Example :
 Returns : A Bio::Tools::Prediction::Gene object.
 Args    :

=cut

sub next_feature {
    my ($self,@args) = @_;
    # even though next_prediction doesn't expect any args (and this method
    # does neither), we pass on args in order to be prepared if this changes
    # ever
    return $self->next_prediction(@args);
}

=head2 next_prediction

 Title   : next_prediction
 Usage   : while($gene = $estscan->next_prediction()) {
                  # do something
           }
 Function: Returns the next gene structure prediction of the ESTScan result
           file. Call this method repeatedly until FALSE is returned.

           So far, this method DOES NOT work for reverse strand predictions,
           even though the code looks like.
 Example :
 Returns : A Bio::Tools::Prediction::Gene object.
 Args    :

=cut

sub next_prediction {
    my ($self) = @_;
    my ($gene, $seq, $predobj);
    my $numins = 0;

    # predictions are in the format of FASTA sequences and can be parsed one
    # at a time
    $seq = $self->_fasta_stream()->next_seq();
    return unless $seq;
    # there is a new prediction
    $gene = Bio::Tools::Prediction::Gene->new('-primary' => "ORFprediction",
                                              '-source' => "ESTScan");
    # score starts the description
    $seq->desc() =~ /^([\d.]+)\s*(.*)/ or
	$self->throw("unexpected format of description: no score in " .
		     $seq->desc());
    $gene->score($1);
    $seq->desc($2);
    # strand may end the description
    if($seq->desc() =~ /(.*)minus strand$/) {
	my $desc = $1;
	$desc =~ s/;\s+$//;
	$seq->desc($desc);
	$gene->strand(-1);
    } else {
	$gene->strand(1);
    }
    # check for the format: default or 'all-in-one' (option -a)
    if($seq->desc() =~ /^(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s*(.*)/) {
	# default format
	$seq->desc($5);
	$predobj = Bio::Tools::Prediction::Exon->new('-source' => "ESTScan",
						     '-start' => $3,
						     '-end' => $4);
	$predobj->strand($gene->strand());
	$predobj->score($gene->score()); # FIXME or $1, or $2 ?
	$predobj->primary_tag("InternalExon");
	$predobj->seqname($seq->display_id());
	# add to gene structure object
	$gene->add_exon($predobj);
	# add predicted CDS
	$gene->predicted_cds($seq);
	if($gene->strand() == -1) {
	    $self->warn("reverse strand ORF, but unable to reverse coordinates!");
	}
    } else {
	#
	# All-in-one format (hopefully). This encodes the following information
	# into the sequence:
	# 1) untranslated regions: stretches of lower-case letters
	# 2) translated regions: stretches of upper-case letters
	# 3) insertions in the translated regions: capital X
	# 4) deletions in the translated regions: a single lower-case letter
	#
	# if reverse strand ORF, save a lot of hassle by reversing the sequence
	if($gene->strand() == -1) {
	    $seq = $seq->revcom();
	}
	my $seqstr = $seq->seq();
	while($seqstr =~ /^([a-z]*)([A-Z].*)$/) {
	    # leading 5'UTR
	    my $utr5 = $1;
	    # exon + 3'UTR
	    my $exonseq = $2;
	    # strip 3'UTR and following exons
	    if($exonseq =~ s/([a-z]{2,}.*)$//) {
		$seqstr = $1;
	    } else {
		$seqstr = "";
	    }
	    # start: take care of yielding the absolute coordinate
	    my $start = CORE::length($utr5) + 1;
	    if($predobj) {
		$start += $predobj->end() + $numins;
	    }
	    # for the end coordinate, we need to subtract the insertions
	    my $cds = $exonseq;
	    $cds =~ s/[X]//g;
	    my $end = $start + CORE::length($cds) - 1;
	    # construct next exon object
	    $predobj = Bio::Tools::Prediction::Exon->new('-start' => $start,
							 '-end' => $end);
	    $predobj->source_tag("ESTScan");
	    $predobj->primary_tag("InternalExon");
	    $predobj->seqname($seq->display_id());
	    $predobj->strand($gene->strand());
	    $predobj->score($gene->score());
	    # add the exon to the gene structure object
	    $gene->add_exon($predobj);
	    # add the predicted CDS
	    $cds = $exonseq;
	    $cds =~ s/[a-z]//g; # remove the deletions, but keep the insertions
	    $cds = Bio::PrimarySeq->new('-seq' => $cds,
					'-display_id' => $seq->display_id(),
					'-desc' => $seq->desc(),
					'-moltype' => "dna");
	    # in case of a reverse strand prediction, we reversed the sequence
	    # initially, so reverse the predicted CDS back
	    $cds = $cds->revcom() if($gene->strand() == -1);
	    # only store the first one in the overall prediction
	    $gene->predicted_cds($cds) unless $gene->predicted_cds();
	    $predobj->predicted_cds($cds);
	    # add the predicted insertions and deletions as subfeatures
	    # of the exon
	    my $fea = undef;
	    while($exonseq =~ /([a-zX])/g) {
		my $indel = $1;
		# start and end: look after previous feature (because of the
		# 1-based indexing of start() as opposed to 0-based indexing
		# of index(), we don't need to add 1 in order to be ahead of
		# the previous match)
		$start = ($fea) ?
		    $fea->start()+$numins : $predobj->start()+$numins;
		$start = index($seq->seq(), $indel, $start) + 1 - $numins;
		$fea = Bio::SeqFeature::Generic->new('-start' => $start,
						     '-end' => $start);
		$fea->source_tag("ESTScan");
		$fea->seqname($seq->display_id());
		$fea->strand($predobj->strand());
		if($indel eq 'X') {
		    # an insertion (depends on viewpoint: to get the 'real'
		    # CDS, a base has to be inserted, i.e., the HMMER model
		    # inserted a base; however, the sequencing process deleted
		    # a base that was there).
		    $fea->primary_tag("insertion");
		    # we need to count insertions because these are left out
		    # of any coordinates saved in the objects (which is correct
		    # because insertions change the original sequence, so
		    # coordinates wouldn't match)
		    $numins++;
		} else {
		    # a deletion (depends on viewpoint: to get the 'real'
		    # CDS, a base has to be deleted, i.e., the HMMER model
		    # deleted a base; however, the sequencing process inserted
		    # a base that wasn't there).
		    $fea->primary_tag("deletion");
		    $fea->add_tag_value('base', $indel);
		}
		$predobj->add_sub_SeqFeature($fea);
	    }
	}
    }
    
    return $gene;
}

=head2 close

 Title   : close
 Usage   : $result->close()
 Function: Closes the file handle associated with this result file.
           Inherited method, overridden.
 Example :
 Returns :
 Args    :

=cut

sub close {
   my ($self, @args) = @_;

   delete($self->{'_fastastream'});
   $self->SUPER::close(@args);
}

=head2 _fasta_stream

 Title   : _fasta_stream
 Usage   : $result->_fasta_stream()
 Function: Gets/Sets the FASTA sequence IO stream for reading the contents of
           the file associated with this MZEF result object.

           If called for the first time, creates the stream from the filehandle
           if necessary.
 Example :
 Returns :
 Args    :

=cut

sub _fasta_stream {
    my ($self, $stream) = @_;
    
    if($stream || (! exists($self->{'_fastastream'}))) {
	if(! $stream) {
	    $stream = Bio::SeqIO->new('-fh' => $self->_filehandle(),
				      '-format' => "fasta");
	}
	$self->{'_fastastream'} = $stream;
    }
    return $self->{'_fastastream'};
}

1;

