#
# BioPerl module for Bio::Tools::Genemark
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Mark Fiers <hlapp@gmx.net>
#
# Copyright Hilmar Lapp, Mark Fiers
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Tools::Genemark - Results of one Genemark run

=head1 SYNOPSIS

   $Genemark = Bio::Tools::Genemark->new(-file => 'result.Genemark');
   # filehandle:
   $Genemark = Bio::Tools::Genemark->new( -fh  => \*INPUT );

   # parse the results
   # note: this class is-a Bio::Tools::AnalysisResult which implements
   # Bio::SeqAnalysisParserI, i.e., $Genemark->next_feature() is the same
   while($gene = $Genemark->next_prediction()) {
       # $gene is an instance of Bio::Tools::Prediction::Gene, which inherits
       # off Bio::SeqFeature::Gene::Transcript.
       #
       # $gene->exons() returns an array of
       # Bio::Tools::Prediction::Exon objects
       # all exons:
       @exon_arr = $gene->exons();

       # initial exons only
       @init_exons = $gene->exons('Initial');
       # internal exons only
       @intrl_exons = $gene->exons('Internal');
       # terminal exons only
       @term_exons = $gene->exons('Terminal');
       # singleton exons:
       ($single_exon) = $gene->exons();
   }

   # essential if you gave a filename at initialization (otherwise the file
   # will stay open)
   $Genemark->close();

=head1 DESCRIPTION

The Genemark module provides a parser for Genemark gene structure
prediction output. It parses one gene prediction into a
Bio::SeqFeature::Gene::Transcript- derived object.

This module has been developed around genemark.hmm for eukaryots v2.2a
and will probably not work with other versions.


This module also implements the Bio::SeqAnalysisParserI interface, and
thus can be used wherever such an object fits. See
L<Bio::SeqAnalysisParserI>.

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

=head1 AUTHOR - Hilmar Lapp, Mark Fiers

Email hlapp@gmx.net
      m.w.e.j.fiers@plant.wag-ur.nl

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Tools::Genemark;
use strict;
use Symbol;

use Bio::Root::Root;
use Bio::Tools::Prediction::Gene;
use Bio::Tools::Prediction::Exon;
use Bio::Seq;
use Bio::Factory::FTLocationFactory;

use base qw(Bio::Tools::AnalysisResult);

=head2 new

 Title   : new
 Usage   : my $obj = Bio::Tools::Genemark->new();
 Function: Builds a new Bio::Tools::Genemark object
 Returns : an instance of Bio::Tools::Genemark
 Args    : seqname


=cut

sub new {
  my($class,@args) = @_;

  my $self = $class->SUPER::new(@args);

  my ($seqname) = $self->_rearrange([qw(SEQNAME)], @args);

  # hardwire seq_id when creating gene and exon objects
  $self->_seqname($seqname) if defined($seqname);

  return $self;
}

sub _initialize_state {
    my ($self,@args) = @_;

    # first call the inherited method!
    $self->SUPER::_initialize_state(@args);

    # our private state variables
    $self->{'_preds_parsed'} = 0;
    $self->{'_has_cds'} = 0;
    # array of pre-parsed predictions
    $self->{'_preds'} = [];
    # seq stack
    $self->{'_seqstack'} = [];
}

=head2 analysis_method

 Usage     : $Genemark->analysis_method();
 Purpose   : Inherited method. Overridden to ensure that the name matches
             /GeneMark.hmm/i.
 Returns   : String
 Argument  : n/a

=cut

#-------------
sub analysis_method {
#-------------
    my ($self, $method) = @_;
    if($method && ($method !~ /Genemark\.hmm/i)) {
	$self->throw("method $method not supported in " . ref($self));
    }
    return $self->SUPER::analysis_method($method);
}

=head2 next_feature

 Title   : next_feature
 Usage   : while($gene = $Genemark->next_feature()) {
                  # do something
           }
 Function: Returns the next gene structure prediction of the Genemark result
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
 Usage   : while($gene = $Genemark->next_prediction()) {
                  # do something
           }
 Function: Returns the next gene structure prediction of the Genemark result
           file. Call this method repeatedly until FALSE is returned.

 Example :
 Returns : A Bio::Tools::Prediction::Gene object.
 Args    :

=cut

sub next_prediction {
    my ($self) = @_;
    my $gene;

    # if the prediction section hasn't been parsed yet, we do this now
    $self->_parse_predictions() unless $self->_predictions_parsed();

    # get next gene structure
    $gene = $self->_prediction();

    return $gene;
}

=head2 _parse_predictions

 Title   : _parse_predictions()
 Usage   : $obj->_parse_predictions()
 Function: Parses the prediction section. Automatically called by
           next_prediction() if not yet done.
 Example :
 Returns :

=cut

sub _parse_predictions {
    my ($self) = @_;
    my %exontags = ('Initial' => 'Initial',
		    'Internal' => 'Internal',
		    'Terminal' => 'Terminal',
		    'Single' => '',
		    '_na_' => '');
    my $exontag;
    my $gene;
    my $exontype;
    my $current_gene_no = -1;

    # The prediction report does not contain a sequence identifier
    # (at least the prokaryotic version doesn't)
    my $seqname = $self->_seqname();
    
    while(defined($_ = $self->_readline())) {

	if( (/^\s*(\d+)\s+(\d+)/) || (/^\s*(\d+)\s+[\+\-]/)) {

	    #  this is an exon, Genemark doesn't predict anything else
	    # $prednr corresponds to geneno.
	    my $prednr = $1;

	    #exon no:
	    my $signalnr = 0;
	    if ($2) { my $signalnr = $2; } # used in tag: exon_no
	
	    # split into fields
	    chomp();
	    my @flds = split(' ', $_);

	    # create the feature (an exon) object
	    my $predobj = Bio::Tools::Prediction::Exon->new();

		
	    # define info depending on it being eu- or prokaryot
	    my ($start, $end, $orientation, $prediction_source);

	    if ($self->analysis_method() =~ /PROKARYOTIC/i) {
	        $prediction_source = "Genemark.hmm.pro";
	       	$orientation = ($flds[1] eq '+') ? 1 : -1;
	        ($start, $end) = @flds[(2,3)];
		$exontag = "_na_";

	    } else {		
	        $prediction_source = "Genemark.hmm.eu";
	       	$orientation = ($flds[2] eq '+') ? 1 : -1;
	        ($start, $end) = @flds[(4,5)];
		$exontag = $flds[3];
	    }

            # instatiate a location object via
            # Bio::Factory::FTLocationFactory (to handle
            # inexact coordinates)
            my $location_string = join('..', $start, $end);
            my $location_factory = Bio::Factory::FTLocationFactory->new();
            my $location_obj = $location_factory->from_string($location_string);
            $predobj->location($location_obj);

            #store the data in the exon object
            $predobj->source_tag($prediction_source);
	    $predobj->strand($orientation);

	    $predobj->primary_tag($exontags{$exontag} . "Exon");

	    $predobj->add_tag_value('exon_no',"$signalnr") if ($signalnr);

    	    $predobj->is_coding(1);

            $predobj->seq_id($seqname) if (defined($seqname) && ($seqname ne 'unknown'));
		
	    # frame calculation as in the genscan module
	    # is to be implemented...
	
	    #If the $prednr is not equal to the current gene, we
	    #need to make a new gene and close the old one
	    if($prednr != $current_gene_no) {
 	        # a new gene, store the old one if it exists
		if (defined ($gene)) {
		    $gene->seq_id($seqname);
		    $gene = undef ;
		}
		#and make a new one
		$gene = Bio::Tools::Prediction::Gene->new
		    (
		     '-primary' => "GenePrediction$prednr",
		     '-source' => $prediction_source);
                $self->_add_prediction($gene);		
		$current_gene_no = $prednr;
                $gene->seq_id($seqname) if (defined($seqname) && ($seqname ne 'unknown'));
	    }
	
	    # Add the exon to the gene
	    $gene->add_exon($predobj, ($exontag eq "_na_" ?
				       undef : $exontags{$exontag}));

	}

	if(/^(Genemark\.hmm\s*[PROKARYOTIC]*)\s+\(Version (.*)\)$/i) {
	    $self->analysis_method($1);

	    my $gm_version = $2;

	    $self->analysis_method_version($gm_version);
	    next;
	}

       #Matrix file for eukaryot version
       if (/^Matrices file:\s+(\S+)?/i)  {
	    $self->analysis_subject($1);
	    # since the line after the matrix file is always the date
	    # (in the output file's I have seen!) extract and store this
	    # here
	     if (defined(my $_date = $self->_readline())) {
	         chomp ($_date);
	     	 $self->analysis_date($_date);
	     }
	}			
	
        #Matrix file for prokaryot version
       if (/^Model file name:\s+(\S+)/) {
	    $self->analysis_subject($1);
	    # since the line after the matrix file is always the date
	    # (in the output file's I have seen!) extract and store this
	    # here
	    my $_date = $self->_readline() ;
	    if (defined($_date = $self->_readline())) {
	         chomp ($_date);
	     	 $self->analysis_date($_date);
	     }
	}
	
	if(/^Sequence[ file]? name:\s+(.+)\s*$/i) {
	    $seqname = $1;
	    #    $self->analysis_subject($seqname);
	    next;
	}
	

	/^>/ && do {		
    	    $self->_pushback($_);

	    # section of predicted aa sequences on recognition
	    # of a fasta start, read all sequences and find the
	    # appropriate gene
            while (1) {
	       my ($aa_id, $seq) = $self->_read_fasta_seq();
	       last unless ($aa_id);

	       #now parse through the predictions to add the pred. protein
	       FINDPRED: foreach my $gene (@{$self->{'_preds'}}) {
	            $gene->primary_tag() =~ /[^0-9]([0-9]+)$/;
		    my $geneno = $1;
		    if ($aa_id =~ /\|gene.$geneno\|/) {
		          #print "x SEQ : \n $seq \nXXXX\n";
  			  my $seqobj = Bio::Seq->new('-seq' => $seq,
	                     		             '-display_id' => $aa_id,
					              '-alphabet' => "protein");
			$gene->predicted_protein($seqobj);
			last FINDPRED;
		    }	

	       }
           }				

 	   last;
	};
    }

    # if the analysis query object contains a ref to a Seq of PrimarySeq
    # object, then extract the predicted sequences and add it to the gene
    # object.
    if (defined $self->analysis_query()) {
        my $orig_seq = $self->analysis_query();
        FINDPREDSEQ: foreach my $gene (@{$self->{'_preds'}}) {
	   my $predseq = "";
	   foreach my $exon ($gene->exons()) {
		#print $exon->start() . " " . $exon->end () . "\n";
		$predseq .= $orig_seq->subseq($exon->start(), $exon->end());
	   }

	   my $seqobj = Bio::PrimarySeq->new('-seq' => $predseq,
	                     		     '-display_id' => "transl");
	   $gene->predicted_cds($seqobj);
	}
    }


    $self->_predictions_parsed(1);
}

=head2 _prediction

 Title   : _prediction()
 Usage   : $gene = $obj->_prediction()
 Function: internal
 Example :
 Returns :

=cut

sub _prediction {
    my ($self) = @_;

    return unless(exists($self->{'_preds'}) && @{$self->{'_preds'}});
    return shift(@{$self->{'_preds'}});
}

=head2 _add_prediction

 Title   : _add_prediction()
 Usage   : $obj->_add_prediction($gene)
 Function: internal
 Example :
 Returns :

=cut

sub _add_prediction {
    my ($self, $gene) = @_;

    if(! exists($self->{'_preds'})) {
	$self->{'_preds'} = [];
    }
    push(@{$self->{'_preds'}}, $gene);
}

=head2 _predictions_parsed

 Title   : _predictions_parsed
 Usage   : $obj->_predictions_parsed
 Function: internal
 Example :
 Returns : TRUE or FALSE

=cut

sub _predictions_parsed {
    my ($self, $val) = @_;

    $self->{'_preds_parsed'} = $val if $val;
    if(! exists($self->{'_preds_parsed'})) {
	$self->{'_preds_parsed'} = 0;
    }
    return $self->{'_preds_parsed'};
}

=head2 _has_cds

 Title   : _has_cds()
 Usage   : $obj->_has_cds()
 Function: Whether or not the result contains the predicted CDSs, too.
 Example :
 Returns : TRUE or FALSE

=cut

sub _has_cds {
    my ($self, $val) = @_;

    $self->{'_has_cds'} = $val if $val;
    if(! exists($self->{'_has_cds'})) {
	$self->{'_has_cds'} = 0;
    }
    return $self->{'_has_cds'};
}

=head2 _read_fasta_seq

 Title   : _read_fasta_seq()
 Usage   : ($id,$seqstr) = $obj->_read_fasta_seq();
 Function: Simple but specialised FASTA format sequence reader. Uses
           $self->_readline() to retrieve input, and is able to strip off
           the traling description lines.
 Example :
 Returns : An array of two elements.

=cut

sub _read_fasta_seq {
    my ($self) = @_;
    my ($id, $seq);
    local $/ = ">";

    return 0 unless (my $entry = $self->_readline());

    $entry =~ s/^>//;
    # complete the entry if the first line came from a pushback buffer
    while(! ($entry =~ />$/)) {
	last unless ($_ = $self->_readline());
	$entry .= $_;
    }

    # delete everything onwards from an new fasta start (>)
    $entry =~ s/\n>.*$//s;
    # id and sequence

    if($entry =~ s/^(.+)\n//) {
	$id = $1;
	$id =~ s/ /_/g;
	$seq = $entry;
	$seq =~ s/\s//g;	
	#print "\n@@ $id \n@@ $seq \n##\n";
    } else {
	$self->throw("Can't parse Genemark predicted sequence entry");
    }
    $seq =~ s/\s//g; # Remove whitespace
    return ($id, $seq);
}

=head2 _seqname

 Title   : _seqname
 Usage   : $obj->_seqname($seqname)
 Function: internal
 Example :
 Returns : String

=cut

sub _seqname {
    my ($self, $val) = @_;

    $self->{'_seqname'} = $val if $val;
    if(! exists($self->{'_seqname'})) {
        $self->{'_seqname'} = 'unknown';
    }
    return $self->{'_seqname'};
}

1;

