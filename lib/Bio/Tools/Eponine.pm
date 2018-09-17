#
# BioPerl module for Bio::Tools::Eponine
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Tania Oh <gisoht@nus.edu.sg>
#
# Copyright Tania Oh 
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Tools::Eponine - Results of one Eponine run

=head1 SYNOPSIS

 use Bio::Tools::Eponine;
 use strict;
    my $seq = "/data/seq.fa";
    my $threshold  = "0.999";
    my @params = ( -seq => $seq,
                   -threshold => $threshold);

   my $factory = Bio::Tools::Run::Eponine->new(@params);
     # run eponine against fasta 
        my $r = $factory->run_eponine($seq);
        my $parser = Bio::Tools::Eponine->new($r);

       while (my $feat = $parser->next_prediction){
                #$feat contains array of SeqFeature
               foreach my $orf($feat) {
                   print $orf->seq_id. "\n";
               }
       }

=head1 DESCRIPTION

Parser for Eponine, a probabilistic transcription start site detector
optimized for mammalian genomic sequence. This module inherits off
Bio::Tools::AnalysisResult and therefore implements 
Bio::SeqAnalysisParserI (see L<Bio::Tools::AnalysisResult> and
L<Bio::SeqAnalysisParserI>).

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

=head1 AUTHOR - Tania Oh 

E<lt>gisoht-at-nus.edu.sgE<gt>

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Tools::Eponine;
use strict;

use Bio::Tools::Prediction::Gene;
use Bio::Tools::Prediction::Exon;

use base qw(Bio::Tools::AnalysisResult);

sub _initialize_state {
    my($self,@args) = @_;

    # first call the inherited method!
    my $make = $self->SUPER::_initialize_state(@args);

    # handle our own parameters

    # our private state variables
    $self->{'_preds_parsed'} = 0;
    #array of Bio::SeqFeatures
    $self->{'_flist'} =[];
}

=head2 analysis_method

 Usage     : $mzef->analysis_method();
 Purpose   : Inherited method. Overridden to ensure that the name matches
             /mzef/i.
 Returns   : String
 Argument  : n/a

=cut

#-------------
sub analysis_method { 
#-------------
    my ($self, $method) = @_;  
    if($method && ($method !~ /epo/i)) {
	$self->throw("method $method not supported in " . ref($self));
    }
    return $self->SUPER::analysis_method($method);
}

=head2 next_feature

 Title   : next_feature
 Usage   : while($gene = $mzef->next_feature()) {
                  # do something
           }
 Function: Returns the next gene structure prediction of the MZEF result
           file. Call this method repeatedly until FALSE is returned.

           The returned object is actually a SeqFeatureI implementing object.
           This method is required for classes implementing the
           SeqAnalysisParserI interface, and is merely an alias for 
           next_prediction() at present.

           Note that with the present version of MZEF there will only be one
           object returned, because MZEF does not predict individual genes
           but just potential internal exons.
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
 Usage   : while($gene = $mzef->next_prediction()) {
                  # do something
           }
 Function: Returns the next gene structure prediction of the MZEF result
           file. Call this method repeatedly until FALSE is returned.

           Note that with the present version of MZEF there will only be one
           object returned, because MZEF does not predict individual genes
           but just potential internal exons.
 Example :
 Returns : A Bio::Tools::Prediction::Gene object.
 Args    :

=cut

sub next_prediction {
    my ($self) = @_;
    my $gene;

    # if the prediction section hasn't been parsed yet, we do this now
    $self->_parse_predictions() unless $self->_predictions_parsed();

    # return the next gene structure (transcript)
    return $self->_prediction();
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

    while(defined($_ = $self->_readline())) {
        if (! /^\#/){ #ignore introductory lines

	      my @element = split;
	      my (%feature);
	      $feature {name} = $element[0];
	      $feature {score} = $element[5];
	      $feature {start} = $element[3];
	      $feature {end} = $element[4];
	      $feature {strand} = $element[6];
	      $feature {source}= 'Eponine';
	      $feature {primary}= 'TSS';
	      $feature {program} = 'eponine-scan';
	      $feature {program_version} = '2';
            
	      $self->create_feature(\%feature);
	            next;

	}
    }
    $self->_predictions_parsed(1);
}

=head2 create_feature

    Title   :   create_feature
    Usage   :   obj->create_feature($feature)
    Function:   Returns an array of features
    Returns :   Returns an array of features
    Args    :   none

=cut

sub create_feature {
    my ($self, $feat) = @_;
     #create and fill Bio::EnsEMBL::Seqfeature object

      my $tss = Bio::SeqFeature::Generic->new
                    (   -seq_id  => $feat->{'name'},
                        -start   => $feat->{'start'},
                        -end     => $feat->{'end'},
                        -strand  => $feat->{'strand'},
            		-score   => $feat->{'score'},
                        -source_tag  => $feat->{'source'},
		        -primary_tag => $feat->{'primary'});

		

  if ($tss) {
         # add to _flist
      push(@{$self->{'_flist'}}, $tss);
   }

   #print $tss->gff_string;
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

    return unless(exists($self->{'_flist'}) && @{$self->{'_flist'}});
    return shift(@{$self->{'_flist'}});
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
    # array of pre-parsed predictions
    if(! exists($self->{'_preds_parsed'})) {
	$self->{'_preds_parsed'} = 0;
    }
    return $self->{'_preds_parsed'};
}


1;
