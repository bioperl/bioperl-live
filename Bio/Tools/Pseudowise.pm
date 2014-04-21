# BioPerl module for Bio::Tools::Pseudowise
#
# 
# Copyright Jason Stajich, Fugu Team 
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Tools::Pseudowise - Results of one Pseudowise run

=head1 SYNOPSIS

  use Bio::Tools::Pseudowise;

  my $parser = Bio::Tools::Pseudowise->new(-file=>"pw.out");
  while(my $feat = $parser->next_result){
      push @feat, $feat;
  }

=head1 DESCRIPTION

Pseudowise is a pseudogene prediction program written by Ewan Birney
as part of the Wise Package. This module is the parser for the output
of the program.

http://www.sanger.ac.uk/software/wise2

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

=head1 AUTHOR - Jason Stajich

Previous committed by the Fugu Team 

Re-written by Jason Stajich jason-at-bioperl-dot-org

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Tools::Pseudowise;
use strict;
use Symbol;

use Bio::Root::Root;
use Bio::SeqFeature::Generic;
use Bio::SeqFeature::Gene::Exon;
use Bio::SeqFeature::FeaturePair;
use Bio::SeqFeature::Gene::Transcript;
use Bio::SeqFeature::Gene::GeneStructure;

use base qw(Bio::Tools::AnalysisResult);

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

 Usage     : $pseudowise->analysis_method();
 Purpose   : Inherited method. Overridden to ensure that the name matches
             /pseudowise/i.
 Returns   : String
 Argument  : n/a

=cut

#-------------
sub analysis_method { 
#-------------
    my ($self, $method) = @_;  
    if($method && ($method !~ /pseudowise/i)) {
	$self->throw("method $method not supported in " . ref($self));
    }
    return $self->SUPER::analysis_method($method);
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

See Also  L<Bio::SeqFeatureI>

=cut

sub next_feature {
    return shift->next_prediction(@_);
}


=head2 next_prediction

 Title   : next_prediction
 Usage   : while($gene = $pseudowise->next_prediction()) {
                  # do something
           }
 Function: Returns the gene of the Pseudowise result
           file. Call this method repeatedly until FALSE is returned.

 Example :
 Returns : a Bio::SeqFeature::Generic 
 Args    : none

See Also L<Bio::SeqFeature::Generic>

=cut

sub next_prediction {
    my ($self) = @_;
    # if the prediction section hasn't been parsed yet, we do this now
    $self->_parse_predictions unless $self->_predictions_parsed;

    # get next gene structure
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
    my $gene;
    my @genes;

    local $/= "\n";
    local($_);
    my %tags;
    while (defined( $_ = $self->_readline)){ 
	if( /^(Total codons|\S+)\s+:\s+(\S+)/ ) {
	    $tags{$1} = $2;
	} elsif(m!^//! ) {
	    if( $gene ) {
		$gene = undef;
		%tags = ();
	    }
	} elsif (/Gene\s+(\d+)\s*$/i) {
	    $gene = Bio::SeqFeature::Generic->new 
		( -primary => 'pseudogene',
		  -source  => 'pseudowise',
		  -tag     => \%tags);
	    push @genes, $gene;
	} elsif( /Gene\s+(\d+)\s+(\d+)/i ) {
	    if( $1 < $2 ) {
		$gene->start($1);
		$gene->end($2);
		$gene->strand(1);
	    } else {
		$gene->start($2);
		$gene->end($1);
		$gene->strand(-1);
	    }
	} elsif (/Exon\s+(\d+)\s+(\d+)\s+phase\s+(\S+)/i) {
	    my ($s,$e,$st) = ($1,$2,1);
	    if( $s > $e) {
		($s,$e,$st)=($e,$s,-1);
	    }
	    my $exon = Bio::SeqFeature::Generic->new 
		( -start   => $s,
		  -end     => $e,
		  -strand  => $st,
		  -primary => 'exon',
		  -source  => 'pseudowise',
		  -tag     => {'frame'  => $3});
	    $gene->add_sub_SeqFeature($exon);
	} 
    }
    $self->_add_prediction(\@genes);
    $self->_predictions_parsed(1);
}

=head1 _prediction

 Title   : _prediction()
 Usage   : $gene = $obj->_prediction()
 Function: internal
 Example :
 Returns : 

=cut

sub _prediction {
    my ($self) = @_;
    return shift(@{$self->{'_preds'} || []});
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
    $self->{'_preds'} ||= [];

    if( ref($gene) =~ /ARRAY/ ) {
	push(@{$self->{'_preds'}}, @$gene);
    } else {
	push(@{$self->{'_preds'}}, $gene);
    }
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

1;
