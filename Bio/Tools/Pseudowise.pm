# BioPerl module for Bio::Tools::Pseudowise
#
# Copyright Fugu Team 
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

Pseudowise is a pseudogene prediction program written by Ewan Birney as part of the 
Wise Package. This module is the parser for the output of the program.

http://www.sanger.ac.uk/software/wise2

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
  http://bugzilla.bioperl.org/

=head1 AUTHOR - Fugu Team 

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Tools::Pseudowise;
use vars qw(@ISA);
use strict;
use Symbol;

use Bio::Root::Root;
use Bio::Tools::AnalysisResult;
use Bio::SeqFeature::Generic;
use Bio::SeqFeature::Gene::Exon;
use Bio::Tools::Run::WrapperBase;
use Bio::SeqFeature::FeaturePair;
use Bio::SeqFeature::Gene::Transcript;
use Bio::SeqFeature::Gene::GeneStructure;

@ISA = qw(Bio::Tools::AnalysisResult);

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

=head2 next_prediction

 Title   : next_prediction
 Usage   : while($gene = $pseudowise->next_prediction()) {
                  # do something
           }
 Function: Returns the gene of the Pseudowise result
           file. Call this method repeatedly until FALSE is returned.

 Example :
 Returns : a Bio::SeqFeature::Generic 
 Args    :

=cut

sub next_prediction {
    my ($self,$filehandle) = @_;
    my $gene;

    # if the prediction section hasn't been parsed yet, we do this now
    $self->_parse_predictions($filehandle) unless $self->_predictions_parsed();

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
    my ($self, $filehandle) = @_;
    my $gene;
    my @genes;
    #The big parsing loop - parses exons and predicted peptides
    while (<$filehandle>)
    {
      if (/Gene/i)
      {
        $gene = new Bio::SeqFeature::Generic (
                    -primary => 'pseudogene',
                    -source => 'pseudowise');
        push @genes, $gene;

        while(<$filehandle>) {
          my @gene_elements = split;
          my $no = scalar(@gene_elements);
          if ((/Gene/i) && $no == 3) {
            my @element = split;
            my $no = scalar(@element);
            my $gene_start = $element[1];
            my $gene_end = $element[2];
            $gene->start($gene_start);
            $gene->end($gene_end);
          }
          elsif (/Exon/i) {
            my @element = split;
            my $no = scalar(@element);
            my $exon_start = $element[1];
            my $exon_end = $element[2];
            my $exon_phase = $element[4];
            my $exon = new Bio::SeqFeature::Generic (
                           -start => $exon_start,
                           -end => $exon_end,
                           -primary => 'exon',
                           -source => 'pseudowise',
                           -frame  => $exon_phase);
            $gene->add_sub_SeqFeature($exon);
          }
          elsif ((/Gene/i) && $no != 3) {
            $gene = new Bio::SeqFeature::Generic (
                        -primary => 'pseudogene',
                        -source => 'pseudowise');
            push @genes, $gene;
          }
        }
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

    return undef unless(exists($self->{'_preds'}) && @{$self->{'_preds'}});
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

1;
