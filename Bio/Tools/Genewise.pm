# BioPerl module for Bio::Tools::Genewise
#
# Copyright Fugu Team 
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Tools::Genewise - Results of one Genewise run

=head1 SYNOPSIS

=head1 DESCRIPTION

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

=head1 AUTHOR - Fugu Team 

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Tools::Genewise;
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

 Usage     : $genewise->analysis_method();
 Purpose   : Inherited method. Overridden to ensure that the name matches
             /genewise/i.
 Returns   : String
 Argument  : n/a

=cut

#-------------
sub analysis_method { 
#-------------
    my ($self, $method) = @_;  
    if($method && ($method !~ /genewise/i)) {
	$self->throw("method $method not supported in " . ref($self));
    }
    return $self->SUPER::analysis_method($method);
}

=head2 next_prediction

 Title   : next_prediction
 Usage   : while($gene = $genewise->next_prediction()) {
                  # do something
           }
 Function: Returns the gene structure prediction of the Genewise result
           file. Call this method repeatedly until FALSE is returned.

 Example :
 Returns : a Bio::SeqFeature::Gene::GeneStructure object 
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

=head2 _get_strand

 Title   : _get_strand
 Usage   : $obj->_get_strand
 Function: takes start and end values, swap them if start>end and returns end
 Example :
 Returns :$start,$end,$strand

=cut

sub _get_strand {
  my ($self,$start,$end) = @_;
  $start || $self->throw("Need a start");
  $end   || $self->throw("Need an end");
  my $strand;
  if ($start > $end) {
    my $tmp = $start;
    $start = $end;
    $end = $tmp;
    $strand = -1;
  }
  else {
    $strand = 1;
  }
  return ($start,$end,$strand);
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
    my $genes = new Bio::SeqFeature::Gene::GeneStructure ;
    my $transcript = new Bio::SeqFeature::Gene::Transcript ;
    $/ = "//";
    my $score;
    while (<$filehandle>) {

        ($score) = $_=~m/Score\s+(\d+[\.][\d]+)/ unless defined $score;

        next unless /Gene\s+\d+\n/;

        #grab exon + supporting feature info
        my @exons =$_=~ m/(Exon .+\s+Supporting .+)/g;
        my $nbr = 1;;

        #loop through each exon-supporting feature pair
        foreach my $e (@exons){
          my ($e_start,$e_end,$phase) = $_=~ m/Exon\s+(\d+)\s+(\d+)\s+phase\s+(\d+)/;
          my $e_strand;
          ($e_start,$e_end,$e_strand) = $self->_get_strand($e_start,$e_end);
          $transcript->strand($e_strand) unless $transcript->strand != 0;

          my $exon = new Bio::SeqFeature::Gene::Exon (-seqname=>"Exon $nbr", -start=>$e_start, -end=>$e_end, -strand=>$e_strand);
          $nbr++;
          $exon->add_tag_value('phase',$phase);

          my ($geno_start,$geno_end,$prot_start,$prot_end) = $_=~m/Supporting\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)/;
          my $prot_strand;
          ($prot_start,$prot_end,$prot_strand) = $self->_get_strand($prot_start,$prot_end);
          my $pf = new Bio::SeqFeature::Generic( -start   => $prot_start,
                                              -end     => $prot_end,
                                              -seqname => 'protein',
                                              -score   => $score,
                                              -strand  => $prot_strand,
                                              -source_tag => 'genewise',
                                              -primary_tag => 'supporting_protein_feature',
                                              );
          my $geno_strand;
          ($geno_start,$geno_end,$geno_strand) = $self->_get_strand($geno_start,$geno_end);
          my $gf = new Bio::SeqFeature::Generic( -start   => $geno_start,
                                              -end     => $geno_end,
                                              -seqname => 'genomic',
                                              -score   => $score,
                                              -strand  => $geno_strand,
                                              -source_tag => 'genewise',
                                              -primary_tag => 'supporting_genomic_feature',
                                              );
          $exon->add_tag_value( 'supporting_protein_feature' => $pf );
          $exon->add_tag_value( 'supporting_genomic_feature' => $gf );
          $transcript->add_exon($exon);
        }

    $genes->add_transcript($transcript);
  }

  $self->_add_prediction($genes);
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
