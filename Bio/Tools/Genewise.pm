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
    my $curr_exon;
    my $score;
    my $seqname;
    my $start;
    my $end;
    my $strand;
    my $phaseno;
    my $pf;
    my $gf;
    #The big parsing loop - parses exons and predicted peptides
    #$/ = "\n//\n";
    while (<$filehandle>) {
        chomp;
        my @f = split;
        if($f[0] eq 'Score'){
          $score = $f[1];
        }
        elsif (($f[0] eq 'Gene') && (scalar(@f)==2)) {
          $seqname = $f[0]." ".$f[1];
        }
        elsif ($f[0] eq 'Exon') {
          $start = $f[1];
          $end = $f[2];
          $phaseno = $f[4];
          $strand = 1;
          if ( $f[1] > $f[2] ) {
              $strand = -1;
              $start = $f[2];
              $end = $f[1];
          }
        }
        elsif ($f[0] eq 'Supporting') {
          my $gstart = $f[1];
          my $gend = $f[2];
          my $gstrand = 1;
          if ($gstart > $gend){
              $gstart = $f[2];
              $gend = $f[1];
              $gstrand = -1;
          }
            if ( $gstrand != $strand ) {
              $self->throw("incompatible strands between exon and supporting feature - cannot add suppfeat\n");
            }

          my $pstart = $f[3];
          my $pend = $f[4];
          my $pstrand = 1;
            if($pstart > $pend){
              $self->warn("Protein start greater than end! Skipping this suppfeat\n");
            }

          $pf = new Bio::SeqFeature::Generic( -start   => $pstart,
              -end     => $pend,
              -seqname => 'protein',
              -score   => $score,
              -strand  => $pstrand,
              -source_tag => 'genewise',
              -primary_tag => 'supporting_protein_feature',
              );
          $pf->source_tag('genewise');
          $pf->primary_tag('supporting_protein_feature');
          $gf  = new Bio::SeqFeature::Generic( -start   => $gstart,
              -end     => $gend,
              -seqname => 'genomic',
              -score   => $score,
              -strand  => $gstrand,
              -source_tag => 'genewise',
              -primary_tag => 'supporting_genomic_feature',
              );
          $gf->source_tag('genewise');
          $gf->primary_tag('supporting_genomic_feature');
        }
         # for listing out elements of the array
         # for ( my $i=0; $i<scalar(@f); $i++) {
         #     print "$i "."$f[$i]\n";
         # }
    }
    $curr_exon = new Bio::SeqFeature::Gene::Exon (-seqname=>$seqname, -start=>$start, -end=>$end, -strand=>$strand);
    $curr_exon->add_tag_value( 'phase' => $phaseno );
    $curr_exon->add_tag_value( 'supporting_protein_feature' => $pf );
    $curr_exon->add_tag_value( 'supporting_genomic_feature' => $gf );

    $transcript->add_exon($curr_exon);
    $genes->add_transcript($transcript);

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
