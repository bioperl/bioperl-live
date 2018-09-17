# Parser module for Signalp Bio::Tools::Signalp
#
# 
# Based on the EnsEMBL module
# Bio::EnsEMBL::Pipeline::Runnable::Protein::Signalp originally
# written by Marc Sohrmann (ms2@sanger.ac.uk) Written in BioPipe by
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Balamurugan Kumarasamy <savikalpa@fugu-sg.org> Cared for by the Fugu
# Informatics team (fuguteam@fugu-sg.org)

# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=head1 NAME

Bio::Tools::Signalp - parser for Signalp output

=head1 SYNOPSIS

 use Bio::Tools::Signalp;

 my $parser = Bio::Tools::Signalp->new(-fh =>$filehandle );

 while( my $sp_feat = $parser->next_result ) {
   if ($sp_feat->score > 0.9) {
      push @likely_sigpep, $sp_feat;
   }
 }

=head1 DESCRIPTION

C<SignalP> predicts the presence and location of signal peptide
cleavage sites in amino acid sequences.

L<Bio::Tools::Signalp> parses the output of C<SignalP> to provide a 
L<Bio::SeqFeature::Generic> object describing the signal peptide
found, if any. It returns a variety of tags extracted from the NN and HMM
analysis. Most importantly, the C<score()> attribute contains the
NN probability of this being a true signal peptide.


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
of the bugs and their resolution. Bug reports can be submitted va the
web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR

# Please direct questions and support issues to I<bioperl-l@bioperl.org> 

Based on the EnsEMBL module Bio::EnsEMBL::Pipeline::Runnable::Protein::Signalp
originally written by Marc Sohrmann (ms2_AT_sanger.ac.uk). Written in BioPipe by
Balamurugan Kumarasamy savikalpa_AT_fugu-sg.org. Cared for by the Fugu
Informatics team (fuguteam_AT_fugu-sg.org)

=head1 CONTRIBUTORS

Torsten Seemann - torsten.seemann AT infotech.monash.edu.au

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

package Bio::Tools::Signalp;
use strict;

use Bio::SeqFeature::Generic;
use base qw(Bio::Root::Root Bio::Root::IO);



=head2 new

 Title   : new
 Usage   : my $obj = Bio::Tools::Signalp->new();
 Function: Builds a new Bio::Tools::Signalp object
 Returns : Bio::Tools::Signalp
 Args    : -fh/-file => $val, # for initing input, see Bio::Root::IO

=cut

sub new {
      my($class,@args) = @_;

      my $self = $class->SUPER::new(@args);
      $self->_initialize_io(@args);

      return $self;
}

=head2 next_result

 Title   : next_result
 Usage   : my $feat = $signalp->next_result
 Function: Get the next result set from parser data
 Returns : Bio::SeqFeature::Generic
 Args    : none

=cut

sub next_result {
        my ($self) = @_;
        
        while (my $line=$self->_readline()) {
           chomp $line;
           
           if ($line=~/^\>(\S+)/) {
              $self->_seqname($1);
           }
           elsif ($line=~/max\.\s+Y\s+(\S+)\s+\S+\s+\S+\s+(\S+)/) {
              $self->_fact1($2);
           }
           elsif ($line=~/mean\s+S\s+(\S+)\s+\S+\s+\S+\s+(\S+)/) {
              my $fact2 = $2;
              
              if ($fact2 eq 'YES' and $self->_fact1 eq 'YES') {
                  
                  my $line = $self->_readline();
                  
                  ###########################################
                  # modification to suit new SignalP output
                  ###########################################
					chomp $line;
					#print STDERR "********** <$line>\n";
					if ($line =~ /\s+D\s+.*/) {
				  		$line = $self->_readline();
					}
					#print STDERR "********** <$line>\n";
		  			my $end;
		  		  ###########################################
		  
              
                  if ($line =~ /Most likely cleavage site between pos\.\s+(\d+)/) {
                      my $end = $1;
                      my (%feature);
                      $feature{seq_id} = $self->_seqname;
                      $feature{start} = 1;
                      $feature{end} = $end;
                      $feature{source_tag} = 'Signalp';
                      $feature{primary}= 'signal_peptide';
                      $self->_parse_hmm_result(\%feature);
                      my $new_feat = $self->_create_feature (\%feature);
                      return $new_feat;
                  }
                  else {
                      $self->throw ("parsing problem in signalp");
                  }
                  
              }
           }
        
        }
}

=head2 _parse_hmm_result

 Title   : _parse_hmm_result
 Usage   : $self->_parse_hmm_result(\%feature)
 Function: Internal (not to be used directly)
 Returns : hash of feature values
 Args    : hash of more feature values

=cut

sub _parse_hmm_result {
    my ($self, $feature_hash) = @_;
    while(my $line = $self->_readline){
        chomp $line;
        if($line =~ /Prediction: (.+)$/){
            $feature_hash->{hmmProdiction} = $1;
        }elsif($line =~ /Signal peptide probability: ([0-9\.]+)/){
            $feature_hash->{peptideProb} = $1;
        }elsif($line =~ /Signal anchor probability: ([0-9\.]+)/){
            $feature_hash->{anchorProb} = $1;
            last;
        }
    }
}

=head2 _create_feature

 Title   : _create_feature
 Usage   : $self->create_feature(\%feature)
 Function: Internal (not to be used directly)
 Returns : hash of feature values
 Args    : hash of more feature values

=cut

sub _create_feature {
    my ($self, $feat) = @_;

    # create feature object
    my $feature = Bio::SeqFeature::Generic->new(
         -seq_id      => $feat->{name},
         -start       => $feat->{start},
         -end         => $feat->{end},
         -score       => $feat->{score},
         -source      => $feat->{source},
         -primary     => $feat->{primary},
         -logic_name  => $feat->{logic_name}, 
    );
           
    $feature->score($feat->{peptideProb});
    $feature->add_tag_value('peptideProb', $feat->{peptideProb});
    $feature->add_tag_value('anchorProb', $feat->{anchorProb});
    $feature->add_tag_value('evalue',$feat->{anchorProb});
    $feature->add_tag_value('percent_id','NULL');
    $feature->add_tag_value("hid",$feat->{primary});
    $feature->add_tag_value('SignalpPrediction', $feat->{hmmProdiction});
    return $feature; 

}

=head2 _seqname

 Title   : _seqname
 Usage   : $self->_seqname($name)
 Function: Internal (not to be used directly)
 Returns :
 Args    :

=cut

sub _seqname{
    my ($self,$seqname)=@_;

    if (defined$seqname){
        $self->{'seqname'}=$seqname;
    }
    return $self->{'seqname'};
}

=head2 _fact1

 Title   : _fact1
 Usage   : $self->fact1($fact1)
 Function: Internal (not to be used directly)
 Returns : 
 Args    :

=cut

sub _fact1{
    my ($self, $fact1)=@_;

    if (defined $fact1){
       $self->{'fact1'}=$fact1;
    }
    return $self->{'fact1'};
}



1;


