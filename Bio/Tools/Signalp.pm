# Parser module for SignalP Bio::Tools::Signalp
#
# Based on the EnsEMBL module Bio::EnsEMBL::Pipeline::Runnable::Protein::Signalp
# originally written by Marc Sohrmann (ms2@sanger.ac.uk)
# Written in BioPipe by Balamurugan Kumarasamy <savikalpa@fugu-sg.org>
# Cared for by the Fugu Informatics team (fuguteam@fugu-sg.org)

# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=head1 NAME

Bio::Tools::SignalP

=head1 SYNOPSIS

 use Bio::Tools::SignalP;
 my $parser = new Bio::Tools::SignalP(-fh =>$filehandle );
 while( my $sp_feat = $parser->next_result ) {
       #do something
       #eg
       push @sp_feat, $sp_feat;
 }

=head1 DESCRIPTION

 Parser for SignalP output

=head1 FEEDBACK

=head2 Mailing Lists

 User feedback is an integral part of the evolution of this and other
 Bioperl modules. Send your comments and suggestions preferably to
 the Bioperl mailing list.  Your participation is much appreciated.

 bioperl-l@bioperl.org              - General discussion
 http://bioperl.org/MailList.shtml  - About the mailing lists

=head2 Reporting Bugs

 Report bugs to the Bioperl bug tracking system to help us keep track
 of the bugs and their resolution. Bug reports can be submitted via
 email or the web:

 bioperl-bugs@bio.perl.org
 http://bugzilla.bioperl.org/

=head1 AUTHOR

 Based on the EnsEMBL module Bio::EnsEMBL::Pipeline::Runnable::Protein::Signalp
 originally written by Marc Sohrmann (ms2@sanger.ac.uk)
 Written in BioPipe by Balamurugan Kumarasamy <savikalpa@fugu-sg.org>
 Cared for by the Fugu Informatics team (fuguteam@fugu-sg.org)

=head1 APPENDIX

 The rest of the documentation details each of the object methods.
 Internal methods are usually preceded with a _


=cut

package Bio::Tools::Signalp;
use vars qw(@ISA);
use strict;

use Bio::Root::Root;
use Bio::SeqFeature::FeaturePair;
use Bio::Root::IO;
use Bio::SeqFeature::Generic;
@ISA = qw(Bio::Root::Root Bio::Root::IO );



=head2 new

 Title   : new
 Usage   : my $obj = new Bio::Tools::SignalP();
 Function: Builds a new Bio::Tools::SignalP object
 Returns : Bio::Tools::SignalP
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
        
        my $line;
        # parse
        my $id;
        my ( $fact1, $fact2, $end);
        while ($_=$self->_readline()) {
           $line = $_;
           chomp $line;
           
           if ($line=~/^\>(\S+)/) {
              $id = $1;
              $self->seqname($id);
              next;
           }
           elsif ($line=~/max\.\s+Y\s+(\S+)\s+\S+\s+\S+\s+(\S+)/) {
              $fact1 = $2;
              $self->fact1($fact1);
              next;
           }
           elsif ($line=~/mean\s+S\s+(\S+)\s+\S+\s+\S+\s+(\S+)/) {
              $fact2 = $2;
                 $fact1 = $self->fact1;
                 $id = $self->seqname;
              
              if ($fact1 eq "YES" && $fact2 eq "YES") {
                  
                  my $line = $self->_readline();
              
                  if ($line =~ /Most likely cleavage site between pos\.\s+(\d+)/) {
                      $end = $1;
                  }
                  else {
                      $self->throw ("parsing problem in signalp");
                  }
                  my (%feature);
                  $feature{name} = $id;
                  $feature{start} = 1;
                  $feature{end} = $end;
                  $feature{source} = 'Signalp';
                  $feature{primary}= 'signal_peptide';
                  $feature{program} = 'Signalp';
                  $feature{logic_name} = 'signal_peptide';
                  
                  my $new_feat = $self->create_feature (\%feature);
                  return $new_feat;
                  
              }
               next;
               
           }
        
           next;
        
        }
        
}

=head2 create_feature

 Title   : create_feature
 Usage   : obj->create_feature(\%feature)
 Function: Internal(not to be used directly)
 Returns :
 Args    :


=cut

sub create_feature {
       my ($self, $feat) = @_;


       # create feature object
       my $feature = Bio::SeqFeature::Generic->new(
                                                 -seq_id=>$feat->{name},
                                                 -start       => $feat->{start},
                                                 -end         => $feat->{end},
                                                 -score       => $feat->{score},
                                                 -source      => $feat->{source},
                                                 -primary     => $feat->{primary},
                                                 -logic_name  => $feat->{logic_name}, 
                                               );
           

          $feature->add_tag_value('evalue',0);
          $feature->add_tag_value('percent_id','NULL');
          $feature->add_tag_value("hid",$feat->{primary});
          
          return $feature; 

}
=head2 seqname

 Title   : seqname
 Usage   : obj->seqname($name)
 Function: Internal(not to be used directly)
 Returns :
 Args    :


=cut

sub seqname{
    my ($self,$seqname)=@_;

    if (defined$seqname){

        $self->{'seqname'}=$seqname;
    }

    return $self->{'seqname'};

}

=head2 fact1

 Title   : fact1
 Usage   : obj->fact1($fact1)
 Function: Internal(not to be used directly)
 Returns :
 Args    :


=cut

sub fact1{
    my ($self,$fact1)=@_;

    if (defined$fact1){

       $self->{'fact1'}=$fact1;
    }

    return $self->{'fact1'};

}



1;


