# Parser module for Coil Bio::Tools::Coil
#
# Based on the EnsEMBL module Bio::EnsEMBL::Pipeline::Runnable::Protein::Coil
# originally written by Marc Sohrmann (ms2@sanger.ac.uk)
# Written in BioPipe by Balamurugan Kumarasamy <savikalpa@fugu-sg.org>
# Cared for by the Fugu Informatics team (fuguteam@fugu-sg.org)

# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=head1 NAME

Bio::Tools::Coil

=head1 SYNOPSIS

 use Bio::Tools::Coil
 my $parser = new Bio::Tools::Coil();
 while( my $sp_feat = $parser->next_result($file) ) {
       #do something
       #eg
       push @sp_feat, $sp_feat;
 }

=head1 DESCRIPTION

 Parser for Coil output

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

 Based on the EnsEMBL module Bio::EnsEMBL::Pipeline::Runnable::Protein::Coil
 originally written by Marc Sohrmann (ms2@sanger.ac.uk)
 Written in BioPipe by Balamurugan Kumarasamy <savikalpa@fugu-sg.org>
 Cared for by the Fugu Informatics team (fuguteam@fugu-sg.org)

=head1 APPENDIX

 The rest of the documentation details each of the object methods.
 Internal methods are usually preceded with a _


=cut

package Bio::Tools::Coil;
use vars qw(@ISA);
use strict;

use Bio::Root::Root;
use Bio::SeqFeature::FeaturePair;
use Bio::Root::IO;
use Bio::SeqFeature::Generic;
@ISA = qw(Bio::Root::Root Bio::Root::IO);



sub new {
      my($class,@args) = @_;

      my $self = $class->SUPER::new(@args);
      $self->_initialize_io(@args);

      return $self;
}

=head2 parse_results

 Title   : parse_results
 Usage   : obj->parse_results
 Function: Parses the coil output. Automatically called by
           next_result() if not yet done.
 Example :
 Returns :

=cut

sub parse_results {
  my ($self,$resfile) = @_;
        my $filehandle = $resfile;
        my %result_hash =_read_fasta($filehandle);#bala no file handle
        my @ids = keys %result_hash;
        my @feats; 
        foreach my $id (keys %result_hash){      
                my $pep = reverse ($result_hash{$id});
                my $count = my $switch = 0;
                my ($start, $end);
                while (my $aa = chop $pep) {
                       $count++;
                       if (!$switch && $aa eq "x") {
                            $start = $count;
                            $switch = 1;
                        }
                       elsif ($switch && $aa ne "x") {
                            $end = $count-1;
                            my (%feature);
                            $feature{name}       = $id;
                            $feature{start}      = $start;
                            $feature{end}        = $end;
                            $feature{source}     = "Coils";
                            $feature{primary}    = 'ncoils';
                           ($feature{program})   = 'ncoils';
                            $feature{logic_name} = 'Coils';
                            my $new_feat = $self->create_feature (\%feature);
                            $self->_add_prediction($new_feat);
                            $switch = 0;
                       }
                 }
         }
        
        $self->_predictions_parsed(1);
        
}
=head2 next_result

 Title   : next_result
 Usage   : while($feat = $coil->next_result($file)) {
                  # do something
           }
 Function: Returns the next protein feature of the coil output file
 Returns : 
 Args    :

=cut

sub next_result{
    
    my ($self,$resfile) = @_;
    my $gene;

    $self->parse_results($resfile) unless $self->_predictions_parsed();

    $gene = $self->_result();

    return $gene;

}

=head2 _result

 Title   : _result
 Usage   : $feat = $obj->_result()
 Function: internal
 Example :
 Returns :

=cut

sub _result{
    my ($self) = @_;

    return undef unless(exists($self->{'_feats'}) && @{$self->{'_feats'}});
    return shift(@{$self->{'_feats'}});

}

=head2 _add_prediction

 Title   : _add_prediction()
 Usage   : $obj->_add_prediction($feat)
 Function: internal
 Example :
 Returns :

=cut

sub _add_prediction {
    my ($self, $gene) = @_;

    if(! exists($self->{'_feats'})) {
        $self->{'_feats'} = [];
    }
    push(@{$self->{'_feats'}}, $gene);
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
       my $feature = Bio::SeqFeature::Generic->new(-seq_id     => $feat->{name},
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

=head2 _read_fasta

 Title   : _read_fasta
 Usage   : obj->_read_fasta($file)
 Function: Internal(not to be used directly)
 Returns :
 Args    :


=cut

sub _read_fasta {
        local (*FILE) = @_;
        my( $id , $seq , %name2seq);#bala
        while (<FILE>) {
       chomp;#bala
             if (/^>(\S+)/) {
              
              my $new_id = $1;
                    if ($id) {
                           $name2seq{$id} = $seq;
                    }
                    $id = $new_id ; $seq = "" ;
             }
             elsif (eof) {
                    if ($id) {
                           $seq .= $_ ;#bala line instead of $_
                           $name2seq{$id} = $seq;
                    }
             }
             else {
                     $seq .= $_
             }
        }
                                                                                                                                                   return %name2seq;
}



1;


