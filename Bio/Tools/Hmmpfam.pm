#BioPerl module for Bio::Tools::Hmmpfam
#
# Cared for by  Balamurugan Kumarasamy
#
# You may distribute this module under the same terms as perl itself
# POD documentation - main docs before the code
#

=head1 NAME

Bio::Tools::Hmmpfam  

=head1 SYNOPSIS

  use Bio::Tools::Hmmpfam;
  my $hmmpfam_parser = new Bio::Tools::Hmmpfam(-fh =>$filehandle );
  while( my $hmmpfam_feat = $hmmpfam_parser->next_result ) {
        push @hmmpfam_feat, $hmmpfam_feat;
  }

=head1 DESCRIPTION

 Parser for Hmmpfam  program

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

 bioperl-bugs@bioperl.org
 http://bugzilla.bioperl.org/

=head1 AUTHOR - Balamurugan Kumarasamy

 Email: fugui@worf.fugu-sg.org

=head1 APPENDIX

 The rest of the documentation details each of the object methods.
 Internal methods are usually preceded with a _


=cut

package Bio::Tools::Hmmpfam;
use vars qw(@ISA);
use strict;

use Bio::Root::Root;
use Bio::SeqFeature::FeaturePair;
use Bio::Root::IO;
use Bio::SeqFeature::Generic;
@ISA = qw(Bio::Root::Root Bio::Root::IO );



=head2 new

 Title   : new
 Usage   : my $obj = new Bio::Tools::Hmmpfam(-fh=>$filehandle);
 Function: Builds a new Bio::Tools::Hmmpfam object
 Returns : Bio::Tools::Hmmpfam
 Args    : -filename
           -fh (filehandle)

=cut

sub new {
      my($class,@args) = @_;

      my $self = $class->SUPER::new(@args);
      $self->_initialize_io(@args);

      return $self;
}


=head2 next_result

 Title   : next_result
 Usage   : my $feat = $hmmpfam_parser->next_result
 Function: Get the next result set from parser data
 Returns : L<Bio::SeqFeature::Generic>
 Args    : none

=cut

sub next_result {
    my ($self) = @_;
    my $filehandle;
    
 my $line;

    my $id;
    while ($_=$self->_readline()) {
         $line = $_;
         chomp $line;
    
        
        last if $line=~m/^Alignments of top-scoring domains/;
        next if ($line=~m/^Model/ || /^\-/ || /^$/);
        
        if ($line=~m/^Query sequence:\s+(\S+)/) {
           $id = $1;
           $self->seqname($id);
        }
       
       if (my ($hid, $start, $end, $hstart, $hend, $score, $evalue) = $line=~m/^(\S+)\s+\S+\s+(\d+)\s+(\d+)\s+\S+\s+(\d+)\s+(\d+)\s+\S+\s+(\S+)\s+(\S+)/) {
            my %feature;
            
            ($feature{name}) = $self->seqname;
            $feature{score} = $score;
            $feature{p_value} = sprintf ("%.3e", $evalue);
            $feature{start} = $start;
            $feature{end} = $end;
            $feature{hname} = $hid;
            $feature{hstart} = $hstart;
            $feature{hend} = $hend;
            ($feature{source}) = 'pfam';
            $feature{primary} = $hid;
            ($feature{program}) = 'pfam';
            ($feature{db}) = 'db1';
            ($feature{logic_name}) = 'hmmpfam';
            my $new_feat = $self->create_feature (\%feature);
            return $new_feat
        
        }
        next;

    }
    return;
}

=head2 create_feature

 Title   : create_feature
 Usage   : my $feat=$hmmpfam_parser->create_feature($feature,$seqname)
 Function: creates a SeqFeature Generic object
 Returns : L<Bio::SeqFeature::Generic>
 Args    :


=cut

sub create_feature {
    my ($self, $feat) = @_;



    my $feature1= Bio::SeqFeature::Generic->new( -seqname    =>$feat->{name},
                                                -start      =>$feat->{start},
                                                -end        =>$feat->{end},
                                                -score      =>$feat->{score},
                                                -source     =>$feat->{source},
                                                -primary    =>$feat->{primary},
                                                   );
    


    my $feature2= Bio::SeqFeature::Generic->new(
                                                 -start      =>$feat->{hstart},
                                                 -end        =>$feat->{hend},
                                                  );




    my $featurepair = Bio::SeqFeature::FeaturePair->new;
    $featurepair->feature1 ($feature1);
    $featurepair->feature2 ($feature2);
   
   $featurepair->add_tag_value('evalue',$feat->{p_value});
   $featurepair->add_tag_value('percent_id','NULL');
   $featurepair->add_tag_value("hid",$feat->{primary});
    return  $featurepair; 
        
}

=head2 seqname

 Title   :   seqname
 Usage   :   obj->seqname($seqname)
 Function:   Internal(not to be used directly)
 Returns :
 Args    :   seqname

=cut

sub seqname{
    my($self,$seqname)=@_;

    if(defined($seqname))
    {
        $self->{'seqname'}=$seqname;
    }

    return $self->{'seqname'};

}

1;


