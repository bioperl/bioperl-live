#
# Cared for by Bala Savikalpa
#
# 
#
# You may distribute this module under the same terms as perl itself
# POD documentation - main docs before the code

=head1 NAME

Bio::Tools::Prints - Parser for FingerPRINTScanII program 

=head1 SYNOPSIS

  use Bio::Tools::Prints;
  my  $parser = Bio::Tools::Prints->new(-fh=>$filehandle);

  while(my $prints_feature = $parser->next_result){
    push @features, $prints_feature;
  }

=head1 DESCRIPTION

FingerPRINTScan II is a PRINTS fingerprint identification algorithm.
Copyright (C) 1998,1999  Phil Scordis

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists.  Your participation is much appreciated.

  bioperl-l@bioperl.org          - General discussion
  http://bio.perl.org/MailList.html             - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
 the bugs and their resolution.  Bug reports can be submitted via
 email or the web:

  bioperl-bugs@bio.perl.org
  http://bio.perl.org/bioperl-bugs/

=head1 AUTHOR - Bala Savikalpa

Email: fugui@worf.fugu-sg.org

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

package Bio::Tools::Prints;
use vars qw(@ISA);
use strict;

use Bio::Root::Root;
use Bio::SeqFeature::FeaturePair;
use Bio::Root::IO;
use Bio::SeqFeature::Generic;
@ISA = qw(Bio::Root::Root Bio::Root::IO Bio::Tools::AnalysisResult);


=head2 new

 Title   : new
 Usage   : my $obj = new Bio::Tools::Prints(-fh=>$filehandler);
 Function: Builds a new Bio::Tools::Prints object
 Returns : Bio::Tools::Prints
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
 Usage   : my $r = $prints_parser->next_result
 Function: Get the next result set from parser data
 Returns : L<Bio::SeqFeature::Generic>
 Args    : none


=cut

sub next_result {
      my ($self) = @_;
     my %printsac;
     my @features;
     my $line;
     my $sequenceId;
     
     while ($_=$self->_readline()) {
      
      
           $line = $_;
           chomp $line;

           if ($line =~ s/^Sn;//) { # We have identified a Sn; line so there should be the following:

              ($sequenceId) = $line =~ /^\s*(\w+)/;
              $self->seqname($sequenceId);
              next;
           }
         
              
           if ($line =~ s/^1TBH//) {
               my  ($id) = $line =~ /^\s*(\w+)/;
               my ($ac) = $line =~ /(PR\w+)\s*$/;
               $printsac{$id} = $ac;
               $self->print_sac(\%printsac);
               next;
           }
             
             
           if ($line =~ s/^3TB//) {
              
                  
              
              
              if ($line =~ s/^[HN]//) {
                   my ($num,$temp,$tot) = "";
                   
                   $line =~ s/^\s+//;

                   my @elements = split /\s+/, $line;

                   my ($fingerprintName,$motifNumber,$temp,$tot,$percentageIdentity,$profileScore,$pvalue,$subsequence,$motifLength,$lowestMotifPosition,$matchPosition,$highestMotifPosition) = @elements;
    
                   my $start = $matchPosition;
                   my $end = $matchPosition + $motifLength - 1;
                   my $print_sac = $self->print_sac;
                   
                   my %printsac =  %{$print_sac};
                   my $print =  $printsac{$fingerprintName};
                   my $seqname=$self->seqname;
                   my $feat = "$print,$start,$end,$percentageIdentity,$profileScore,$pvalue";
                   my $new_feat =  $self->create_feature($feat,$seqname);
                   return $new_feat;
               }
               if ($line =~ s/^F//) {
                   return;  
               }
                   next;                                                                                                                               }
            next;         
 
      }


        
}

=head2 create_feature

 Title   : create_feature
 Usage   : my $r = $prints_parser->create_feature
 Function: Creates a Seqfeature Generic object
 Returns : L<Bio::SeqFeature::Generic>
 Args    :  

=cut

sub create_feature {
       my ($self, $feat,$sequenceId) = @_;

       my @f = split (/,/,$feat);
       # create feature object
        my $feature= Bio::SeqFeature::Generic->new(-seqname    =>$sequenceId,
                                                   -start=>$f[1],
                                                   -end  => $f[2],
                                                   -score      => $f[4],
                                                   -source     => "PRINTS",
                                                   -primary    =>$f[0],
                                                   -logic_name => "PRINTS",
                                                   );
        $feature->add_tag_value('evalue',$f[5]);
        $feature->add_tag_value('percent_id',$f[3]);
        

     
    return  $feature; 
        
}

=head2 print_sac

 Title   : print_sac
 Usage   : $prints_parser->print_sac($print_sac)
 Function: get/set for pritn_sac 
 Returns : 
 Args    : 

=cut


sub print_sac{
    my($self,$printsac)=@_;
 
   if(defined($printsac))
   {
       $self->{'print_sac'}=$printsac;
   }
    return $self->{'print_sac'};

}

=head2 seqname

 Title   : seqname
 Usage   : $prints_parser->seqname($seqname)
 Function: get/set for seqname
 Returns :
 Args    :

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


