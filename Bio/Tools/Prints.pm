package Bio::Tools::Prints;
use vars qw(@ISA);
use strict;

use Bio::Root::Root;
use Bio::SeqFeature::FeaturePair;
use Bio::Root::IO;
use Bio::SeqFeature::Generic;
@ISA = qw(Bio::Root::Root Bio::Root::IO Bio::Tools::AnalysisResult);



sub new {
      my($class,@args) = @_;

      my $self = $class->SUPER::new(@args);
      $self->_initialize_io(@args);

      return $self;
}


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

sub print_sac{
    my($self,$printsac)=@_;
 
   if(defined($printsac))
   {
       $self->{'print_sac'}=$printsac;
   }
    return $self->{'print_sac'};

}

sub seqname{
    my($self,$seqname)=@_;

    if(defined($seqname))
    {
        $self->{'seqname'}=$seqname;
    }

    return $self->{'seqname'};

}

1;


