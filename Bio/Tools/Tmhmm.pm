package Bio::Tools::Tmhmm;
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
        
        my $line;
       
        # parse
        my $id;
        while ($_=$self->_readline()) { 
           $line = $_;
           chomp $line;
           
           
           next if /^$/;
           if ($line=~/^#\s+(\S+)/) {
               
                    $id = $1;
                    my ($junk, $values) = split /:/;
                   $self->seqname($id);
                    next;
           }
             
           elsif ($line=~/^(\S+)\s+(\S+)\s+(\w+)\s+(\d+)\s+(\d+)/) {
                    my $orien = $3;
                    my $start = $4;
                    my $end = $5;
                    $orien = uc ($orien);
                  
                    if ($orien eq "TMHELIX") {
                         my (%feature);
                         $feature{name} = $self->seqname;
                         $feature{start} = $start;
                         $feature{end} = $end;
                         $feature{source} ='tmhmm';
                         $feature{primary}= 'transmembrane';
                         $feature{program} ='tmhmm';
                         $feature{logic_name} = 'TMHelix';
                         my $new_feat= $self->create_feature(\%feature);
                         return $new_feat;
                    }
                    next;
           }
           next;
        }
}

sub create_feature {
       my ($self, $feat) = @_;


       # create feature object
       my $feature = Bio::SeqFeature::Generic->new(-seqname     => $feat->{name},
                                                 -start       => $feat->{start},
                                                 -end         => $feat->{end},
                                                 -score       => $feat->{score},
                                                 -source      => $feat->{source},
                                                 -primary     => $feat->{primary},
                                                 -logic_name  => $feat->{logic_name}, 
                                               );
 
       $feature->add_tag_value('evalue',0);
       $feature->add_tag_value('percent_id','NULL');
                        
       return $feature;                 
                        
                                                                                                                                           }
sub seqname{
    my ($self,$seqname)=@_;

    if (defined$seqname){

        $self->{'seqname'}=$seqname;
    }

    return $self->{'seqname'};

}




1;


