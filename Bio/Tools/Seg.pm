package Bio::Tools::Seg;
use vars qw(@ISA);
use strict;

use Bio::Root::Root;
use Bio::SeqFeature::FeaturePair;
use Bio::Root::IO;
use Bio::SeqFeature::Generic;
@ISA = qw(Bio::Root::Root Bio::Root::IO Bio::Tools::AnalysisResult);



my @feats;
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
           if ($line=~/^\>/) {
               $line=~/^\>\s*(\S+)\s*\((\d+)\-(\d+)\)\s*complexity=(\S+)/;
               my $id = $1;
               my $start = $2;
               my $end = $3;
               my $score = $4;
               my (%feature);
               $feature{name} = $id;
               $feature{score} = $score;
               $feature{start} = $start;
               $feature{end} = $end;
               $feature{source} = "Seg";
               $feature{primary} = 'low_complexity';
               $feature{program} = "Seg";
               $feature{logic_name} = 'low_complexity';
               my $new_feat =  $self->create_feature (\%feature);
               return $new_feat;
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
1;


