package Bio::Tools::Prints;
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


sub parse_results {
        my ($self,$resfile) = @_;
        my $filehandle;
        
        $filehandle = $resfile;
        my %printsac;
        my $line;
        my @features;
        # parse
        my $sequenceId;
        while (<$filehandle>) {
           $line = $_;
           chomp $line;
           # Pattern match the Sn; field which should contain the SequenceId and Accession

           if ($line =~ s/^Sn;//) { # We have identified a Sn; line so there should be the following:

              #ENSP00000003603 Gene:ENSG00000000003 Clone:AL035608 Contig:AL035608.00001 Chr:chrX basepair:97227305
              ($sequenceId) = $line =~ /^\s*(\w+)/;
           }
         
              
           if ($line =~ s/^1TBH//) {
               my  ($id) = $line =~ /^\s*(\w+)/;
               my ($ac) = $line =~ /(PR\w+)\s*$/;
               $printsac{$id} = $ac;
           }
 
           if ($line =~ s/^3TB//) {
               if ($line =~ s/^[HN]//) {
                   my ($num,$temp,$tot) = "";
                   # Grab these lines
                   #   1433ZETA     1  of  6  88.19   1328    1.00e-16  ELTVEERNLLSVAYKNVIGARRASWRIITS  30   35   36   48
                   # split line on space, hence strip off all leading spaces first.
                   $line =~ s/^\s+//;

                   # Place all elements of list into an array
                   my @elements = split /\s+/, $line;

                   # Name each of the elements in the array
                   my ($fingerprintName,$motifNumber,$temp,$tot,$percentageIdentity,$profileScore,$pvalue,$subsequence,$motifLength,$lowestMotifPosition,$matchPosition,$highestMotifPosition) = @elements;
    
                   my $start = $matchPosition;
                   my $end = $matchPosition + $motifLength - 1;
                   my $print =  $printsac{$fingerprintName};

                   my $feat = "$print,$start,$end,$percentageIdentity,$profileScore,$pvalue";

                   push (@features,$feat);
               }
               if ($line =~ s/^F//) {
 
                   foreach my $feats (@features) {
                      $self->create_feature($feats,$sequenceId);
                   }
                   @features = ();
                                                                                                                                                                                                                                                                                                     }
                                                                                                                                                                                                                                                                                                 }
                                                                                                                                                                                                                                                                                               }


        close $filehandle;
        return @feats;
}

sub create_feature {
       my ($self, $feat,$sequenceId) = @_;

       my @f = split (/,/,$feat);
       # create feature object
       my $feat1 = Bio::SeqFeature::Generic->new(-seqname     =>$sequenceId,
                                                 -start       => $f[1],
                                                 -end         => $f[2],
                                                 -score       => $f[4],
                                                 -source      => $feat->{source},
                                                 -primary     => $feat->{primary},
                                                 -logic_name  => $feat->{logic_name}, 
                                                 -percent_id  => $f[3],
                                                 -p_value     => $f[5]
                                               );
     
     my $feat2 = Bio::SeqFeature::Generic->new(-start =>0,
                                               -end => 0,
                                               -seqname => $f[0]);
     
     my $feature =  Bio::SeqFeature::FeaturePair->new(-feature1 => $feat1,
                                                 -feature2 => $feat2);

     
     if ($feature) {
         
         push(@feats ,$feature);
         
       }

        
}
1;


