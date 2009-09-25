# -*-Perl-*- Test Harness script for Bioperl
# $Id$


use strict;

BEGIN {
   use lib '.';
   use Bio::Root::Test;
   
   test_begin(-tests => 14);
   
   use_ok('Bio::Matrix::PSM::IO');
}

# Test psiblast reading functionality.
my $psmIO =  Bio::Matrix::PSM::IO->new(-format => 'psiblast', 
			      -file   => test_input_file('atp1.matrix'));
ok $psmIO;

my $psm = $psmIO->next_psm;
ok $psm;

# Verify that getting IUPAC sequence is functional
my $IUPAC = 'MEMSINPSEISSIIKEQIENYDTKAEVSEVGTVLSVGDGIARVYGLDNVMAGEMVEFPSGVKGMALNLEEDNVGVVLLGDDTGIKEGDLVKRTGKIVEVPVGEALLGRVVDPLGNPIDAKGPIKTDERRPVEVKAPGIIPRKSVHEPLQTGLKAIDSLVPIGRGQRELIIGDRQTGKTAIAIDTIINQKRINDESTDEGKKVYCIYVAIGQKRSTVAQVVQTLREAGALEYTIIVAATAAAPAPAQYLSAYAGCAIGEAFADNGAAACIIHDDLSRQAVAYAIISLLLRRPPGREAYPGDVFYLHSRLLERAAKLSDELGGGSLTALPIIETQAGDVSAYIPTNVISITDGQIFLETDLFNSGIRPAINVGLSVSRVGSAAQIKAMKKVAGSLKLELAQYRELAAFAQFGSDLDAATQAQLNRGARLTELLKQPQYSPLPVEEQVVILYAGVNGYLDDIPVEDIRDFEKELLEYLKSNHPEILESIRTGKLSDEIEKALKEAIKEFV';
is $psm->IUPAC, $IUPAC;

## Lets try to compress and uncompress the log odds and the
## frequencies, see if there is no considerable loss of data.
SKIP: {
   skip('TODO: Module incomplete',10); 
   my $fA=$psm->get_compressed_freq('A');
   my @check=Bio::Matrix::PSM::SiteMatrix::_uncompress_string($fA,1,1);
   my @A=$psm->get_array('A');
   my ($var,$max) = (0,0);
   
   for (my $i = 0; $i<@check;$i++) {
     my $diff=abs(abs($check[$i])-abs($A[$i]));
     $var += $diff;
     $max=$diff if ($diff>$max);
   }
   my $avg=$var/@check;
   cmp_ok $avg,'<',0.01; #Loss of data under 1 percent
   is $psm->sequence_match_weight('CAGAAAAATAAAATGGCCACCACCC'),2015;
   
   my $lA=$psm->get_compressed_logs('A');
   @check=Bio::Matrix::PSM::SiteMatrix::_uncompress_string($lA,1000,2);
   @A=$psm->get_logs_array('A');
   ($var,$max) = (0,0);
   for (my $i = 0;$i<@check;$i++) {
     my $diff=abs(abs($check[$i])-abs($A[$i]));
     $var += $diff;
     $max=$diff if ($diff>$max);
   }
   $avg=$var/@check;
   cmp_ok $avg,'<',10; #Loss of data under 1 percent
   
   my $matrix=$psm->matrix;
   ok $matrix;
   my $psm2=$psm;
   $psm2->matrix($matrix);
   is $psm,$psm2;
   
   is $IUPAC,'CAGAAAAATWVAATYCCCACCHCCC';
   is $IUPAC,$psm2->IUPAC;
   is $IUPAC,$matrix->IUPAC;
   
   my $instances=$psm->instances;
   ok $instances;
   
   foreach my $instance (@{$instances}) {
     my $id=$instance->primary_id;
     is $instance->strand,1;
     last if (ok $id);
   }
}
