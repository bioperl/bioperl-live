# $Id$
#---------------------------------------------------------

# tests for Bio::Matrix::PSM::ProtPsm
# written by James Thompson <tex@biosysadmin.com>

use strict;

BEGIN {
   # to handle systems with no installed Test module
   # we include the t dir (where a copy of Test.pm is located)
   # as a fallback
   eval { require Test; };
   if( $@ ) {
	use lib 't';
   }
   use Test;

   plan tests => 1;
}

use Bio::Matrix::PSM::ProtPsm;
ok 1;
use Bio::Matrix::PSM::IO;
ok 2;

# Test psiblast reading functionality.
my $psmIO =  new Bio::Matrix::PSM::IO(-format => 'psiblast', 
				      -file   => Bio::Root::IO->catfile(qw(data atp1.matrix)));
ok $psmIO;

my $psm = $psmIO->next_psm;
ok $psm;

# Verify that getting IUPAC sequence is functional
my $IUPAC = 'MEMSINPSEISSIIKEQIENYDTKAEVSEVGTVLSVGDGIARVYGLDNVMAGEMVEFPSGVKGMALNLEEDNVGVVLLGDDTGIKEGDLVKRTGKIVEVPVGEALLGRVVDPLGNPIDAKGPIKTDERRPVEVKAPGIIPRKSVHEPLQTGLKAIDSLVPIGRGQRELIIGDRQTGKTAIAIDTIINQKRINDESTDEGKKVYCIYVAIGQKRSTVAQVVQTLREAGALEYTIIVAATAAAPAPAQYLSAYAGCAIGEAFADNGAAACIIHDDLSRQAVAYAIISLLLRRPPGREAYPGDVFYLHSRLLERAAKLSDELGGGSLTALPIIETQAGDVSAYIPTNVISITDGQIFLETDLFNSGIRPAINVGLSVSRVGSAAQIKAMKKVAGSLKLELAQYRELAAFAQFGSDLDAATQAQLNRGARLTELLKQPQYSPLPVEEQVVILYAGVNGYLDDIPVEDIRDFEKELLEYLKSNHPEILESIRTGKLSDEIEKALKEAIKEFV';
ok $psm->IUPAC, $IUPAC;

##Lets try to compress and uncompress the log odds and the frequencies, see if there is no
##considerable loss of data.
#my $fA=$psm->get_compressed_freq('A');
#my @check=Bio::Matrix::PSM::SiteMatrix::_uncompress_string($fA,1,1);
#my @A=$psm->get_array('A');
#my ($var,$max) = (0,0);
#for (my $i = 0; $i<@check;$i++) {
#  my $diff=abs(abs($check[$i])-abs($A[$i]));
#  $var += $diff;
#  $max=$diff if ($diff>$max);
#}
#my $avg=$var/@check;
#ok $avg<0.01; #Loss of data under 1 percent
#ok $psm->sequence_match_weight('CAGAAAAATAAAATGGCCACCACCC'),2015;
#
#my $lA=$psm->get_compressed_logs('A');
#@check=Bio::Matrix::PSM::SiteMatrix::_uncompress_string($lA,1000,2);
#@A=$psm->get_logs_array('A');
#($var,$max) = (0,0);
#for (my $i = 0;$i<@check;$i++) {
#  my $diff=abs(abs($check[$i])-abs($A[$i]));
#  $var += $diff;
#  $max=$diff if ($diff>$max);
#}
#$avg=$var/@check;
#ok $avg<10; #Loss of data under 1 percent
#
#my $matrix=$psm->matrix;
#ok $matrix;
#my $psm2=$psm;
#$psm2->matrix($matrix);
#ok $psm,$psm2;


#ok $IUPAC,'CAGAAAAATWVAATYCCCACCHCCC';
#ok $IUPAC,$psm2->IUPAC;
#ok $IUPAC,$matrix->IUPAC;
#
#my $instances=$psm->instances;
#ok $instances;
#
#foreach my $instance (@{$instances}) {
#  my $id=$instance->primary_id;
#  ok $instance->strand,1;
#  last if (ok $id);
#}
