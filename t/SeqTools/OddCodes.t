# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN { 
    use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 11);
	
	use_ok('Bio::PrimarySeq');
	use_ok('Bio::Tools::OddCodes');
}

my ($seqobj, $oddcode_obj);

$seqobj = Bio::PrimarySeq->new('-seq'=>'ABCDEFGHIJKLMNOPQRSTUVWXYZ',
			       '-alphabet'=>'protein', 
			       '-id'=>'test');
$oddcode_obj  =  Bio::Tools::OddCodes->new('-seq' => $seqobj);

isa_ok $oddcode_obj, 'Bio::Tools::OddCodes';

is ${$oddcode_obj->structural()}, 'ABAEEIAEIJEIIEOAEEAAUIAXAZ';
is ${$oddcode_obj->chemical()}, 'LBSAARLCLJCLSMOIMCHHULRXRZ';
is ${$oddcode_obj->functional()}, 'HBPAAHPCHJCHHPOHPCPPUHHXPZ';
is ${$oddcode_obj->charge()}, 'NBNAANNCNJCNNNONNCNNUNNXNZ';
is ${$oddcode_obj->hydrophobic()}, 'IBOOOIOOIJOIIOOIOOOOUIIXOZ';
is ${$oddcode_obj->Dayhoff()}, 'CBADDGCEFJEFFDOCDECCUFGXGZ';
is ${$oddcode_obj->Sneath()}, 'CBEFFHCHAJGADDOCDGEEUAHXHZ';
is ${$oddcode_obj->Stanfel()}, 'ABACCDAEAJEAACOACEAAUADXDZ';
