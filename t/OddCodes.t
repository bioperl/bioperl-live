# -*-Perl-*-
## Bioperl Test Harness Script for Modules
##$Id$

use Test;
use strict;

BEGIN { plan tests => 10}

use Bio::PrimarySeq;
use Bio::Tools::OddCodes;
ok 1;

my ($seqobj, $oddcode_obj, $s);

$seqobj = Bio::PrimarySeq->new('-seq'=>'ABCDEFGHIIJKLMNOPQRSTUVWXYZ',
			       '-moltype'=>'protein', 
			       '-id'=>'test');
$oddcode_obj  =  Bio::Tools::OddCodes->new($seqobj);

ok defined($oddcode_obj) && ref($oddcode_obj) && 
    $oddcode_obj->isa('Bio::Tools::OddCodes');

ok ${$oddcode_obj->structural()}, 'ABAEEIAEIIJEIIEOAEEAAUIAXAZ';
ok ${$oddcode_obj->chemical()}, 'LBSAARLCLLJCLSMOIMCHHULRXRZ';
ok ${$oddcode_obj->functional()}, 'HBPAAHPCHHJCHHPOHPCPPUHHXPZ';
ok ${$oddcode_obj->charge()}, 'NBNAANNCNNJCNNNONNCNNUNNXNZ';
ok ${$oddcode_obj->hydrophobic()}, 'IBOOOIOOIIJOIIOOIOOOOUIIXOZ';
ok ${$oddcode_obj->Dayhoff()}, 'CBADDGCEFFJEFFDOCDECCUFGXGZ';
ok ${$oddcode_obj->Sneath()}, 'CBEFFHCHAAJGADDOCDGEEUAHXHZ';
ok ${$oddcode_obj->Stanfel()}, 'ABACCDAEAAJEAACOACEAAUADXDZ';
