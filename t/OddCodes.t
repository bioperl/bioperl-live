# -*-Perl-*-
## Bioperl Test Harness Script for Modules
##$Id$

use Test;
use strict;

BEGIN { plan tests => 10}

use Bio::PrimarySeq;
use Bio::Tools::OddCodes;
ok 1;

my ($seqobj, $oddcode_obj);

$seqobj = Bio::PrimarySeq->new('-seq'=>'ABCDEFGHIJKLMNOPQRSTUVWXYZ',
			       '-moltype'=>'protein', 
			       '-id'=>'test');
$oddcode_obj  =  Bio::Tools::OddCodes->new($seqobj);

ok defined($oddcode_obj) && ref($oddcode_obj) && 
    $oddcode_obj->isa('Bio::Tools::OddCodes');

ok ${$oddcode_obj->structural()}, 'ABAEEIAEIJEIIEOAEEAAUIAXAZ';
ok ${$oddcode_obj->chemical()}, 'LBSAARLCLJCLSMOIMCHHULRXRZ';
ok ${$oddcode_obj->functional()}, 'HBPAAHPCHJCHHPOHPCPPUHHXPZ';
ok ${$oddcode_obj->charge()}, 'NBNAANNCNJCNNNONNCNNUNNXNZ';
ok ${$oddcode_obj->hydrophobic()}, 'IBOOOIOOIJOIIOOIOOOOUIIXOZ';
ok ${$oddcode_obj->Dayhoff()}, 'CBADDGCEFJEFFDOCDECCUFGXGZ';
ok ${$oddcode_obj->Sneath()}, 'CBEFFHCHAJGADDOCDGEEUAHXHZ';
ok ${$oddcode_obj->Stanfel()}, 'ABACCDAEAJEAACOACEAAUADXDZ';
