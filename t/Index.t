# -*-Perl-*-
## Bioperl Test Harness Script for Modules
##
# CVS Version
# $Id$


# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

#-----------------------------------------------------------------------
## perl test harness expects the following output syntax only!
## 1..3
## ok 1  [not ok 1 (if test fails)]
## 2..3
## ok 2  [not ok 2 (if test fails)]
## 3..3
## ok 3  [not ok 3 (if test fails)]
##
## etc. etc. etc. (continue on for each tested function in the .t file)
#-----------------------------------------------------------------------


## We start with some black magic to print on failure.
BEGIN { $| = 1; print "1..8\n"; 
	use vars qw($loaded); }

END {print "not ok 1\n" unless $loaded;}

use strict;
use Bio::Index::Fasta;
use Bio::Index::SwissPfam;
use Bio::Index::EMBL;
use Bio::Index::GenBank;
use Bio::Index::Swissprot;

my $dir;
eval {
    use Cwd;
    $dir = cwd;
};
if( $@) {
    # CWD not installed, revert to unix behavior, best we can do
    $dir = `pwd`;
}

$loaded = 1;
print "ok 1\n";    # 1st test passes.

sub test ($$;$) {
    my($num, $true,$msg) = @_;
    print($true ? "ok $num\n" : "not ok $num $msg\n");
}

chomp( $dir );
{
    my $ind = Bio::Index::Fasta->new(-filename => 'Wibbl', -write_flag => 1,
				     -verbose => 0);
    $ind->make_index("$dir/t/multifa.seq");
    $ind->make_index("$dir/t/seqs.fas");
    $ind->make_index("$dir/t/multi_1.fa");
}

test 2, ( -e "Wibbl" );

# Test that the sequences we've indexed
# are all retrievable, and the correct length
{
    my %t_seq = (
        HSEARLOBE               => 321,
        HSMETOO                 => 134,
        MMWHISK                 => 62,
        'gi|238775|bbs|65126'   => 70,
    );

    # Tests opening index of unknown type
    my $ind = Bio::Index::Abstract->new(-FILENAME => 'Wibbl');

    my $ok_3 = 1;
    while (my($name, $length) = each %t_seq) {
        my $seq = $ind->fetch($name);
        if( defined($seq) and $seq->isa('Bio::SeqI') ) {
            my $r_length = $seq->length;
	    unless ($r_length == $length) {
                warn "$name - retrieved length '$r_length' doesn't match known length '$length'\n";
                $ok_3 = 0;
            }
        } else {
            warn "Didn't get sequence '$name' from index\n";
            $ok_3 = 0;
        }
    }
    test 3, $ok_3;

    my $stream = $ind->get_PrimarySeq_stream();
    my $ok_4 = 1;
    while( my $seq2 = $stream->next_primary_seq ) {
	unless ($seq2->isa('Bio::PrimarySeqI')) {
	    $ok_4 = 0;
	    last; # no point continuing...
	}
    }
    test 4, $ok_4;
}

{
    my $ind = Bio::Index::SwissPfam->new(-filename=>'Wibbl2', 
					 -write_flag=>1);
    $ind->make_index("$dir/t/swisspfam.data");
    test 5, ( -e "Wibbl2" );
}

{
    my $ind = Bio::Index::EMBL->new(-filename=>'Wibbl3', 
				    -write_flag=>1);
    $ind->make_index("$dir/t/test.embl");
    test 6, ( -e "Wibbl3" ) && $ind->fetch('AL031232')->length == 4870;
}

{
    my $ind = Bio::Index::Swissprot->new(-filename=>'Wibbl4', 
				    -write_flag=>1);
    $ind->make_index("$dir/t/roa1.swiss");
    test 7, ( -e "Wibbl4" ) && ($ind->fetch('P09651')->display_id() eq 'ROA1_HUMAN');
}

{
    my $ind = Bio::Index::GenBank->new(-filename=>'Wibbl5', 
				       -write_flag=>1, 
				       -verbose => 0);
    $ind->make_index("$dir/t/roa1.genbank");
    test 8, ( -e "Wibbl5" ) && $ind->fetch('AI129902')->length == 37;
}

END { 
    unlink qw( Wibbl Wibbl2 Wibbl3 Wibbl4 Wibbl5);
}



