# -*-Perl-*- Test Harness script for Bioperl
# $Id$

#####
#
# this script simply tests parsing ace* files
# - it cares nothing about the chromat_dir,phd_dir,edit_dir types of things
#
#####

use strict;

BEGIN {
    use lib '.';
	use Bio::Root::Test;
	
	test_begin(-tests => 15);
    
	use_ok('Bio::Tools::Alignment::Consed');
}

my $DEBUG = test_debug();

# scope some variables
my($o_consed,@singlets,@singletons,@pairs,@doublets,@multiplets,$invoker);

# instantiate a new object
my $passed_in_acefile = test_input_file('acefile.ace.1');
$o_consed = Bio::Tools::Alignment::Consed->new(-acefile => $passed_in_acefile);
ok defined $o_consed, 'new CSM::Consed object was created';

$o_consed->verbose($DEBUG);

isa_ok($o_consed,'Bio::Tools::Alignment::Consed');

isnt($o_consed->set_singlets(), 1,  'singlets can be successfully set');
	
@singlets = $o_consed->get_singlets();
is (scalar(@singlets), 65, 'singlets can be retrieved');

isnt ($o_consed->set_doublets(), 1, 'doublets can be set');

ok @doublets = $o_consed->get_doublets(), 'doublets can be retreived';

print(scalar(@doublets)." doublets were found\n") if ($DEBUG > 0);
is (scalar(@doublets), 45, 'doublets can be retrieved');

@pairs = $o_consed->get_pairs();
is (scalar(@pairs),1, 'pairs can be retrieved');

@multiplets = $o_consed->get_multiplets();
is (scalar(@multiplets), 4, 'multiplets can be retrieved');

@singletons = $o_consed->get_singletons();
is (scalar(@singletons), 3, 'singletons can be retrieved');
my($total_object_sequences, $total_grep_sequences);
is($total_grep_sequences = $o_consed->count_sequences_with_grep(), 179, 'how many sequences there are in the acefile _and_ in the singlets file');

is($total_object_sequences = $o_consed->sum_lets("total_only"),179, 'statistics from the Bio::Tools::Alignment::Consed object to compare the total number of sequences accounted for there to the number of sequences found via grep');
print("Match?\n") if($DEBUG > 0) ;
is ($total_object_sequences, $total_grep_sequences);

print("These are the statistics. Look right? ".$o_consed->sum_lets()."\n") if($DEBUG > 0);
is($o_consed->sum_lets(),'Singt/singn/doub/pair/mult/total : 65,3,45(90),1(2),4(19),179');

print("Dumping out the hash in a compact way...\n")if($DEBUG > 0)  ;
$o_consed->dump_hash_compact() if($DEBUG > 0)  ;

# print("Dumping out the hash in an ugly way...\n");
# $o_consed->dump_hash();

sub allele_script {
	my($a,$trunc,$rev);
	ok defined $a,
	isa_ok $a, 'Bio::Variation::Allele';
	
	is $a->accession_number(), 'X677667';
	is $a->seq(), 'ACTGACTGACTG';
	is $a->display_id(),'new-id' ;
	is $a->desc, 'Sample Bio::Seq object';
	is $a->moltype(), 'dna';

	ok defined($trunc = $a->trunc(1,4));
	is $trunc->seq(), 'ACTG', "Expecting ACTG. Got ". $trunc->seq();

	ok defined($rev = $a->revcom());
	is $rev->seq(), 'CAGTCAGTCAGT';

	$a->is_reference(1);
	ok $a->is_reference;

	$a->repeat_unit('ACTG');
	is $a->repeat_unit, 'ACTG';
	
	$a->repeat_count(3);
	is $a->repeat_count, 3;
}
