# -*-Perl-*-
# Bioperl Test Harness Script for Modules
# $Id$
#
# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

my $error = 0;

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
	plan tests => 7;
}

if( $error == 1 ) {
	exit(0);
}

use Bio::SeqIO;
use Bio::Root::IO;

my $verbose = $ENV{'BIOPERLDEBUG'};
ok(1);

my $t_file = Bio::Root::IO->catfile("t","data","test.ace");
my( $before );
{
	local $/ = undef;
	local *BEFORE;
	open BEFORE, $t_file;
	$before = <BEFORE>;
	close BEFORE;
}

my $a_in = Bio::SeqIO->new( -FILE => $t_file,
									 -verbose => $verbose,
									 -FORMAT => 'ace');
my( @a_seq );
while (my $a = $a_in->next_seq) {
	push(@a_seq, $a);
}

ok @a_seq, 3, 'wrong number of sequence objects';

my $esc_name = $a_seq[1]->display_id;
ok( $esc_name , 'Name; 4% strewn with \ various / escaped characters',
	 "bad unescaping of characters, $esc_name");

ok $a_seq[0]->alphabet, 'protein', 'alphabets incorrectly detected';
ok $a_seq[1]->alphabet, 'dna', 'alphabets incorrectly detected';

my $o_file = Bio::Root::IO->catfile("t","data","test.out.ace");
my $a_out = Bio::SeqIO->new(-FILE => "> $o_file",
									 -verbose => $verbose,
									 -FORMAT => 'ace');
my $a_out_ok = 1;
foreach my $a (@a_seq) {
	$a_out->write_seq($a) or $a_out_ok = 0;
}
undef($a_out);  # Flush to disk
ok $a_out_ok,1,'error writing sequence';

my( $after );
{
	local $/ = undef;
	local *AFTER;
	open AFTER, $o_file;
	$after = <AFTER>;
	close AFTER;
}
unlink($o_file);

ok( ($before and $after and ($before eq $after)),1,
	 'test output file differs from input');

