# -*-Perl-*-
# Bioperl Test Harness Script for Modules
# $Id$
#
# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

my $error = 0;

use strict;
BEGIN {
	eval { require Test::More; };
	if( $@ ) {
		use lib 't/lib';
	}
	use Test::More;
	plan tests => 8;
	use_ok('Bio::SeqIO');
	use_ok('Bio::Root::IO');
}

my $verbose = $ENV{'BIOPERLDEBUG'};

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

is @a_seq, 3, 'number of sequence objects';

my $esc_name = $a_seq[1]->display_id;
is( $esc_name , 'Name; 4% strewn with \ various / escaped characters',
	 "unescaping of characters, $esc_name");

is $a_seq[0]->alphabet, 'protein', 'alphabets detected';
is $a_seq[1]->alphabet, 'dna', 'alphabets detected';

my $o_file = Bio::Root::IO->catfile("t","data","test.out.ace");
my $a_out = Bio::SeqIO->new(-FILE => "> $o_file",
									 -verbose => $verbose,
									 -FORMAT => 'ace');
my $a_out_ok = 1;
foreach my $a (@a_seq) {
	$a_out->write_seq($a) or $a_out_ok = 0;
}
undef($a_out);  # Flush to disk
is $a_out_ok,1,'writing sequence';

my( $after );
{
	local $/ = undef;
	local *AFTER;
	open AFTER, $o_file;
	$after = <AFTER>;
	close AFTER;
}
unlink($o_file);

is( ($before and $after and ($before eq $after)),1,
	 'test output');

