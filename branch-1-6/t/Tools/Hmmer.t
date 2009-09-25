# -*-Perl-*- Test Harness script for Bioperl
# $Id: Hmmer.t 14989 2008-11-11 19:52:02Z cjfields $

use strict;

BEGIN {     
    use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 29);
	
	use_ok('Bio::Tools::HMMER::Domain');
	use_ok('Bio::Tools::HMMER::Set');
	use_ok('Bio::Tools::HMMER::Results');
}

my ($domain,$set,$homol,$rev,$res,$dom,@doms);
$domain = Bio::Tools::HMMER::Domain->new(-verbose=>1);

is ref($domain), 'Bio::Tools::HMMER::Domain';

$domain->start(50);
$domain->end(200);
$domain->hstart(10);
$domain->hend(100);
$domain->seqbits(50);
$domain->bits(20);
$domain->evalue(0.0001);
$domain->seq_id('silly');


# test that we can get out forward and reverse homol_SeqFeatures
$homol = $domain->feature2();
is $homol->start(), 10;

$rev = $domain;

is $rev->start(), 50;

$set = Bio::Tools::HMMER::Set->new();
$set->add_Domain($domain);

@doms = $set->each_Domain();
$dom = shift @doms;

is $dom->start(), 50;

$set->bits(300);
$set->evalue(0.0001);
$set->name('sillyname');
$set->desc('a desc');
$set->accession('fakeaccesssion');
is $set->bits(), 300;
is $set->evalue(), 0.0001;
is $set->name(), 'sillyname';
is $set->desc, 'a desc';
is $set->accession, 'fakeaccesssion';

$res = Bio::Tools::HMMER::Results->new( -file => test_input_file('hmmsearch.out') , -type => 'hmmsearch');
my $seen =0;
is $res->hmmfile, "HMM";
is $res->seqfile, "HMM.dbtemp.29591";

my $first = 0;
foreach $set ( $res->each_Set) {
    foreach $domain ( $set->each_Domain ) {
    #print STDERR "Got domain ",$domain->seq_id," start ",$domain->start," end ",$domain->end,"\n";
    # do nothing for the moment
      $seen = 1;
  }
}
is $seen, 1;

is $res->number, 1215;

$res = Bio::Tools::HMMER::Results->new( -file => test_input_file('hmmpfam.out') , 
					-type => 'hmmpfam');

is ($res->number, 2);

# parse HMM 2.2 files

$res = Bio::Tools::HMMER::Results->new( -file => test_input_file('L77119.hmmer'),
					-type => 'hmmpfam');
$seen =0;
is $res->hmmfile, 'Pfam';
is $res->seqfile, 'L77119.faa';
foreach $set ( $res->each_Set) {
    # only one set anyways

    is($set->name, 'gi|1522636|gb|AAC37060.1|');
    is($set->desc, 'M. jannaschii predicted coding region MJECS02 [Methanococcus jannaschii]');
    is($set->accession, '[none]');
    foreach $domain ( $set->each_Domain ) {
	#print STDERR "Got domain ",$domain->seq_id," start ",$domain->start," end ",$domain->end,"\n";
    # do nothing for the moment
	is($domain->start, 280);
	is($domain->end, 481);
	is($domain->bits, -105.2);
	is($domain->evalue, 0.0022 );
    }
}
is ($res->number, 1);

# test for bugs #(1189,1034,1172)
$res = Bio::Tools::HMMER::Results->new( -file => test_input_file('hmmsearch.out') , 
					-type => 'hmmsearch');
my $res2 = $res->filter_on_cutoff(100,50);
ok($res2);
is($res2->number, 604);
