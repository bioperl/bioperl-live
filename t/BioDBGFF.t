#-*-Perl-*-
## Bioperl Test Harness Script for Modules

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

use strict;
use ExtUtils::MakeMaker;
use Bio::Root::IO;
use constant TEST_COUNT => 113;
use constant FASTA_FILES => Bio::Root::IO->catfile('t','data','dbfa');
use constant GFF_FILE    => Bio::Root::IO->catfile('t','data',
						   'biodbgff','test.gff');

BEGIN {
    # to handle systems with no installed Test module
    # we include the t dir (where a copy of Test.pm is located)
    # as a fallback
    eval { require Test; };
    if( $@ ) {
        use lib 't';
    }
    use Test;
    plan test => TEST_COUNT;
}


sub bail ($;$);
sub user_prompt ($;$);
sub fail ($);
use lib './blib/lib';
use Bio::DB::GFF;

my %EXCLUDE_DRIVERS    = ('ExampleP' =>1, 'Multiplex'=>1, 'Proxy'=>1);
my @PREFERRED_DRIVERS  = qw(mysql Pg Oracle Sybase mSQL Solid Informix Illustra CSV);
my $idx = 1;
my %PREFERRED_DRIVERS  = map {$_=>$idx++} @PREFERRED_DRIVERS;

my $adaptor = shift || 'dbi::mysqlopt';
my @args;
if ($adaptor =~ /^dbi/) {
  eval { require DBI } or bail(TEST_COUNT,'DBI driver is missing');
  my $cfg = ChooseDrivers() or bail(TEST_COUNT,'Skipped by user');
  QueryDb($cfg,{'prompt'=>1,'verbose'=>1});
  @args = ( '-adaptor'  => $adaptor,
	     '-dsn'     => $cfg->{test_dsn},
	     );
  push @args,('-user' => $cfg->{test_user}) if $cfg->{test_user};
  push @args,('-pass' => $cfg->{test_pass}) if $cfg->{test_pass};
} else {
  @args = ('-adaptor' => $adaptor);
}

my $db = eval { Bio::DB::GFF->new(@args) };
warn $@ if $@;
ok($db);
fail(TEST_COUNT - 1) unless $db;

# exercise the loader
ok($db->initialize(1));
ok($db->load_gff(GFF_FILE));
ok($db->load_fasta(FASTA_FILES));

# exercise db->types
my @types = sort $db->types;
ok(scalar @types,8);
ok($types[0],'CDS:confirmed');
ok($types[-1],'transposon:tc1');
my %types = $db->types('-enumerate'=>1);
ok($types{'transposon:tc1'},2);

# exercise segment
my $segment1 = $db->segment('Contig1');
ok($segment1);
ok($segment1->length,37450);
ok($segment1->start,1);
ok($segment1->end,37450);
ok($segment1->strand,1);

my $segment2  = $db->segment('Contig1',1=>1000);
ok($segment2->length,1000);
ok($segment2->start,1);
ok($segment2->end,1000);
ok($segment2->strand,1);

my $segment3 = $db->segment('Contig1',10=>1);
ok($segment3->start,10);
ok($segment3->end,1);
ok($segment3->strand,-1);

# exercise attribute fetching
my @t = $db->fetch_feature_by_name(Transcript => 'trans-1');
my ($t) = grep {$_->type eq 'transcript:confirmed'} @t;
ok($t->attributes('Note'),'function unknown');
ok(join(' ',sort $t->attributes('Gene')),'abc-1 xyz-2');
my $att = $t->attributes;
ok(scalar @{$att->{Gene}},2);
@t = $db->fetch_feature_by_attribute('Gene'=>'abc-1');
ok(@t>0);
ok($t[0] eq $t);
my $seg = $db->segment('Contig1');
@t = $seg->features(-attributes=>{'Gene'=>'abc-1'});
ok(@t>0);
@t = $seg->features(-attributes=>{'Gene'=>'xyz-2',Note=>'Terribly interesting'});
  ok(@t==1);

# exercise dna() a bit
my $dna = $segment2->dna;
ok(length $dna,1000);
ok(substr($dna,0,10),'gcctaagcct');
ok($segment3->dna,'aggcttaggc');

# exercise ref()
my $segment4 = $db->segment('-name'=>'c128.1','-class'=>'Transposon');
ok($segment4->length,1000);
ok($segment4->start,1);
ok($segment4->end,1000);
ok($segment4->ref,'c128.1');
ok($segment4->strand,1);
ok(!$segment4->absolute);

$segment4->absolute(1);
ok($segment4->absolute);
ok($segment4->ref,'Contig1');
ok($segment4->start,5001);
$segment4->absolute(0);
my $tmp = $db->segment('Contig1',5001=>6000);
ok($segment4->dna,$tmp->dna);

$segment4->ref('Contig1');
ok($segment4->ref,'Contig1');
ok($segment4->start,5001);
ok($segment4->end,6000);

my $segment5 = $db->segment('-name'=>'c128.2','-class'=>'Transposon');
ok($segment5->length,1000);
ok($segment5->start,1);
ok($segment5->end,1000);
ok($segment5->ref,'c128.2');
ok($segment5->strand,1);

$tmp = $db->segment('Contig1',9000,8001);
ok($segment5->dna,$tmp->dna);
$segment5->absolute(1);
ok($segment5->strand,-1);

# rel/rel addressing
# first two positive strand features
$segment4 = $db->segment('-name'=>'c128.1','-class'=>'Transposon');
my $start4 = $segment4->abs_start;
$segment5  = $db->segment('Transcript' => 'trans-1');
my $start5 = $segment5->abs_start;
$segment4->ref($segment5);
ok($segment4->strand,1);
ok($segment4->start,$start4-$start5+1);
ok($segment4->stop,$start4-$start5+$segment4->length);

$segment4->ref('Transposon' => 'c128.1');
$segment5->ref('Transcript' => 'trans-1');
$segment5->ref($segment4);
ok($segment5->start,$start5-$start4+1);

# now a positive on a negative strand feature
my $segment6 = $db->segment('Transcript'=>'trans-2');
my $start6 = $segment6->abs_start;
ok($segment6->strand,1);
ok($segment6->abs_strand,-1);
$segment6->ref($segment4);
ok($segment6->start,$start6-$start4+1);
ok($segment6->strand,-1);

$segment4->ref($segment6);
ok($segment4->start,$start6-$start4+1);
ok($segment4->strand,-1);
ok($segment4->ref eq $segment6);

# the reference sequence shouldn't affect the dna
$segment6 = $db->segment('Transcript'=>'trans-2');
$dna = $segment6->dna;
$segment6->ref($segment4);
ok($segment6->dna,$dna);

# segments should refuse to accept a reference sequence on a foreign segment
undef $@;
my $result = eval { $segment6->ref('Contig2') };
ok(!$result);
ok("$@" =~ /are on different sequence segments/);

# types across a segment
$segment1 = $db->segment('Contig1');
@types = sort $segment1->types;
ok(scalar @types,6);
ok($types[0],'CDS:confirmed');
ok($types[-1],'transposon:tc1');
%types = $segment1->types('-enumerate'=>1);
ok($types{'similarity:est'},3);

# features across a segment
my @features = $segment1->features('-automerge'=>0);
ok(scalar @features,17);
my %types_seen;
foreach (@features) {
  $types_seen{$_->type}++;
}
my $inconsistency = 0;
foreach (keys %types,keys %types_seen) {
  $inconsistency++ unless $types_seen{$_} == $types{$_};
}
ok(!$inconsistency);

@features = sort {$a->start<=>$b->start} @features;
ok($features[0]->type,'Component:reference');
ok($features[-1]->type,'exon:confirmed');

# make sure that we can use features to get at dna
ok($features[0]->dna,$db->segment('Contig1',$features[0]->start,$features[0]->end)->dna);

# check three forward features and three reverse features
# (This depends on the test.gff data)
for (1..3,-3..-1) {
  $segment2 = $db->segment($features[$_],50,100);
  if ($features[$_]->strand >= 0) {
    ok($segment2->dna,$db->segment('Contig1',
				   $features[$_]->start+50-1,
				   $features[$_]->start+100-1)->dna)
  } else {
    ok($segment2->dna,$db->segment('Contig1',
				   $features[$_]->start-50+1,
				   $features[$_]->start-100+1)->dna)
  }
}

# exercise the aggregator
my $aggregator = Bio::DB::GFF::Aggregator->new('-method'      => 'aggregated_transcript',
					       '-main_method' => 'transcript',
					       '-sub_parts'   => ['exon','CDS']);
$db->add_aggregator($aggregator);
$segment1 = $db->segment('Contig1');
@features = $segment1->features('aggregated_transcript');
ok(scalar @features,2);
ok($features[0]->Exon > 0);
ok($features[0]->Cds > 0);

# Test that sorting is correct.  The way that test.gff is set up, the lower one is
# on the + strand and the higher is on the -.
@features = sort {$a->start <=> $b->start} @features;
ok($features[0]->strand,1);
ok($features[1]->strand,-1);

my $last = 0;
$inconsistency = 0;
foreach ($features[0]->Exon) {
  $inconsistency++ if $_->start > $_->end;
  $inconsistency++ if $last && $_->start < $last;
  $last = $_->start;
}
ok(!$inconsistency);

$inconsistency = $last = 0;
foreach ($features[1]->Exon) {
  $inconsistency++ if $_->start < $_->end;
  $inconsistency++ if $last && $_->start > $last;
  $last = $_->start;
}
ok(!$inconsistency);

# relative addressing in aggregated features
my $transcript1 = $db->segment($features[0]);
$transcript1->ref($features[0]);
my @overlap     = sort {$a->start <=> $b->start } $transcript1->features;
ok(scalar(@overlap),11);
ok($overlap[0]->start,-999);

$transcript1 = $db->segment('Transcript' => 'trans-1');
@overlap     = sort {$a->start <=> $b->start } $transcript1->features;
ok($overlap[0]->start,-999);

# test strandedness of features
$segment1 = $db->segment('-class' => 'Transcript',
			 '-name'  => 'trans-3',
			 '-start' => 1,
			 '-stop'  => 6000);
ok($segment1->strand,1);
@overlap  = sort {$a->start <=> $b->start} $segment1->features('transcript');
ok(scalar(@overlap),2);
ok($overlap[0]->name,'trans-3');
ok($overlap[1]->name,'trans-4');
ok($overlap[0]->strand,1);
ok($overlap[1]->strand,-1);

# testing feature id and group_id
my $tf = $overlap[0];
ok(defined $tf->id);
ok(defined $tf->group_id);
my $t1 = $db->fetch_feature_by_id($tf->id);
ok($t1->id,$tf->id);
my $t2 = $db->fetch_feature_by_gid($tf->group_id);
ok($t2->group_id,$tf->group_id);
ok($t2->group_id,$t1->group_id);

$segment1 = $db->segment('-class' => 'Transcript',
			 '-name'  => 'trans-4',
			 '-start' => 1,
			 '-stop'  => 6000);
ok($segment1->strand,1);
@overlap = sort {$a->start <=> $b->start} $segment1->features('transcript');
ok($overlap[0]->name,'trans-4');
ok($overlap[1]->name,'trans-3');
ok($overlap[0]->strand,1);
ok($overlap[1]->strand,-1);

@overlap = sort {$a->start <=> $b->start} $segment1->features('Component');
ok($overlap[0]->strand,0);

# test iterator across a segment
$segment1 = $db->segment('Contig1');
my $i = $segment1->features('-automerge'=>0,'-iterator'=>1);
my %strand;
while (my $s = $i->next_feature) {
  $strand{$s->strand}++;
}
ok(keys %strand == 3);

# test iterator across entire database
$i = $db->features('-automerge'=>0,'-iterator'=>1);
%strand = ();
while (my $s = $i->next_feature) {
  $strand{$s->strand}++;
}
ok(keys %strand == 3);

# test iterator across a segment, limited by an attribute
$i = $seg->get_feature_stream(-attributes=>{'Gene'=>'abc-1',Note=>'function unknown'});
my $count = 0;
while ($i->next_seq) {
  $count++;
}
ok($count,2);

# test that aliases work
my $st1 = $db->segment(Transcript => 'trans-3');
ok($st1);
my $st2 = $db->segment(Transcript => 'trans-18');  # this is an alias!
ok($st2);
ok($st1 eq $st2);
my @transcripts = $st1->features('transcript');
ok(($transcripts[0]->aliases)[0] eq 'trans-18');

END {
  unlink FASTA_FILES."/directory.index";
}

############################################################################
#
#   Name:    ChooseDrivers
#
#   Purpose: Choose drivers for test suite
#
#   Inputs:  none
#
#   Returns: configuration hash
#
############################################################################
sub ChooseDrivers {
  my $cfg = shift || {};
  my @drivers = grep {!$EXCLUDE_DRIVERS{$_}} DBI->available_drivers;
  my @preferred_drivers = sort {($PREFERRED_DRIVERS{$a}||99) <=> ($PREFERRED_DRIVERS{$b}||99)} @drivers;
  print STDERR <<END;

*** Bio::DB::GFF TEST SETUP ***
    Bio::DB::GFF requires a running database server (preferably MySQL).
    To run this test suite you must provide the name and host of a test database
    for which you have insert privileges.
END
;
  my $proceed = user_prompt("Do you want to run the tests? y/n",'n');
  return unless $proceed =~/^[yY]/;
  my $choice = user_prompt("Which database driver do you wish to use (@preferred_drivers)",$preferred_drivers[0]);
  $cfg->{dbd_driver} = $choice;
  $cfg;
}

############################################################################
#
#   Name:    QueryDb
#
#   Purpose: Query settings for running the test suite
#
#   Inputs:  $cfg - Config hash ref
#
#   Returns: Nothing; creates
#            $cfg->{$driver}->{test_(db|host|dsn|user|pass)}
#
############################################################################

sub QueryDb {
  my $cfg     = shift;
  my $options = shift || {};
  my $db     = $cfg->{'description'} || $cfg->{'dbd_driver'};
  my $driver = $cfg->{'dbd_driver'};

  my $d = lc $driver;
  my $prompt = $options->{'prompt'};

  my $test_db =  exists($options->{"$d-test-db"}) ?
    $options->{"$d-test-db"} : ($cfg->{'test_db'} || 'test');
  $test_db = user_prompt
    ("Which database should I use for testing the $db drivers?",
     $test_db) if $prompt;

  my $test_host = exists($options->{"$d-test-host"}) ?
    $options->{"$d-test-host"} : ($cfg->{'test_host'} || 'localhost');
  $test_host = user_prompt
    ("On which host is database $test_db running (hostname, ip address or host:port)", $test_host) if $prompt;

  my($test_user, $test_pass);

  $test_user = exists($options->{"$d-test-user"}) ?
    $options->{"$d-test-user"} : ($cfg->{'test_user'} || "undef");
  $test_user = user_prompt
    ("User name for connecting to database $test_db?", $test_user)
      if $prompt;
  $test_user = undef if $test_user eq 'undef';

  $test_pass = exists($options->{"$d-test-pass"}) ?
    $options->{"$d-test-pass"} : ($cfg->{'test_pass'} || "undef");
  $test_pass = user_prompt
    ("Password for connecting to database $test_db?", $test_pass)
      if $prompt;
  $test_pass = undef if $test_pass eq 'undef';

  $cfg->{'test_db'}   = $test_db;
  $cfg->{'test_host'} = $test_host;
  $cfg->{'test_user'} = $test_user;
  $cfg->{'test_pass'} = $test_pass;
  if ($test_host eq 'undef'  ||  $test_host eq 'localhost') {
    $test_host = '';
  }
  my $test_dsn = "DBI:$driver:database=$test_db";
  $cfg->{test_dsn} = $test_dsn;
  if ($test_host) {
    $cfg->{'test_dsn'} .= ";host=$test_host";
  }


  if ($options->{'verbose'}) {
    local $^W=0;  # some uninitialized variable warning coming through
    print("Driver $driver is using the following settings for tests:\n",
	  "    Database $test_db\n",
	  "    Host     $test_host\n",
	  "    DSN      $test_dsn\n",
	  "    User     $test_user\n",
	  "    Password $test_pass\n");
  }
}

sub bail ($;$) {
  my $count = shift;
  my $explanation = shift;
  for (1..$count) {
    skip($explanation,1);
  }
  exit 0;
}

sub fail ($) {
  my $count = shift;
  for (1..$count) {
    ok(0);
  }
  exit 0;
}

sub user_prompt ($;$) {
    my($mess,$def)=@_;
    Carp::confess("prompt function called without an argument") unless defined $mess;
    my $dispdef = defined $def ? "[$def] " : " ";
    $def = defined $def ? $def : "";
    my $ans;
    local $|=1;
    print STDERR "$mess $dispdef";
    chomp($ans = <STDIN>);
    return ($ans ne '') ? $ans : $def;
}
