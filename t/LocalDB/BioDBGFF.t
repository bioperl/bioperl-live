# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;
use Module::Build;
use Data::Dumper;

BEGIN {
    use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 275);
	
	use_ok('Bio::DB::GFF');
}

my $fasta_files = test_input_file('dbfa');
my $gff_file1   = test_input_file('biodbgff', 'test.gff');
my $gff_file2   = test_input_file('biodbgff', 'test.gff3');

my $build = Module::Build->current;
my $test_dsn = $build->notes('test_dsn');

my $adaptor = $test_dsn ? $test_dsn : 'memory';
$adaptor    = shift if @ARGV;

if ($adaptor =~ /sqlite/i) {
    $adaptor = 'memory';
}

my @args;
if ($adaptor =~ /^dbi/) {
  my $cfg = {};
  $cfg->{dbd_driver} = $build->notes('dbd_driver');
  $cfg->{test_db} = $build->notes('test_db');
  $cfg->{test_host} = $build->notes('test_host');
  $cfg->{test_user} = $build->notes('test_user');
  $cfg->{test_pass} = $build->notes('test_pass');
  $cfg->{test_dsn} = $build->notes('test_dsn');
  
  $adaptor = "dbi::$cfg->{dbd_driver}" if $cfg->{dbd_driver};
  @args = ( '-adaptor'  => $adaptor,
	    '-dsn'     => $cfg->{test_dsn},
	  );
  push @args,('-user' => $cfg->{test_user}) if $cfg->{test_user};
  push @args,('-pass' => $cfg->{test_pass}) if $cfg->{test_pass};
} else {
  @args = ('-adaptor' => $adaptor,
	   '-create'  => 1);
}

push @args,('-aggregators' => ['transcript','processed_transcript']);

SKIP: {
for my $FILE ($gff_file1,$gff_file2) {

  my $db = eval { Bio::DB::GFF->new(@args) };
  skip "DB load failed? Skipping all! $@", 278 if $@;
  ok($db);

  $db->debug(0);
  $db->gff3_name_munging(1);

  # set the preferred groups
  $db->preferred_groups( [ 'transcript', 'gene', 'mRNA' ] );
  my @pg = $db->preferred_groups;
  is(scalar(@pg), 3);
  is($pg[1], 'gene'); 

  # exercise the loader
  ok($db->initialize(1));
  ok($db->load_gff($FILE));
  ok($db->load_fasta($fasta_files));

  # exercise db->types
  my @types = sort $db->types;
  is(scalar @types,11);
  is($types[0],'CDS:confirmed');
  is($types[-1],'transposon:tc1');
  my %types = $db->types('-enumerate'=>1);
  is($types{'transposon:tc1'},2);

  # exercise segment
  my $segment1 = $db->segment('Contig1');

  ok($segment1);
  is($segment1->length,37450);
  is($segment1->start,1);
  is($segment1->end,37450);
  is($segment1->strand,1);
  
  my $segment2  = $db->segment('Contig1',1=>1000);
  is($segment2->length,1000);
  is($segment2->start,1);
  is($segment2->end,1000);
  is($segment2->strand,1);
  
  my $segment3 = $db->segment('Contig1',10=>1);
  is($segment3->start,10);
  is($segment3->end,1);
  is($segment3->strand,-1);

  # exercise attribute fetching
  my @t = $db->fetch_feature_by_name(Transcript => 'trans-1');
  my ($t) = grep {$_->type eq 'transcript:confirmed'} @t;
  is($t->attributes('Note'),'function unknown');
  is(join(' ',sort $t->attributes('Gene')),'abc-1 xyz-2');
  my $att = $t->attributes;
  is(scalar @{$att->{Gene}},2);
  @t = sort {$a->display_name cmp $b->display_name} $db->fetch_feature_by_attribute('Gene'=>'abc-1');
  cmp_ok(@t,'>',0);
  is($t[0], $t);
  my $seg = $db->segment('Contig1');
  @t = $seg->features(-attributes=>{'Gene'=>'abc-1'});
  cmp_ok(@t,'>',0);
  is($seg->feature_count, 17);
  @t = $seg->features(-attributes=>{'Gene'=>'xyz-2',Note=>'Terribly interesting'});
  is(@t,1);

  # exercise dna() a bit
  my $dna = $segment2->dna;
  is(length $dna,1000);
  is(substr($dna,0,10),'gcctaagcct');
  is($segment3->dna,'aggcttaggc');
  is($segment1->dna, $db->dna($segment1->ref));

  # exercise ref()
  my $segment4 = $db->segment('-name'=>'c128.1','-class'=>'Transposon');
  is($segment4->length,1000);
  is($segment4->start,1);
  is($segment4->end,1000);
  is($segment4->ref,'c128.1');
  is($segment4->strand,1);
  ok(!$segment4->absolute);

  $segment4->absolute(1);
  ok($segment4->absolute);
  is($segment4->ref,'Contig1');
  is($segment4->start,5001);
  $segment4->absolute(0);
  my $tmp = $db->segment('Contig1',5001=>6000);
  is($segment4->dna,$tmp->dna);

  $segment4->ref('Contig1');
  is($segment4->ref,'Contig1');
  is($segment4->start,5001);
  is($segment4->end,6000);

  my $segment5 = $db->segment('-name'=>'c128.2','-class'=>'Transposon');
  is($segment5->length,1000);
  is($segment5->start,1);
  is($segment5->end,1000);
  is($segment5->ref,'c128.2');
  is($segment5->strand,1);

  $tmp = $db->segment('Contig1',9000,8001);
  is($segment5->dna,$tmp->dna);
  $segment5->absolute(1);
  is($segment5->strand,-1);

  # rel/rel addressing
  # first two positive strand features
  $segment4 = $db->segment('-name'=>'c128.1','-class'=>'Transposon');
  my $start4 = $segment4->abs_start;
  $segment5  = $db->segment('Transcript' => 'trans-1');
  my $start5 = $segment5->abs_start;
  $segment4->ref($segment5);
  is($segment4->strand,1);
  is($segment4->start,$start4-$start5+1);
  is($segment4->stop,$start4-$start5+$segment4->length);

  $segment4->ref('Transposon' => 'c128.1');
  $segment5->ref('Transcript' => 'trans-1');
  $segment5->ref($segment4);
  is($segment5->start,$start5-$start4+1);

  # now a positive on a negative strand feature
  my $segment6 = $db->segment('Transcript'=>'trans-2');
  my $start6 = $segment6->abs_start;
  is($segment6->strand,1);
  is($segment6->abs_strand,-1);
  $segment6->ref($segment4);
  is($segment6->start,$start6-$start4+1);
  is($segment6->strand,-1);

  $segment4->ref($segment6);
  is($segment4->start,$start6-$start4+1);
  is($segment4->strand,-1);
  is($segment4->ref,$segment6);

  # the reference sequence shouldn't affect the dna
  $segment6 = $db->segment('Transcript'=>'trans-2');
  $dna = $segment6->dna;
  $segment6->ref($segment4);
  is($segment6->dna,$dna);

  # segments should refuse to accept a reference sequence on a foreign segment
  undef $@;
  my $result = eval { $segment6->ref('Contig2') };
  ok(!$result);
  like($@, qr/are on different sequence segments/);

  # types across a segment
  $segment1 = $db->segment('Contig1');
  @types = sort $segment1->types;
  is(scalar @types,6);
  is($types[0],'CDS:confirmed');
  is($types[-1],'transposon:tc1');
  %types = $segment1->types('-enumerate'=>1);
  is($types{'similarity:est'},3);

  # features across a segment
  my @features = $segment1->features('-automerge'=>0);
  is(scalar @features,17);
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

  # make sure that we can use features to get at dna
  is($features[0]->dna,$db->segment('Contig1',$features[0]->start,$features[0]->end)->dna);

  # check three forward features and three reverse features
  # (This depends on the test.gff data)
  for (1..3,-3..-1) {
    $segment2 = $db->segment($features[$_],50,100);
    if ($features[$_]->strand >= 0) {
      is($segment2->dna,$db->segment('Contig1',
				     $features[$_]->start+50-1,
				     $features[$_]->start+100-1)->dna)
    } else {
      is($segment2->dna,$db->segment('Contig1',
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
  @features = sort $segment1->features('aggregated_transcript');   # sort so that trans-1 comes first
  is(scalar @features,2);
  cmp_ok($features[0]->Exon, '>', 0);
  cmp_ok($features[0]->Cds,'>', 0);

  # Test that sorting is correct.  The way that test.gff is set up, the lower one is
  # on the + strand and the higher is on the -.
  @features = sort {$a->start <=> $b->start} @features;
  is($features[0]->strand,1);
  is($features[1]->strand,-1);

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
  is(scalar(@overlap),5);
  is($overlap[0]->start,-999);

  $transcript1 = $db->segment('Transcript' => 'trans-1');
  @overlap     = sort {$a->start <=> $b->start } $transcript1->features;
  is($overlap[0]->start,-999);

  # test strandedness of features
  $segment1 = $db->segment('-class' => 'Transcript',
			   '-name'  => 'trans-3',
			   '-start' => 1,
			   '-stop'  => 6000);
  is($segment1->strand,1);
  @overlap  = sort {$a->start <=> $b->start} $segment1->features('transcript');
  is(scalar(@overlap),2);
  is($overlap[0]->name,'trans-3');
  is($overlap[1]->name,'trans-4');
  is($overlap[0]->strand,1);
  is($overlap[1]->strand,-1);

  # testing feature id and group_id
  my $tf = $overlap[0];
  ok(defined $tf->id);
  my $t1 = $db->fetch_feature_by_id($tf->id);
  is($t1->id,$tf->id);

  SKIP: {
    if (defined $tf->group_id) {
      my $t2 = $db->fetch_feature_by_gid($tf->group_id);
      is($t2->group_id,$tf->group_id);
      is($t2->group_id,$t1->group_id);
    } else {
      skip("fetch_feature_by_gid() not implemented by this adaptor",2);
    }
  }

  $segment1 = $db->segment('-class' => 'Transcript',
			   '-name'  => 'trans-4',
			   '-start' => 1,
			   '-stop'  => 6000);
  is($segment1->strand,1);
  @overlap = sort {$a->start <=> $b->start} $segment1->features('transcript');
  is($overlap[0]->name,'trans-4');
  is($overlap[1]->name,'trans-3');
  is($overlap[0]->strand,1);
  is($overlap[1]->strand,-1);

  @overlap = sort {$a->start <=> $b->start} $segment1->features('Component');
  is($overlap[0]->strand,0);

SKIP: {
  # test preferred group assignments
  if ($FILE =~ /\.gff$/) {
    my @gene = $db->get_feature_by_name( gene => 'gene-9' );
    my @mrna = $db->get_feature_by_name( mRNA => 'trans-9' );
    is($gene[0]->ref, 'Contig4');
    is(scalar(@gene), 2);
    is(scalar(@mrna), 1);
  } else {
    skip('preferred groups are not supported by gff3',3);
  }
}

  # test iterator across a segment
  $segment1 = $db->segment('Contig1');
  my $i = $segment1->features('-automerge'=>0,'-iterator'=>1);
  my %strand;
  while (my $s = $i->next_feature) {
    $strand{$s->strand}++;
  }
  is(keys %strand, 3);

  # test iterator across entire database
  $i = $db->features('-automerge'=>0,'-iterator'=>1);
  %strand = ();
  while (my $s = $i->next_feature) {
    $strand{$s->strand}++;
  }
  is(keys %strand, 3);

  # test iterator across a segment, limited by an attribute
  $i = $seg->get_feature_stream(-attributes=>{'Gene'=>'abc-1',Note=>'function unknown'});
  my $count = 0;
  while ($i->next_seq) {
    $count++;
  }
  is($count,2);

  # test that aliases work
  my $st1 = $db->segment(Transcript => 'trans-3');
  ok($st1);
  my $st2 = $db->segment(Transcript => 'trans-18');  # this is an alias!
  ok($st2);
  is($st1,$st2);
  my @transcripts = $st1->features('transcript');
  is(($transcripts[0]->aliases)[0],'trans-18');

  # test truncation
  $db->strict_bounds_checking(1);
  my $tseg = $db->segment(-name=>'trans-1',-class=>'Transcript',-start=>1,-stop=>500);
  ok(!$tseg->truncated);
  $tseg    = $db->segment(-name=>'trans-1',-class=>'Transcript',-start=>1,-stop=>50000);
  ok($tseg->truncated);
  $db->strict_bounds_checking(0);
  $tseg    = $db->segment(-name=>'trans-1',-class=>'Transcript',-start=>1,-stop=>50000);
  ok(!$tseg->truncated);

  # test the processed_transcript aggregator
  $db->clear_aggregators;
  $db->add_aggregator('processed_transcript');
  my @f = $db->fetch_feature_by_name(mRNA => 'trans-8');
  is(scalar @f,1);
  is($f[0]->length,35000-32000+1);
  is(scalar $f[0]->CDS,3);
  is(scalar $f[0]->UTR,2);

  # test deletions
  # segment delete() method
  my $clone = $db->segment(Clone=>'M7.3');
  my $overlapping_feature_count = $clone->features(-range_type =>'overlaps');
  my $contained_feature_count   = $clone->features(-range_type =>'contains');
  is(scalar $clone->delete(-range_type=>'contains'),$contained_feature_count);
  is(scalar $clone->features,$overlapping_feature_count - $contained_feature_count);

  # database delete() method
  is($db->delete(-type=>['mRNA:confirmed','transposon:tc1']),4);
  is($db->delete(-type=>'UTR',-ref=>'Contig29'),undef);
  is($db->delete(-type=>'CDS',-ref=>'AL12345.2',-class=>'Clone'),3);
  is($db->delete_features(1,2,3),3);

  SKIP: {
    $result = eval {
      is($db->delete_groups(1,2,3,4,5),5);
      my @features = $db->get_feature_by_name(Sequence => 'Contig2');
      is($db->delete_groups(@features),1);
      1;
    };
    if (!$result && $@ =~ /not implemented/i) {
      skip("delete_groups() not implemented by this adaptor",2);
    }
  }
  
  SKIP: {
	test_skip(-tests => 1, -excludes_os => 'mswin');
	
	# test ability to pass adaptors across a fork
	if (my $child = open(F,"-|")) { # parent reads from child
		ok(scalar <F>);
		close F;
	}
	else { # in child
		$db->clone;
		my @f = $db->features();
		print @f>0;
		exit 0;
	}
  }

  ok(!defined eval{$db->delete()});
  ok($db->delete(-force=>1));
  is(scalar $db->features,0);
  ok(!$db->segment('Contig1'));

}

}



END {
  unlink $fasta_files."/directory.index";
}
