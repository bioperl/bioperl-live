# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
    use lib '.';
    use Bio::Root::Test;

    test_begin(-tests => 1427,
               -requires_module => 'Bio::ASN1::EntrezGene');

    use_ok('Bio::SeqIO::entrezgene');
}

my @species=('Homo sapiens','Mus musculus', 'Caenorhabditis elegans');
my @pubmed=qw(15461460 15221005 14702039 12477932 8889549 3610142 3458201 2591067);

my %pmed=(1=>8, 
            2=>55,
            3=>1,
            4=>0,
            5=>0,
            6=>0,
            7=>0,
            8=>1,
            9=>32,
            10=>58,
            11=>1,
            12=>76,
            13=>7,
            14=>5,
            15=>13,
            9996=>0,
            11286=>0,
            11287=>5,
            11288=>0,
            11289=>0,
            11293=>0,
            11294=>0,
            11295=>0,
            11296=>0,
            11297=>0,
            11298=>3,
            11299=>0,
            11300=>0,
            11301=>0,
            11302=>9,
            11303=>54,
            11304=>11,
            11305=>3,
            11306=>9,
            171590=>0,
            171591=>0,
            171592=>0,
            171593=>0,
            171594=>0);
            
my %asym=(1=>['A1B', 'ABG', 'GAB', 'HYST2477', 'DKFZp686F0970'],
            2=>['FWP007','S863-7','DKFZp779B086'], 4=>['A12M1'], 5=>['A12M2'],6=>['A12M3'],7=>['A12M4'],
            9=>['AAC1'],10=>['AAC2'],11=>['NATP'],
            12=>['ACT','AACT','MGC88254'],13=>['DAC'],15=>['SNAT','AA-NAT'],
            14=>[''],
            11287=>['A1m','A2m','MAM'],
            11298=>['Nat4','SNAT','Nat-2'],
            11302=>['AATYK','mKIAA0641'],11303=>['Abc1'],
            11304=>['RmP','Abcr','Abc10','D430003I15Rik'],
            11305=>['Abc2','mKIAA1062','D2H0S1474E'],
            11306=>['Abc7'],
            171590=>['Y74C9A.3','CELK05052'],
            171591=>['Y74C9A.2','CELK01753'],
            171592=>['Y74C9A.4a','Y74C9A.4b','CELK08126'],
            171593=>['Y74C9A.5','CELK09643'],
            171594=>['Y48G1C.4','CELK05819']);
            
my @ids=qw(1
2
3
4
5
6
7
8
9
10
11
12
13
14
15
9996
11286
11287
11288
11289
11293
11294
11295
11296
11297
11298
11299
11300
11301
11302
11303
11304
11305
11306
171590
171591
171592
171593
171594);

my @loop_counts = ([1,1,1,1,5,1,12,1,1,1,14,8,16],
				   [1,1,1,1,3,1,40,1,1,1,14,31],
				   [1,1,1,1,0,1,4,1,10,5],
				   [1,1,0,1,1,1,1,7,0],
				   [1,1,0,1,1,1,1,7,0],
				   [1,1,0,1,1,1,1,7,0],
				   [1,1,0,1,1,1,1,7,0],
				   [1,0,1,1,0,1,4,1,1,1,13,0],
				   [1,1,1,1,1,1,33,1,1,1,14,41],
				   [1,1,1,1,1,1,51,1,1,1,14,51],
				   [1,0,1,1,1,1,1,1,10,1],
				   [1,1,1,1,3,1,28,1,1,1,14,33],
				   [1,1,1,1,1,1,17,1,1,1,14,10],
				   [1,1,1,1,0,1,11,1,1,1,13,20],
				   [1,1,1,1,2,1,16,1,1,1,14,23],
				   [1,0,0,0,0,0,0,2,0],
				   [1,0,0,0,0,0,0,3,0],
				   [1,1,1,1,3,1,10,1,13,10],
				   [1,0,0,0,0,0,0,3,0],
				   [1,0,0,0,0,0,0,2,0],
				   [1,0,0,0,0,0,0,3,0],
				   [1,0,0,0,0,0,0,3,0],
				   [1,0,0,0,0,0,0,3,0],
				   [1,0,0,0,0,0,0,2,0],
				   [1,0,0,0,0,0,0,2,0],
				   [1,1,1,1,3,1,9,1,13,5,16],
				   [1,0,0,0,0,0,0,3,0],
				   [1,0,0,0,0,0,0,3,0],
				   [1,0,0,0,0,0,0,2,0],
				   [1,1,1,1,2,1,10,1,12,19],
				   [1,1,1,1,1,1,50,1,13,14],
				   [1,1,1,1,4,1,9,1,13,12],
				   [1,1,1,1,3,1,9,1,13,8,8],
				   [1,1,1,1,1,1,11,1,13,12],
				   [1,1,0,1,2,0,0,1,1,1,1,9,4],
				   [1,1,0,1,2,0,0,1,1,1,1,9,4],
				   [1,1,0,1,3,0,0,1,1,1,1,9,8],
				   [1,1,0,1,2,0,0,1,1,1,1,9,4],
				   [1,1,0,1,2,0,0,1,1,1,1,9,4]);

my $fs='!';
my @revkeys=('Entrez Gene Status','RefSeq status','Official Full Name','chromosome','cyto','Reference','dblink',
'ALIAS_SYMBOL','OntologyTerm','Index terms','Official Symbol','cM','Property');

ok my $eio=Bio::SeqIO->new(-file=>test_input_file('entrezgene.dat'), -format=>'entrezgene', -debug=>'on',-service_record=>'yes');

my ($seq,$struct,$uncapt);
my $num_of_seqs = 0;
while (1) {
	my $seq;
	($seq,$struct,$uncapt)=$eio->next_seq;
	last unless ($seq);
	my @lc = @{$loop_counts[$num_of_seqs]};
	$num_of_seqs++;
	
	#T0: GENERAL TESTS
	ok $seq;
	is ref($struct),'Bio::Cluster::SequenceFamily';
	my $acc=$seq->accession_number;
	
	#T1: ORGANISM
	my $org=$seq->species->binomial;
	is grep(/\b$org\b/,@species),1;
	
	#T2: SUMMARY test
	ok $seq->desc if ($acc eq '1')||($acc eq '2')||($acc eq '11304');
	ok !defined $seq->desc if ($acc eq '171592')||($acc eq '11306');
	
	#Are we supposed to have this in our test?
	ok grep(/\b$acc\b/,@ids);
	
	my $ann=$seq->annotation();
	my $tcount;
	
	#T3: ENTREZGENE STATUS TESTS
	my @egstatus=$ann->get_Annotations('Entrez Gene Status');
	my $loop_count = 0;
	foreach my $status (@egstatus) {
		$loop_count++;
		STATUS: {
			if ($acc==1) {is $status->value,'live'; last STATUS;}
			if ($acc==2) {is $status->value,'live'; last STATUS;}
			if ($acc==4) {is $status->value,'discontinued'; last STATUS;}
			if ($acc==6) {is $status->value,'discontinued'; last STATUS;}
			if ($acc==11288) {is $status->value,'secondary'; last STATUS;}
			if ($acc==11293) {is $status->value,'secondary'; last STATUS;} 
			if ($acc==171594) {is $status->value,'live'; last STATUS;} 
		}
	}
	is $loop_count, shift @lc, "correct number of loops for T3";
	$loop_count = 0;
	
	#T4: REFSEQ STATUS TESTS
	my @refstatus=$ann->get_Annotations('RefSeq status');
	foreach my $status (@refstatus) {
		$loop_count++;
		STATUS: {
			if ($acc==1) {is $status->value,'REVIEWED'; last STATUS;}
			if ($acc==2) {is $status->value,'REVIEWED'; last STATUS;}
			if ($acc==3) {is $status->value,'PROVISIONAL'; last STATUS;}
			if ($acc==4) {is $status->value,'WITHDRAWN'; last STATUS;}
			if ($acc==9) {is $status->value,'VALIDATED'; last STATUS;}
			if ($acc==11300) {is $status->value,''; last STATUS;}
			if ($acc==11306) {is $status->value,'MODEL'; last STATUS;}
			if ($acc==11293) {is $status->value,'secondary'; last STATUS;} 
			if ($acc==171594) {is $status->value,'Reviewed'; last STATUS;} 
		}
	}
	is $loop_count, shift @lc, "correct number of loops for T4";
	$loop_count = 0;
	
	#T5: GENE NAME TESTS
	my @ofname=$ann->get_Annotations('Official Full Name');
	foreach my $name (@ofname) {
		$loop_count++;
		STATUS: {
			if ($acc==10) {is $name->value,'N-acetyltransferase 2 (arylamine N-acetyltransferase)'; last STATUS;}
			if ($acc==13) {is $name->value,'arylacetamide deacetylase (esterase)'; last STATUS;}
			if ($acc==14) {is $name->value,'angio-associated, migratory cell protein'; last STATUS;}
			if ($acc==11287) {is $name->value,'pregnancy zone protein'; last STATUS;}
			if ($acc==11298) {is $name->value,'arylalkylamine N-acetyltransferase'; last STATUS;}
			if ($acc==11304) {is $name->value,'ATP-binding cassette, sub-family A (ABC1), member 4'; last STATUS;}
			if ($acc==11306) {is $name->value,'ATP-binding cassette, sub-family B (MDR/TAP), member 7'; last STATUS;} 
		}
	}
	is $loop_count, shift @lc, "correct number of loops for T5";
	$loop_count = 0;
	
	#T6: CHROMOSOME TESTS
	my @chr=$ann->get_Annotations('chromosome');
	foreach my $chr (@chr) {
		$loop_count++;
		STATUS: {
			if ($acc==5) {is $chr->value,1; last STATUS;}
			if ($acc==6) {is $chr->value,1; last STATUS;}
			if ($acc==7) {is $chr->value,17; last STATUS;}
			if ($acc==11306) {is $chr->value,'X'; last STATUS;}
			if ($acc==11304) {is $chr->value,3; last STATUS;}
			if ($acc==171590) {is $chr->value,'I'; last STATUS;}
			if ($acc==171592) {is $chr->value,'I'; last STATUS;} 
		}
	}
	is $loop_count, shift @lc, "correct number of loops for T6";
	$loop_count = 0;
	
	#T7: GENE SYMBOL ALIAS TESTS
	my @sym=$ann->get_Annotations('ALIAS_SYMBOL');
	foreach my $sym (@sym) {
		$loop_count++;
        my $val = $sym->display_text;
		next if (($val eq '')||!defined($val));
		is grep(/\b$val\b/,@{$asym{$acc}}),1;
	}
	is $loop_count, shift @lc, "correct number of loops for T7";
	$loop_count = 0;
	
	#T8: CYTO LOCATION TESTS
	my @map=$ann->get_Annotations('cyto');
	foreach my $map (@map) {
		$loop_count++;
		STATUS: {
			if ($acc==10) {is $map->value,'8p22'; last STATUS;}
			if ($acc==11) {is $map->value,'8p22'; last STATUS;}
			if ($acc==13) {is $map->value,'3q21.3-q25.2'; last STATUS;}
			if ($acc==11306) {is $map->value,'X C-D'; last STATUS;}
			if ($acc==11305) {is $map->value,'2 A2-B'; last STATUS;}
			if ($acc==11304) {is $map->value,'3 G1'; last STATUS;}
			if ($acc==11303) {is $map->value,'4 A5-B3'; last STATUS;} 
		}
	}
	is $loop_count, shift @lc, "correct number of loops for T8";
	$loop_count = 0;
	
	#T9: REFERENCE NUMBER TEST
	my @refs=$ann->get_Annotations('Reference');
	my $refs=$#refs+1||0;
	is $pmed{$acc},$refs;
	
	my @dblinks=$ann->get_Annotations('dblink');
	my @keys=$ann->get_all_annotation_keys;
	
	#T10: GENERIF AND OTHER DBLINK TESTS
	my @url=qw(HGMD Ensembl KEGG Homologene);#Only validate the URL
	foreach my $dblink (@dblinks) {
		$loop_count++;
		my $dbname=$dblink->database||'';
		DB: {
			if ( $dbname eq 'generif') {#Should have ID and text
				ok $dblink->primary_id;
				ok $dblink->comment->text;
				last DB;
			}
			if ($acc==2) {
				if (($dbname eq 'MIM')&&($dblink->authority)&&($dblink->authority eq 'phenotype')) {
					ok $dblink->optional_id;
					last DB;
				}
				if ($dbname eq 'Evidence viewer') {
					ok $dblink->url; #We may even validate the urls?
					is $dblink->primary_id,2;
					last DB;
				}
				if ($dbname eq 'Model maker') {
					ok $dblink->url; #We may even validate the urls?
					is $dblink->primary_id,2;
					last DB;
				}
				if ($dbname eq 'AceView') {
					ok $dblink->url; #We may even validate the urls?
					is $dblink->primary_id,2;
					last DB; 
				}
				if (grep(/$dbname/,@url)) {
					ok $dblink->url; #We may even validate the urls?
					last DB;
				}
				if ($dbname eq 'GDB') {
					is $dblink->primary_id,'GDB:119639'; #We may even validate the urls?
					last DB;
				}
				if ($dbname eq 'UniGene') {
					ok $dblink->url; #We may even validate the urls?
					is $dblink->primary_id,'Hs.212838';
					last DB;
				}
				if ($dbname eq 'PharmGKB') {
					is $dblink->primary_id,'PA24357';
					last DB;
				}
				if ($dbname eq 'MGC') {
					ok $dblink->url; #We may even validate the urls?
					is $dblink->primary_id,'BC040071';
					last DB;
				}
			}
		}
	}
	is $loop_count, shift @lc, "correct number of loops for T10";
	$loop_count = 0;
	
	#T11: SOME EXTERNAL DATABASE IDS TESTS
	foreach my $key (@keys) {
		$loop_count++;
		next if grep(/\b$key\b/, @revkeys);
		my @all=$ann->get_Annotations($key);
		#Checking xref to some databases- OMIM, Wormbase and HGNC, others later
		my $loop_count_internal = 0;
		foreach my $pid (@all) {
			$loop_count_internal++;
			DBID: {
				if (($acc==8)&&($key eq 'MIM')) {is $pid->value,'108985'; last DBID;}
				if (($acc==9)&&($key eq 'HGNC')) {is $pid->value,'7645'; last DBID;}
				if (($acc==11298)&&($key eq 'MGI')) {is $pid->value,'1328365'; last DBID;}
				if (($acc==171593)&&($key eq 'AceView/WormGenes')) {is $pid->value,'1A502'; last DBID;} 
				if (($acc==171594)&&($key eq 'WormBase')) {is $pid->value,'Y48G1C.4'; last DBID;} 
			}
		}
		is $loop_count_internal, shift @lc, "correct number of loops for T11a";
	}
	is $loop_count, shift @lc, "correct number of loops for T11";
	$loop_count = 0;
	
	#T12: REFERENCE RECORD TEST
	if ($acc==1) {
		foreach my $ref (@refs) {
			$loop_count++;
			my $pmed=$ref->medline;
			is grep(/\b$pmed\b/,@pubmed),1;
		}
		is $loop_count, shift @lc, "correct number of loops for T12";
		$loop_count = 0;
	}
	
	#T13/14: STS Markers and Gene Ontology
	my @syn=('MGI:707739','MPC786');
	my @evid=qw(IEA TAS ISS);
	my (%pmeds,%go);
	$go{11305}=['5524', '16887', '5215', '8203', '6810', '16021' ,'5765'];
	$go{11298}=['8080', '8415', '4060', '16740'];
	$pmeds{11305}=['12466851']; 
	my @types=qw(Function Component Process);
	if (($acc==11305)||($acc==11298)) { #Let's check just this two...
		foreach my $ot ($ann->get_Annotations('OntologyTerm')) {
			$loop_count++;
			if (($ot->term->authority)&&($ot->term->authority eq 'STS marker')) {
				if ($acc==11305) {
					is $ot->name,'AI413825';
					is $ot->term->namespace,'UniSTS';
					is $ot->identifier,158928;
				}
				else {
					is $ot->name,'D11Mit102';
					is $ot->term->namespace,'UniSTS';
					is $ot->identifier,126289;
					foreach my $syn ($ot->get_synonyms) {
						is grep(/\b$syn\b/,@syn),1;
					}
				}
				next;
			}
			my $evid=$ot->comment;
			$evid=~s/evidence: //i;
			my $type=$ot->ontology->name;
			my @ref=$ot->term->get_references;
			my $id=$ot->identifier;
			my $thispmed=$ref[0]->medline if (@ref);
			is grep(/\b$type\b/,@types),1;
			is grep(/\b$id\b/,@{$go{$acc}}),1;
			is grep(/\b$thispmed\b/,@{$pmeds{$acc}}),1 if ($thispmed);
			ok $ot->name;
		}
		is $loop_count, shift @lc, "correct number of loops for T13/14";
		$loop_count = 0;
	}
	
	#T15/16/17: GENOMIC LOCATION TESTS/SEQUENCE TYPES TESTS/CONSERVED DOMAINS TESTS
	my @gffs=('SEQ	entrezgene	gene location	63548355	63556668	.	+	.',
				 'SEQ	entrezgene	genestructure	63548355	63556668	.	+	.',
				 'SEQ	entrezgene	gene location	31124733	31133046	.	+	.',
				 'SEQ	entrezgene	genestructure	31124733	31133046	.	+	.',
				 'SEQ	entrezgene	gene location	8163589	8172398	.	+	.',
				 'SEQ	entrezgene	genestructure	8163589	8172398	.	+	.');
	my @contigs=$struct->get_members;
	my @auth=('mrna','genomic','product','mrna sequence','protein','peptide');#Known types....
	foreach my $contig (@contigs) {
		$loop_count++;
		my $stype=$contig->authority;
		is grep(/^$stype$/i,@auth),1;
		if ($acc==1) {#Do just 1?
			if (($contig->authority eq 'genomic')||($contig->authority eq 'Genomic')) {
				foreach my $sf ($contig->get_SeqFeatures) {
					$sf->source_tag('entrezgene');
					my $gff=$sf->gff_string;
					$gff=~s/[\t\s]+$//g;
					foreach my $gffstr (@gffs) {
						if ($gffstr eq $gff) {
							ok(1);
							last;
						}
					}
				}
			}
			if ($contig->authority eq 'Product') {
				is $contig->id,'NP_570602';
				is $contig->accession_number,21071030;
				foreach my $sf ($contig->get_SeqFeatures) {
				foreach my $dblink ($sf->annotation->get_Annotations('dblink')) {
						my $key=$dblink->{_anchor}?$dblink->{_anchor}:$dblink->optional_id;
						my $db=$dblink->database;
						next unless (($db =~/cdd/i)||($sf->primary_tag=~ /conserved/i));
						my $desc;
						if ($key =~ /:/) {
							($key,$desc)=split(/:/,$key);
						}
						$desc=~s/^\s+//;#THIS SHOULD GO IN entrezgene.pm!!!
						is $desc,'IGc2; Immunoglobulin C-2 Type';
						is $key,'smart00408';
						is $sf->score,103;
						is $db,'CDD';
						is $sf->start,223;
						is $sf->end,282;
				}
				}
			}
		}
	}
	cmp_ok( $loop_count,'>=', shift @lc, "correct number of loops for T15");
	$loop_count = 0;
}
is $num_of_seqs, 39, 'looped through correct number of sequences';


#, -locuslink=>'convert');
#See if we can convert to locuslink
#T18: BACKCOMPATIBILITY TESTS
my @llsp =('OFFICIAL_GENE_NAME','CHR','MAP','OFFICIAL_SYMBOL');
ok my $eio_b=Bio::SeqIO->new(-file=>test_input_file('entrezgene.dat'),-format=>'entrezgene', -debug=>'on',-service_record=>'yes',-locuslink=>'convert');
my $loop_count = 0;
while (my $seq=$eio_b->next_seq) {
    $loop_count++;
    ok $seq;
    my $acc=$seq->accession_number;
    is grep(/\b$acc\b/,@ids),1;
    my $ann=$seq->annotation;
    last if ($acc==4);#3 is enough? and 4 does not have gene name, so....
    foreach my $key (@llsp) {
        my @vals=$ann->get_Annotations($key);
        ok @vals;
    }
}
is $loop_count, 4, "correct number of loops for T18";

# Test for Bug #3453
ok my $eio_c = Bio::SeqIO->new(-format => 'entrezgene',
                               -file   => test_input_file('entrezgene_bug3453.dat') );
my $entry = 0;
while ( my ( $gene, $genestructure, $uncaptured ) = $eio_c->next_seq ) {
    $entry++;
    if ($entry == 1) {
        is $gene->accession_number, 3581;
        is scalar @{ $uncaptured }, 55;
    }
    elsif ($entry == 2) {
        is $gene->accession_number, 56111;
        is scalar @{ $uncaptured }, 32;
    }
}
