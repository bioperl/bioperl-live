#!/usr/bin/perl

use strict;
use Bio::Root::IO;
use Data::Dumper;
use vars qw($DEBUG $NUMTESTS $ASNOK);
$DEBUG = $ENV{'BIOPERLDEBUG'} || 0;

BEGIN {
	eval { require Test; };
	if( $@ ) {
		use lib 't';
	}
	use Test;
	eval {
		require Bio::ASN1::EntrezGene;
		$ASNOK=1;
	};
	if ($@) {
		$ASNOK = 0;
		warn "Bio::ASN1::EntrezGene not installed, skipping tests\n";
	}
    plan tests => ($NUMTESTS = 1003);
}

END {
	foreach ( $Test::ntest..$NUMTESTS) {
		skip('Cannot complete entrezgene tests',1);
	}
}

exit(0) unless $ASNOK;

use Bio::SeqIO; 

my @species=('Homo sapiens','Mus musculus', 'Caenorhabditis elegans');
my @pubmed=qw(15461460
15221005
14702039
12477932
8889549
3610142
3458201
2591067);

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
ok(1);

my $fs='!';
my @revkeys=('Entrez Gene Status','RefSeq status','Official Full Name','chromosome','cyto','Reference','dblink',
'ALIAS_SYMBOL','OntologyTerm','Index terms','Official Symbol','cM','Property');


my $eio=new Bio::SeqIO(-file=>Bio::Root::IO->catfile("t","data",
							 "entrezgene.dat"),-format=>'entrezgene', -debug=>'on',-service_record=>'yes');
ok $eio;
my ($seq,$struct,$uncapt);
while (1) {
my $seq;
($seq,$struct,$uncapt)=$eio->next_seq;
last unless ($seq);

#T0: GENERAL TESTS
ok $seq;
ok ref($struct),'Bio::Cluster::SequenceFamily';
my $acc=$seq->accession_number;

#T1: ORGANISM
my $org=$seq->species->binomial;
ok grep(/\b$org\b/,@species),1,$org;

#T2: SUMMARY test
ok $seq->desc if ($acc eq '1')||($acc eq '2')||($acc eq '11304');
ok !defined $seq->desc if ($acc eq '171592')||($acc eq '11306');

#Are we supposed to have this in our test?
ok grep(/\b$acc\b/,@ids),1;

my $ann=$seq->annotation();
my $tcount;

#T3: ENTREZGENE STATUS TESTS
my @egstatus=$ann->get_Annotations('Entrez Gene Status');
foreach my $status (@egstatus) {
 STATUS: {
		if ($acc==1) {ok $status->value,'live'; last STATUS;}
		if ($acc==2) {ok $status->value,'live'; last STATUS;}
		if ($acc==4) {ok $status->value,'discontinued'; last STATUS;}
		if ($acc==6) {ok $status->value,'discontinued'; last STATUS;}
		if ($acc==11288) {ok $status->value,'secondary'; last STATUS;}
		if ($acc==11293) {ok $status->value,'secondary'; last STATUS;} 
		if ($acc==171594) {ok $status->value,'live'; last STATUS;} 
	}
}

#T4: REFSEQ STATUS TESTS
my @refstatus=$ann->get_Annotations('RefSeq status');
foreach my $status (@refstatus) {
 STATUS: {
		if ($acc==1) {ok $status->value,'REVIEWED'; last STATUS;}
		if ($acc==2) {ok $status->value,'REVIEWED'; last STATUS;}
		if ($acc==3) {ok $status->value,'PROVISIONAL'; last STATUS;}
		if ($acc==4) {ok $status->value,'WITHDRAWN'; last STATUS;}
		if ($acc==9) {ok $status->value,'VALIDATED'; last STATUS;}
		if ($acc==11300) {ok $status->value,''; last STATUS;}
		if ($acc==11306) {ok $status->value,'MODEL'; last STATUS;}
		if ($acc==11293) {ok $status->value,'secondary'; last STATUS;} 
		if ($acc==171594) {ok $status->value,'Reviewed'; last STATUS;} 
	}
}

#T5: GENE NAME TESTS
my @ofname=$ann->get_Annotations('Official Full Name');
foreach my $name (@ofname) {
 STATUS: {
		if ($acc==10) {ok $name->value,'N-acetyltransferase 2 (arylamine N-acetyltransferase)'; last STATUS;}
		if ($acc==13) {ok $name->value,'arylacetamide deacetylase (esterase)'; last STATUS;}
		if ($acc==14) {ok $name->value,'angio-associated, migratory cell protein'; last STATUS;}
		if ($acc==11287) {ok $name->value,'pregnancy zone protein'; last STATUS;}
            if ($acc==11298) {ok $name->value,'arylalkylamine N-acetyltransferase'; last STATUS;}
		if ($acc==11304) {ok $name->value,'ATP-binding cassette, sub-family A (ABC1), member 4'; last STATUS;}
		if ($acc==11306) {ok $name->value,'ATP-binding cassette, sub-family B (MDR/TAP), member 7'; last STATUS;} 
	}
}

#T6: CHROMOSOME TESTS
my @chr=$ann->get_Annotations('chromosome');
foreach my $chr (@chr) {
 STATUS: {
		if ($acc==5) {ok $chr->value,1; last STATUS;}
		if ($acc==6) {ok $chr->value,1; last STATUS;}
		if ($acc==7) {ok $chr->value,17; last STATUS;}
		if ($acc==11306) {ok $chr->value,'X'; last STATUS;}
		if ($acc==11304) {ok $chr->value,3; last STATUS;}
		if ($acc==171590) {ok $chr->value,'I'; last STATUS;}
		if ($acc==171592) {ok $chr->value,'I'; last STATUS;} 
	}
}

#T7: GENE SYMBOL ALIAS TESTS
my @sym=$ann->get_Annotations('ALIAS_SYMBOL');
foreach my $sym (@sym) {
    next if (($sym eq '')||!defined($sym));
    ok grep(/\b$sym\b/,@{$asym{$acc}}),1;
}

#T8: CYTO LOCATION TESTS
my @map=$ann->get_Annotations('cyto');
foreach my $map (@map) {

  STATUS: {
		 if ($acc==10) {ok $map->value,'8p22'; last STATUS;}
		 if ($acc==11) {ok $map->value,'8p22'; last STATUS;}
		 if ($acc==13) {ok $map->value,'3q21.3-q25.2'; last STATUS;}
		 if ($acc==11306) {ok $map->value,'X C-D'; last STATUS;}
		 if ($acc==11305) {ok $map->value,'2 A2-B'; last STATUS;}
		 if ($acc==11304) {ok $map->value,'3 G1'; last STATUS;}
		 if ($acc==11303) {ok $map->value,'4 A5-B3'; last STATUS;} 
	 }
 }

#T9: REFERENCE NUMBER TEST
my @refs=$ann->get_Annotations('Reference');
my $refs=$#refs+1||0;
ok $pmed{$acc},$refs;


my @dblinks=$ann->get_Annotations('dblink');
my @keys=$ann->get_all_annotation_keys;

#T10: GENERIF AND OTHER DBLINK TESTS
my @url=qw(HGMD Ensembl KEGG Homologene);#Only validate the URL
foreach my $dblink (@dblinks) {
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
            ok $dblink->primary_id,2;
            last DB;
        }
        if ($dbname eq 'Model maker') {
            ok $dblink->url; #We may even validate the urls?
            ok $dblink->primary_id,2;
            last DB;
        }
        if ($dbname eq 'AceView') {
            ok $dblink->url; #We may even validate the urls?
            ok $dblink->primary_id,2;
            last DB; 
        }
        if (grep(/$dbname/,@url)) {
            ok $dblink->url; #We may even validate the urls?
            last DB;
        }
        if ($dbname eq 'GDB') {
            ok $dblink->primary_id,'GDB:119639'; #We may even validate the urls?
            last DB;
        }
        if ($dbname eq 'UniGene') {
            ok $dblink->url; #We may even validate the urls?
            ok $dblink->primary_id,'Hs.212838';
            last DB;
        }
        if ($dbname eq 'PharmGKB') {
            ok $dblink->primary_id,'PA24357';
            last DB;
        }
        if ($dbname eq 'MGC') {
            ok $dblink->url; #We may even validate the urls?
            ok $dblink->primary_id,'BC040071';
            last DB;
        }
    }
}
}

#T11: SOME EXTERNAL DATABASE IDS TESTS
foreach my $key (@keys) {
	next if grep(/\b$key\b/, @revkeys);
	my @all=$ann->get_Annotations($key);
	#Checking xref to some databases- OMIM, Wormbase and HGNC, others later
	foreach my $pid (@all) {
	 DBID: {
			if (($acc==8)&&($key eq 'MIM')) {ok $pid->value,'108985'; last DBID;}
			if (($acc==9)&&($key eq 'HGNC')) {ok $pid->value,'7645'; last DBID;}
			if (($acc==11298)&&($key eq 'MGI')) {ok $pid->value,'1328365'; last DBID;}
			if (($acc==171593)&&($key eq 'AceView/WormGenes')) {ok $pid->value,'1A502'; last DBID;} 
			if (($acc==171594)&&($key eq 'WormBase')) {ok $pid->value,'Y48G1C.4'; last DBID;} 
		}
	}
}

#T12: REFERENCE RECORD TEST
if ($acc==1) {
    foreach my $ref (@refs) {
        my $pmed=$ref->medline;
        ok grep(/\b$pmed\b/,@pubmed),1;
    }
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
		if (($ot->term->authority)&&($ot->term->authority eq 'STS marker')) {
			if ($acc==11305) {
				ok $ot->name,'AI413825';
				ok $ot->term->namespace,'UniSTS';
				ok $ot->identifier,158928;
			}
			else {
				ok $ot->name,'D11Mit102';
				ok $ot->term->namespace,'UniSTS';
				ok $ot->identifier,126289;
				foreach my $syn ($ot->get_synonyms) {
					ok grep(/\b$syn\b/,@syn),1;
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
		ok grep(/\b$type\b/,@types),1;
		ok grep(/\b$id\b/,@{$go{$acc}}),1;
		ok grep(/\b$thispmed\b/,@{$pmeds{$acc}}),1 if ($thispmed);
		ok $ot->name;
	}
}

#T15/16/17: GENOMIC LOCATION TESTS/SEQUENCE TYPES TESTS/CONSERVED DOMAINS TESTS
my @gffs=('SEQ	entrezgene	gene location	63548355	63556668	.	+	.',
			 'SEQ	entrezgene	genestructure	63548355	63556668	.	+	.',
			 'SEQ	entrezgene	gene location	31124733	31133046	.	+	.',
			 'SEQ	entrezgene	genestructure	31124733	31133046	.	+	.',
			 'SEQ	entrezgene	gene location	8163589	8172398	.	+	.',
			 'SEQ	entrezgene	genestructure	8163589	8172398	.	+	.');
my @contigs=$struct->get_members;
my @auth=('mrna','genomic','product','mrna sequence','protein');#Known types....
foreach my $contig (@contigs) {
	my $stype=$contig->authority;
	ok grep(/^$stype$/i,@auth),1;
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
			ok $contig->id,'NP_570602';
			ok $contig->accession_number,21071030;
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
					ok $desc,'IGc2; Immunoglobulin C-2 Type';
					ok $key,'smart00408';
					ok $sf->score,103;
					ok $db,'CDD';
					ok $sf->start,223;
					ok $sf->end,282;
            }
			}
		}
	}
}
}


#, -locuslink=>'convert');
#See if we can convert to locuslink
#T18: BACKCOMPATIBILITY TESTS
my @llsp =('OFFICIAL_GENE_NAME','CHR','MAP','OFFICIAL_SYMBOL');
my $eio_b=new Bio::SeqIO(-file=>Bio::Root::IO->catfile("t","data",
							 "entrezgene.dat"),-format=>'entrezgene', -debug=>'on',-service_record=>'yes',-locuslink=>'convert');

while (my $seq=$eio_b->next_seq) {
    ok $seq;
    my $acc=$seq->accession_number;
    ok grep(/\b$acc\b/,@ids),1;
    my $ann=$seq->annotation;
    last if ($acc==4);#3 is enough? and 4 does not have gene name, so....
    foreach my $key (@llsp) {
        my @vals=$ann->get_Annotations($key);
        ok @vals;
    }
}
