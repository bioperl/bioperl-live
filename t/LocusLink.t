# -*-Perl-*- mode (to keep my emacs happy)
# $Id$

use strict;
use vars qw($DEBUG $NUMTESTS $HAVEGRAPHDIRECTED);
$DEBUG = $ENV{'BIOPERLDEBUG'} || 0;

BEGIN {
	eval { require Test::More; };
	if( $@ ) {
		use lib 't/lib';
	}
	use Test::More;
	eval {
		require Graph::Directed; 
	};
	if ($@) {
		plan skip_all => "Graph::Directed not installed, skipping tests";
	} else {
		plan tests => ($NUMTESTS = 26);
	}
	use_ok('Bio::SeqIO');
	use_ok('Bio::Root::IO');
	use_ok('Bio::SeqFeature::Generic');
	use_ok('Bio::SeqFeature::AnnotationAdaptor');
}

END {
	unlink("locuslink-test.out.embl") if -e "locuslink-test.out.embl";
}

my $seqin = Bio::SeqIO->new(-file => Bio::Root::IO->catfile("t","data",
							 "LL-sample.seq"),
			    -format => 'locuslink');
ok $seqin;
my $seqout = Bio::SeqIO->new(-file => ">locuslink-test.out.embl",
			     -format => 'embl');

# process and write to output
my @seqs = ();

while(my $seq = $seqin->next_seq()) {
    push(@seqs, $seq);
    
    # create an artificial feature to stick the annotation on
    my $fea = Bio::SeqFeature::Generic->new(-start => 1, -end => 9999,
					    -strand => 1,
					    -primary => 'annotation');
    my $ac = Bio::SeqFeature::AnnotationAdaptor->new(-feature => $fea);
    foreach my $k ($seq->annotation->get_all_annotation_keys()) {
	foreach my $ann ($seq->annotation->get_Annotations($k)) {
	    next unless $ann->isa("Bio::Annotation::SimpleValue");
	    $ac->add_Annotation($ann);
	}
    }
    $seq->add_SeqFeature($fea);
    $seqout->write_seq($seq);
}

is (scalar(@seqs), 2);

is ($seqs[0]->desc,
    "amiloride binding protein 1 (amine oxidase (copper-containing))");
is ($seqs[0]->accession, "26");
is ($seqs[0]->display_id, "ABP1");
is ($seqs[0]->species->binomial, "Homo sapiens");


my @dblinks = $seqs[0]->annotation->get_Annotations('dblink');
my %counts = map { ($_->database(),0) } @dblinks;
foreach (@dblinks) { $counts{$_->database()}++; }

is ($counts{GenBank}, 11);
is ($counts{RefSeq}, 4);
is ($counts{UniGene}, 1);
is ($counts{Pfam}, 1);
is ($counts{STS}, 2);
is ($counts{MIM}, 1);
is ($counts{PUBMED}, 6);
is (scalar(@dblinks), 27);

is ($seqs[1]->desc, "v-abl Abelson murine leukemia viral oncogene homolog 2 (arg, Abelson-related gene)");
is ($seqs[1]->display_id, "ABL2");

my $ac = $seqs[1]->annotation;
my @keys = $ac->get_all_annotation_keys();
is (scalar(@keys), 19);

my ($cmt) = $ac->get_Annotations('comment');
is (length($cmt->text), 403);

my @isoforms = qw(a b);
foreach ($ac->get_Annotations('PRODUCT')) {
    is ($_->value,
	"v-abl Abelson murine leukemia viral oncogene homolog 2 isoform ".
	shift(@isoforms));
}

my @goann = ();
foreach my $k (@keys) {
    foreach my $ann ($ac->get_Annotations($k)) {
	next unless $ann->isa("Bio::Ontology::TermI");
	push(@goann, $ann);
    }
}
is (scalar(@goann), 4);
@goann = sort { $a->as_text() cmp $b->as_text() } @goann;
is ($goann[2]->as_text, "cellular component|cytoplasm|");


