# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
	use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 26,
               -requires_module => 'Graph::Directed');
	
	use_ok('Bio::SeqIO::locuslink');
	use_ok('Bio::SeqFeature::Generic');
	use_ok('Bio::SeqFeature::AnnotationAdaptor');
}

my $seqin = Bio::SeqIO->new(-file => test_input_file('test.locuslink'),
			    -format => 'locuslink');
ok $seqin;
isa_ok($seqin, 'Bio::SeqIO');
my $seqout = Bio::SeqIO->new(-file => ">".test_output_file(),
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
