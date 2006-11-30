# -*-Perl-*- mode (to keep my emacs happy)
# $GNF: projects/gi/symgene/src/perl/seqproc/t/table.t,v 1.3 2006/01/19 04:20:36 hlapp Exp $

use strict;
use vars qw($DEBUG $ERROR);
use constant NUMTESTS => 449;
use constant NONEXCELTESTS => 337;

BEGIN {     
    # to handle systems with no installed Test module
    # we include the t dir (where a copy of Test.pm is located)
    # as a fallback
    eval { require Test; };
    if( $@ ) { 
        use lib 't';
    }
    use Test;
	
	# we seem to need IO::Scalar for this
	eval {
		require IO::Scalar;
	};
	if ($@) {
		$ERROR = 1;
	}
	
    plan tests => NUMTESTS;
}

if ($ERROR) {
	foreach (1..NUMTESTS) { 
        skip ("IO::Scalar not installed, skipping all tests",1,1); 
    }
    exit(0);
}

use Bio::Tools::CodonTable;
use Bio::SeqIO;
use Bio::Root::IO;

ok(1); # if this fails already we're in trouble

my @names = qw(A6
               A6r
               A6ps1
               A6ps2
               CaMK2d
               CaMKK2
               AMPKa1
               AMPKa2
               MARK3
               MARK2);
my @accs = qw(SK001
              SK512
              SK752
              SK766
              SK703
              SK482
              SK032
              SK033
              SK096
              SK120);
my @num_anns = (5, 5, 5, 5, 6, 7, 7, 7, 7, 7);
my @psg = (0, 0, 1, 1, 0, 0, 0, 0, 0, 0);
my @rs = (0, 0, 0, 0, 1, 1, 1, 1, 1, 1);

my $io = "Bio::Root::IO";
my $seqin = Bio::SeqIO->new(-file => $io->catfile("t","data","kinases.tsv"),
			    -format  => 'table',
                            -species => "Homo sapiens",
                            -delim   => "\t",
                            -header  => 1,
                            -display_id => 1,
                            -accession_number => 2,
                            -seq => 7,
                            -annotation => 1,
                            -trim => 1);
ok $seqin;
run_tests([@names],[@accs],[@num_anns],[@psg],[@rs]);

$seqin->close();

$seqin = Bio::SeqIO->new(-file => $io->catfile("t","data","kinases.tsv"),
                         -format  => 'table',
                         -species => "Homo sapiens",
                         -delim   => "\t",
                         -header  => 1,
                         -display_id => 1,
                         -accession_number => 2,
                         -seq => 7,
                         -colnames => "[Family,Subfamily,Pseudogene?,Protein,Novelty]",
                         -trim => 1);
ok $seqin;
run_tests([@names],[@accs],[4,4,4,4,4,5,5,5,5,5],[@psg],[@rs]);

$seqin->close();

$seqin = Bio::SeqIO->new(-file => $io->catfile("t","data","kinases.tsv"),
                         -format  => 'table',
                         -species => "Homo sapiens",
                         -delim   => "\t",
                         -header  => 1,
                         -display_id => 1,
                         -accession_number => 2,
                         -seq => 7,
                         -annotation => "[4,5,6,8,10]",
                         -trim => 1);
ok $seqin;
run_tests([@names],[@accs],[4,4,4,4,4,5,5,5,5,5],[@psg],[@rs]);

$seqin->close();

# need Spreadsheet::ParseExcel installed for testing Excel format
eval {
    require Spreadsheet::ParseExcel;
};
if ($@) {
    foreach ((NONEXCELTESTS+1)..NUMTESTS) { 
        skip ("Skip Excel format test because Spreadsheet::ParseExcel not installed",1,1); 
    }
    exit(0);
}

$seqin = Bio::SeqIO->new(-file => $io->catfile("t","data","kinases.xls"),
                         -format  => 'excel',
                         -species => "Homo sapiens",
                         -header  => 1,
                         -display_id => 1,
                         -accession_number => 2,
                         -seq => 7,
                         -annotation => 1,
                         -trim => 1);
ok $seqin;
run_tests([@names],[@accs],[@num_anns],[@psg],[@rs]);

$seqin->close();

sub run_tests {
    my ($names_,$accs_,$num_anns_,$psg_,$rs_) = @_;

    my @names = @$names_;
    my @accs = @$accs_;
    my @num_anns = @$num_anns_;
    my @psg = @$psg_;
    my @rs = @$rs_;

    my $n = 0;
    my $translator = Bio::Tools::CodonTable->new(-id => 1);
    while (my $seq = $seqin->next_seq()) {
        $n++;
        ok ($seq->display_id, shift(@names));
        ok ($seq->accession_number, shift(@accs));
        ok ($seq->species);
        ok ($seq->species->binomial, "Homo sapiens");
        my @anns = $seq->annotation->get_Annotations();
        ok (scalar(@anns), shift(@num_anns));
        @anns = grep { $_->value eq "Y"; 
                     } $seq->annotation->get_Annotations("Pseudogene?");
        ok (scalar(@anns), shift(@psg));
        
        # check sequences and that they translate to what we expect
        if (($n >= 5) && ($seq->display_id ne "MARK3")) {
            my $dna = $seq->seq;
            my $protein = "";
            my $frame = 0;
            while ($frame <= 2) {
                my $inframe = substr($dna,$frame);
                # translate to protein
                my $protseq = $translator->translate($inframe);
                # chop off everything after the stop and before the first Met
                while ($protseq =~ /(M[^\*]+)/g) {
                    $protein = $1 if length($1) > length($protein);
                }
                $frame++;
            }
            # retrieve expected result from annotation and compare
            my ($protann) = $seq->annotation->get_Annotations("Protein");
            ok ($protann);
            ok ($protein, $protann->value);
        }
        
        @anns = grep { $_->value eq "Known - Refseq"; 
                     } $seq->annotation->get_Annotations("Novelty");
        ok (scalar(@anns), shift(@rs));
        @anns = $seq->annotation->get_Annotations("Subfamily");
        ok (scalar(@anns), ($n <= 5) ? 0 : 1);
        @anns = $seq->annotation->get_Annotations("Family");
        ok (scalar(@anns), 1);
        ok (substr($anns[0]->value,0,4), ($n <= 4) ? "A6" : "CAMK");    
    }
    
    ok ($n, 10);
}

