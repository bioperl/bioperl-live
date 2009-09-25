# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {     
    use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 450,
			   -requires_module => 'IO::Scalar');
	
	use_ok('Bio::Tools::CodonTable');
	use_ok('Bio::SeqIO::table');
}

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

ok my $seqin = Bio::SeqIO->new(-file => test_input_file("test.tsv"),
			    -format  => 'table',
                            -species => "Homo sapiens",
                            -delim   => "\t",
                            -header  => 1,
                            -display_id => 1,
                            -accession_number => 2,
                            -seq => 7,
                            -annotation => 1,
                            -trim => 1);
run_tests([@names],[@accs],[@num_anns],[@psg],[@rs]);

$seqin->close();

ok $seqin = Bio::SeqIO->new(-file => test_input_file("test.tsv"),
                         -format  => 'table',
                         -species => "Homo sapiens",
                         -delim   => "\t",
                         -header  => 1,
                         -display_id => 1,
                         -accession_number => 2,
                         -seq => 7,
                         -colnames => "[Family,Subfamily,Pseudogene?,Protein,Novelty]",
                         -trim => 1);
run_tests([@names],[@accs],[4,4,4,4,4,5,5,5,5,5],[@psg],[@rs]);

$seqin->close();

ok $seqin = Bio::SeqIO->new(-file => test_input_file("test.tsv"),
                         -format  => 'table',
                         -species => "Homo sapiens",
                         -delim   => "\t",
                         -header  => 1,
                         -display_id => 1,
                         -accession_number => 2,
                         -seq => 7,
                         -annotation => "[4,5,6,8,10]",
                         -trim => 1);
run_tests([@names],[@accs],[4,4,4,4,4,5,5,5,5,5],[@psg],[@rs]);

$seqin->close();

# need Spreadsheet::ParseExcel installed for testing Excel format
SKIP: {
	test_skip(-tests => 112, -requires_module => 'Spreadsheet::ParseExcel');

	ok $seqin = Bio::SeqIO->new(-file => test_input_file("test.xls"),
							 -format  => 'excel',
							 -species => "Homo sapiens",
							 -header  => 1,
							 -display_id => 1,
							 -accession_number => 2,
							 -seq => 7,
							 -annotation => 1,
							 -trim => 1);
	run_tests([@names],[@accs],[@num_anns],[@psg],[@rs]);
	
	$seqin->close();
}

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
        is ($seq->display_id, shift(@names));
        is ($seq->accession_number, shift(@accs));
        ok ($seq->species);
        is ($seq->species->binomial, "Homo sapiens");
        my @anns = $seq->annotation->get_Annotations();
        is (scalar(@anns), shift(@num_anns));
        @anns = grep { $_->value eq "Y"; 
                     } $seq->annotation->get_Annotations("Pseudogene?");
        is (scalar(@anns), shift(@psg));
        
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
            ok (defined $protann);
            is ($protein, $protann->value);
        }
        
        @anns = grep { $_->value eq "Known - Refseq"; 
                     } $seq->annotation->get_Annotations("Novelty");
        is (scalar(@anns), shift(@rs));
        @anns = $seq->annotation->get_Annotations("Subfamily");
        is (scalar(@anns), ($n <= 5) ? 0 : 1);
        @anns = $seq->annotation->get_Annotations("Family");
        is (scalar(@anns), 1);
        is (substr($anns[0]->value,0,4), ($n <= 4) ? "A6" : "CAMK");    
    }
    
    is ($n, 10);
}
