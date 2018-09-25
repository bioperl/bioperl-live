# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
    use Bio::Root::Test;

    test_begin(-tests => 27);

    use_ok('Bio::Seq');
    use_ok('Bio::SeqIO');
}

ok my $str = Bio::SeqIO->new(-file   => test_input_file('U58726.gb'),
                             -format => 'GenBank');
my $seq;
ok ( $seq = $str->next_seq() );

# Here is a cute way to verify the sequence by seeing if the
# the translation matches what is annotated in the file -js
foreach my $ft ( grep { $_->primary_tag eq 'CDS'}
                 $seq->top_SeqFeatures ) {
    if( $ft->has_tag('translation') ) {
        my ($translation) = $ft->get_tag_values('translation');
        my $t = $ft->spliced_seq(-nosort => 1);
        my $pepseq = $t->translate()->seq();
        chop($pepseq); # chop is to remove stop codon
        is($translation, $pepseq);
    }
}

my $stream = Bio::SeqIO->new(-file   => test_input_file('M12730.gb'),
                             -format => 'genbank');
# Jump down to M12730 which lists CDS join(1959..2355,1..92)
while ($seq->accession ne "M12730") {
    $seq = $stream->next_seq;
}
ok(my @features = $seq->get_SeqFeatures(), "get_SeqFeatures()");
my $feat;
foreach my $feat2 ( @features ) {
    next unless ($feat2->primary_tag eq "CDS");
    my @db_xrefs = $feat2->get_tag_values("db_xref");
    if (grep { $_ eq "GI:150830" } @db_xrefs) {
       $feat = $feat2;
       last;
    }
}
my ($protein_seq) = $feat->get_tag_values("translation");
like($protein_seq, qr(^MKERYGTVYKGSQRLIDE.*ANEKQENALYLIIILSRTSIT$),
                   "protein sequence");
my ($nucleotide_seq) = $feat->spliced_seq(-nosort => 1)->seq;
like($nucleotide_seq, qr(^ATGAAAGAAAGATATGGA.*TCAAGGACTAGTATAACATAA$),
                     "nucleotide sequence - correct CDS range");
is(length($nucleotide_seq), 489, "nucleotide length");

#  Test for Fix spliced seq #72
my $str2 = Bio::SeqIO->new(-file   => test_input_file('AF032047.gbk'),
                           -format => 'GenBank');
my @feats = $str2-> next_seq -> get_SeqFeatures;
# feat[1] has 2 exons from remote sequence AF032048.1
my $len_nodb;
warnings_like { $len_nodb = length($feats[1]->spliced_seq()->seq); }
              [ {carped => qr/cannot get remote location for/},
                {carped => qr/cannot get remote location for/}
               ],
              "appropriate warning if db not provided for remote sequence";
ok($len_nodb == 374, "correct number of Ns added if remote sequence not provided");

# Test for cut by origin features
my $seq_obj = Bio::Seq->new(-display_id => 'NC_008309',
                            -seq        => 'AAAAACCCCCGGGGGTTTTT');
$seq_obj->is_circular(1);
my $loc_obj = Bio::Factory::FTLocationFactory->from_string('join(16..20,1..2)');
my $cut_feat = Bio::SeqFeature::Generic->new(-primary_tag => 'CDS',
                                             -location    => $loc_obj,
                                             -tag => { locus_tag  => 'HS_1792',
                                                       product    => 'hypothetical protein',
                                                       protein_id => 'YP_718205.1',
                                                      } );
$seq_obj->add_SeqFeature($cut_feat);
is $cut_feat->seq->seq,         'TTTTTAA', 'cut by origin sequence using $feat->seq';
is $cut_feat->spliced_seq->seq, 'TTTTTAA', 'cut by origin sequence using $feat->spliced_seq';
is $cut_feat->start,             16,       'cut by origin start using $feat->start';
is $cut_feat->end,               2,        'cut by origin end using $feat->end';
is $cut_feat->location->start,   16,       'cut by origin start using $feat->location->start';
is $cut_feat->location->end,     2,        'cut by origin end using $feat->location->end';

SKIP: {
    test_skip(-tests => 3,
              -requires_modules    => [qw(Bio::DB::GenBank
                                          LWP::UserAgent )],
              -requires_networking => 1);
    my $db_in;
    eval {
        ok $db_in = Bio::DB::GenBank->new();
        my $seq_obj = $db_in->get_Seq_by_id('AF032048.1');
    };
    if ($@) {
        print "$@\n";
        skip  "Warning: Problem accessing GenBank entry AF032048.1 "
            . "to test spliced_seq on remote DBs", 2;
    }

    my $len_w_db;
    warning_is { $len_w_db = length($feats[1]->spliced_seq(-db => $db_in)->seq) }
               [],
               "no warnings if GenBank db provided for remote sequence";
    ok($len_w_db == 374, "correct length if remote sequence is provided")
}
