# -*-Perl-*-
## Bioperl Test Harness Script for Modules
## $Id$

use strict;

BEGIN {
	# to handle systems with no installed Test module
	# we include the t dir (where a copy of Test.pm is located)
	# as a fallback
    eval { require Test; };
    if( $@ ) {
        use lib 't';
    }
    use Test;
    plan tests => 38;
}


my $DEBUG = $ENV{'BIOPERLDEBUG'};


use Bio::Seq::Meta;
ok(1);

ok my $seq = Bio::Seq::Meta->new
    ( -seq => "",
      -meta => "",
      -alphabet => 'dna',
      -id => 'myid'
    );

# create a sequence object
ok $seq = Bio::Seq::Meta->new( -seq => "AT-CGATCGA",
                               -id => 'test',
                               -verbose => 2
                             );

ok $seq->meta, "          ";

# create some random meta values, but gap in the wrong place
my $metastring = "a-abb  bb ";
eval {
    $seq->meta($metastring);
};
ok 1 if $@ =~ 'column [2]';
$seq->verbose(1);

# create some random meta values, but not for the last residue
$metastring = "aa-bb  bb";
ok $seq->meta($metastring), $metastring. " ";

# truncate the sequence by assignment
$seq->seq('AT-CGA');
$seq->alphabet('dna');
ok $seq->meta, 'aa-bb ';
ok $seq->meta_text, 'aa-bb ';

# truncate the sequence with trunc()
ok $seq->strand(-1), -1;
ok $seq = $seq->trunc(1,5);
ok $seq->seq, 'AT-CG';
ok $seq->meta, 'aa-bb';
ok $seq->strand, -1;

# revcom
ok $seq = $seq->revcom;
ok $seq->seq, 'CG-AT';
ok $seq->meta, 'bb-aa';
ok $seq->strand, 1;

# submeta
ok $seq->subseq(2,4), 'G-A';
ok $seq->submeta(2,4), 'b-a';
ok $seq->submeta(2,undef, 'c-c'), 'c-c';
ok $seq->submeta(2,4), 'c-c';
ok $seq->meta, 'bc-ca';

ok $seq->meta(''), '     ';
ok $seq->submeta(2,undef, 'c-c'), 'c-c';
ok $seq->meta, ' c-c ';

# add named meta annotations

my $first = '11-22';
ok $seq->named_meta('first', $first), $first;
ok $seq->named_meta('first'), $first;

my $second = '[[-]]';
ok $seq->named_meta('second', $second);

# undefined range arguments
ok $seq->named_submeta('second', 3, 4), '-]';
ok $seq->named_submeta('second', 3), '-]]';
ok $seq->named_submeta('second'), '[[-]]';

my @names =  $seq->meta_names;
ok @names, 3;
ok $names[0], 'DEFAULT';



#
# IO tests
#

sub diff {
    my ($infile, $outfile) = @_;
    my ($in, $out);
    open FH, $infile;
    $in .= $_ while (<FH>);
    close FH;

    open FH, $outfile;
    $out .= $_ while (<FH>);
    close FH;
    print "|$in||$out|" if $DEBUG;
    ok $in, $out;

}

use Bio::SeqIO;
use Bio::AlignIO;
use Bio::Root::IO;

# SeqIO
my $str = Bio::SeqIO->new
    ( '-file'=> Bio::Root::IO->catfile("t","data","test.metafasta"),
      '-format' => 'metafasta');
ok  $seq = $str->next_seq;

my $strout = Bio::SeqIO->new
    ('-file'=> ">". Bio::Root::IO->catfile("t","data","test.metafasta.out"),
     '-format' => 'metafasta');
ok $strout->write_seq($seq);

diff (Bio::Root::IO->catfile("t","data","test.metafasta"),
      Bio::Root::IO->catfile("t","data","test.metafasta.out")
     );


# AlignIO

$str = Bio::AlignIO->new
    ( '-file'=> Bio::Root::IO->catfile("t","data","testaln.metafasta"),
      '-format' => 'metafasta');
ok my $aln = $str->next_aln;

$strout = Bio::AlignIO->new
    ('-file'=> ">". Bio::Root::IO->catfile("t","data","testaln.metafasta.out"),
     '-format' => 'metafasta');
ok $strout->write_aln($aln);

diff (Bio::Root::IO->catfile("t","data","testaln.metafasta"),
      Bio::Root::IO->catfile("t","data","testaln.metafasta.out")
     );


END {
    unlink(Bio::Root::IO->catfile("t","data","test.metafasta.out"));
    unlink(Bio::Root::IO->catfile("t","data","testaln.metafasta.out"));
}
