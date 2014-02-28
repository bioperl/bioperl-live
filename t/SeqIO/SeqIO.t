# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {     
    use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 58);
    
    use_ok 'Bio::SeqIO';
}

my $verbose = test_debug();

my @formats = qw(gcg fasta raw pir tab ace );
# The following files or formats are failing: swiss genbank interpro embl

for my $format (@formats) {
    print "======== $format ========\n" if $verbose;
    my $seq;
    my $str = Bio::SeqIO->new( -file   => test_input_file("test.$format"),
                               -format => $format                          );

    is $str->format(), $format;

    ok $seq = $str->next_seq();
    print "Sequence 1 of 2 from $format stream:\n", $seq->seq, "\n\n" if  $verbose;
    unless ($format eq 'raw') {
        is $seq->id, 'roa1_drome',"ID for format $format";
        is $seq->length, 358;
    }
    
    unless ($format eq 'gcg') { # GCG file can contain only one sequence
        ok $seq = $str->next_seq();
        print "Sequence 2 of 2 from $format stream:\n", $seq->seq, $seq->seq, "\n" if $verbose;
    }
    
    my $outfile = test_output_file();
    my $out = Bio::SeqIO->new( -file   => ">$outfile",
                               -format => $format      );
    ok $out->write_seq($seq);
    if ($format eq 'fasta') {
        my $id_type;
        ok($id_type = $out->preferred_id_type('accession.version'),
            'accession.version');
    }
    
    ok -s $outfile;
}



# from testformats.pl
SKIP: {
    test_skip(-tests => 6, -requires_modules => [qw(Algorithm::Diff
                                                    IO::ScalarArray
                                                    IO::String)]);
    use_ok 'Algorithm::Diff';
    eval "use Algorithm::Diff qw(diff LCS);";
    use_ok 'IO::ScalarArray';
    use_ok 'IO::String';
    
    my %files = ( 
             #'test.embl'      => 'embl',
             #'test.ace'       => 'ace',
              'test.fasta'     => 'fasta',
             #'test.game'      => 'game',
              'test.gcg'       => 'gcg',
             #'test.genbank'   => 'genbank',
              'test.raw'       => 'raw',
             #'test_badlf.gcg' => 'gcg'
              );
    
    while( my ($file, $type) = each %files ) {
        my $filename = test_input_file($file);
        print "processing file $filename\n" if $verbose;
        open my $FILE, '<', $filename or die "Could not read file '$filename': $!\n";
        my @datain = <$FILE>;
        close $FILE;

        my $in = IO::String->new( join('', @datain) );
        my $seqin = Bio::SeqIO->new( -fh     => $in,
                                     -format => $type );
        my $out = IO::String->new();
        my $seqout = Bio::SeqIO->new( -fh     => $out,
                                      -format => $type );
        my $seq;
        while( defined($seq = $seqin->next_seq) ) {
            $seqout->write_seq($seq);
        }
        $seqout->close();
        $seqin->close();
        my $strref = $out->string_ref;
        my @dataout = map { $_."\n"} split(/\n/, $$strref );
        my @diffs = &diff( \@datain, \@dataout);
        is @diffs, 0;
        
        if(@diffs && $verbose) {
            foreach my $d ( @diffs ) {
                foreach my $diff ( @$d ) {
                    chomp($diff->[2]);
                    print $diff->[0], $diff->[1], "\n>", $diff->[2], "\n";
                }
            }
            print "in is \n", join('', @datain), "\n";
            print "out is \n", join('',@dataout), "\n";
        }
    }
}

# simple tests specific to Bio::SeqIO interface (applicable to all SeqIO
# modules)

##############################################
# test format() and variant() in Bio::RootIO
##############################################

my $in = Bio::SeqIO->new(
   -file    => test_input_file('bug2901.fa'),
   -format  => "fasta",
);
is $in->format, 'fasta';
is $in->variant, undef;

$in = Bio::SeqIO->new(
   -file    => test_input_file('fastq', 'illumina_faked.fastq'),
   -format  => "fastq",
   -variant => 'illumina',
);
is $in->format, 'fastq';
is $in->variant, 'illumina';


######################################################
# test format detection from different inputs
######################################################

$in = Bio::SeqIO->new( -file => test_input_file('test.fastq') );
is $in->format, 'fastq';

open my $fh, '<', test_input_file('test.genbank') or die "Could not read file 'test.genbank': $!\n";
$in = Bio::SeqIO->new( -fh => $fh );
is $in->format, 'genbank';
close $fh;

my $string = ">seq\nACGATCG\n";
$in = Bio::SeqIO->new( -string => $string );
is $in->format, 'fasta';


############ EXCEPTION HANDLING ############

TODO: {
    local $TODO = 'file/fh-based tests should be in Bio::Root::IO, see issue #3204';
    throws_ok {
        Bio::SeqIO->new();
    } qr/No file, fh, or string argument provided/, 'Must pass a file or file handle';
}

throws_ok {
    Bio::SeqIO->new(-fh => undef);
} qr/fh argument provided, but with an undefined value/,
    'Must pass a file or file handle';

throws_ok {
    Bio::SeqIO->new(-file => undef);
} qr/file argument provided, but with an undefined value/,
    'Must pass a file or file handle';

throws_ok {
    Bio::SeqIO->new(-file => 'foo.bar');
} qr/Could not read file 'foo.bar':/,
    'Must pass a real file';
