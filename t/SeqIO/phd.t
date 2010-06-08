# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
    use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 21);
	
	use_ok('Bio::SeqIO::phd');
}

my $DEBUG = test_debug();

print("Checking to see if Bio::Seq::Quality objects can be created from a file...\n") if ($DEBUG);
my $in_phd  = Bio::SeqIO->new('-file' => test_input_file('test.phd'),
			      '-format'  => 'phd',
			      '-verbose' => $DEBUG);
isa_ok($in_phd,'Bio::SeqIO::phd');



my $phd = $in_phd->next_seq();
is($phd->quality_levels,'99',"Did you get the 'QUALITY_LEVELS' comment?");
isa_ok($phd,"Bio::Seq::Quality");


if( $DEBUG ) {
    my $position = 6;
    print("I saw these in phredfile.phd:\n\n");
    print $_->tagname,": ",$_->display_text || 0," \n"
        for ($phd->annotation->get_Annotations('header'));

    print("What is the base at position $position (using subseq)?\n");
    print($phd->subseq($position,$position)."\n");
    print("What is the base at position $position (using baseat)?\n");
    print($phd->baseat($position)."\n");
    print("What is the quality at $position? (using subqual)\n");

my @qualsretr = @{$phd->subqual($position,$position)};
    print($qualsretr[0]."\n");
    print("What is the quality at $position? (using qualat)\n");
    print($phd->qualat($position)."\n");
    print("What is the trace at $position? (using trace_index_at)\n");
    print($phd->trace_index_at($position)."\n");
    print("What is the trace at $position? (using subtrace)\n");
    my @tracesretr = @{$phd->subtrace($position,$position)};
    print($tracesretr[0]."\n");
}

print("OK. Now testing write_phd...\n") if($DEBUG);

my $outfile = test_output_file();
my $out_phd = Bio::SeqIO->new(-file => ">$outfile",
			      '-format' => 'phd');
isa_ok($out_phd,"Bio::SeqIO::phd");

$out_phd->write_seq($phd);

ok( -s $outfile);

# Bug 2120

my @qual = q(9 9 12 12 8 8 9 8 8 8 9);
my @trace = q(113 121 130 145 153 169 177 203 210 218 234);

$in_phd  = Bio::SeqIO->new('-file' => test_input_file('bug2120.phd'),
			      '-format'  => 'phd',
			      '-verbose' => $DEBUG);

my $seq = $in_phd->next_seq;
is($seq->subseq(10,20),'gggggccttat','$seq->subseq()');
my @seq_qual =$seq->subqual_text(10,20);
is_deeply(\@seq_qual,\@qual,'$seq->subqual_tex()');
my @seq_trace = $seq->subtrace_text(10,20);
is_deeply(\@seq_trace,\@trace,'$seq->subqual_tex()');

if($DEBUG) {
    print "\nDefault header ... \n\n";
    use Bio::Seq::Quality;
    my $seq = Bio::Seq::Quality->new('-seq' => 'GAATTC');
    $out_phd->_fh(\*STDOUT);
    $out_phd->write_header($seq);
    print "Complete output\n\n";
    $out_phd->write_seq($seq);
}

print("Testing the header manipulation\n") if($DEBUG);
is($phd->chromat_file(),'ML4924R','$phd->chromat_file()');
$phd->chromat_file('ML4924R.esd');
is($phd->chromat_file(), 'ML4924R.esd','$phd->chromat_file()');
$phd->touch();

# Commented out 1/17/09.
# This isn't exactly a stable regression test as the comparison tests
# localtime() called from two different timepoints. They can differ if the calls
# occurred before and after a change in seconds, for example.

#my $localtime = localtime();
#is($phd->time, "$localtime", $phd->time.':'.$localtime);

if ($DEBUG){
    print "Testing the sequence ...\n";
    print ">",$phd->id," ",$phd->desc,"\n",$phd->seq,"\n";
    my $revcom = $phd->revcom;
    print ">revcom\n",$revcom->seq,"\n";
    print ">revcom_qual at 6\n",$revcom->qualat(6),"\n";
    print ">revcom_trace at 6 !!\n",$revcom->trace_index_at(6),"\n";
    my $trunc = $phd->trunc(10,20);
    print ">TRUNC 10,20\n",$trunc->seq,"\n>qual\n@{$trunc->qual}\n>trace\n@{$trunc->trace}\n";
}

# Multiple seqs in one file

$in_phd  = Bio::SeqIO->new('-file' => test_input_file('multi.phd'),
			      '-format'  => 'phd',
			      '-verbose' => $DEBUG);

@qual = qq(9 9 15 17 17 22 22 25 25 22 22);
@trace = qq(98 105 119 128 143 148 162 173 185 197 202);

$seq = $in_phd->next_seq;
is($seq->id, 'ML4924F');
is($seq->subseq(10,20),'tctcgagggta','$seq->subseq()');
@seq_qual =$seq->subqual_text(10,20);
is_deeply(\@seq_qual,\@qual,'$seq->subqual_tex()');
@seq_trace = $seq->subtrace_text(10,20);
is_deeply(\@seq_trace,\@trace,'$seq->subqual_tex()');

@qual = qq(11 9 6 6 9 19 20 32 34 34 39);
@trace = qq(98 104 122 128 140 147 159 167 178 190 200);

$seq = $in_phd->next_seq;
is($seq->id, 'ML4924R');
is($seq->subseq(10,20),'gcctgcaggta','$seq->subseq()');
@seq_qual =$seq->subqual_text(10,20);
is_deeply(\@seq_qual,\@qual,'$seq->subqual_tex()');
@seq_trace = $seq->subtrace_text(10,20);
is_deeply(\@seq_trace,\@trace,'$seq->subqual_tex()');

#if($DEBUG) {
#    print "\nDefault header ... \n\n";
#    use Bio::Seq::Quality;
#    my $seq = Bio::Seq::Quality->new('-seq' => 'GAATTC');
#    $out_phd->_fh(\*STDOUT);
#    $out_phd->write_header($seq);
#    print "Complete output\n\n";
#    $out_phd->write_seq($seq);
#}

##print("Testing the header manipulation\n") if($DEBUG);
#is($phd->chromat_file(),'ML4924R','$phd->chromat_file()');
#$phd->chromat_file('ML4924R.esd');
#is($phd->chromat_file(), 'ML4924R.esd','$phd->chromat_file()');
#$phd->touch();
#my $localtime = localtime();
#is($phd->time, "$localtime");
#if ($DEBUG){
#    print "Testing the sequence ...\n";
#    print ">",$phd->id," ",$phd->desc,"\n",$phd->seq,"\n";
#    my $revcom = $phd->revcom;
#    print ">revcom\n",$revcom->seq,"\n";
#    print ">revcom_qual at 6\n",$revcom->qualat(6),"\n";
#    print ">revcom_trace at 6 !!\n",$revcom->trace_index_at(6),"\n";
#    my $trunc = $phd->trunc(10,20);
#    print ">TRUNC 10,20\n",$trunc->seq,"\n>qual\n@{$trunc->qual}\n>trace\n@{$trunc->trace}\n";
#}
#

# Whole-read tags in the file
$in_phd  = Bio::SeqIO->new('-file' => test_input_file('multiseq_tags.phd'),
			      '-format'  => 'phd',
			      '-verbose' => $DEBUG);
isa_ok($in_phd,'Bio::SeqIO::phd');
my @seqs = ();
while (my $seq = $in_phd->next_seq){
  push @seqs, $seq;
}
is( scalar @seqs, 2 );
