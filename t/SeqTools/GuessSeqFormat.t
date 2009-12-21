# -*-Perl-*- Test Harness script for Bioperl
# $Id$


use strict;

BEGIN {
   use lib '.';
   use Bio::Root::Test;
   
   test_begin(-tests => 52);
   
   use_ok('Bio::SeqIO');
   use_ok('Bio::AlignIO');
   use_ok('Bio::Tools::GuessSeqFormat');
}

my @seqformats = qw{ ace embl fasta fastq gcg genbank mase
                        pfam pir raw swiss tab game};

my $format;
my $verbose = test_debug();
#
# Seqio formats
#

#not tested:  waba

my %no_seqio_module = map {$_=>1} qw {gcgblast gcgfasta mase pfam};

my $guessed_format = Bio::Tools::GuessSeqFormat->new
        (-file => test_input_file('test.waba'))->guess;
is $guessed_format, undef ;

my $seq;

eval {
   my $input = Bio::SeqIO->new
       (-file=>test_input_file('test.waba'));
   $seq = $input->next_seq();
};

ok !$seq;

$@ ? ok 1 : ok 0;

foreach $format (@seqformats) {
   SKIP: {
      if ($format eq 'game') {
         test_skip(-tests => 2, -requires_modules => [qw(XML::Writer XML::Parser::PerlSAX)]);
      }
      
      my $guessed_format = Bio::Tools::GuessSeqFormat->new
          (-file => test_input_file("test.$format"),
           #-verbose=> $verbose;
          )->guess;
      $format =~ s/\..*$//;
      is $guessed_format, $format, "Guessed:$format";
      next if $no_seqio_module{$format};
     
      eval {
          my $input = Bio::SeqIO->new
              (-file=>test_input_file("test.$format"));
          $seq = $input->next_seq();
      };
      
      my $implemented = $format eq 'ace' ? 'Bio::PrimarySeqI' : 'Bio::SeqI';
      
      isa_ok $seq, $implemented;
      
      is 0, 1, $@ if $@;
   }
}

#
# AlignIO formats
#

@seqformats = qw{ aln:clustalw fasta fastq mase msf nexus pfam phylip
                  prodom stockholm}; # not selex (same as pfam, mainly)

my %no_alignio_module = map {$_=>1} qw { fastq };

foreach my $ext (@seqformats) {
    my $format;
    ($ext, $format) = split /:/, $ext;
    my $guesser = Bio::Tools::GuessSeqFormat->new
        (-file => test_input_file("testaln.$ext"));
    $format ||= $ext;
    ok $guesser->guess(), $format;

    next if $no_alignio_module{$format};

    eval {
        my $input = Bio::AlignIO->new
            (-file=>test_input_file("testaln.$ext"));
        $seq = $input->next_aln();
    };
    
    isa_ok $seq, 'Bio::Align::AlignI';
    #ok 0, 1, $@ if $@;
}


#
# File handle tests
#
SKIP: {
   test_skip(-tests => 3, -requires_modules => [qw(IO::String)]);

    my $string = ">test1 no comment
agtgctagctagctagctagct
>test2 no comment
gtagttatgc
";

    my $stringfh = new IO::String($string);
    
    my $seqio = Bio::SeqIO->new(-fh => $stringfh);
    while( my $seq = $seqio->next_seq ) {
	ok $seq->id =~ /test/;
    }
    
#
# text guessing
#

    ok new Bio::Tools::GuessSeqFormat( -text => $string )->guess, 'fasta';
}
