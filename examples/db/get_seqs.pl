#!/usr/local/bin/perl
use strict;
use vars qw($USAGE);
use Carp;
use Getopt::Long;
use Bio::SeqIO;

$USAGE = "get_seqs.pl\t[--db=DBNAME] [--format=FORMAT] \n\t\t[--output=FILENAME] [--proxy=PROXY] accession1, accession2, ...\n Defaults: db=GenBank format=fasta output=STDOUT proxy=none\n See LWP::UserAgent for more information on proxy syntax";
my %dbs = (
  'genbank'   => 'Bio::DB::GenBank',
  'embl'      => 'Bio::DB::EMBL',
  'swissprot' => 'Bio::DB::SwissProt', 
);

my ($db,$format,$file,$proxy,$help) = ( 'genbank', 'fasta' );

&GetOptions
    (
     'db:s'       => \$db,
     'f|format:s' => \$format,
     "file|out|output:s" => \$file,
     'proxy:s'           => \$proxy,
     "h|\?|help"  => \$help ,     
     );

if( $help ) { print $USAGE, "\n";exit; }

if( $db =~ /gb|gen|genbank/i ) {
    $db = 'genbank';
} elsif( $db =~ /embl|em|e/i ) {
    $db = 'embl';
} elsif( $db =~ /swiss|sp/i ) {
    $db = 'swissprot';
} else { 
    croak("Unknown db parameter '$db' valid parameters are (" . join(',', keys %dbs) . ")");
}

my %params = ( '-format' => $format );

if( defined $file ) {
    $params{'-file'} = ">$file";
} else { 
    $params{'-fh'} = \*STDOUT;
}

my $seqio = new Bio::SeqIO(%params);

my $remotedb;

eval {
    my $filename = "$dbs{$db}.pm";
    $filename =~ s!::!/!g;
    require $filename;
    $remotedb = "$dbs{$db}"->new();
};

die($@) unless ! $@;

if( defined $proxy ) { $remotedb->proxy($proxy); }

my $stream;

if( $remotedb->can('get_Stream_by_batch') ) {
    $stream = $remotedb->get_Stream_by_batch(@ARGV);
} else {
    $stream = $remotedb->get_Stream_by_acc(\@ARGV);
}

while( my $seq = $stream->next_seq ) {
    $seqio->write_seq($seq);
}
