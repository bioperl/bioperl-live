#!/usr/local/bin/perl -w
use Bio::SearchIO;
use Bio::DB::GenBank;
use Bio::DB::GenPept;
use Bio::SeqIO;
use Getopt::Long;

use strict;
use vars qw($usage $alphabet $file $db $blastfile $outformat $ext);
$usage = "usage: blast_fetch [ -alphabet protein|nuc ] blastfile\n";

$alphabet = 'DNA';
my $result = GetOptions( 'help' => sub { print STDERR "$usage"; exit; },
			 'alphabet=s' => sub { $alphabet = &assign_type(pop @_)},
			 );
( $blastfile) = shift @ARGV;
unless ( defined $blastfile && -e $blastfile ) { print STDERR "$usage"; exit }

if( $alphabet eq 'DNA' ) {
    $db = new Bio::DB::GenBank;
    $outformat = 'genbank';
    $ext = 'gb';
} elsif( $alphabet eq 'PROT' ) {
    $db = new Bio::DB::GenPept;
    $outformat = 'swiss';
    $ext = 'sp';
} else { 
    die $usage;
}

my $blast = Bio::SearchIO->new(-file   =>$blastfile,
			       -format => 'blast'
			       );

my ($an,$seq);
my $thres ='1e-5';

while( my $result = $blast->next_result ) {
  HIT: while( my $hit = $result->next_hit ) {
      while( my $hsp = $hit->next_hsp ) {
	  if( $hsp->evalue > $thres ) { next HIT; }
      }     
      eval {
	  ($an) = split(/\s+/,$hit->name);
	  $seq = $db->get_Seq_by_id($an);
	  if( defined $seq ) {
	      my ($subname) = split(/\s+/,$seq->display_id);
	      my $seqio = new Bio::SeqIO(-format=>$outformat, 
					 -file=>">$subname.$ext");
	      $seqio->write_seq($seq);	    
	  } 
      };
      if( $@ ) {
	  print STDERR $@, "\n", " for name, an $an ", $hit->name, "\n";
      }
  }
}

sub assign_type {
    my $type = shift @_;
    if( $type =~ /[dr]na/i || $type =~ /nuc/i  ) {
	return 'DNA';
    } elsif( $type =~ /prot/i ) {
	return 'PROT';
    } else {
	return '';
    }
}
