#!/usr/local/bin/perl -w
use Bio::Tools::Blast;
use Bio::DB::GenBank;
use Bio::DB::GenPept;
use Bio::SeqIO;
use Getopt::Long;

use strict;
use vars qw($usage $moltype $file $db $blastfile );
$usage = "usage: blast_fetch [ -moltype protein|nuc ] blastfile\n";

$moltype = 'DNA';
my $result = GetOptions( 'help' => sub { print STDERR "$usage"; exit; },
			 'moltype=s' => sub { $moltype = &assign_type(pop @_)},
			 );
( $blastfile) = @ARGV;

if( $moltype eq 'DNA' ) {
    $db = new Bio::DB::GenBank;
} elsif( $moltype eq 'PROT' ) {
    $db = new Bio::DB::GenPept;
} else { 
    die $usage;
}

my $blast = Bio::Tools::Blast->new(-file   =>$blastfile,
					-signif => 1e-5,
					-parse  => 1,
					-stats  => 1,
					-check_all_hits => 1,
					);
my ($an,$seq);
foreach my $hit ( $blast->hits) {
    (undef,$an) = split(/\|/, $hit->name);
    if( !defined $an || $an eq '' ) {
	print STDERR "could not find accession number in line ", $hit->name(), ", skipping...\n";
	next;
    }
    eval { 
	$seq = $db->get_Seq_by_id($an);
	if( defined $seq ) {
	    my $seqio = new Bio::SeqIO(-format=>'swiss', -file=>">$an.sp");
	    $seqio->write_seq($seq);	    
	} 
    };
    if( $@ ) {
	print STDERR $@, "\n", " for name, an $an ", $hit->name, "\n";	
    }
}

sub assign_type {
    my $type = shift @_;
    if( $type =~ /[dr]na/i ) {
	return 'DNA';
    } elsif( $type =~ /prot/i ) {
	return 'PROT';
    } else {
	return '';
    }
}
