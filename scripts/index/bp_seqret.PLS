#!/usr/bin/perl -w
# $Id: bp_seqret.PLS,v 1.2 2005-08-27 19:58:23 jason Exp $
# -*-Perl-*- mode (for emacs)

=head1 NAME

bp_seqret - bioperl implementation of sequence fetch from local db (like EMBOSS seqret)

=head1 USAGE

bp_seqret [-f/--format outputformat] [-o/--out/--outfile outfile] [-d/--db dbname] [-i/--id/-s/--seqname seqname1]

Example usage:

   bp_seqret -f fasta -db db.fa -i seq1 -i seq2 > output.fas
   bp_seqret db.fa:seq1 output.fas
   bp_seqret db.fa:seq1 -o output.fas
   bp_seqret -db db.fa -o output.fas seq1 seq2 seq3
   bp_seqret -db db.fa seq1 seq2 seq3 output.fas
   bp_seqret -db db.fa seq1 seq2 seq3 - > output.fas  

The DB is expected to be a Fasta formatted sequence file with multiple
sequences.

Output format is Fasta by default.

If no output filename is provided then output is written to STDOUT.
Providing '-' as the output filename will accomplish the same thing.


=head1 AUTHOR

Jason Stajich jason_AT_bioperl-dot-org

=cut

use strict;
use Bio::DB::Fasta;
use Bio::SeqIO;
use Getopt::Long;

my $dbname;
my @names;
my $format = 'fasta';
my $outfile;
GetOptions(
	   'f|format:s'   => \$format,
	   'o|out|outfile:s' => \$outfile,
	'd|db|dbname:s'   => \$dbname,
	'i|id|s|seqname:s' => \@names);


if( ! $dbname ) {
    die "need a dbname\n" unless @ARGV;
    $dbname = shift @ARGV;	
    if( $dbname =~ s/^([^:]+):// ) {
	push @names, $dbname;
	$dbname = $1;
    }				
}

my $db = Bio::DB::Fasta->new($dbname, -glob => "*.{fa,fas,fsa,fasta,pep,aa,seq,cds,peps}");
if( ! $outfile ) {
    $outfile = pop @ARGV;
}
my $out;
if( $outfile ) {
    $out = Bio::SeqIO->new(-format => $format,
			   -file   => ">$outfile");
} else {
    $out = Bio::SeqIO->new(-format => $format);
}
for my $nm ( @names ) {    
    if( my $seq = $db->get_Seq_by_acc($nm) ) {
	my ($id,$desc) = split(/\s+/,$db->header($nm),2);
	$out->write_seq(Bio::PrimarySeq->new(-display_id => $id,
					     -description => $desc,
					     -seq => $seq->seq));
    } else {
	warn("$nm not found\n");
    }
}

