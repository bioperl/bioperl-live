#!/usr/bin/perl
=pod

=head1 NAME

extract_genes.pl - extract genomic sequences from NCBI files
using BioPerl

=head1 DESCRIPTION

This script is a simple solution to the problem of
extracting genomic regions corresponding to genes. There are other 
solutions, this particular approach uses genomic sequence 
files from NCBI and gene coordinates from Entrez Gene.

The first time this script is run it will be slow as it will
extract species-specific data from the gene2accession file and create
a storable hash (retrieving the positional data from this hash is
significantly faster than reading gene2accession each time the script
runs). The subsequent runs should be fast.

=head1 INSTALLATION

=head2

Install BioPerl, full instructions at http://bioperl.org.

=head2 Download gene2accession.gz

Download this file from ftp://ftp.ncbi.nlm.nih.gov/gene/DATA into 
your working directory and gunzip it.

=head2 Download sequence files

Create one or more species directories in the working directory, the
directory names do not have to match those at NCBI (e.g. "Sc", "Hs").

Download the nucleotide fasta files for a given species from its CHR*
directories at ftp://ftp.ncbi.nlm.nih.gov/genomes and put these files into a 
species directory. The sequence files will have the suffix ".fna" or 
"fa.gz", gunzip if necessary.

=head2 Determine Taxon id

Determine the taxon id for the given species. This id is the first column
in the gene2accession file. Modify the %species hash in this script
such that name of your species directory is a key and the taxon id is the 
value.

=head2 Command-line options

  -i   Gene id
  -s   Name of species directory
  -h   Help

Example:

  extract_genes.pl -i 850302 -s Sc

=cut

use strict;
use Bio::DB::Fasta;
use Getopt::Long;
use Storable;

my %species = ( "Sc" => 4932,  # Saccharomyces cerevisiae
				     "Ec" => 83333, # Escherichia coli K12
					  "Hs" => 9606   # H. sapiens
				   );

my ($help,$id,$name);

GetOptions( "s=s"  =>  \$name,
            "i=i"  =>  \$id,
				"h"    =>  \$help );

usage() if ($help || !$id || !$name);

my $storedHash = $name . ".dump";

# create index for a directory of fasta files
my $db = Bio::DB::Fasta->new($name, -makeid => \&make_my_id);

# extract species-specific data from gene2accession
unless (-e $storedHash) {
	my $ref;
	# extract species-specific information from gene2accession
	my $file = 'gene2accession';
	open my $MYIN, '<', $file or die "Could not read file '$file': $!\n";
	while (my $line = <$MYIN>) {
		my @arr = split "\t", $line;
		if ($arr[0] == $species{$name} && $arr[9] =~ /\d+/ && $arr[10] =~ /\d+/) {
			($ref->{$arr[1]}->{"start"}, $ref->{$arr[1]}->{"end"}, 
			 $ref->{$arr[1]}->{"strand"}, $ref->{$arr[1]}->{"id"}) =	
				($arr[9], $arr[10], $arr[11], $arr[7]);
		}
	}
	close $MYIN;
	# save species-specific information using Storable
	store $ref, $storedHash;
} 

# retrieve the species-specific data from a stored hash
my $ref = retrieve($storedHash);

# retrieve sequence and sub-sequence
if (defined $ref->{$id}) {
	my $chr = $db->get_Seq_by_id($ref->{$id}->{"id"});
	my $seq = $chr->trunc($ref->{$id}->{"start"},$ref->{$id}->{"end"});
	$seq = $seq->revcom if ($ref->{$id}->{"strand"} eq "-");

	# Insert SeqIO options here...
	print $seq->seq,"\n";
} else {
	print "Cannot find id: $id\n";
}

sub make_my_id {
	my $line = shift;
	$line =~ /ref\|([^|]+)/;
	$1;
}

sub usage {
	system "perldoc $0";
	exit;
}

__END__
