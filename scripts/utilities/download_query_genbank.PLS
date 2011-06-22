#!/usr/bin/perl
#$Id$

=head1 NAME

download_query_genbank - script to query Genbank and retrieve records

=head1 USAGE

 download_query_genbank --query "Neurospora[ORGN]" --db nucest -o Ncrassa_ESTs.fa --format fasta

 download_query_genbank --queryfile 'filewithquery' --db nucest -o Ncrassa_ESTs.fa --format fasta 

=head2 Other options

 Provide ONE of:

  --query query string OR
  --queryfile profile file with query OR
  --gi --gifile file with list of GIs to download

 Database type:

 --db database (nucest, protein, nucleotide)

 -f --format sequence file output format
 -v --verbose debugging output

=head2 Query options

 --maxids maximum number of IDs to retrieve in a set
 --reldate 
 --maxdate maxdate for a record
 --mindate minimum date for record
 --datetype edat or mdat (entered or modified)

=head1 AUTHOR Jason Stajich

Jason Stajich, jason-AT-bioperl.org

=cut

use strict;
use Bio::DB::GenBank;
use Bio::DB::GenPept;
use Bio::DB::Query::GenBank;
use Bio::SeqIO;
use Getopt::Long;

my ($queryfile,$outfile,$format,$debug,%options);

$format = 'fasta';

$options{'-maxids'} = '100';
$options{'-db'} = 'nucleotide'; # can be nucleotide, nucest, protein 
my $gifile;
GetOptions(
		   'h|help' => sub { exec('perldoc', $0); 
									exit(0);
								},
			  'v|verbose'       => \$debug,
			  'f|format:s'      => \$format,
			  'queryfile:s'     => \$queryfile,
			  'o|out|outfile:s' => \$outfile,
			  'gi|gifile|gis:s' => \$gifile,
			  # DB::Query options	   
			  'd|db:s'     => \$options{'-db'},
			  'mindate:s'  => \$options{'-mindate'},
			  'maxdate:s'  => \$options{'-maxdate'},
			  'reldate:s'  => \$options{'-reldate'}, 
			  'datetype:s' => \$options{'-datetype'}, # edat or mdat
			  'maxids:i'   => \$options{'-maxids'},
			  'q|query:s'  => \$options{'-query'},
			 );

my $out;

if( $outfile ) {
	$out = Bio::SeqIO->new(-format => $format,
								  -file   => ">$outfile");
} else {
	$out = Bio::SeqIO->new(-format => $format); # write to STDOUT
}

my $dbh;
if( $options{'-db'} eq 'protein' ) {
	$dbh = Bio::DB::GenPept->new(-verbose => $debug);
} else {
	$dbh = Bio::DB::GenBank->new(-verbose => $debug);
}
my $query;
if( $gifile ) {
	my @ids;
	open( my $fh => $gifile ) || die $!;
	while(<$fh>) {
		push @ids, split;
	}
	close($fh);	
	while( @ids ) {
		my @mini_ids = splice(@ids, 0, $options{'-maxids'});
		$query = Bio::DB::Query::GenBank->new(%options,
														  -ids => \@mini_ids,
														 );
		my $stream = $dbh->get_Stream_by_query($query);
		while( my $seq = $stream->next_seq ) {
			$out->write_seq($seq);
		}
	}
	exit;
} elsif( $options{'-query'}) {
	$query = Bio::DB::Query::GenBank->new(%options);
} elsif( $queryfile ) {
	open(my $fh => $queryfile) || die $!;
	while(<$queryfile>) {
		chomp;
		$options{'-query'} .= $_;
	}
	$query = Bio::DB::Query::GenBank->new(%options);
	close($fh);
} else {
	die("no query string or gifile\n");
}
my $stream = $dbh->get_Stream_by_query($query);
while( my $seq = $stream->next_seq ) {
	$out->write_seq($seq);
}
