use Bio::Tools::Blast;
use Bio::DB::GenBank;
use Bio::SeqIO;

die("usage: blast_fetch blastfile outputfile") if( @ARGV != 2 );
my $db = new Bio::DB::GenBank;
my $blast = Bio::Tools::Blast->new(-file   =>@ARGV[0],
					-signif => 1e-5,
					-parse  => 1,
					-stats  => 1,
					-check_all_hits => 1,
					);
my $seqio = new Bio::SeqIO(-format=>'Fasta', -file=>">$ARGV[1]");
my ($an,$seq);
foreach my $hit ( $blast->hits) {
    (undef,$an) = split(/\|/, $hit->name);
    $seq = $db->get_Seq_by_id($an);
    if( defined $seq ) {
	$seqio->write_seq($seq);
    } 
}
