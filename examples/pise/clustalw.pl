
use Bio::Factory::Pise;
use Bio::AlignIO;

$in = Bio::AlignIO->new('-file' => $ARGV[0]);
$aln = $in->next_aln();

$email = ''; # your email
$factory = Bio::Factory::Pise->new(-email => $email);
my $clustalw = $factory->program('clustalw');

my $job = $clustalw->run( -infile => $aln);
if ($job->error) {
	print $job->error_message, "\n";
}

print STDERR "jobid: ", $job->jobid, "\n";



