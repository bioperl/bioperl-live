
use Bio::Factory::Pise;
use Bio::SeqIO;
use Bio::Tools::Genscan;

$email = $ENV{USER} . "\@pasteur.fr";

$_in_seq = Bio::SeqIO->newFh (-file => $ARGV[0],
			      -format => "fasta");
my $seq = <$_in_seq>; 

my $factory = Bio::Factory::Pise->new(-email => $email);

my $genscan = $factory->program('genscan',
				-parameter_file => "Arabidopsis.smat"
				);
$genscan->seq($seq);
my $job = $genscan->run();

if ($job->error) {
    print ".............error: ",$job->error_message,".............\n";
    exit;
}

#print STDERR $job->content('genscan.out');

my $parser_file = $job->save('genscan.out');
$parser = Bio::Tools::Genscan->new(-file => $parser_file);
while(my $gene = $parser->next_prediction()) {
    my $prot = $gene->predicted_protein;
    print "protein (", ref($prot), ") :\n", $prot->seq, "\n\n";
}

print STDERR "jobid: ", $job->jobid, "\n";

