
use Bio::Factory::Pise;

$email = $ENV{USER} . "\@pasteur.fr";

my $factory = Bio::Factory::Pise->new(-email => $email);

my $dnadist = $factory->program('dnadist');
my $job = $dnadist->run(-infile => $ARGV[0]);

if ($job->error) {
    print ".............error: ",$job->error_message,".............\n";
    exit;
}

print STDERR "jobid: ", $job->jobid, "\n";
print $job->content('outfile');

my $neighbor = $factory->program('neighbor',
				 -infile => $job->fh('outfile'));
my $job = $neighbor->run;
if ($job->error) {
    print ".............error: ",$job->error_message,".............\n";
    exit;
}

print $job->content('outtree');


