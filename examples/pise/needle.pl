use Bio::Factory::Pise;

my $factory = Bio::Factory::Pise->new(-email => ''); # put your email

my $needle = $factory->program('needle');

my $job = $needle->run(-sequencea => $ARGV[0],
		       -seqall => $ARGV[1],
		       -gapopen => 5,
		       -gapextend => 1);

if ($job->error) {
    print ".............error: ",$job->error_message,".............\n";
    exit;
}

print STDERR "jobid: ", $job->jobid, "\n";
print $job->content('outfile.align');


