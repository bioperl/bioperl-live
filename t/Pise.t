use lib '/home1/sis/letondal/packages/bioperl-live';

use strict;
BEGIN {
    use vars qw($NTESTS);
    # to handle systems with no installed Test module
    # we include the t dir (where a copy of Test.pm is located)
    # as a fallback
    eval { require Test; };
    if( $@ ) {
	use lib 't';
    }
    use Test;
    $NTESTS = 8;
    plan tests => $NTESTS }

use Bio::Factory::Pise;
use Bio::Tools::Genscan;
use Bio::SeqIO;

END { 

}

my $verbose = $ENV{'BIOPERLDEBUG'} || -1;
ok(1);

my $factory = Bio::Factory::Pise->new(-email => 'letondal@pasteur.fr');
ok($factory);
my $golden = $factory->program('golden', 
			       -db => 'swissprot', -query => '100K_RAT');
ok($golden->isa('Bio::Tools::Run::PiseApplication::golden'));
my $job = $golden->run();
ok($job->isa('Bio::Tools::Run::PiseJob'));
if ($job->error) {
    print STDERR "Error: ", $job->error_message, "\n";
}
ok(! $job->error);

# testing a program with a Bio::Seq + bioperl parsing of output
my $in = Bio::SeqIO->new ( -file   => 't/data/Genscan.FastA',
			   -format => 'fasta');

my $seq = $in->next_seq();
my $genscan = $factory->program('genscan',
				-parameter_file => "HumanIso.smat"
				);
ok($genscan);
$genscan->seq($seq);
my $job = $genscan->run();
ok($job);
my $parser = Bio::Tools::Genscan->new(-fh => $job->fh('genscan.out'));
ok($parser->isa('Bio::Tools::Genscan'));


