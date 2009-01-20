# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN { 
	use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 18,
			   -requires_modules => [qw(IO::String
									    LWP::UserAgent
										HTTP::Request::Common)],
			   -requires_networking => 1);
	
	use_ok('Bio::DB::GenBank');
	use_ok('Bio::DB::GenPept');
	use_ok('Bio::DB::SwissProt');
	use_ok('Bio::DB::RefSeq');
	use_ok('Bio::DB::EMBL');
	use_ok('Bio::DB::BioFetch');
}

my $verbose = test_debug();

sub fetch {
    my ($id, $class) = @_;
    print "###################### $class  ####################################\n" if $verbose;
    my $seq;
    ok defined( my $gb = $class->new('-verbose'=>$verbose,
									 '-delay'=>0,
									 '-retrievaltype' => 'tempfile') ), "defined for $class";
    eval { $seq = $gb->get_Seq_by_id($id) };
    if ($@ || !defined $seq) {
		ok 1, "error or undef for $class";
		return;
    }
    ok 0, "failure for $class";
}

my @classes = qw( Bio::DB::BioFetch Bio::DB::GenBank Bio::DB::GenPept
                  Bio::DB::SwissProt Bio::DB::RefSeq Bio::DB::EMBL );

my $id = 'XXX111';  # nonsense id

for (@classes) {
    fetch($id, $_);
}
