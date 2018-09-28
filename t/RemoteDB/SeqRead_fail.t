# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
    use Bio::Root::Test;

    test_begin(-tests => 12,
               -requires_modules => [qw(IO::String
                                        LWP::UserAgent
                                        HTTP::Request::Common)],
               -requires_networking => 1);
}

my $verbose = test_debug();

sub fetch {
    my ($id, $class) = @_;
    print "###################### $class  ####################################\n" if $verbose;
    my $seq;
    ok defined( my $gb = $class->new('-verbose'       => $verbose,
                                     '-delay'         => 0,
                                     '-retrievaltype' => 'tempfile') ), "defined for $class";

    if ($class eq 'Bio::DB::SwissProt') {
        test_skip(-tests => 1, -requires_module => 'Data::Stag');
        next if $@;
    }

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

## This is really stupid since many of this modules are not longer
## part of this distribution.  However, they are split over many
## distributions and we don't want to have this test code duplicated
## all over the place.  We should instead have this a Bio::Test module
## but that's work.  See bioperl-live issue #290
for my $class (@classes) {
    SKIP: {
        eval "require $class";
        skip "failed to use $class (guessing it's not available)", 2 if $@;
        fetch($id, $class);
    }
}
