#-*-Perl-*- mode
# $Id$

BEGIN {
    eval { require Test; };
    if( $@ ) {
        use lib 't';
    }
    use Test;
    use vars qw($NUMTESTS);
    $NUMTESTS = 20;
    plan tests => $NUMTESTS;
    @ARGV = (-1);
    require 'bptutorial.pl';
}

END {
    unlink 'bptutorial.out';
}

# run the first 21 tests
for my $test ( 1..6, 8..21 ) {
    ok(&run_examples($test));
}


