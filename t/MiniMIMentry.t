# -*-Perl-*-
## Bioperl Test Harness Script for Modules
## # $Id$

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

use strict;
BEGIN {
    # to handle systems with no installed Test module
    # we include the t dir (where a copy of Test.pm is located)
    # as a fallback
    eval { require Test::More; };
    if( $@ ) {
        use lib 't/lib';
    }
    use Test::More;
    plan tests => 15;
    use_ok('Bio::Phenotype::OMIM::MiniMIMentry');
}
  
my  $mm = Bio::Phenotype::OMIM::MiniMIMentry->new( -description  => "The central form of ...",
                                                   -created      => "Victor A. McKusick: 6/4/1986",
                                                   -contributors => "Kelly A. Przylepa - revised: 03/18/2002",
                                                   -edited       => "alopez: 06/03/1997" );

isa_ok( $mm, "Bio::Phenotype::OMIM::MiniMIMentry");

ok( $mm->to_string() );

is( $mm->description(), "The central form of ..." );
is( $mm->created(), "Victor A. McKusick: 6/4/1986" );
is( $mm->contributors(), "Kelly A. Przylepa - revised: 03/18/2002" );
is( $mm->edited(), "alopez: 06/03/1997" );

$mm->init();

is( $mm->description(), "" );
is( $mm->created(), "" );
is( $mm->contributors(), "" );
is( $mm->edited(), "" );


is( $mm->description( "A" ), "A" );
is( $mm->created( "B" ), "B" );
is( $mm->contributors( "C" ), "C" );
is( $mm->edited( "D" ), "D" );






