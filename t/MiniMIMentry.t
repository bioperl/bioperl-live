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
    eval { require Test; };
    if( $@ ) {
        use lib 't';
    }
    use Test;
    plan tests => 14;
}

use Bio::Phenotype::OMIM::MiniMIMentry;
  
my  $mm = Bio::Phenotype::OMIM::MiniMIMentry->new( -description  => "The central form of ...",
                                                   -created      => "Victor A. McKusick: 6/4/1986",
                                                   -contributors => "Kelly A. Przylepa - revised: 03/18/2002",
                                                   -edited       => "alopez: 06/03/1997" );

ok( $mm->isa( "Bio::Phenotype::OMIM::MiniMIMentry" ) );

ok( $mm->to_string() );


ok( $mm->description(), "The central form of ..." );
ok( $mm->created(), "Victor A. McKusick: 6/4/1986" );
ok( $mm->contributors(), "Kelly A. Przylepa - revised: 03/18/2002" );
ok( $mm->edited(), "alopez: 06/03/1997" );

$mm->init();

ok( $mm->description(), "" );
ok( $mm->created(), "" );
ok( $mm->contributors(), "" );
ok( $mm->edited(), "" );


ok( $mm->description( "A" ), "A" );
ok( $mm->created( "B" ), "B" );
ok( $mm->contributors( "C" ), "C" );
ok( $mm->edited( "D" ), "D" );






