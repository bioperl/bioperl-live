# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
    use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 15);
	
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
