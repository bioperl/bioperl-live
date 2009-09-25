# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
    use lib '.';
	use Bio::Root::Test;
	
	test_begin(-tests => 27);
	
    use_ok('Bio::Tools::ECnumber');
}



my $EC1 = Bio::Tools::ECnumber->new( -ec_string => " EC  01. 02.03.00022  ",
                                     -comment   => "is 1.2.3.22" );


my $EC2 = Bio::Tools::ECnumber->new( -ec_string => "ec:1.2.3.-",
                                     -comment   => "is 1.2.3.-" );


my $EC3 = $EC1->copy();

isa_ok( $EC1,"Bio::Tools::ECnumber" );

isa_ok( $EC3,"Bio::Tools::ECnumber");

is( $EC1->EC_string(), "1.2.3.22" );

is( $EC1->EC_string(), "1.2.3.22" );

is( $EC1->to_string(), "1.2.3.22" );

is( $EC1->comment(),   "is 1.2.3.22" );

is( $EC1->enzyme_class(), "1" );

is( $EC1->sub_class(), "2" );

is( $EC1->sub_sub_class(), "3" );

is( $EC1->serial_number(), "22" );

ok( $EC3->is_equal( $EC1 ) );

ok( $EC3->is_equal( "1.2.3.22" ) );

ok( ! $EC3->is_equal( "1.2.3.-" ) );

ok( ! $EC3->is_equal( "1.2.3.23" ) );

ok( $EC1->is_member( $EC2 ) );

ok( $EC1->is_member( "1.2.3.-" ) );

$EC1->init();

ok( $EC2->is_member( $EC1 ) );

is( $EC1->to_string(), "-.-.-.-" );

$EC1->enzyme_class( 44 );

$EC1->sub_class( "033" );

$EC1->sub_sub_class( 22 );

$EC1->serial_number( "-" );

is( $EC1->to_string(), "44.33.22.-" );

ok( ! $EC1->is_member( "44.33.23.-" ) );

ok( ! $EC1->is_member( "44.33.22.1" ) );

ok( $EC1->is_member( "-.-.-.-" ) );

ok( $EC1->is_member( "44.-.-.-" ) );

ok( $EC1->is_member( "44.33.-.-" ) );

ok( $EC1->is_member( "EC 44.33.22.-" ) );

ok( ! $EC1->is_member( "45.33.22.-" ) );
