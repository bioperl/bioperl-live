# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
    use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 49,
			   -requires_modules => [qw(XML::SAX XML::SAX::Writer)]);
	
	use_ok('Bio::SeqIO::tigrxml');
}

my $verbose = test_debug();

my $ast = Bio::SeqIO->new(-format => 'tigrxml' ,
			  -verbose => $verbose,
			  -file => test_input_file('test.tigrxml'));
isa_ok($ast, 'Bio::SeqIO');
$ast->verbose($verbose);
ok my $as = $ast->next_seq();
is($as->display_id, 'chr9');

my $first = 1;
for my $f ( sort { $a->start * $a->strand <=> $b->start * $b->strand } $as->get_SeqFeatures ) {
    ok($f);

    my ($name);
    for my $tag ( qw(Parent ID) ) {
	if( $f->has_tag($tag) ) {
	    ($name) = $f->get_tag_values($tag);
	    last;
	}
    }
    if( $name eq '162.t00500' || $name eq '162.m02638' ) {
	if( $f->primary_tag eq 'gene' ) {
	    is($f->start, 185408);
	    is($f->end, 187155);
	    # warn($f->gff_string, "\n");
	} elsif( $f->primary_tag eq 'mRNA' ) { 
	    is($f->start, 185408); # the values list for COORD are start/end of CDS not whole transcript
	    is($f->end, 187155);    
	    is($f->strand, 1);
	} elsif( $f->primary_tag eq "five_prime_UTR" ) {
	    my ($id) = $f->get_tag_values('ID');
	    if( $id =~ /UTR1$/ ) {
		is($f->start, 185408);
		is($f->end,   185433);
	    } elsif( $id =~ /UTR2$/ ) {
		is($f->start, 185487);
		is($f->end,   185793);
	    } else {
		ok(0, 'expected only two UTRS');
	    }	    
	} elsif( $f->primary_tag eq "three_prime_UTR" ) {
	    is($f->start, 187042);
	    is($f->end, 187155);
	} elsif( $f->primary_tag eq 'CDS' ) {
	    is($f->start, 185794);
	    is($f->end, 187041);
	}
    } elsif ( $name eq '162.t00448' || $name eq '162.m02967' ) {
	if( $f->primary_tag eq 'gene' ) {
	    is($f->start, 59343);
	    is($f->end, 61061);
	} elsif( $f->primary_tag eq 'mRNA' ) { 
	    is($f->start, 59343); # the values list for COORD are start/end of CDS not whole transcript
	    is($f->end, 61061);    
	    is($f->strand, -1);
	} elsif( $f->primary_tag eq "five_prime_UTR" ) {
	    my ($id) = $f->get_tag_values('ID');
	    is($f->start, 60834);
	    is($f->end, 61061);
	    is($f->strand, -1);
	} elsif( $f->primary_tag eq "three_prime_UTR" ) {
	    is($f->start, 59343);
	    is($f->end,   59632);
	    is($f->strand, -1);
	} elsif( $f->primary_tag eq 'CDS' ) {
	    if( $first ) { 
		is($f->start, 60801);
		is($f->end,   60833);
		is($f->strand, -1);
		$first = 0;
	    }
	}
    } else { 
	ok(0, "unexpected name '$name'\n");
    }
}
