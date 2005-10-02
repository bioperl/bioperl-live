# -*-Perl-*- mode (to keep my emacs happy)
# $Id$

use strict;

BEGIN {
    eval { require Test; };
    if ( $@ ) {
	use lib 't';
    }
    use Test;
    plan tests => 48;
}

use Bio::SeqIO;
use Bio::Root::IO;

ok(1);

my $verbose = $ENV{'BIOPERLDEBUG'} || 0;

my $ast = Bio::SeqIO->new(-format => 'tigrxml' ,
			  -verbose => $verbose,
			  -file => Bio::Root::IO->catfile
			  (qw(t data test.tigrxml)));
$ast->verbose($verbose);
my $as = $ast->next_seq();
ok($as);
ok($as->display_id, 'chr9');

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
	    ok($f->start, 185408);
	    ok($f->end, 187155);
	    # warn($f->gff_string, "\n");
	} elsif( $f->primary_tag eq 'mRNA' ) { 
	    ok($f->start, 185408); # the values list for COORD are start/end of CDS not whole transcript
	    ok($f->end, 187155);    
	    ok($f->strand, 1);
	} elsif( $f->primary_tag eq "five_prime_UTR" ) {
	    my ($id) = $f->get_tag_values('ID');
	    if( $id =~ /UTR1$/ ) {
		ok($f->start, 185408);
		ok($f->end,   185433);
	    } elsif( $id =~ /UTR2$/ ) {
		ok($f->start, 185487);
		ok($f->end,   185793);
	    } else {
		ok(0, , 'expected only two UTRS');
	    }	    
	} elsif( $f->primary_tag eq "three_prime_UTR" ) {
	    ok($f->start, 187042);
	    ok($f->end, 187155);
	} elsif( $f->primary_tag eq 'CDS' ) {
	    ok($f->start, 185794);
	    ok($f->end, 187041);
	}
    } elsif ( $name eq '162.t00448' || $name eq '162.m02967' ) {
	if( $f->primary_tag eq 'gene' ) {
	    ok($f->start, 59343);
	    ok($f->end, 61061);
	} elsif( $f->primary_tag eq 'mRNA' ) { 
	    ok($f->start, 59343); # the values list for COORD are start/end of CDS not whole transcript
	    ok($f->end, 61061);    
	    ok($f->strand, -1);
	} elsif( $f->primary_tag eq "five_prime_UTR" ) {
	    my ($id) = $f->get_tag_values('ID');
	    ok($f->start, 60834);
	    ok($f->end, 61061);
	    ok($f->strand, -1);
	} elsif( $f->primary_tag eq "three_prime_UTR" ) {
	    ok($f->start, 59343);
	    ok($f->end,   59632);
	    ok($f->strand, -1);
	} elsif( $f->primary_tag eq 'CDS' ) {
	    if( $first ) { 
		ok($f->start, 60801);
		ok($f->end,   60833);
		ok($f->strand, -1);
		$first = 0;
	    }
	}
    } else { 
	warn("name is $name\n");
    }
}
