#!/usr/bin/perl

use lib '../blib/lib';
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
    plan tests => 31 }

use Bio::Tools::SwissProtParser;
ok(1);					#1

open (DATA,'./data/SwissProt.dat') or die "no such file!";
my $parser = Bio::Tools::SwissProtParser->new(\*DATA);
ok defined $parser;			#2

my $report;
ok $report = $parser->nextSbjct;	#3
ok !ref $report->ID;			#4
ok !ref $report->AC;			#5
ok my($t) = $report->AC;		#6
ok !ref $report->DT;			#7
ok !ref $report->DE;			#8
ok my($t) = $report->DE;		#9
ok !ref $report->GN;			#10
ok my($t) = $report->GN;		#11
ok !ref $report->OS;			#12
ok my($t) = $report->OG;		#13
ok !ref $report->OG;			#14
ok my($t) = $report->OC;		#15
ok !ref $report->OC;			#16
ok my($t) = $report->OX;		#17
ok !ref $report->OX;			#18
ok my($t) = $report->KW;		#19
ok !ref $report->KW;			#20
ok !ref $report->SS;			#21
ok my($t) = $report->FT;		#22
ok !ref $report->FT;			#23
ok my($t) = $report->SQ;		#24
ok !ref $report->SQ;			#25
ok my($t) = $report->DR;		#26
ok !ref $report->DR;			#27
ok my($t) = $report->RR;		#28
ok !ref $report->RR;			#29
ok my($t) = $report->CC;		#30
ok !ref $report->CC;			#31
