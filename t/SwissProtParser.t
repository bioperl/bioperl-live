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

open (DATA,'./t/data/SwissProt.dat') or die "no such file!";
my $parser = Bio::Tools::SwissProtParser->new(\*DATA);
ok defined $parser;			#2

my $report;
my $t;
ok $report = $parser->nextSbjct;	#3
ok $report->ID;			#4
ok $report->AC;			#5
ok $report->AC;		#6
ok $report->DT;			#7
ok $report->DE;			#8
ok $report->DE;		#9
ok $report->GN;			#10
ok $report->GN;		#11
ok $report->OS;			#12
ok $report->OG;		#13
ok $report->OG;			#14
ok $report->OC;		#15
ok $report->OC;			#16
ok $report->OX;		#17
ok $report->OX;			#18
ok $report->KW;		#19
ok $report->KW;			#20
ok $report->SS;			#21
ok $report->FT;		#22
ok $report->FT;			#23
ok $report->SQ;		#24
ok $report->SQ;			#25
ok $report->DR;		#26
ok $report->DR;			#27
ok $report->RR;		#28
ok $report->RR;			#29
ok $report->CC;		#30
ok $report->CC;			#31





