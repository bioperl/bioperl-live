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
ok !ref $report->ID;			#4
ok !ref $report->AC;			#5
ok @_ = $report->AC;		#6
ok !ref $report->DT;			#7
ok !ref $report->DE;			#8
ok @_ = $report->DE;		#9
ok !ref $report->GN;			#10
ok @_ = $report->GN;		#11
ok !ref $report->OS;			#12
ok @_ = $report->OG;		#13
ok !ref $report->OG;			#14
ok @_ = $report->OC;		#15
ok !ref $report->OC;			#16
ok @_ = $report->OX;		#17
ok !ref $report->OX;			#18
ok @_ = $report->KW;		#19
ok !ref $report->KW;			#20
ok !ref $report->SS;			#21
ok @_ = $report->FT;		#22
ok !ref $report->FT;			#23
ok @_ = $report->SQ;		#24
ok !ref $report->SQ;			#25
ok @_ = $report->DR;		#26
ok !ref $report->DR;			#27
ok @_ = $report->RR;		#28
ok !ref $report->RR;			#29
ok @_ = $report->CC;		#30
ok !ref $report->CC;			#31
