## Test for a suspected bug and tests for debugging.

use strict;
use warnings;
use Data::Dumper;

BEGIN {
  use Bio::Root::Test;

  test_begin( -tests => 174 );

  use_ok('Bio::Location::Simple');
  use_ok('Bio::Coordinate::Pair');
}


## Set up two coordinate systems for the same sequence

## The contig
ok my $ctg = Bio::Location::Simple->
  new( -seq_id => 'ctg',
       -start  =>    1,
       -end    => 1001,
       -strand =>   +1,
     );

isa_ok $ctg, 'Bio::Location::Simple';


## The contig in the chromosome
ok my $ctg_on_chr_f = Bio::Location::Simple->
  new( -seq_id => 'ctg on chr f',
       -start  =>           5001,
       -end    =>           6001,
       -strand =>             +1,
     );

isa_ok $ctg_on_chr_f, 'Bio::Location::Simple';


## The contig in the chromosome (again)
ok my $ctg_on_chr_r = Bio::Location::Simple->
  new( -seq_id => 'ctg on chr r',
       -start  =>           5001,
       -end    =>           6001,
       -strand =>             -1,
     );

isa_ok $ctg_on_chr_r, 'Bio::Location::Simple';





## Set up the mapping between them

ok my $agp_f = Bio::Coordinate::Pair->
  new( -in  => $ctg,
       -out => $ctg_on_chr_f
     );

isa_ok $agp_f, 'Bio::Coordinate::Pair';


ok my $agp_r = Bio::Coordinate::Pair->
  new( -in  => $ctg,
       -out => $ctg_on_chr_r
     );

isa_ok $agp_r, 'Bio::Coordinate::Pair';





## Perform some very basic sanity testing on the resulting map objects

## f

ok $agp_f->test;

is $agp_f->in->seq_id, 'ctg';
is $agp_f->in->start,      1;
is $agp_f->in->end,     1001;
is $agp_f->in->strand,    +1;

is $agp_f->out->seq_id, 'ctg on chr f';
is $agp_f->out->start,            5001;
is $agp_f->out->end,              6001;
is $agp_f->out->strand,             +1;


## r

ok $agp_r->test;

is $agp_r->in->seq_id, 'ctg';
is $agp_r->in->start,      1;
is $agp_r->in->end,     1001;
is $agp_r->in->strand,    +1;

is $agp_r->out->seq_id, 'ctg on chr r';
is $agp_r->out->start,            5001;
is $agp_r->out->end,              6001;
is $agp_r->out->strand,             -1;





##
## Map a particular match through both map objects
##

## Define the match 1
ok my $match_on_ctg_1 = Bio::Location::Simple->
  new( -seq_id => 'hit 1',
       -start  =>      25,
       -end    =>     225,
       -strand =>      +1,
     );

isa_ok $match_on_ctg_1, 'Bio::LocationI';


# Convert the match from contig into chromosomal coordinates

ok my $match_on_chr_1_f =
  $agp_f->map( $match_on_ctg_1 );

isa_ok $match_on_chr_1_f, 'Bio::Coordinate::Result';


ok my $match_on_chr_1_r =
  $agp_r->map( $match_on_ctg_1 );

isa_ok $match_on_chr_1_r, 'Bio::Coordinate::Result';



## Perform some very basic sanity testing on the match objects

is $match_on_ctg_1->seq_id, 'hit 1';
is $match_on_ctg_1->start,       25;
is $match_on_ctg_1->end,        225;
is $match_on_ctg_1->strand,      +1;

is $match_on_chr_1_f->seq_id, 'ctg on chr f';
is $match_on_chr_1_f->start,            5025;
is $match_on_chr_1_f->end,              5225;
is $match_on_chr_1_f->strand,             +1;

is $match_on_chr_1_r->seq_id, 'ctg on chr r';
is $match_on_chr_1_r->start,            5777;
is $match_on_chr_1_r->end,              5977;
is $match_on_chr_1_r->strand,             -1;



## Define the match 2
ok my $match_on_ctg_2 = Bio::Location::Simple->
  new( -seq_id => 'hit 2',
       -start  =>      25,
       -end    =>     225,
       -strand =>      -1,
     );

isa_ok $match_on_ctg_2, 'Bio::LocationI';


# Convert the match from contig into chromosomal coordinates

ok my $match_on_chr_2_f =
  $agp_f->map( $match_on_ctg_2 );

isa_ok $match_on_chr_2_f, 'Bio::Coordinate::Result';


ok my $match_on_chr_2_r =
  $agp_r->map( $match_on_ctg_2 );

isa_ok $match_on_chr_2_r, 'Bio::Coordinate::Result';




## Perform some very basic sanity testing on the match objects

is $match_on_ctg_2->seq_id, 'hit 2';
is $match_on_ctg_2->start,       25;
is $match_on_ctg_2->end,        225;
is $match_on_ctg_2->strand,      -1;

is $match_on_chr_2_f->seq_id, 'ctg on chr f';
is $match_on_chr_2_f->start,            5025;
is $match_on_chr_2_f->end,              5225;
is $match_on_chr_2_f->strand,             -1;

is $match_on_chr_2_r->seq_id, 'ctg on chr r';
is $match_on_chr_2_r->start,            5777;
is $match_on_chr_2_r->end,              5977;
is $match_on_chr_2_r->strand,             +1;







## OK, now we can get down to some debugging...



## TEST ONE

## Create a match that goes off the end of the contig

## Define the match 3
ok my $match_on_ctg_3 = Bio::Location::Simple->
  new( -seq_id => 'hit 3',
       -start  =>     925,
       -end    =>    1125,
       -strand =>      +1,
     );

isa_ok $match_on_ctg_3, 'Bio::LocationI';


# Convert the match from contig into chromosomal coordinates

ok my $match_on_chr_3_f =
  $agp_f->map( $match_on_ctg_3 );

isa_ok $match_on_chr_3_f, 'Bio::Coordinate::Result';


ok my $match_on_chr_3_r =
  $agp_r->map( $match_on_ctg_3 );

isa_ok $match_on_chr_3_r, 'Bio::Coordinate::Result';



## Perform some very basic sanity testing on the match objects

is $match_on_ctg_3->seq_id, 'hit 3';
is $match_on_ctg_3->start,      925;
is $match_on_ctg_3->end,       1125;
is $match_on_ctg_3->strand,      +1;

is $match_on_chr_3_f->seq_id, 'ctg on chr f';
is $match_on_chr_3_f->start,            5925;
isnt $match_on_chr_3_f->end,            6125; # Gets truncated to maximum!
is $match_on_chr_3_f->end,              6001; # Gets truncated to maximum!
is $match_on_chr_3_f->strand,             +1;

#print Dumper $match_on_ctg_3;
#print Dumper $match_on_chr_3_f;

is $match_on_chr_3_r->seq_id, 'ctg on chr r';
isnt $match_on_chr_3_r->start,          4877; # Gets truncated to minimum!
is $match_on_chr_3_r->start,            5001; # Gets truncated to minimum!
is $match_on_chr_3_r->end,              5077;
#is $match_on_chr_3_r->strand,             -1; # FAIL
is $match_on_chr_3_r->strand,          undef; # See Bio::Location::Split

#print Dumper $match_on_ctg_3;
#print Dumper $match_on_chr_3_r;



## Define the match 4
ok my $match_on_ctg_4 = Bio::Location::Simple->
  new( -seq_id => 'hit 4',
       -start  =>     925,
       -end    =>    1125,
       -strand =>      -1,
     );

isa_ok $match_on_ctg_4, 'Bio::LocationI';


# Convert the match from contig into chromosomal coordinates

ok my $match_on_chr_4_f =
  $agp_f->map( $match_on_ctg_4 );

isa_ok $match_on_chr_4_f, 'Bio::Coordinate::Result';


ok my $match_on_chr_4_r =
  $agp_r->map( $match_on_ctg_4 );

isa_ok $match_on_chr_4_r, 'Bio::Coordinate::Result';



## Perform some very basic sanity testing on the match objects

is $match_on_ctg_4->seq_id, 'hit 4';
is $match_on_ctg_4->start,      925;
is $match_on_ctg_4->end,       1125;
is $match_on_ctg_4->strand,      -1;

is $match_on_chr_4_f->seq_id, 'ctg on chr f';
is $match_on_chr_4_f->start,            5925;
isnt $match_on_chr_4_f->end,            6125; # Gets truncated to maximum!
is $match_on_chr_4_f->end,              6001; # Gets truncated to maximum!
is $match_on_chr_4_f->strand,             -1;

#print Dumper $match_on_ctg_4;
#print Dumper $match_on_chr_4_f;

is $match_on_chr_4_r->seq_id, 'ctg on chr r';
isnt $match_on_chr_4_r->start,          4877; # Gets truncated to minimum!
is $match_on_chr_4_r->start,            5001; # Gets truncated to minimum!
is $match_on_chr_4_r->end,              5077;
#is $match_on_chr_4_r->strand,             +1; # FAIL
is $match_on_chr_4_r->strand,          undef; # See Bio::Location::Split

#print Dumper $match_on_ctg_4;
#print Dumper $match_on_chr_4_r;







###
### NOW! NONE OF THE ABOVE SHOULD BE AFFECTED BY LEAVING OFF seq_id
### NOW SHOULD IT?!
###

## Try commenting out the three -seq_id lines below to observe strange
## interactions!

## The contig
ok my $ctg_x = Bio::Location::Simple->
  new( -seq_id => 'ctg',
       -start  =>    1,
       -end    => 1001,
       -strand =>   +1,
     );

isa_ok $ctg_x, 'Bio::Location::Simple';

## The contig in the chromosome
ok my $ctg_on_chr_f_x = Bio::Location::Simple->
  new( -seq_id => 'ctg on chr f',
       -start  =>           5001,
       -end    =>           6001,
       -strand =>             +1,
     );

isa_ok $ctg_on_chr_f_x, 'Bio::Location::Simple';

## The contig in the chromosome (again)
ok my $ctg_on_chr_r_x = Bio::Location::Simple->
  new( -seq_id => 'ctg on chr r',
       -start  =>           5001,
       -end    =>           6001,
       -strand =>             -1,
     );

isa_ok $ctg_on_chr_r_x, 'Bio::Location::Simple';



## Set up the mapping between them

ok my $agp_xf = Bio::Coordinate::Pair->
  new( -in  => $ctg_x,
       -out => $ctg_on_chr_f_x
     );

isa_ok $agp_xf, 'Bio::Coordinate::Pair';


ok my $agp_xr = Bio::Coordinate::Pair->
  new( -in  => $ctg_x,
       -out => $ctg_on_chr_r_x
     );

isa_ok $agp_xr, 'Bio::Coordinate::Pair';





## Perform some very basic sanity testing on the resulting map objects

## f

ok $agp_xf->test;

is $agp_xf->in->start,      1;
is $agp_xf->in->end,     1001;
is $agp_xf->in->strand,    +1;

is $agp_xf->out->start,            5001;
is $agp_xf->out->end,              6001;
is $agp_xf->out->strand,             +1;


## r

ok $agp_r->test;

is $agp_xr->in->start,      1;
is $agp_xr->in->end,     1001;
is $agp_xr->in->strand,    +1;

is $agp_xr->out->start,            5001;
is $agp_xr->out->end,              6001;
is $agp_xr->out->strand,             -1;





##
## Map a particular match through both map objects
##

# Convert the match from contig into chromosomal coordinates

ok my $match_on_chr_1_xf =
  $agp_xf->map( $match_on_ctg_1 );

isa_ok $match_on_chr_1_xf, 'Bio::Coordinate::Result';


ok my $match_on_chr_1_xr =
  $agp_xr->map( $match_on_ctg_1 );

isa_ok $match_on_chr_1_xr, 'Bio::Coordinate::Result';

## Perform some very basic sanity testing on the match objects

is $match_on_chr_1_xf->start,            5025;
is $match_on_chr_1_xf->end,              5225;
is $match_on_chr_1_xf->strand,             +1;

is $match_on_chr_1_xr->start,            5777;
is $match_on_chr_1_xr->end,              5977;
is $match_on_chr_1_xr->strand,             -1;



# Convert the match from contig into chromosomal coordinates

ok my $match_on_chr_2_xf =
  $agp_xf->map( $match_on_ctg_2 );

isa_ok $match_on_chr_2_xf, 'Bio::Coordinate::Result';


ok my $match_on_chr_2_xr =
  $agp_xr->map( $match_on_ctg_2 );

isa_ok $match_on_chr_2_xr, 'Bio::Coordinate::Result';

## Perform some very basic sanity testing on the match objects

is $match_on_chr_2_xf->start,            5025;
is $match_on_chr_2_xf->end,              5225;
is $match_on_chr_2_xf->strand,             -1;

is $match_on_chr_2_xr->start,            5777;
is $match_on_chr_2_xr->end,              5977;
is $match_on_chr_2_xr->strand,             +1;



# Convert the match from contig into chromosomal coordinates

ok my $match_on_chr_3_xf =
  $agp_xf->map( $match_on_ctg_3 );

isa_ok $match_on_chr_3_xf, 'Bio::Coordinate::Result';


ok my $match_on_chr_3_xr =
  $agp_xr->map( $match_on_ctg_3 );

isa_ok $match_on_chr_3_xr, 'Bio::Coordinate::Result';

## Perform some very basic sanity testing on the match objects

is $match_on_chr_3_xf->start,            5925;
isnt $match_on_chr_3_xf->end,            6125; # Gets truncated to maximum!
is $match_on_chr_3_xf->end,              6001; # Gets truncated to maximum!
is $match_on_chr_3_xf->strand,             +1;

isnt $match_on_chr_3_xr->start,          4877; # Gets truncated to minimum!
is $match_on_chr_3_xr->start,            5001; # Gets truncated to minimum!
is $match_on_chr_3_xr->end,              5077;
#is $match_on_chr_3_xr->strand,             -1; # FAIL
is $match_on_chr_3_xr->strand,          undef; # See Bio::Location::Split


# Convert the match from contig into chromosomal coordinates

ok my $match_on_chr_4_xf =
  $agp_xf->map( $match_on_ctg_4 );

isa_ok $match_on_chr_4_xf, 'Bio::Coordinate::Result';


ok my $match_on_chr_4_xr =
  $agp_xr->map( $match_on_ctg_4 );

isa_ok $match_on_chr_4_xr, 'Bio::Coordinate::Result';

## Perform some very basic sanity testing on the match objects

is $match_on_chr_4_xf->start,            5925;
isnt $match_on_chr_4_xf->end,            6125; # Gets truncated to maximum!
is $match_on_chr_4_xf->end,              6001; # Gets truncated to maximum!
is $match_on_chr_4_xf->strand,             -1;

isnt $match_on_chr_4_xr->start,          4877; # Gets truncated to minimum!
is $match_on_chr_4_xr->start,            5001; # Gets truncated to minimum!
is $match_on_chr_4_xr->end,              5077;
#is $match_on_chr_4_xr->strand,             +1; # FAIL
is $match_on_chr_4_xr->strand,          undef; # See Bio::Location::Split
