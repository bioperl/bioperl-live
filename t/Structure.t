# -*-Perl-*-
## Bioperl Test Harness Script for Modules
## $Id$

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
    plan tests => 69;
}
use Bio::Structure::Entry;
use Bio::Structure::Model;
use Bio::Structure::Chain;
use Bio::Structure::Residue;
use Bio::Structure::Atom;

ok(1);

my $entry = Bio::Structure::Entry->new;
ok(1);
ok defined  $entry;
ok ref($entry), 'Bio::Structure::Entry';

my $model = Bio::Structure::Model->new;
ok(1);
ok defined $model;
ok ref($model), 'Bio::Structure::Model';

my $chain = Bio::Structure::Chain->new;
ok(1);
ok defined $chain;
ok ref($chain), 'Bio::Structure::Chain';

my $residue = Bio::Structure::Residue->new;
ok(1);
ok defined $residue;
ok ref($residue), 'Bio::Structure::Residue';

my $atom = Bio::Structure::Atom->new;
ok(1);
ok defined $atom;
ok ref($atom), 'Bio::Structure::Atom';

# does entry has a Model and a Chain
ok $entry->model, 1;
my ($m) = $entry->model; 
ok $m->chain, 1;

# does model has a Chain
ok $model->chain, 1;

# adding/removing to Entry
$entry->model($model);
ok $entry->model,1;
my $m2 = Bio::Structure::Model->new;
$entry->add_model($m2);
ok $entry->model, 2;

$entry->model($m2);
ok $entry->model, 1;
my @m = ($model, $m2);
$entry->model(\@m);
ok $entry->model, 2;

ok ref($model->entry), 'Bio::Structure::Entry';
ok ref($m2->entry), 'Bio::Structure::Entry';
ok $model->entry, $m2->entry;

# does $m2 gest orphaned
$entry->model($model);
ok $m2->entry, undef;


# adding/removing to Model 
$model->chain($chain);
ok $model->chain,1;
my $c2 = Bio::Structure::Chain->new;
$model->add_chain($c2);
ok $model->chain, 2;

$model->chain($c2);
ok $model->chain, 1;
my @c = ($chain, $c2);
$model->chain(\@c);
ok $model->chain, 2;

ok ref($chain->model), 'Bio::Structure::Model';
ok ref($c2->model), 'Bio::Structure::Model';
ok $chain->model, $c2->model;

# does $c2 gest orphaned
$model->chain($chain);
ok $c2->model, undef;


# adding/removing to Chain 
$chain->residue($residue);
ok $chain->residue,1;
my $r2 = Bio::Structure::Residue->new;
$chain->add_residue($r2);
ok $chain->residue, 2;

$chain->residue($r2);
ok $chain->residue, 1;
my @r = ($residue, $r2);
$chain->residue(\@r);
ok $chain->residue, 2;

ok ref($residue->chain), 'Bio::Structure::Chain';
ok ref($r2->chain), 'Bio::Structure::Chain';
ok $residue->chain, $r2->chain;

# does $r2 gest orphaned
$chain->residue($residue);
ok $r2->chain, undef;


# adding/removing to Residue 
$residue->atom($atom);
ok $residue->atom,1;
my $a2 = Bio::Structure::Atom->new;
$residue->add_atom($a2);
ok $residue->atom, 2;

$residue->atom($a2);
ok $residue->atom, 1;
my @a = ($atom, $a2);
$residue->atom(\@a);
ok $residue->atom, 2;

ok ref($atom->residue), 'Bio::Structure::Residue';
ok ref($a2->residue), 'Bio::Structure::Residue';
ok $atom->residue, $a2->residue;

# does $a2 gest orphaned
$residue->atom($atom);
ok $a2->residue, undef;


$atom->x(10.234);
ok $atom->x, 10.234;
my $y = 12.345;
$atom->y($y);
ok $atom->y, $y;
my $z = $atom->x - $y;
$atom->z($z);
ok $atom->z, -2.111;
ok ($atom->xyz), 3;
my @xyz = $atom->xyz;
ok $xyz[0], 10.234;
ok $xyz[1], 12.345;
ok $xyz[2], -2.111;

my $e2 = Bio::Structure::Entry->new(-id => "Entry 2");
ok (1);
ok $e2->id, "Entry 2";
my $m3 = Bio::Structure::Model->new(-id => "Model 2");
ok (1);
ok $m3->id, "Model 2";
my $c3 = Bio::Structure::Chain->new(-id => "Chain 2");
ok (1);
ok $c3->id, "Chain 2";
my $r3 = Bio::Structure::Residue->new(-id => "Residue 2");
ok (1);
ok $r3->id, "Residue 2";
my $a2 = Bio::Structure::Atom->new(-id => "Atom 2");
ok (1);
ok $a2->id, "Atom 2";



$r3->add_atom($atom);
$r3->add_atom($a2);
$c3->add_residue($r3);
ok $atom->residue->chain->id, "Chain 2";
