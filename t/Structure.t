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
    plan tests => 52;
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

# adding/removing to Entry
my $m1 = Bio::Structure::Model->new;
$entry->model($m1);
ok $entry->model,1;
my $m2 = Bio::Structure::Model->new;
$entry->add_model($m2);
ok $entry->get_models, 2;

$entry->model($m2);
ok $entry->model, 1;
my @m = ($m1, $m2);
$entry->model(\@m);
ok $entry->model, 2;

# does $m2 gest orphaned
$entry->model($m1);
ok $entry->parent($m2), undef;


# adding/removing to Model 
my $c1 = Bio::Structure::Chain->new;
$entry->add_chain($model,$c1);
ok $entry->get_chains($model),1;
my $c2 = Bio::Structure::Chain->new;
$entry->add_chain($model,$c2);
ok $entry->get_chains($model), 2;
ok ref($entry->parent($c1)), 'Bio::Structure::Model';
ok $entry->parent($c1), $entry->parent($c2);


# adding/removing to Chain 
my $r1 = Bio::Structure::Residue->new;
$entry->add_residue($chain,$r1);
ok $entry->get_residues($chain),1;
my $r2 = Bio::Structure::Residue->new;
$entry->add_residue($chain,$r2);
ok $entry->get_residues($chain), 2;

ok ref($entry->parent($r2)), 'Bio::Structure::Chain';
ok $entry->parent($r1), $entry->parent($r2);

# adding/removing to Residue 
$entry->add_atom($residue,$atom);
ok $entry->get_atoms($residue),1;
my $a2 = Bio::Structure::Atom->new;
my $a3 = Bio::Structure::Atom->new;
my $a4 = Bio::Structure::Atom->new;
my $a5 = Bio::Structure::Atom->new;
my $a6 = Bio::Structure::Atom->new;
$entry->add_atom($residue,$a2);
ok $entry->get_atoms($residue), 2;

my @a = ($a3, $a4, $a5);
$entry->add_atom($r2,\@a);
ok $entry->get_atoms($r2), 3;

ok ref($entry->parent($a2)), 'Bio::Structure::Residue';
ok $entry->parent($a3), $entry->parent($a5);



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

$entry->add_atom($r3,$a6);
$entry->add_residue($c3,$r3);
ok $entry->parent( $entry->parent($a6) )->id , "Chain 2";
