# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN { 
    use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 51);
	
	use_ok('Bio::Structure::Entry');
	use_ok('Bio::Structure::Model');
	use_ok('Bio::Structure::Chain');
	use_ok('Bio::Structure::Residue');
	use_ok('Bio::Structure::Atom');
}

ok my $entry = Bio::Structure::Entry->new;
isa_ok $entry, 'Bio::Structure::Entry';

ok my $model = Bio::Structure::Model->new;
isa_ok $model, 'Bio::Structure::Model';

ok my $chain = Bio::Structure::Chain->new;
isa_ok $chain, 'Bio::Structure::Chain';

ok my $residue = Bio::Structure::Residue->new;
isa_ok $residue, 'Bio::Structure::Residue';

ok my $atom = Bio::Structure::Atom->new;
isa_ok $atom, 'Bio::Structure::Atom';

# adding/removing to Entry
my $m1 = Bio::Structure::Model->new;
$entry->model($m1);
is $entry->model,1;
my $m2 = Bio::Structure::Model->new;
$entry->add_model($m2);
is $entry->get_models, 2;

$entry->model($m2);
is $entry->model, 1;
my @m = ($m1, $m2);
$entry->model(\@m);
is $entry->model, 2;

# does $m2 gest orphaned
$entry->model($m1);
is $entry->parent($m2), undef;


# adding/removing to Model 
my $c1 = Bio::Structure::Chain->new;
$entry->add_chain($model,$c1);
is $entry->get_chains($model),1;
my $c2 = Bio::Structure::Chain->new;
$entry->add_chain($model,$c2);
is $entry->get_chains($model), 2;
isa_ok $entry->parent($c1), 'Bio::Structure::Model';
is $entry->parent($c1), $entry->parent($c2);


# adding/removing to Chain 
my $r1 = Bio::Structure::Residue->new;
$entry->add_residue($chain,$r1);
is $entry->get_residues($chain),1;
my $r2 = Bio::Structure::Residue->new;
$entry->add_residue($chain,$r2);
is $entry->get_residues($chain), 2;

isa_ok $entry->parent($r2), 'Bio::Structure::Chain';
is $entry->parent($r1), $entry->parent($r2);

# adding/removing to Residue 
$entry->add_atom($residue,$atom);
is $entry->get_atoms($residue),1;
my $a2 = Bio::Structure::Atom->new;
my $a3 = Bio::Structure::Atom->new;
my $a4 = Bio::Structure::Atom->new;
my $a5 = Bio::Structure::Atom->new;
my $a6 = Bio::Structure::Atom->new;
$entry->add_atom($residue,$a2);
is $entry->get_atoms($residue), 2;

my @a = ($a3, $a4, $a5);
$entry->add_atom($r2,\@a);
is $entry->get_atoms($r2), 3;

isa_ok $entry->parent($a2), 'Bio::Structure::Residue';
is $entry->parent($a3), $entry->parent($a5);


$atom->x(10.234);
is $atom->x, 10.234;
my $y = 12.345;
$atom->y($y);
is $atom->y, $y;
my $z = $atom->x - $y;
$atom->z($z);
is $atom->z, -2.111;
my @xyz = $atom->xyz;
is (scalar (@xyz), 3);

is $xyz[0], 10.234;
is $xyz[1], 12.345;
is $xyz[2], -2.111;

ok my $e2 = Bio::Structure::Entry->new(-id => "Entry 2");
is $e2->id, "Entry 2";
ok my $m3 = Bio::Structure::Model->new(-id => "Model 2");
is $m3->id, "Model 2";
ok my $c3 = Bio::Structure::Chain->new(-id => "Chain 2");
is $c3->id, "Chain 2";
ok my $r3 = Bio::Structure::Residue->new(-id => "Residue 2");
is $r3->id, "Residue 2";
ok $a2 = Bio::Structure::Atom->new(-id => "Atom 2");
is $a2->id, "Atom 2";

$entry->add_atom($r3,$a6);
$entry->add_residue($c3,$r3);
is $entry->parent( $entry->parent($a6) )->id , "Chain 2";
