# -*-Perl-*- Test Harness script for Bioperl
# $Id: TreeIO.t 14580 2008-03-01 17:01:30Z cjfields $

use strict;

BEGIN {
    use lib '.';
    use Bio::Root::Test;
    use File::Temp qw(tempfile);
    
    test_begin(-tests => 19);
    use_ok('Bio::TreeIO');
}

my $verbose = 0; #test_debug();

my $treeio = Bio::TreeIO->new(
  -format => 'nhx',
  -verbose => $verbose,
  -file   => test_input_file('test.nhx'),
  );
my $tree;
ok($treeio);
$tree = $treeio->next_tree;
isa_ok($tree, 'Bio::Tree::TreeI');

my @nodes = $tree->get_nodes;
is(@nodes, 12, "Total Nodes");
#print STDERR "TREE: ".$tree->as_text('nhx')."\n";

my $adhy = $tree->find_node('ADHY');
is($adhy->branch_length, 0.1);
is(($adhy->get_tag_values('S'))[0], 'nematode');
is(($adhy->get_tag_values('E'))[0], '1.1.1.1');

test_roundtrip('((a,b),c);','simple newick');
test_roundtrip('((x:0.05,y:0.06),a:0.1[&&NHX:G=dummy]);','bug 1471 test');
test_roundtrip('((x:0.05[&&NHX:label=x],y:0.06)[&&NHX:label=int_node],a:0.1[&&NHX:label=a]);','different combinations of label, NHX, and branch length');

test_roundtrip('(a:1,b:2,c:3,d:4)TEST:1.2345;','doot node branch length');
test_roundtrip('(A:0.1,B:0.2,(C:0.3,D:0.4)E:0.5)F;','Example from Wikipedia');

test_roundtrip('(((ADH2:0.1[&&NHX:E=1.1.1.1:S=human],ADH1:0.11[&&NHX:E=1.1.1.1:S=human]):0.05[&&NHX:B=100:D=Y:E=1.1.1.1:S=Primates],ADHY:0.1[&&NHX:E=1.1.1.1:S=nematode],ADHX:0.12[&&NHX:E=1.1.1.1:S=insect]):0.1[&&NHX:D=N:E=1.1.1.1:S=Metazoa],(ADH4:0.09[&&NHX:E=1.1.1.1:S=yeast],ADH3:0.13[&&NHX:E=1.1.1.1:S=yeast],ADH2:0.12[&&NHX:E=1.1.1.1:S=yeast],ADH1:0.11[&&NHX:E=1.1.1.1:S=yeast]):0.1[&&NHX:S=Fungi])[&&NHX:D=N:E=1.1.1.1];','ADH NHX tree');
test_roundtrip('(gene1_Hu[&&NHX:S=Hu_Homo_sapiens],(gene2_Hu[&&NHX:S=Hu_Homo_sapiens],gene2_Mu[&&NHX:S=Mu_Mus_musculus]));','notung nhx example http://www.cs.cmu.edu/~aiton/split/Manual-2.6.master014.html');
test_roundtrip('(cow_gene1,(mouse_gene2,cow_gene2)[&&NHX:B=100]);','notung nhx bootstrap http://www.cs.cmu.edu/~aiton/split/Manual-2.6.master014.html');

# Read in some larger trees from data files...
test_roundtrip(read_file(test_input_file('nhx-bacteria.nhx')),'r-sig-phylo mailing list http://www.mail-archive.com/r-sig-phylo@r-project.org/msg00516.html');
test_roundtrip(read_file(test_input_file('ex1.nucl.nhx')),'treebest example nhx');
# Note: these files aren't reproduced exactly in their online form. We need to round-trip them once
# before including them in the test, because the ordering of annotation keys is not a well-defined
# part of the NHX format. Since nhx.pm sorts the keys before output, once they've been through
# one time, the ordering becomes stable.
test_roundtrip(read_file(test_input_file('wellcome_tol.nhx')),'Wellcome Trust ToL (from http://iphylo.blogspot.com/2009/02/thoughts-on-wellcome-interactive-tree.html)');

# Uncomment to run (takes a long time!!)
#test_roundtrip(read_file(test_input_file('tol-2010-02-18.nhx')),'Tolweb.org converted to NHX');

test_roundtrip(read_file(test_input_file('biorecipe.nhx')),'Biorecipes NHX file (http://www.biorecipes.com/Orthologues/StatusPage/pics/TreeEukaryota.nt)');

sub test_roundtrip {
  my $string = shift;
  my $desc = shift;

  my $in = Bio::TreeIO->new(-format => 'nhx',
                            -string => $string,
                            -verbose => $verbose
                            );
  
  my $t = $in->next_tree;
  my $out;
  if (defined $t) {
    $out = $t->as_text('nhx');
  }

  $desc = "Roundtrip: $desc";
  return is($out,$string,$desc);
}

sub read_file {
  my $file = shift;
  local $/=undef;
  my $string;
  open my $IN, '<', $file or die "Could not read file '$file': $!\n";
  binmode $IN;
  $string = <$IN>;
  close $IN;
  $string =~ s/\n//g;
  $string =~ s/\r//g; # For files with Windows line-endings
  #print STDERR "STR: $string\n";
  return $string;
}
