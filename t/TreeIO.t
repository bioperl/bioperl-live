# -*-Perl-*-
## Bioperl Test Harness Script for Modules
## $Id$

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

my $error = 0;

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
    plan tests => 10; 

#    eval { require XML::Parser::PerlSAX; };
#    if( $@ ) {
#	print STDERR "XML::Parser::PerlSAX not loaded. This means TreeIO::phyloxml test cannot be executed. Skipping\n";
#	foreach ( 1..43 ) {
#	    skip(1,1);
#	}
#       $error = 1;
#	
#    } 

}

if( $error == 1 ) {
    exit(0);
}

use vars qw($FILE1 $FILE2);

$FILE1= 'testnewick.phylip';
$FILE2= 'testlarge.phy';

END { 
	unlink $FILE1;
	unlink $FILE2;
}
use Bio::TreeIO;
use Bio::Root::IO;
my $verbose = $ENV{'BIOPERLDEBUG'};
ok(1);

my $treeio = new Bio::TreeIO(-verbose => $verbose,
			     -format => 'newick',
			     -file   => Bio::Root::IO->catfile('t','data', 
							       'cysprot1b.newick'));

ok($treeio);
my $tree = $treeio->next_tree;

ok(ref($tree) && $tree->isa('Bio::Tree::TreeI'));

my @nodes = $tree->get_nodes;
ok(@nodes, 6);

if($verbose ) { 
    foreach my $node ( $tree->get_root_node()->each_Descendent() ) {
	print "node: ", $node->to_string(), "\n";
	my @ch = $node->each_Descendent();
	if( @ch ) {
	    print "\tchildren are: \n";
	    foreach my $node ( $node->each_Descendent() ) {
		print "\t\t ", $node->to_string(), "\n";
	    }
	}
    }
}
$treeio = new Bio::TreeIO(-verbose => $verbose,
			  -format => 'newick', 
			  -file   => ">$FILE1");
$treeio->write_tree($tree);
undef $treeio;
ok( -s $FILE1 );
$treeio = new Bio::TreeIO(-verbose => $verbose,
			  -format => 'newick',
			  -file   => Bio::Root::IO->catfile('t','data', 
							    'LOAD_Ccd1.dnd'));
ok($treeio);
$tree = $treeio->next_tree;

ok(ref($tree) && $tree->isa('Bio::Tree::TreeI'));

@nodes = $tree->get_nodes;
ok(@nodes, 52);

if($verbose ) { 
    foreach my $node ( @nodes ) {    
	print "node: ", $node->to_string(), "\n";
	my @ch = $node->each_Descendent();
	if( @ch ) { 
	    print "\tchildren are: \n";
	    foreach my $node ( $node->each_Descendent() ) {
		print "\t\t ", $node->to_string(), "\n";
	    }
	}
    }
}

ok($tree->total_branch_length, 7.12148);

$treeio = new Bio::TreeIO(-verbose => $verbose,
			  -format => 'newick', 
			  -file   => ">$FILE2");
$treeio->write_tree($tree);
undef $treeio;
ok(-s $FILE2);
