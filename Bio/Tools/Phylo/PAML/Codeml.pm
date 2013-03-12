#
# BioPerl module for Bio::Tools::Phylo::PAML::Codeml
#
# Cared for by Jason Stajich <jason@bioperl.org>
#
# Copyright Jason Stajich, Aaron J Mackey
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Tools::Phylo::PAML::Codeml - Parses output from the PAML program codeml. 

=head1 SYNOPSIS

  #!/usr/bin/perl -Tw
  use strict;

  use Bio::Tools::Phylo::PAML::Codeml;

  # need to specify the output file name (or a fh) (defaults to
  # -file => "codeml.mlc"); also, optionally, the directory in which
  # the other result files (rst, 2ML.dS, etc) may be found (defaults
  # to "./")
  my $parser = new Bio::Tools::Phylo::PAML::Codeml::Parser
    (-file => "./results/mlc", -dir => "./results/");

  # get the first/next result; a Bio::[...]::Codeml::Result object
  my $result = $parser->next_result();

  # get the sequences used in the analysis; returns Bio::PrimarySeq
  # objects (OTU = Operational Taxonomic Unit).
  my @otus = $result->get_seqs();

  # codon summary: codon usage of each sequence [ arrayref of {
  # hashref of counts for each codon } for each sequence and the
  # overall sum ], and positional nucleotide distribution [ arrayref
  # of { hashref of frequencies for each nucleotide } for each
  # sequence and overall frequencies ].

  my ($codonusage, $ntdist) = $result->get_codon_summary();

  # example manipulations of $codonusage and $ntdist:
  printf "There were %d '%s' codons in the first seq (%s)\n",
    $codonusage->[0]->{AAA}, 'AAA', $otus[0]->id();
  printf "There were %d '%s' codons used in all the sequences\n",
    $codonusage->[$#{$codonusage}]->{AAA}, 'AAA';
  printf "Nucleotide '%c' was present %g of the time in seq %s\n",
    'A', $ntdist->[1]->{A}, $otus[1]->id();

  # get Nei & Gojobori dN/dS matrix:
  my $NGmatrix = $result->get_NGmatrix();

  # get ML-estimated dN/dS matrix, if calculated; this corresponds to
  # the runmode = -2, pairwise comparison usage of codeml
  my $MLmatrix = $result->get_MLmatrix();

  # These matrices are length(@otu) x length(@otu) "strict lower
  # triangle" 2D-matrices, which means that the diagonal and
  # everything above it is undefined.  Each of the defined cells is a
  # hashref of estimates for "dN", "dS", "omega" (dN/dS ratio), "t",
  # "S" and "N".  If a ML matrix, "lnL" will also be defined.  Any
  # additional ML parameters estimated by the model will be in an
  # array ref under "params"; it's up to the user to know which
  # position corresponds to which parameter (since PAML doesn't label
  # them, and we can't guess very well yet (a TODO I guess).

  printf "The omega ratio for sequences %s vs %s was: %g\n",
    $otus[0]->id, $otus[1]->id, $MLmatrix->[0]->[1]->{omega};

  # with a little work, these matrices could also be passed to
  # Bio::Tools::Run::Phylip::Neighbor, or other similar tree-building
  # method that accepts a matrix of "distances" (using the LOWTRI
  # option):
  my $distmat = [ map { [ map { $$_{omega} } @$_ ] } @$MLmatrix ];

  # for runmode's other than -2, get tree topology with estimated
  # branch lengths; returns a Bio::Tree::TreeI-based tree object with
  # added PAML parameters at each node
  my $tree = $result->get_tree();
  for my $node ($tree->get_nodes()) {
     # inspect the tree: the "t" (time) parameter is available via
     # $node->branch_length(); all other branch-specific parameters
     # ("omega", "dN", etc.) are available via $node->param('omega');
  }

  # get any general model parameters: kappa (the
  # transition/transversion ratio), NSsites model parameters ("p0",
  # "p1", "w0", "w1", etc.), etc.
  my $params = $result->get_model_params();
  printf "M1 params: p0 = %g\tp1 = %g\n", $params->{p0}, $params->{p1};

  # for NSsites models, obtain posterior probabilities for membership
  # in each class for every position; probabilities correspond to
  # classes w0, w1, ... etc.
  my @probs = $result->get_posteriors();

  # find, say, positively selected sites!
  if ($params->{w2} > 1) {
    for (my $i = 0; $i < @probs ; $i++) {
      if ($probs[$i]->[2] > 0.5) {
         # assumes model M1: three w's, w0, w1 and w2 (positive selection)
         printf "position %d: (%g prob, %g omega, %g mean w)\n",
           $i, $probs[$i]->[2], $params->{w2}, $probs[$i]->[3];
      }
    }
  } else { print "No positive selection found!\n"; }

=head1 DESCRIPTION

This module is used to parse the output from the PAML program codeml.
You can use the Bio::Tools::Run::Phylo::Phylo::PAML::Codeml module to
actually run codeml; this module is only useful to parse the output.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

  bioperl-l@bioperl.org              - General discussion
  http://bioperl.org/MailList.shtml  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution. Bug reports can be submitted via
email or the web:

  bioperl-bugs@bioperl.org
  http://bioperl.org/bioperl-bugs/

=head1 AUTHOR - Jason Stajich, Aaron Mackey

Email jason@bioperl.org
Email amackey@virginia.edu

=head1 TODO

This module should also be able to handle "codemlsites" batch
output...

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Tools::Phylo::PAML::Codeml;
use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::Root::Root

use Bio::Root::Root;
use Bio::Root::IO;
use Bio::TreeIO;
use IO::String;

@ISA = qw(Bio::Root::Root Bio::Root::IO );

=head2 new

 Title   : new
 Usage   : my $obj = new Bio::Tools::Phylo::PAML::Codeml();
 Function: Builds a new Bio::Tools::Phylo::PAML::Codeml object 
 Returns : Bio::Tools::Phylo::PAML::Codeml
 Args    :


=cut

sub new {
  my($class,@args) = @_;

  my $self = $class->SUPER::new(@args);
  $self->_initialize_io(@args);
  $self->_parse_mlc();
  return $self;
}

=head2 get_trees

 Title   : get_trees
 Usage   : my @trees = $codemlparser->get_trees();
 Function: Returns a list of trees (if any) are in the output file
 Returns : List of L<Bio::Tree::TreeI> objects
 Args    : none


=cut

sub get_trees{
   my ($self) = @_;


}

=head2 get_statistics

 Title   : get_statistics
 Usage   : my $data = $codemlparser->get_statistics
 Function: Retrieves the set of pairwise comparisons 
 Returns : Hash Reference keyed as 'seqname' -> 'seqname' -> 'datatype'
 Args    : none


=cut

sub get_statistics {
   my ($self) = @_;
   

}


# parse the mlc file

sub _parse_mlc {
    my ($self) = @_;
    my %data;
    while( defined( $_ = $self->_readline) ) {
	print;
	# Aaron this is where the parsing should begin

	# I'll do the Tree objects if you like - 
	# I'd do it by building an IO::String for the
	# the tree data 
	# or does it make more sense to parse this out of a collection of 
	# files?
	if( /^TREE/ ) {
	    # ...
	    while( defined($_ = $self->_readline) ) {
		if( /^\(/) {
		    my $treestr = IO::String->new($_);
		    my $treeio = Bio::TreeIO->new(-fh => $treestr,
						 -format => 'newick');
		    # this is very tenative here!!
		    push @{$self->{'_trees'}}, $treeio->next_tree;
		}
	    }
	}
    }
}

1;
