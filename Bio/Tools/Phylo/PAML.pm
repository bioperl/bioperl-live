# $Id$
#
# BioPerl module for Bio::Tools::Phylo::PAML
#
# Cared for by Jason Stajich <jason@bioperl.org>
#
# Copyright Jason Stajich, Aaron J Mackey
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Tools::Phylo::PAML - Parses output from the PAML programs codeml,
baseml, basemlg, codemlsites and yn00

=head1 SYNOPSIS

  #!/usr/bin/perl -Tw
  use strict;

  use Bio::Tools::Phylo::PAML;

  # need to specify the output file name (or a fh) (defaults to
  # -file => "codeml.mlc"); also, optionally, the directory in which
  # the other result files (rst, 2ML.dS, etc) may be found (defaults
  # to "./")
  my $parser = new Bio::Tools::Phylo::PAML
    (-file => "./results/mlc", -dir => "./results/");

  # get the first/next result; a Bio::Tools::Phylo::PAML::Result object,
  # which isa Bio::SeqAnalysisResultI object.
  my $result = $parser->next_result();

  # get the sequences used in the analysis; returns Bio::PrimarySeq
  # objects (OTU = Operational Taxonomic Unit).
  my @otus = $result->get_seqs();

  # codon summary: codon usage of each sequence [ arrayref of {
  # hashref of counts for each codon } for each sequence and the
  # overall sum ], and positional nucleotide distribution [ arrayref
  # of { hashref of frequencies for each nucleotide } for each
  # sequence and overall frequencies ]:
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
  # "S" and "N".  If a ML matrix, "lnL" will also be defined.
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

  # for NSsites models, obtain arrayrefs of posterior probabilities
  # for membership in each class for every position; probabilities
  # correspond to classes w0, w1, ... etc.
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

This module is used to parse the output from the PAML programs codeml,
baseml, basemlg, codemlsites and yn00.  You can use the
Bio::Tools::Run::Phylo::PAML::* modules to actually run some of the
PAML programs, but this module is only useful to parse the output.

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

check output from pre 1.12

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Tools::Phylo::PAML;
use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::Root::Root

use Bio::Root::Root;
use Bio::SeqAnalysisParserI;
use Bio::Root::IO;
@ISA = qw(Bio::Root::Root Bio::Root::IO Bio::AnalysisParserI);

# other objects used:

use Bio::TreeIO;
use IO::String;
use Bio::Tools::Phylo::PAML::Result;

=head2 new

 Title   : new
 Usage   : my $obj = new Bio::Tools::Phylo::PAML(%args);
 Function: Builds a new Bio::Tools::Phylo::PAML object
 Returns : Bio::Tools::Phylo::PAML
 Args    : Hash of options: -file, -fh, -dir
           -file (or -fh) should contain the contents of the PAML
           outfile; -dir is the (optional) name of the directory in
           which the PAML program was run (and includes other
           PAML-generated files from which we can try to gather data)

=cut

sub new {

  my ($class, @args) = @_;

  my $self = $class->SUPER::new(@args);

  my ($dir) = $self->_rearrange([qw(DIR)], @args);
  $self->{_dir} = $dir if defined $dir;

  return $self;
}

=head2 Implement Bio::AnalysisParserI interface

=cut

=head2 next_result

 Title   : next_result
 Usage   : $result = $obj->next_result();
 Function: Returns the next result available from the input, or
           undef if there are no more results.
 Example :
 Returns : a Bio::Tools::Phylo::PAML::Result object
 Args    : none

=cut

sub next_result {

    my ($self) = @_;

    my %data;

    # get the various codon and other sequence summary data, if necessary:
    $self->_parse_summary
	unless ($self->{_summary} && !$self->{_summary}->{multidata});


    # OK, depending on seqtype and runmode now, one of a few things can happen:
    my $seqtype = $self->{_summary}->{seqtype};
    if ($seqtype eq 'CODONML' || $seqtype eq 'AAML') {
	while ($_ = $self->_readline) {

	    if ($seqtype eq 'CODONML' && m/^pairwise comparison, codon frequencies:/o) {

		# runmode = -2, CODONML
		$self->_pushback($_);
		%data = $self->_parse_PairwiseCodon;
		last;

	    } elsif ($seqtype eq 'AAML' && m/^ML distances of aa seqs\.$/o) {

		# runmode = -2, AAML
		$self->throw( -class => 'Bio::Root::NotYetImplemented',
			      -text  => "Pairwise AA not yet implemented!"
			    );

		# $self->_pushback($_);
		# %data = $self->_parse_PairwiseAA;
		# last;

	    } elsif (m/^Model \d+: /o) {

		# NSSitesBatch
		$self->throw( -class => 'Bio::Root::NotYetImplemented',
			      -text  => "NSsitesBatch not yet implemented!"
			    );

		# $self->_pushback($_);
		# %data = $self->_parse_NSsitesBatch;
		# last;

	    } elsif (m/TREE \d+/) {

		# runmode = 0
		$self->_pushback($_);
		%data = $self->_parse_Forestry;

	    } elsif (m/Heuristic tree search by stepwise addition$/o) {

		# runmode = 3
		$self->throw( -class => 'Bio::Root::NotYetImplemented',
			      -text  => "StepwiseAddition not yet implemented!"
			    );

		# $self->_pushback($_);
		# %data = $self->_parse_StepwiseAddition;

	    } elsif (m/Heuristic tree search by NNI perturbation$/o) {

		# runmode = 4
		$self->throw( -class => 'Bio::Root::NotYetImplemented',
			      -text  => "NNI Perturbation not yet implemented!"
			    );

		# $self->_pushback($_);
		# %data = $self->_parse_Perturbation;

	    } elsif (m/^stage 0:/o) {

		# runmode = (1 or 2)
		$self->throw( -class => 'Bio::Root::NotYetImplemented',
			      -text  => "StarDecomposition not yet implemented!"
			    );

		# $self->_pushback($_);
		# %data = $self->_parse_StarDecomposition;

	    }
	}
    } elsif ($seqtype eq 'BASEML') {
    } elsif ($seqtype eq 'YN00') {
    }


    if (%data) {
	return new Bio::Tools::Phylo::PAML::Result %data;
    } else {
	return undef;
    }
}


sub _parse_summary {

    my ($self) = @_;

    # Depending on whether verbose > 0 or not, and whether the result
    # set comes from a multi-data run, the first few lines could be
    # various things; we're going to throw away any sequence data
    # here, since we'll get it later anyways

    # multidata ? : \n\nData set 1\n
    # verbose ? : cleandata ? : \nBefore deleting alignment gaps. \d sites\n
    #                           [ sequence printout ]
    #                           \nAfter deleting gaps. \d sites\n"
    #           : [ sequence printout ]
    # CODONML (in paml 3.12 February 2002)  <<-- what we want to see!

    my $SEQTYPES = qr( (?: (?: CODON | AA | BASE | CODON2AA ) ML ) | YN00 )x;
    while ($_ = $self->_readline) {
	if ( m/^($SEQTYPES) \s+                      # seqtype: CODONML, AAML, BASEML, CODON2AAML, YN00, etc
	       (?: \(in \s+ ([^\)]+?) \s* \) \s* )?  # version: "paml 3.12 February 2002"; not present < 3.1 or YN00
	       (\S+) \s*                             # tree filename
	       (?: (.+?) )?                          # model description (not there in YN00)
	       \s* $                                 # trim any trailing space
	       /ox
	   ) {

	    @{$self->{_summary}}{qw(seqtype version treefile model)} = ($1, $2, $3, $4);
	    last;

	} elsif (m/^Data set \d$/o) {
	    $self->{_summary} = {};
	    $self->{_summary}->{multidata}++;
	}
    }

    unless (defined $self->{_summary}->{seqtype}) {
	$self->throw( -class => 'Bio::Root::NotImplemented',
		      -text => 'Unknown format of PAML output');
    }


    my $seqtype = $self->{_summary}->{seqtype};

    if ($seqtype == "CODEML") {

	$self->_parse_inputparams(); # settings from the .ctl file that get printed
	$self->_parse_patterns();    # codon patterns - not very interesting
	$self->_parse_seqs();        # the sequences data used for analysis
	$self->_parse_codoncts();    # counts and distributions of codon/nt usage
	$self->_parse_distmat();     # NG distance matrices

    } elsif ($seqtype == "AAML") {
	$self->throw( -class => 'Bio::Root::NotImplemented',
		      -text => 'AAML parsing not yet implemented!');
    } elsif ($seqtype == "CODON2AAML") {
	$self->throw( -class => 'Bio::Root::NotImplemented',
		      -text => 'CODON2AAML parsing not yet implemented!');
    } elsif ($seqtype == "BASEML") {
	$self->throw( -class => 'Bio::Root::NotImplemented',
		      -text => 'BASEML parsing not yet implemented!');
    } elsif ($seqtype == "YN00") {
	$self->throw( -class => 'Bio::Root::NotImplemented',
		      -text => 'YN00 parsing not yet implemented!');
    } else {
	$self->throw( -class => 'Bio::Root::NotImplemented',
		      -text => 'Unknown seqtype, not yet implemented!',
		      -value => $seqtype
		    );
    }

}


sub _parse_inputparams { }
sub _parse_patterns { }
sub _parse_seqs { }
sub _parse_codoncts { }
sub _parse_distmat { }


sub _parse_Forestry {

    my ($self) = @_;
    my %data;


    return %data
};

# parse the mlc file

sub _parse_mlc {
    my ($self) = @_;
    my %data;
    while( defined( $_ = $self->_readline) ) {
	print;
	# Aaron this is where the parsing should begin

	# I'll do the Tree objects if you like - I'd do it by building
	# an IO::String for the the tree data or does it make more
	# sense to parse this out of a collection of files?
	if( /^TREE/ ) {
	    # ...
	    while( defined($_ = $self->_readline) ) {
		if( /^\(/) {
		    my $treestr = new IO::String($_);
		    my $treeio = new Bio::TreeIO(-fh => $treestr,
						 -format => 'newick');
		    # this is very tenative here!!
		    push @{$self->{'_trees'}}, $treeio->next_tree;
		}
	    }
	}
    }
}

1;
