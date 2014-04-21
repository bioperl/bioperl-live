#
# BioPerl module for Bio::Tools::Phylo::PAML
#
# Please direct questions and support issues to <bioperl-l@bioperl.org>
#
# Cared for by Jason Stajich <jason-at-bioperl.org>
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
  my $parser = Bio::Tools::Phylo::PAML->new
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
  printf "There were %d %s codons in the first seq (%s)\n",
    $codonusage->[0]->{AAA}, 'AAA', $otus[0]->id();
  printf "There were %d %s codons used in all the sequences\n",
    $codonusage->[$#{$codonusage}]->{AAA}, 'AAA';
  printf "Nucleotide %c was present %g of the time in seq %s\n",
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
  # "S" and "N".  If a ML matrix, "lnL" and "kappa" will also be defined.
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
  my ($tree) = $result->get_trees();
  for my $node ($tree->get_nodes()) {
     # inspect the tree: the "t" (time) parameter is available via
     # $node->branch_length(); all other branch-specific parameters
     # ("omega", "dN", etc.) are available via
     # ($omega) = $node->get_tag_values('omega');
  }

  # if you are using model based Codeml then trees are stored in each
  # modelresult object
  for my $modelresult ( $result->get_NSSite_results ) {
    # model M0, M1, etc
    print "model is ", $modelresult->model_num, "\n";
    my ($tree) = $modelresult->get_trees();
    for my $node ($tree->get_nodes()) {
     # inspect the tree: the "t" (time) parameter is available via
     # $node->branch_length(); all other branch-specific parameters
     # ("omega", "dN", etc.) are available via
     # ($omega) = $node->get_tag_values('omega');
   }
  }

  # get any general model parameters: kappa (the
  # transition/transversion ratio), NSsites model parameters ("p0",
  # "p1", "w0", "w1", etc.), etc.
  my $params = $result->get_model_params();
  printf "M1 params: p0 = %g\tp1 = %g\n", $params->{p0}, $params->{p1};

  # parse AAML result files
  my $aamat = $result->get_AADistMatrix();
  my $aaMLmat = $result->get_AAMLDistMatrix();

=head1 DESCRIPTION

This module is used to parse the output from the PAML programs codeml,
baseml, basemlg, codemlsites and yn00.  You can use the
Bio::Tools::Run::Phylo::PAML::* modules to actually run some of the
PAML programs, but this module is only useful to parse the output.

This module has fledgling support for PAML version 4.3a (October 2009).
Please report any problems to the mailing list (see FEEDBACK below).

=head1 TO DO

Implement get_posteriors(). For NSsites models, obtain arrayrefs of
posterior probabilities for membership in each class for every
position; probabilities correspond to classes w0, w1, ... etc.

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

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

  bioperl-l@bioperl.org                  - General discussion
  http://bioperl.org/wiki/Mailing_lists  - About the mailing lists

=head2 Support

Please direct usage questions or support issues to the mailing list:

I<bioperl-l@bioperl.org>

rather than to the module maintainer directly. Many experienced and
reponsive experts will be able look at the problem and quickly
address it. Please include a thorough description of the problem
with code and data examples if at all possible.

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution. Bug reports can be submitted via the
web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Jason Stajich, Aaron Mackey

Email jason-at-bioperl.org
Email amackey-at-virginia.edu

=head1 CONTRIBUTORS

Albert Vilella avilella-AT-gmail-DOT-com
Sendu Bala     bix@sendu.me.uk
Dave Messina   dmessina@cpan.org

=head1 TODO

RST parsing -- done, Avilella contributions bug#1506, added by jason 1.29
            -- still need to parse in joint probability and non-syn changes
               at site table

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::Tools::Phylo::PAML;
use vars qw($RSTFILENAME);
use strict;

# Object preamble - inherits from Bio::Root::Root

use base qw(Bio::Root::Root Bio::Root::IO Bio::AnalysisParserI);

BEGIN {
    $RSTFILENAME = 'rst';    # where to get the RST data from
}

# other objects used:
use IO::String;
use File::Spec;
use Bio::TreeIO;
use Bio::Tools::Phylo::PAML::Result;
use Bio::LocatableSeq;
use Bio::PrimarySeq;
use Bio::Matrix::PhylipDist;
use Bio::Tools::Phylo::PAML::ModelResult;

=head2 new

 Title   : new
 Usage   : my $obj = Bio::Tools::Phylo::PAML->new(%args);
 Function: Builds a new Bio::Tools::Phylo::PAML object
 Returns : Bio::Tools::Phylo::PAML
 Args    : Hash of options: -file, -fh, -dir
           -file (or -fh) should contain the contents of the PAML
                 outfile;
           -dir is the (optional) name of the directory in
                which the PAML program was run (and includes other
                PAML-generated files from which we can try to gather data)

=cut

sub new {

    my ( $class, @args ) = @_;

    my $self = $class->SUPER::new(@args);
    $self->_initialize_io(@args);
    my ($dir) = $self->_rearrange( [qw(DIR)], @args );
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

    # parse the RST file, if it doesn't exist or if dir is not set
    # this will just skip the parsing
    $self->_parse_rst();
    my $idlookup;    # a hashreference to SEQID (number) ==> 'SEQUENCENAME'
        # get the various codon and other sequence summary data, if necessary:
    $self->_parse_summary
      unless ( $self->{'_summary'} && !$self->{'_summary'}->{'multidata'} );

    # OK, depending on seqtype and runmode now, one of a few things can happen:
    my $seqtype = $self->{'_summary'}->{'seqtype'};
    if ( $seqtype eq 'CODONML' || $seqtype eq 'AAML' ) {
        my $has_model_line = 0;
        while ( defined( $_ = $self->_readline ) ) {
            if ( $seqtype eq 'CODONML'
                && m/^pairwise comparison, codon frequencies:/ )
            {

                # runmode = -2, CODONML
                $self->debug("pairwise Ka/Ks\n");
                $self->_pushback($_);
                %data = $self->_parse_PairwiseCodon;
                last;
            }
            elsif ( $seqtype eq 'AAML' && m/^ML distances of aa seqs\.$/ ) {
                $self->_pushback($_);

                # get AA distances
                %data = ( '-AAMLdistmat' => $self->_parse_aa_dists() );

                # $self->_pushback($_);
                # %data = $self->_parse_PairwiseAA;
                # last;
            }
            elsif (
                m/^Model\s+(\d+)/
                || ( ( !$has_model_line && m/^TREE/ )
                    && $seqtype eq 'CODONML'
                    && ($self->{'_summary'}->{'version'} !~ /4/))
                    # last bit to keep PAML >= 4 from being caught here
                    # bug 2482. Not sure this is the right fix, but tests
                    # pass and the bug's test case passes.
              )
            {
                $self->_pushback($_);
                my $model = $self->_parse_NSsitesBatch;
                push @{ $data{'-NSsitesresults'} }, $model;
                $has_model_line = 1;
            }
            elsif (m/for each branch/) {
                my %branch_dnds = $self->_parse_branch_dnds;
                if ( !defined $data{'-trees'} ) {
                    $self->warn(
                        "No trees have been loaded, can't do anything\n");
                    next;
                }
                my ($tree) = @{ $data{'-trees'} };
                if (   !$tree
                    || !ref($tree)
                    || !$tree->isa('Bio::Tree::Tree') )
                {
                    $self->warn("no tree object already stored!\n");
                    next;
                }

                # These need to be added to the Node/branches
                while ( my ( $k, $v ) = each %branch_dnds ) {

                    # we can probably do better by caching at some point
                    my @nodes;
                    for my $id ( split( /\.\./, $k ) ) {
                        my @nodes_L =
                          map { $tree->find_node( -id => $_ ) }
                          @{ $idlookup->{$id} };
                        my $n =
                          @nodes_L < 2
                          ? shift(@nodes_L)
                          : $tree->get_lca(@nodes_L);
                        if ( !$n ) {
                            $self->warn("no node for $n\n");
                        }
                        unless ( $n->is_Leaf && $n->id ) {
                            $n->id($id);
                        }
                        push @nodes, $n;
                    }
                    my ( $parent, $child ) = @nodes;
                    while ( my ( $kk, $vv ) = each %$v ) {
                        $child->add_tag_value( $kk, $vv );
                    }
                }
            }
            elsif (m/^TREE/) {

                # runmode = 0
                $self->_pushback($_);
                ( $data{'-trees'}, $idlookup ) = $self->_parse_Forestry;

                #last;
            }
            elsif (m/Heuristic tree search by stepwise addition$/) {

                # runmode = 3
                $self->throw(
                    -class => 'Bio::Root::NotImplemented',
                    -text  => "StepwiseAddition not yet implemented!"
                );

                # $self->_pushback($_);
                # %data = $self->_parse_StepwiseAddition;
                # last;

            }
            elsif (m/Heuristic tree search by NNI perturbation$/) {

                # runmode = 4
                $self->throw(
                    -class => 'Bio::Root::NotImplemented',
                    -text  => "NNI Perturbation not yet implemented!"
                );

                # $self->_pushback($_);
                # %data = $self->_parse_Perturbation;
                # last;

            }
            elsif (m/^stage 0:/) {

                # runmode = (1 or 2)
                $self->throw(
                    -class => 'Bio::Root::NotImplemented',
                    -text  => "StarDecomposition not yet implemented!"
                );

                $self->_pushback($_);
                %data = $self->_parse_StarDecomposition;
                last;
            }
        }
    }
    elsif ( $seqtype eq 'BASEML' ) {
        while ( defined( $_ = $self->_readline ) ) {
            if (/^Distances:/) {
                $self->_pushback($_);
                my ( $kappa, $alpha ) = $self->_parse_nt_dists();
                %data = (
                    '-kappa_distmat' => $kappa,
                    '-alpha_distmat' => $alpha
                );
            }
            elsif (/^TREE/) {
                $self->_pushback($_);
                ( $data{'-trees'}, $idlookup ) = $self->_parse_Forestry;
            }
        }
    }
    elsif ( $seqtype eq 'YN00' ) {
        while ( $_ = $self->_readline ) {
            if (
m/^Estimation by the method|\(B\) Yang & Nielsen \(2000\) method/
              )
            {
                $self->_pushback($_);
                %data = $self->_parse_YN_Pairwise;
                last;
            }
        }
    }
    if (%data) {
        $data{'-version'}      = $self->{'_summary'}->{'version'};
        $data{'-seqs'}         = $self->{'_summary'}->{'seqs'};
        $data{'-patterns'}     = $self->{'_summary'}->{'patterns'};
        $data{'-ngmatrix'}     = $self->{'_summary'}->{'ngmatrix'};
        $data{'-codonpos'}     = $self->{'_summary'}->{'codonposition'};
        $data{'-codonfreq'}    = $self->{'_summary'}->{'codonfreqs'};
        $data{'-model'}        = $self->{'_summary'}->{'model'};
        $data{'-seqfile'}      = $self->{'_summary'}->{'seqfile'};
        $data{'-aadistmat'}    = $self->{'_summary'}->{'aadistmat'};
        $data{'-stats'}        = $self->{'_summary'}->{'stats'};
        $data{'-aafreq'}       = $self->{'_summary'}->{'aafreqs'};
        $data{'-ntfreq'}       = $self->{'_summary'}->{'ntfreqs'};
        $data{'-input_params'} = $self->{'_summary'}->{'inputparams'};
        $data{'-rst'}          = $self->{'_rst'}->{'rctrted_seqs'};
        $data{'-rst_persite'}  = $self->{'_rst'}->{'persite'};
        $data{'-rst_trees'}    = $self->{'_rst'}->{'trees'};
        return Bio::Tools::Phylo::PAML::Result->new(%data);
    }
    else {
        return;
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
    my $line;
    $self->{'_already_parsed_seqs'} = $self->{'_already_parsed_seqs'} ? 1 : 0;
    my @lines;
    while ( $_ = $self->_readline ) {
        push @lines, $_;
        if (m/^($SEQTYPES) \s+                 # seqtype: CODONML, AAML, BASEML, CODON2AAML, YN00, etc
	     (?: \(in \s+ ([^\)]+?) \s* \) \s* )?  # version: "paml 3.12 February 2002"; not present < 3.1 or YN00
	     (\S+) \s*                             # tree filename
	     (?: (.+?) )?                          # model description (not there in YN00)
	     \s* $                                 # trim any trailing space
	    /ox
          )
        {
            @{ $self->{'_summary'} }{qw(seqtype version seqfile model)} =
              ( $1, $2, $3, $4 );

            # in 4.3, the model is on its own line
            if ( !defined $self->{'_summary'}->{'model'} ) {
                my $model_line = $self->_readline;
                chomp $model_line;
                if ($model_line =~ /^Model:/) {
                    $self->{'_summary'}->{'model'} = $model_line;
                }
            }

            defined $self->{'_summary'}->{'model'}
              && $self->{'_summary'}->{'model'} =~ s/Model:\s+//;
            $self->_pushback($_)
              if $self->{'_summary'}->{'seqtype'} eq 'AAMODEL';
            last;
        }
        elsif ((m/\s+?\d+?\s+?\d+?/) && ( $self->{'_already_parsed_seqs'} != 1 )) {
			$self->_parse_seqs;
		}
        elsif (m/^Data set \d$/) {
            $self->{'_summary'} = {};
            $self->{'_summary'}->{'multidata'}++;
        }
        elsif (m/^Before\s+deleting\s+alignment\s+gaps/) {    #Gap
            my ($phylip_header) = $self->_readline;
            $self->_parse_seqs;
        }
        elsif ( ( @lines >= 3 ) && ( $self->{'_already_parsed_seqs'} != 1 ) )
        {                                                     #No gap
            # don't start parsing seqs yet if we're on a blank line
            # (gives another opportunity to match one of the other regexes)
            unless (/^\n$/) {
                $self->_parse_seqs;
            }
        }
        elsif ( (/Printing out site pattern counts/)
            && ( $self->{'_already_parsed_seqs'} != 1 ) ) {
            $self->_parse_patterns;
        }
    }

    unless ( defined $self->{'_summary'}->{'seqtype'} ) {
        $self->throw(
            -class => 'Bio::Root::NotImplemented',
            -text  => 'Unknown format of PAML output did not see seqtype'
        );
    }
    my $seqtype = $self->{'_summary'}->{'seqtype'};
    if ( $seqtype eq "CODONML" ) {
        $self->_parse_inputparams();    # settings from the .ctl file
                                        # that get printed
        $self->_parse_patterns();       # codon patterns - not very interesting
        $self->_parse_seqs();           # the sequences data used for analysis
        $self->_parse_codoncts();       # counts and distributions of codon/nt
                                        # usage
        $self->_parse_codon_freqs();    # codon frequencies
        $self->_parse_distmat();        # NG distance matrices
    }
    elsif ( $seqtype eq "AAML" ) {
        $self->_parse_inputparams;
        $self->_parse_patterns();
        $self->_parse_seqs();           # the sequences data used for analysis
        $self->_parse_aa_freqs();       # AA frequencies
                                        # get AA distances
        $self->{'_summary'}->{'aadistmat'} = $self->_parse_aa_dists();

    }
    elsif ( $seqtype eq "CODON2AAML" ) {
        $self->throw(
            -class => 'Bio::Root::NotImplemented',
            -text  => 'CODON2AAML parsing not yet implemented!'
        );
    }
    elsif ( $seqtype eq "BASEML" ) {
        $self->_parse_patterns();
        $self->_parse_seqs();
        $self->_parse_nt_freqs();

    }
    elsif ( $seqtype eq "YN00" ) {
        $self->_parse_codon_freqs();
        $self->_parse_codoncts();
        $self->_parse_distmat();    # NG distance matrices
    }
    else {
        $self->throw(
            -class => 'Bio::Root::NotImplemented',
            -text  => 'Unknown seqtype, not yet implemented!',
            -value => $seqtype
        );
    }

}

sub _parse_inputparams {
    my ($self) = @_;
    while ( defined( $_ = $self->_readline ) ) {
        if (/^((?:Codon frequencies)|(?:Site-class models))\s*:\s+(.+)/) {
            my ( $param, $val ) = ( $1, $2 );
            $self->{'_summary'}->{'inputparams'}->{$param} = $val;
        }
        elsif (/^\s+$/) {
            next;
        }
        elsif ( /^ns\s+=\s+/ || /^Frequencies/ ) {
            $self->_pushback($_);
            last;
        }
    }
}

sub _parse_codon_freqs {
    my ($self) = @_;
    my ( $okay, $done ) = ( 0, 0 );

    while ( defined( $_ = $self->_readline ) ) {
        if (/^Nei|\(A\) Nei/) { $self->_pushback($_); last }
        last if ($done);
        next if (/^\s+/);
        next
          unless ( $okay || /^Codon position x base \(3x4\) table\, overall/ );
        $okay = 1;
        if (s/^position\s+(\d+):\s+//) {
            my $pos = $1;
            s/\s+$//;
            my @bases = split;
            foreach my $str (@bases) {
                my ( $base, $freq ) = split( /:/, $str, 2 );
                $self->{'_summary'}->{'codonposition'}->[ $pos - 1 ]->{$base} =
                  $freq;
            }
            $done = 1 if $pos == 3;
        }
    }
    $done = 0;
    while ( defined( $_ = $self->_readline ) ) {
        if (/^Nei\s\&\sGojobori|\(A\)\sNei-Gojobori/) {
            $self->_pushback($_);
            last;
        }
        last if ($done);
        if (/^Codon frequencies under model, for use in evolver/) {
            while ( defined( $_ = $self->_readline ) ) {
                last if (/^\s+$/);
                s/^\s+//;
                s/\s+$//;
                push @{ $self->{'_summary'}->{'codonfreqs'} }, [split];
            }
            $done = 1;
        }
    }
}

sub _parse_aa_freqs {
    my ($self) = @_;
    my ( $okay, $done, $header ) = ( 0, 0, 0 );
    my (@bases);
    my $numseqs = scalar @{ $self->{'_summary'}->{'seqs'} || [] };
    while ( defined( $_ = $self->_readline ) ) {
        if ( /^TREE/ || /^AA distances/ ) {
            $self->_pushback($_);
            last;
        }
        last if ($done);
        next if ( /^\s+$/ || /^\(Ambiguity/ );
        if (/^Frequencies\./) {
            $okay = 1;
        }
        elsif ( !$okay ) {    # skip till we see 'Frequencies.
            next;
        }
        elsif ( !$header ) {
            s/^\s+//;         # remove leading whitespace
            @bases  = split;  # get an array of the all the aa names
            $header = 1;
            $self->{'_summary'}->{'aafreqs'} = {};    # reset/clear values
            next;
        }
        elsif (
            /^\#\s+constant\s+sites\:\s+
		 (\d+)\s+ # constant sites
		 \(\s*([\d\.]+)\s*\%\s*\)/x
          )
        {
            $self->{'_summary'}->{'stats'}->{'constant_sites'}            = $1;
            $self->{'_summary'}->{'stats'}->{'constant_sites_percentage'} = $2;
        }
        elsif (/^ln\s+Lmax\s+\(unconstrained\)\s+\=\s+(\S+)/x) {
            $self->{'_summary'}->{'stats'}->{'loglikelihood'} = $1;
            $done = 1;    # done for sure
        }
        else {
            my ( $seqname, @freqs ) = split;
            my $basect = 0;
            foreach my $f (@freqs) {

                # this will also store 'Average'
                $self->{'_summary'}->{'aafreqs'}->{$seqname}
                  ->{ $bases[ $basect++ ] } = $f;
            }
        }
    }
}

# This is for parsing the automatic tree output

sub _parse_StarDecomposition {
    my ($self) = @_;
    my %data;

    return %data;
}

sub _parse_aa_dists {
    my ($self) = @_;
    my ( $okay, $seen, $done ) = ( 0, 0, 0 );
    my ( %matrix, @names, @values );
    my $numseqs = scalar @{ $self->{'_summary'}->{'seqs'} || [] };
    my $type = '';
    while ( defined( $_ = $self->_readline ) ) {
        last if $done;
        if (/^TREE/) { $self->_pushback($_); last; }
        if (/^\s+$/) {
            last if ($seen);
            next;
        }
        if (/^(AA|ML) distances/) {
            $okay = 1;
            $type = $1;
            next;
        }
        s/\s+$//g;    # remove trailing space
        if ($okay) {
            my ( $seqname, @vl ) = split;
            $seen = 1;
            my $i = 0;

            # hacky workaround to problem with 3.14 aaml
            if (
                   $type eq 'ML'
                && !@names
                &&    # first entry
                @vl
              )
            {         # not empty
                push @names, $self->{'_summary'}->{'seqs'}->[0]->display_id;
            }
            for my $s (@names) {
                last unless @vl;
                $matrix{$seqname}->{$s} = $matrix{$s}->{$seqname} = shift @vl;
            }
            push @names, $seqname;

            $matrix{$seqname}->{$seqname} = 0;
        }
        $done = 1 if ( scalar @names == $numseqs );
    }
    my %dist;
    my $i = 0;
    @values = ();
    foreach my $lname (@names) {
        my @row;
        my $j = 0;
        foreach my $rname (@names) {
            my $v = $matrix{$lname}->{$rname};
            $v = $matrix{$rname}->{$lname} unless defined $v;
            push @row, $v;
            $dist{$lname}{$rname} = [ $i, $j++ ];
        }
        $i++;
        push @values, \@row;
    }
    return new Bio::Matrix::PhylipDist(
        -program => $self->{'_summary'}->{'seqtype'},
        -matrix  => \%dist,
        -names   => \@names,
        -values  => \@values
    );
}

sub _parse_patterns {
    my ($self) = @_;
    my ( $patternct, @patterns, $ns, $ls );
    return if exists $self->{'_summary'}->{'patterns'};

    while ( defined( $_ = $self->_readline ) ) {
        if ( /^Codon\s+(usage|position)/ || /Model/ ) {
            $self->_pushback($_);
            last;
        }
        elsif ($patternct) {

            #	    last unless ( @patterns == $patternct );
            last if (/^\s+$/);
            s/^\s+//;
            push @patterns, split;
        }
        elsif (/^ns\s+\=\s*(\d+)\s+ls\s+\=\s*(\d+)/) {
            ( $ns, $ls ) = ( $1, $2 );
        }
        elsif (/^\# site patterns \=\s*(\d+)/) {
            $patternct = $1;
        }
        else {

            #	    $self->debug("Unknown line: $_");
        }
    }
    $self->{'_summary'}->{'patterns'} = {
        -patterns => \@patterns,
        -ns       => $ns,
        -ls       => $ls
    };
}

sub _parse_seqs {

    # this should in fact be packed into a Bio::SimpleAlign object instead of
    # an array but we'll stay with this for now
    my ($self) = @_;

    # Use this flag to deal with paml 4 vs 3 differences
    # In PAML 4 the sequences precede the CODONML|BASEML|AAML
    # while in PAML3 the files start off with this
    return 1 if $self->{'_already_parsed_seqs'};
    my ( @firstseq, @seqs );
    while ( defined( $_ = $self->_readline ) ) {
        if (/^(Printing|After|TREE|Codon)/) {
            $self->_pushback($_);
            last;
        }
        last if ( /^\s+$/ && @seqs > 0 );
        next if (/^\s+$/);
        next if (/^\d+\s+$/);

        # we are reading PHYLIP format
        my ( $name, $seqstr ) = split( /\s+/, $_, 2 );
        $seqstr =~ s/\s+//g;    # remove whitespace
        unless (@firstseq) {
            @firstseq = split( //, $seqstr );
            push @seqs,
              Bio::LocatableSeq->new(
                -display_id => $name,
                -seq        => $seqstr
              );
        }
        else {

            my $i = 0;
            my $v;
            while ( ( $v = index( $seqstr, '.', $i ) ) >= $i ) {

                # replace the '.' with the correct seq from the
                substr( $seqstr, $v, 1, $firstseq[$v] );
                $i = $v;
            }
            push @seqs,
              Bio::LocatableSeq->new(
                -display_id => $name,
                -seq        => $seqstr
              );
        }
    }
    if ( @seqs > 0 ) {
        $self->{'_summary'}->{'seqs'} = \@seqs;
        $self->{'_already_parsed_seqs'} = 1;
    }
    1;
}

sub _parse_codoncts { }

sub _parse_distmat {
    my ($self) = @_;
    my @results;
    my $ver = 3.14;
	my $firstseq, my $secondseq;

    while ( defined( $_ = $self->_readline ) ) {
        next if /^\s+$/;

        # We need to get the names of the sequences if this is from YN00:
        if (/^\(A\)\sNei-Gojobori\s\(1986\)\smethod/) {
            $ver = 3.15;
			while ( defined( $_ = $self->_readline ) ) {
				if ($_ =~ m/.*\d+?\.\d+?\s*\(.*/) {
					$secondseq = $_;
					last;
				}
				$firstseq = $_;
			}
        }
        last;
    }

    #return unless (/^Nei\s*\&\s*Gojobori/);

    # skip the next 3 lines
    if ( $self->{'_summary'}->{'seqtype'} eq 'CODONML' ) {
        $self->_readline;
        $self->_readline;
        $self->_readline;
    }
    my $seqct = 0;
    my @seqs;
	if ( $self->{'_summary'}->{'seqtype'} eq 'YN00' ) {
		if ($firstseq) {
			$firstseq =~ s/(.+?)\s+.*/$1/;
			$secondseq =~ s/(.+?)\s+.*/$1/;
			chomp $firstseq;
			chomp $secondseq;
			push @seqs, Bio::PrimarySeq->new( -display_id => $firstseq );
			push @seqs, Bio::PrimarySeq->new( -display_id => $secondseq );
        }
    }
    while ( defined( $_ = $self->_readline ) ) {
        last if ( /^\s+$/ && exists $self->{'_summary'}->{'ngmatrix'} );
        next if ( /^\s+$/ || /^NOTE:/i );
        chomp;

        my ( $seq, $rest );
        if ( $self->{'_summary'}->{'seqtype'} eq 'YN00' ) {
             ( $seq, $rest ) = split( /\s+/, $_, 2 );
        }
        else {
            $_ =~ m/(.+?)\s*(-*\d+?\.\d+?.*)/;
 	        $seq = $1;
 		$rest = $2;
	}
        $rest = '' unless defined $rest;    # get rid of empty messages
        my $j = 0;
        if ( $self->{'_summary'}->{'seqtype'} eq 'YN00' ) {
            push @seqs, Bio::PrimarySeq->new( -display_id => $seq );
        }
        while ($rest
            && $rest =~
            /(\-?\d+(\.\d+)?)\s*\(\-?(\d+(\.\d+)?)\s+(\-?\d+(\.\d+)?)\)/g )
        {
            $self->{'_summary'}->{'ngmatrix'}->[ $j++ ]->[$seqct] = {
                'omega' => $1,
                'dN'    => $3,
                'dS'    => $5
            };
        }
        $seqct++;
    }
    if ( $self->{'_summary'}->{'seqtype'} eq 'YN00' && @seqs ) {
        $self->{'_summary'}->{'seqs'} = \@seqs;
    }

    1;
}

sub _parse_PairwiseCodon {
    my ($self) = @_;
    my @result;
    my ( $a, $b, $log, $model, $t, $kappa, $omega, $fixedkappa );
    # check to see if we have a fixed kappa:
    if ( $self->{'_summary'}->{'model'} =~ /kappa = (\d+?\.\d+?) fixed/) {
		$fixedkappa = $1;
	}
    while ( defined( $_ = $self->_readline ) ) {
        if (/^pairwise comparison, codon frequencies\:\s*(\S+)\./) {
            $model = $1;
        }
        # 1st line of a pair block, e.g.
        # 2 (all_c7259) ... 1 (all_s57600)
        elsif (/^(\d+)\s+\((\S+)\)\s+\.\.\.\s+(\d+)\s+\((\S+)\)/) {
            ( $a, $b ) = ( $1, $3 );
        }
        # 2nd line of a pair block, e.g.
        # lnL = -126.880601
        elsif (/^lnL\s+\=\s*(\-?\d+(\.\d+)?)/) {
            $log = $1;
            if ( defined( $_ = $self->_readline ) ) {
                # 3rd line of a pair block, e.g.
                # 0.19045  2.92330  0.10941
                s/^\s+//;
                ( $t, $kappa, $omega ) = split;
                # if there was a fixed kappa, there will only be two values here ($t, $omega) and $kappa = $fixedkappa.
                if ($omega eq "") {
                	$omega = $kappa;
                	$kappa = $fixedkappa;
                }
            }
        }
        # 5th line of a pair block, e.g.
        # t= 0.1904  S=     5.8  N=   135.2  dN/dS= 0.1094  dN= 0.0476  dS= 0.4353
        # OR lines like (note last field; this includes a fix for bug #3040)
        # t= 0.0439  S=     0.0  N=   141.0  dN/dS= 0.1626  dN= 0.0146  dS=    nan
		elsif (m/^t\=\s*(\d+(\.\d+)?)\s+/)
        {
        	# Breaking out each piece individually so that you can see
        	# what each regexp actually looks for
        	my $parse_string = $_;
        	$parse_string =~ m/.*t\s*\=\s*(\d+?\.\d+?)\s/;
        	my $temp_t = $1;
        	$parse_string =~ m/\sS\s*\=\s*(\d+?\.\d+?)\s/;
        	my $temp_S = $1;
         	$parse_string =~ m/\sN\s*\=\s*(\d+?\.\d+?)\s/;
        	my $temp_N = $1;
         	$parse_string =~ m/\sdN\/dS\s*\=\s*(\d+?\.\d+?)\s/;
        	my $temp_omega = $1;
         	$parse_string =~ m/\sdN\s*\=\s*(\d+?\.\d+?)\s/;
        	my $temp_dN = $1;
         	$parse_string =~ m/\sdS\s*\=\s*(.+)\s/;
        	my $temp_dS = $1;
            $result[ $b - 1 ]->[ $a - 1 ] = {
                'lnL'   => $log,
                't'     => defined $t && length($t) ? $t : $temp_t,
                'S'     => $temp_S,
                'N'     => $temp_N,
                'kappa' => $kappa,
                'omega' => defined $omega && length($omega) ? $omega : $temp_omega,
                'dN'    => $temp_dN,
                'dS'    => $temp_dS
            };
        }
        # 4th line of a pair block (which is blank)
        elsif (/^\s+$/) {
            next;
        }
        elsif (/^\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)/) {
        }
        else {
            $self->debug("unknown line: $_");
        }
    }
    return ( -mlmatrix => \@result );
}

sub _parse_YN_Pairwise {
    my ($self) = @_;
    my @result;
    while ( defined( $_ = $self->_readline ) ) {
        last if (/^seq\.\s+seq\./);
    }
    while ( defined( $_ = $self->_readline ) ) {
        if (
            m/^\s+(\d+)\s+  # seq #
	    (\d+)\s+        # seq #
	    (\d+(\.\d+))\s+ # S
	    (\d+(\.\d+))\s+ # N
	    (\d+(\.\d+))\s+ # t
	    (\d+(\.\d+))\s+ # kappa
	    (\d+(\.\d+))\s+ # omega
	    \-??(\d+(\.\d+))\s+ # dN
	    \+\-\s+
	    \-??(\d+(\.\d+))\s+ # dN SE
	    \-??(\d+(\.\d+))\s+ # dS
	    \+\-\s+
	    \-??(\d+(\.\d+))\s+ # dS SE
	    /ox
          )
        {

            $result[ $2 - 1 ]->[ $1 - 1 ] = {
                'S'     => $3,
                'N'     => $5,
                't'     => $7,
                'kappa' => $9,
                'omega' => $11,
                'dN'    => $13,
                'dN_SE' => $15,
                'dS'    => $17,
                'dS_SE' => $19,
            };
        }
        elsif (/^\s+$/) {
            next;
        }
        elsif (/^\(C\) LWL85, LPB93 & LWLm methods/) {
            $self->_pushback($_);
            last;
        }
    }
    return ( -mlmatrix => \@result );
}

sub _parse_Forestry {
    my ($self) = @_;
    my ( $instancecount, $num_param, $loglikelihood, $score, $done,
        $treelength ) = ( 0, 0, 0, 0, 0, 0 );
    my $okay = 0;
    my ( @ids, %match, @branches, @trees );
    while ( defined( $_ = $self->_readline ) ) {
        last if $done;
        if (s/^TREE\s+\#\s*\d+:\s+//) {
            ($score) = (s/MP\s+score\:\s+(\S+)\s+$//);
            @ids = /(\d+)[\,\)]/g;
        }
        elsif (/^Node\s+\&/
            || /^\s+N37/
            || /^(CODONML|AAML|YN00|BASEML)/
            || /^\*\*/
            || /^Detailed output identifying parameters/ )
        {
            $self->_pushback($_);
            $done = 1;
            last;
        }
        elsif (/^tree\s+length\s+\=\s+(\S+)/) {
            $treelength = $1;    # not going to store this for now
                                 # as it is directly calculated from
                                 # $tree->total_branch_length;
        }
        elsif (/^\s*lnL\(.+np\:\s*(\d+)\)\:\s+(\S+)/) {

            # elsif( /^\s*lnL\(.+\)\:\s+(\S+)/ ) {
            ( $num_param, $loglikelihood ) = ( $1, $2 );
        }
        elsif (/^\(/) {
            s/([\,:])\s+/$1/g;
            my $treestr = IO::String->new($_);
            my $treeio  = Bio::TreeIO->new(
                -fh     => $treestr,
                -format => 'newick'
            );
            my $tree = $treeio->next_tree;
            if ($tree) {
                $tree->score($loglikelihood);
                $tree->id("num_param:$num_param");
                if ( $okay > 0 ) {

                    # we don't save the trees with the number labels
                    if ( !%match && @ids ) {
                        my $i = 0;
                        for my $m (/([^():,]+):/g) {
                            $match{ shift @ids } = [$m];
                        }
                        my %grp;
                        while ( my $br = shift @branches ) {
                            my ( $parent, $child ) = @$br;
                            if ( $match{$child} ) {
                                push @{ $match{$parent} }, @{ $match{$child} };
                            }
                            else {
                                push @branches, $br;
                            }
                        }
                        if ( $self->verbose > 1 ) {
                            for my $k ( sort { $a <=> $b } keys %match ) {
                                $self->debug( "$k -> ",
                                    join( ",", @{ $match{$k} } ), "\n" );
                            }
                        }
                    }

                    # Associate SEs to nodes using tags
                    if ( defined( $self->{_SEs} ) ) {
                        my @SEs = split( " ", $self->{_SEs} );
                        my $i = 0;
                        foreach my $parent_id ( map { /\d+\.\.(\d+)/ }
                            split( " ", $self->{_branch_ids} ) )
                        {
                            my @nodes;
                            my @node_ids = @{ $match{$parent_id} };
                            my @nodes_L =
                              map { $tree->find_node( -id => $_ ) } @node_ids;
                            my $n =
                              @nodes_L < 2
                              ? shift(@nodes_L)
                              : $tree->get_lca(@nodes_L);
                            if ( !$n ) {
                                $self->warn(
"no node could be found for node in SE assignation (no lca?)"
                                );
                            }
                            $n->add_tag_value( 'SE', $SEs[$i] );
                            $i++;
                        }
                    }
                    push @trees, $tree;
                }
            }
            $okay++;
        }
        elsif (/^SEs for parameters/) {
            my $se_line = $self->_readline;
            $se_line =~ s/\n//;
            $self->{_SEs} = $se_line;
        }
        elsif (/^\s*\d+\.\.\d+/) {
            push @branches, map { [ split( /\.\./, $_ ) ] } split;
            my $ids = $_;
            $ids =~ s/\n//;
            $self->{_branch_ids} = $ids;
        }
    }
    return \@trees, \%match;
}

sub _parse_NSsitesBatch {
    my $self = shift;
    my ( %data, $idlookup );
    my ( $okay, $done ) = ( 0, 0 );
    while ( defined( $_ = $self->_readline ) ) {
        last if $done;
        next if /^\s+$/;
        next unless ( $okay || /^Model\s+\d+/ || /^TREE/ );

        if (/^Model\s+(\d+)/) {
            if ($okay) {

                # this only happens if $okay was already 1 and
                # we hit a Model line
                $self->_pushback($_);
                $done = 1;
            }
            else {
                chomp;
                $data{'-model_num'} = $1;
                ( $data{'-model_description'} ) = (/\:\s+(.+)/);
                $okay = 1;
            }
        }
        elsif (/^Time used\:\s+(\S+)/) {
            $data{'-time_used'} = $1;
            $done = 1;
        }
        elsif (/^kappa\s+\(ts\/tv\)\s+\=\s+(\S+)/) {
            $data{'-kappa'} = $1;
        }
        elsif (/^TREE/) {
            $self->_pushback($_);
            ( $data{'-trees'}, $idlookup ) = $self->_parse_Forestry;
            if ( defined $data{'-trees'}
                && scalar @{ $data{'-trees'} } )
            {
                $data{'-likelihood'} = $data{'-trees'}->[0]->score;
            }
            $okay = 1;
        }
        elsif (/^omega\s+\(dn\/ds\)\s+\=\s+(\S+)/i) {

            # for M0 (single ratio for the entire tree)
            # explicitly put '1.00000' rather than '1', because \d+\.\d{5}
            # is reported in all other cases.
            my @p = (q/1.00000/);    # since there is only one class,
            my @w = $1;
            $data{'-dnds_site_classes'} = {
                'p' => \@p,
                'w' => \@w
            };

            # since no K=X is provided, put 1 here
            $data{q/-num_site_classes/} = 1;
        }
        elsif (
/^(Naive Empirical Bayes)|(Bayes Empirical Bayes)|(Positively\sselected\ssites)/i
          )
        {
            $self->_pushback($_);
            my ( $sites, $neb, $beb ) = $self->_parse_Pos_selected_sites;
            $data{'-pos_sites'} = $sites;
            $data{'-neb_sites'} = $neb;
            $data{'-beb_sites'} = $beb;
        }
        elsif (/^dN/i) {
            if (/K\=(\d+)/) {
                $data{'-num_site_classes'} = $1;
                while ( $_ = $self->_readline ) {
                    unless ( $_ =~ /^\s+$/ ) {
                        $self->_pushback($_);
                        last;
                    }
                }
                if (/^site class/) {
                    $self->_readline;
                    my $tmp = $self->_readline;
                    my @p = $tmp =~ /(\d+\.\d{5})/g;
                    $tmp = $self->_readline;
                    my @b_w = $tmp =~ /(\d+\.\d{5})/g;
                    $tmp = $self->_readline;
                    my @f_w = $tmp =~ /(\d+\.\d{5})/g;
                    my @w;

                    foreach my $i ( 0 .. $#b_w ) {
                        push @w,
                          {
                            q/background/ => $b_w[$i],
                            q/foreground/ => $f_w[$i]
                          };
                    }
                    $data{'-dnds_site_classes'} = {
                        q/p/ => \@p,
                        q/w/ => \@w
                    };
                }
                else {
                    my $tmp = $self->_readline;
                    my @p = $tmp =~ /(\d+\.\d{5})/g;
                    $tmp = $self->_readline;
                    my @w = $tmp =~ /(\d+\.\d{5})/g;
                    $data{'-dnds_site_classes'} = {
                        'p' => \@p,
                        'w' => \@w
                    };
                }
            }
            elsif (/for each branch/) {
                my %branch_dnds = $self->_parse_branch_dnds;
                if ( !defined $data{'-trees'} ) {
                    $self->warn(
                        "No trees have been loaded, can't do anything\n");
                    next;
                }
                my ($tree) = @{ $data{'-trees'} };
                if (   !$tree
                    || !ref($tree)
                    || !$tree->isa('Bio::Tree::Tree') )
                {
                    $self->warn("no tree object already stored!\n");
                    next;
                }

                # These need to be added to the Node/branches
                while ( my ( $k, $v ) = each %branch_dnds ) {

                    # we can probably do better by caching at some point
                    my @nodes;
                    for my $id ( split( /\.\./, $k ) ) {
                        my @nodes_L =
                          map { $tree->find_node( -id => $_ ) }
                          @{ $idlookup->{$id} };
                        my $n =
                          @nodes_L < 2
                          ? shift(@nodes_L)
                          : $tree->get_lca(@nodes_L);
                        if ( !$n ) {
                            $self->warn(
                                "no node could be found for $id (no lca?)");
                        }
                        unless ( $n->is_Leaf && $n->id ) {
                            $n->id($id);
                        }
                        push @nodes, $n;
                    }
                    my ( $parent, $child ) = @nodes;
                    while ( my ( $kk, $vv ) = each %$v ) {
                        $child->add_tag_value( $kk, $vv );
                    }
                }
            }
        }
        elsif (/^Parameters in beta:/) {
            $_ = $self->_readline;    # need the next line
            if (/p\=\s+(\S+)\s+q\=\s+(\S+)/) {
                $data{'-shape_params'} = {
                    'shape' => 'beta',
                    'p'     => $1,
                    'q'     => $2
                };
            }
            else {
                $self->warn("unparseable beta parameters: $_");
            }
        }
        elsif (/^Parameters in beta\&w\>1:/) {

            # Parameters in beta&w>1:
            #   p0=  1.00000  p=  0.07642 q=  0.85550
            #  (p1=  0.00000) w=  1.00000
            $_ = $self->_readline;    # need the next line
            my ( $p0, $p, $q, $p1, $w );
            if (/p0\=\s+(\S+)\s+p\=\s+(\S+)\s+q\=\s+(\S+)/) {
                $p0 = $1;
                $p  = $2;
                $q  = $3;
            }
            else {
                $self->warn("unparseable beta parameters: $_");
            }
            $_ = $self->_readline;    # need the next line
            if (/\(p1\=\s+(\S+)\)\s+w\=\s*(\S+)/) {
                $p1                    = $1;
                $w                     = $2;
                $data{'-shape_params'} = {
                    'shape' => 'beta',
                    'p0'    => $p0,
                    'p'     => $p,
                    'q'     => $q,
                    'p1'    => $p1,
                    'w'     => $w
                };
            }
            else {
                $self->warn("unparseable beta parameters: $_");
            }
        }
        elsif (/^alpha\s+\(gamma\)\s+\=\s+(\S+)/) {
            my $gamma = $1;
            $_ = $self->_readline;
            my ( @r, @f );
            if (s/^r\s+\(\s*\d+\)\:\s+//) {
                @r = split;
            }
            $_ = $self->_readline;
            if (s/^f\s*\:\s+//) {
                @f = split;
            }
            $data{'-shape_params'} = {
                'shape' => 'alpha',
                'gamma' => $gamma,
                'r'     => \@r,
                'f'     => \@f
            };
        }
    }
    return new Bio::Tools::Phylo::PAML::ModelResult(%data);
}

sub _parse_Pos_selected_sites {
    my $self    = shift;
    my $okay    = 0;
    my (%sites) = (
        'default' => [],
        'neb'     => [],
        'beb'     => []
    );
    my $type = 'default';
    while ( defined( $_ = $self->_readline ) ) {
        next if ( /^\s+$/ || /^\s+Pr\(w\>1\)/ );
        if ( /^Time used/ || /^TREE/ ) {
            $self->_pushback($_);
            last;
        }
        if (/^Naive Empirical Bayes/i) {
            $type = 'neb';
        }
        elsif (/^Bayes Empirical Bayes/i) {
            $type = 'beb';
        }
        elsif (/^Positively selected sites/) {
            $okay = 1;
        }
        elsif ( $okay
            && /^\s+(\d+)\s+(\S+)\s+(\-?\d+(?:\.\d+)?)(\**)\s+(\-?\d+(?:\.\d+)?)\s+\+\-\s+(\-?\d+(?:\.\d+)?)/
          )
        {
            my $signif = $4;
            $signif = '' unless defined $signif;
            push @{ $sites{$type} }, [ $1, $2, $3, $signif, $5, $6 ];
        }
        elsif ( $okay
            && /^\s+(\d+)\s+(\S+)\s+(\-?\d*(?:.\d+))(\**)\s+(\-?\d+(?:\.\d+)?)/
          )
        {
            my $signif = $4;
            $signif = '' unless defined $signif;
            push @{ $sites{$type} }, [ $1, $2, $3, $signif, $5 ];
        }
        elsif ( $okay && /^\s+(\d+)\s+(\S)\s+([\d\.\-\+]+)(\**)/ ) {
            my $signif = $4;
            $signif = '' unless defined $signif;
            push @{ $sites{$type} }, [ $1, $2, $3, $signif ];
        }
    }
    return ( $sites{'default'}, $sites{'neb'}, $sites{'beb'} );
}

sub _parse_branch_dnds {
    my $self = shift;
    my ($okay) = (0);
    my %branch_dnds;
    my @header;
    while ( defined( $_ = $self->_readline ) ) {
        next if (/^\s+$/);
        next unless ( $okay || /^\s+branch\s+t/ );
        if (/^\s+branch\s+(.+)/) {
            s/^\s+//;
            @header = split( /\s+/, $_ );
            $okay = 1;
        }
        elsif (/^\s*(\d+\.\.\d+)/) {
            my $branch = $1;
            s/^\s+//;
            my $i = 0;

            # fancyness just maps the header names like 't' or 'dN'
            # into the hash so we get at the end of the day
            # 't' => 0.067
            # 'dN'=> 0.001
            $branch_dnds{$branch} = { map { $header[ $i++ ] => $_ } split };
        }
        else {
            $self->_pushback($_);
            last;
        }
    }
    return %branch_dnds;
}

#baseml stuff
sub _parse_nt_freqs {
    my ($self) = @_;
    my ( $okay, $done, $header ) = ( 0, 0, 0 );
    my (@bases);
    my $numseqs = scalar @{ $self->{'_summary'}->{'seqs'} || [] };
    while ( defined( $_ = $self->_readline ) ) {
        if ( /^TREE/ || /^Distances/ ) { $self->_pushback($_); last }
        last if ($done);
        next if ( /^\s+$/ || /^\(Ambiguity/ );
        if (/^Frequencies\./) {
            $okay = 1;
        }
        elsif ( !$okay ) {    # skip till we see 'Frequencies.
            next;
        }
        elsif ( !$header ) {
            s/^\s+//;         # remove leading whitespace
            @bases  = split;  # get an array of the all the aa names
            $header = 1;
            $self->{'_summary'}->{'ntfreqs'} = {};    # reset/clear values
            next;
        }
        elsif (
            /^\#\s+constant\s+sites\:\s+
		 (\d+)\s+	# constant sites
		 \(\s*([\d\.]+)\s*\%\s*\)/ox
          )
        {
            $self->{'_summary'}->{'stats'}->{'constant_sites'}            = $1;
            $self->{'_summary'}->{'stats'}->{'constant_sites_percentage'} = $2;
        }
        elsif (/^ln\s+Lmax\s+\(unconstrained\)\s+\=\s+(\S+)/ox) {
            $self->{'_summary'}->{'stats'}->{'loglikelihood'} = $1;
            $done = 1;    # done for sure
        }
        else {
            my ( $seqname, @freqs ) = split;
            my $basect = 0;
            foreach my $f (@freqs) {

                # this will also store 'Average'
                $self->{'_summary'}->{'ntfreqs'}->{$seqname}
                  ->{ $bases[ $basect++ ] } = $f;
            }
        }
    }
}

sub _parse_nt_dists {
    my ($self) = @_;
    my ( $okay, $seen, $done ) = ( 0, 0, 0 );
    my ( %matrix, @names );
    my $numseqs = scalar @{ $self->{'_summary'}->{'seqs'} || [] };
    my $type = '';
    while ( defined( $_ = $self->_readline ) ) {
        if (/^TREE/) { $self->_pushback($_); last; }
        last if $done;
        next if (/^This matrix is not used in later/);
        if (/^\s+$/) {
            last if ($seen);
            next;
        }
        if (/^Distances:(\S+)\s+\(([^\)]+)\)\s+\(alpha set at (\-?\d+\.\d+)\)/)
        {
            $okay = 1;
            $type = $1;
            next;
        }
        s/\s+$//g;    # remove trailing space
        if ($okay) {
            my ( $seqname, $vl ) = split( /\s+/, $_, 2 );
            $seen = 1;
            my $i = 0;
            if ( defined $vl ) {
                while ( $vl =~ /(\-?\d+\.\d+)\s*\(\s*(\-?\d+\.\d+)\s*\)\s*/g ) {
                    my ( $kappa, $alpha ) = ( $1, $2 );
                    $matrix{$seqname}{ $names[$i] } =
                      $matrix{ $names[$i] }{$seqname} = [ $kappa, $alpha ];

                    $i++;
                }
                unless ($i) {
                    $self->warn("no matches for $vl\n");
                }
            }
            push @names, $seqname;
            $matrix{$seqname}->{$seqname} = [ 0, 0 ];
        }
        $done = 1 if ( scalar @names == $numseqs );
    }
    my %dist;
    my $i = 0;
    my ( @kvalues, @avalues );
    foreach my $lname (@names) {
        my ( @arow, @krow );
        my $j = 0;
        foreach my $rname (@names) {
            my $v = $matrix{$lname}{$rname};

            push @krow, $v->[0];    # kappa values
            push @arow, $v->[1];    # alpha
            $dist{$lname}{$rname} = [ $i, $j++ ];
        }
        $i++;
        push @kvalues, \@krow;
        push @avalues, \@arow;
    }
    return (
        Bio::Matrix::PhylipDist->new(
            -program => $self->{'_summary'}->{'seqtype'},
            -matrix  => \%dist,
            -names   => \@names,
            -values  => \@kvalues
        ),
        Bio::Matrix::PhylipDist->new(
            -program => $self->{'_summary'}->{'seqtype'},
            -matrix  => \%dist,
            -names   => \@names,
            -values  => \@avalues
        )
    );
}

# BASEML
sub _parse_rate_parametes {
    my $self = shift;
    my (%rate_parameters);
    while ( defined( $_ = $self->_readline ) ) {
        if (/^Rate\s+parameters:\s+/) {
            s/\s+$//;
            $rate_parameters{'rate_parameters'} = [ split( /\s+/, $_ ) ];
        }
        elsif (/^Base\s+frequencies:\s+/) {
            s/\s+$//;
            $rate_parameters{'base_frequencies'} = [ split( /\s+/, $_ ) ];
        }
        elsif (
            m/^Rate\s+matrix\s+Q,\s+Average\s+Ts\/Tv\s+(\([^\)+]+\))?\s*\=\s+
		 (\-?\d+\.\d+)/x
          )
        {
            $rate_parameters{'average_TsTv'} = $1;
            while ( defined( $_ = $self->_readline ) ) {

                # short circuit
                last if (/^\s+$/);
                if (/^alpha/) {
                    $self->_pushback($_);
                    last;
                }
                s/^\s+//;
                s/\s+$//;
                push @{ $rate_parameters{'rate_matrix_Q'} }, [split];
            }
        }
        elsif (/^alpha\s+\(gamma,\s+K=\s*(\d+)\s*\)\s*\=\s*(\-?\d+\.\d+)/) {
            $rate_parameters{'K'}     = $1;
            $rate_parameters{'alpha'} = $2;
        }
        elsif (s/^(r|f):\s+//) {
            my ($p) = $1;
            s/\s+$//;
            $rate_parameters{$p} = [split];
        }
    }
}

# RST parsing
sub _parse_rst {
    my ($self) = @_;
    return unless $self->{'_dir'} && -d $self->{'_dir'} && -r $self->{'_dir'};

    my $rstfile = File::Spec->catfile( $self->{'_dir'}, $RSTFILENAME );
    return unless -e $rstfile && !-z $rstfile;
    my $rstio = Bio::Root::IO->new( -file => $rstfile );

    # define whatever data structures you need to store the data
    # key points are to reuse existing bioperl objs (like Bio::Seq)
    # where appropriate
    my ( @firstseq, @seqs, @trees, @per_site_prob );
    my $count;
    while ( defined( $_ = $rstio->_readline ) ) {

        # implement the parsing here
        if (/^TREE\s+\#\s+(\d+)/) {
            while ( defined( $_ = $rstio->_readline ) ) {
                if (/tree\s+with\s+node\s+labels\s+for/) {
                    my $tree = Bio::TreeIO->new(
                        -noclose => 1,
                        -fh      => $rstio->_fh,
                        -format  => 'newick'
                    )->next_tree;

                    # cleanup leading/trailing whitespace
                    for my $n ( $tree->get_nodes ) {
                        my $id = $n->id;
                        $id =~ s/^\s+//;
                        $id =~ s/\s+$//;
                        $n->id($id);
                        if ( defined( my $blen = $n->branch_length ) ) {
                            $blen =~ s/^\s+//;
                            $blen =~ s/\s+$//;
                            $n->branch_length($blen);
                        }
                    }
                    push @trees, $tree;
                    last;
                }
            }
        }
        elsif (/^Prob\sof\sbest\scharacter\sat\seach\snode,\slisted\sby\ssite/)
        {
            $self->{'_rst'}->{'persite'} = [];
            while ( defined( $_ = $rstio->_readline ) ) {
                next if ( /^Site/i || /^\s+$/ );
                if (s/^\s+(\d+)\s+(\d+)\s+([^:]+)\s*:\s*(.+)//) {
                    my ( $sitenum, $freq, $extant, $ancestral ) =
                      ( $1, $2, $3, $4 );
                    my ( @anc_site, @extant_site );
                    @extant_site = {};
                    while ( $extant =~ s/^([A-Z\-]{3})\s+\(([A-Z*])\)\s+//g ) {
                        push @extant_site, { 'codon' => $1, 'aa' => $2 };
                    }
                    while (
                        $ancestral =~ s/^([A-Z\-]{3})\s+([A-Z*])\s+ # codon AA
			                (\S+)\s+                   # Prob
			                \(([A-Z*])\s+(\S+)\)\s*//xg    # AA Prob
                      )
                    {
                        push @anc_site,
                          {
                            'codon'          => $1,
                            'aa'             => $2,
                            'prob'           => $3,
                            'Yang95_aa'      => $4,
                            'Yang95_aa_prob' => $5
                          };
                    }

                    # saving persite
                    $self->{'_rst'}->{'persite'}->[$sitenum] =
                      [ @extant_site, @anc_site ];
                }
                elsif (/^Summary\sof\schanges\salong\sbranches\./) {
                    last;
                }
            }
        }
        elsif (/^Check\sroot\sfor\sdirections\sof\schange\./
            || /^Summary\sof\schanges\salong\sbranches\./ )
        {
            my ( @branches, @branch2node, $branch, $node );
            my $tree = $trees[-1];
            if ( !$tree ) {
                $self->warn("No tree built before parsing Branch changes\n");
                last;
            }
            my @nodes = (
                map  { $_->[0] }
                sort { $a->[1] <=> $b->[1] }
                map  { [ $_, $_->id =~ /^(\d+)\_?/ ] } $tree->get_nodes
            );
            unshift @nodes,
              undef;    # fake first node so that index will match nodeid
            while ( defined( $_ = $rstio->_readline ) ) {
                next if /^\s+$/;
                if (m/^List\sof\sextant\sand\sreconstructed\ssequences/) {
                    $rstio->_pushback($_);
                    last;
                }
                elsif (/^Branch\s+(\d+):\s+(\d+)\.\.(\d+)\s+/) {
                    my ( $left, $right );
                    ( $branch, $left, $right ) = ( $1, $2, $3 );
                    ($node) = $nodes[$right];
                    if ( !$node ) {
                        $self->warn(
"cannot find $right in $tree ($branch $left..$right)\n"
                        );
                        last;
                    }
                    if (/\(n=\s*(\S+)\s+s=\s*(\S+)\)/) {
                        $node->add_tag_value( 'n', $1 );
                        $node->add_tag_value( 's', $2 );
                    }
                    $branch2node[$branch] = $right;
                }
                elsif (
                    /^\s+(\d+)\s+([A-Z*])\s+(\S+)\s+\-\>\s+([A-Z*])\s+(\S+)?/)
                {
                    my ( $site, $anc, $aprob, $derived, $dprob ) =
                      ( $1, $2, $3, $4, $5 );
                    if ( !$node ) {
                        $self->warn("no branch line was previously parsed!");
                        next;
                    }
                    my %c = (
                        'site'       => $site,
                        'anc_aa'     => $anc,
                        'anc_prob'   => $aprob,
                        'derived_aa' => $derived,
                    );
                    $c{'derived_prob'} = $dprob if defined $dprob;
                    $node->add_tag_value( 'changes', \%c );
                }
            }
        }
        elsif (
            /^Overall\s+accuracy\s+of\s+the\s+(\d+)\s+ancestral\s+sequences:/)
        {
            my $line = $rstio->_readline;
            $line =~ s/^\s+//;
            $line =~ s/\s+$//;
            my @overall_site = split( /\s+/, $line );

            # skip next 2 lines, want the third
            for ( 1 .. 3 ) {
                $line = $rstio->_readline;
            }
            $line =~ s/^\s+//;
            $line =~ s/\s+$//;
            my @overall_seq = split( /\s+/, $line );
            if (   @overall_seq != @overall_site
                || @overall_seq != @seqs )
            {
                $self->warn(
                    "out of sync somehow seqs, site scores don't match\n");
                $self->warn("@seqs @overall_seq @overall_site\n");
            }
            for (@seqs) {
                $_->description(
                    sprintf(
                        "overall_accuracy_site=%s overall_accuracy_seq=%s",
                        shift @overall_site,
                        shift @overall_seq
                    )
                );
            }
        }
        elsif (m/^List of extant and reconstructed sequences/o) {
            my $seqcount = 0;
            while ( defined( $_ = $rstio->_readline ) ) {
                last if (/^Overall accuracy of the/);
                if (/^\s+$/) {
                    last if $seqcount && $seqcount == @seqs;
                    next;
                }
                if (/^\s*(\d+)\s+(\d+)\s+$/) { $seqcount = $1; next }

                # runmode = (0)
                # this should in fact be packed into a Bio::SimpleAlign object
                # instead of an array but we'll stay with this for now
                if (/^node/) {
                    my ( $name, $num, $seqstr ) = split( /\s+/, $_, 3 );
                    $name .= $num;
                    $seqstr =~ s/\s+//g;    # remove whitespace
                    unless (@firstseq) {
                        @firstseq = split( //, $seqstr );
                        push @seqs,
                          Bio::LocatableSeq->new(
                            -display_id => $name,
                            -seq        => $seqstr
                          );
                    }
                    else {
                        my $i = 0;
                        my $v;
                        while ( ( $v = index( $seqstr, '.', $i ) ) >= $i ) {

                            # replace the '.' with the correct seq from the
                            substr( $seqstr, $v, 1, $firstseq[$v] );
                            $i = $v;
                        }
                        $self->debug("adding seq $seqstr\n");
                        push @seqs,
                          Bio::LocatableSeq->new(
                            -display_id => $name,
                            -seq        => $seqstr
                          );
                    }
                }
            }
            $self->{'_rst'}->{'rctrted_seqs'} = \@seqs;
        }
        else {
        }
    }
    $self->{'_rst'}->{'trees'} = \@trees;
    return;
}

1;
