# PAML.pm,v 1.3 2002/06/20 18:50:37 amackey Exp
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
  http://bugzilla.bioperl.org/

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
use Bio::AnalysisParserI;
use Bio::Root::IO;

@ISA = qw(Bio::Root::Root Bio::Root::IO Bio::AnalysisParserI);

# other objects used:
use IO::String;
use Bio::TreeIO;
use Bio::Tools::Phylo::PAML::Result;
use Bio::PrimarySeq;
use Bio::Matrix::PhylipDist;
use Bio::Tools::Phylo::PAML::ModelResult;

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
  $self->_initialize_io(@args);

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
	unless ($self->{'_summary'} && !$self->{'_summary'}->{'multidata'});
    
    # OK, depending on seqtype and runmode now, one of a few things can happen:
    my $seqtype = $self->{'_summary'}->{'seqtype'};
    if ($seqtype eq 'CODONML' || $seqtype eq 'AAML') {
	while (defined ($_ = $self->_readline)) {
	    if ($seqtype eq 'CODONML' && 
		m/^pairwise comparison, codon frequencies:/o) {

		# runmode = -2, CODONML
		$self->_pushback($_);
		%data = $self->_parse_PairwiseCodon;
		last;
	    } elsif ($seqtype eq 'AAML' && m/^ML distances of aa seqs\.$/o) {
		$self->_pushback($_);
		# get AA distances
		%data = ( '-AAMLdistmat' => $self->_parse_aa_dists());

		# $self->_pushback($_);
		# %data = $self->_parse_PairwiseAA;
		# last;	    
	    } elsif (m/^Model\s+(\d+)/o) {
		$self->_pushback($_);
		my $model = $self->_parse_NSsitesBatch;
		push @{$data{'-NSsitesresults'}}, $model;
	    } elsif (m/^TREE/) {
		# runmode = 0
		$self->_pushback($_);
		$data{'-trees'} = [$self->_parse_Forestry];
		last;

	    } elsif (m/Heuristic tree search by stepwise addition$/o) {
		
		# runmode = 3
		$self->throw( -class => 'Bio::Root::NotImplemented',
			      -text  => "StepwiseAddition not yet implemented!"
			      );

		# $self->_pushback($_);
		# %data = $self->_parse_StepwiseAddition;
		# last;

	    } elsif (m/Heuristic tree search by NNI perturbation$/o) {

		# runmode = 4
		$self->throw( -class => 'Bio::Root::NotImplemented',
			      -text  => "NNI Perturbation not yet implemented!"
			      );

		# $self->_pushback($_);
		# %data = $self->_parse_Perturbation;
		# last;

	    } elsif (m/^stage 0:/o) {

		# runmode = (1 or 2)
		$self->throw( -class => 'Bio::Root::NotImplemented',
			      -text  => "StarDecomposition not yet implemented!"
			      );

		$self->_pushback($_);
		%data = $self->_parse_StarDecomposition;
		last;

	    }
	}
    } elsif ($seqtype eq 'BASEML') {
    } elsif ($seqtype eq 'YN00') {
	while ($_ = $self->_readline) {
	    if( m/^Estimation by the method/ ) {
		$self->_pushback($_);
		%data = $self->_parse_YN_Pairwise;
		last;
	    }
	}
    }
    if (%data) {
	$data{'-version'}   = $self->{'_summary'}->{'version'};
	$data{'-seqs'}      = $self->{'_summary'}->{'seqs'};
	$data{'-patterns'}  = $self->{'_summary'}->{'patterns'};
	$data{'-ngmatrix'}  = $self->{'_summary'}->{'ngmatrix'};
	$data{'-codonpos'}  = $self->{'_summary'}->{'codonposition'};
	$data{'-codonfreq'} = $self->{'_summary'}->{'codonfreqs'};
	$data{'-model'}     = $self->{'_summary'}->{'model'};
	$data{'-aadistmat'} = $self->{'_summary'}->{'aadistmat'};
	$data{'-stats'}     = $self->{'_summary'}->{'stats'};
	$data{'-aafreq'}    = $self->{'_summary'}->{'aafreqs'};
	
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

	    @{$self->{'_summary'}}{qw(seqtype version treefile model)} = ($1, 
									$2,
									$3,
									$4);
	    defined $self->{'_summary'}->{'model'} &&
		$self->{'_summary'}->{'model'} =~ s/Model:\s+//;
	    last;
	    
	} elsif (m/^Data set \d$/o) {
	    $self->{'_summary'} = {};
	    $self->{'_summary'}->{'multidata'}++;
	}
    }

    unless (defined $self->{'_summary'}->{'seqtype'}) {	
	$self->throw( -class => 'Bio::Root::NotImplemented',
		      -text => 'Unknown format of PAML output');
    }
    my $seqtype = $self->{'_summary'}->{'seqtype'};
    $self->debug( "seqtype is $seqtype\n");
    if ($seqtype eq "CODONML") {

	$self->_parse_inputparams(); # settings from the .ctl file 
	                             # that get printed
	$self->_parse_patterns();    # codon patterns - not very interesting
	$self->_parse_seqs();        # the sequences data used for analysis
	$self->_parse_codoncts();    # counts and distributions of codon/nt
	                             # usage
	$self->_parse_codon_freqs(); # codon frequencies
	$self->_parse_distmat();     # NG distance matrices
	

    } elsif ($seqtype eq "AAML") {

	$self->_parse_inputparams;

	$self->_parse_patterns();
        $self->_parse_seqs();     # the sequences data used for analysis
	$self->_parse_aa_freqs(); # codon frequencies

	# get AA distances
	$self->{'_summary'}->{'aadistmat'} = $self->_parse_aa_dists();

    } elsif ($seqtype eq "CODON2AAML") {
	$self->throw( -class => 'Bio::Root::NotImplemented',
		      -text => 'CODON2AAML parsing not yet implemented!');
    } elsif ($seqtype eq "BASEML") {
	$self->throw( -class => 'Bio::Root::NotImplemented',
		      -text => 'BASEML parsing not yet implemented!');
    } elsif ($seqtype eq "YN00") {
	$self->_parse_codon_freqs();
	$self->_parse_codoncts();	
	$self->_parse_distmat(); # NG distance matrices
	
    } else {
	$self->throw( -class => 'Bio::Root::NotImplemented',
		      -text => 'Unknown seqtype, not yet implemented!',
		      -value => $seqtype
		    );
    }

}


sub _parse_inputparams { 
    my ($self) = @_;
    
}

sub _parse_codon_freqs {
    my ($self) = @_;
    my ($okay,$done) = (0,0);
    
    while( defined($_ = $self->_readline ) ) {
	if( /^Nei/ ) { $self->_pushback($_); last }
	last if( $done);
	next if ( /^\s+/);
	next unless($okay || /^Codon position x base \(3x4\) table\, overall/ );
	$okay = 1;
	if( s/^position\s+(\d+):\s+// ) {
	    my $pos = $1;
	    s/\s+$//;
	    my @bases = split;
	    foreach my $str ( @bases ) {
		my ( $base,$freq) = split(/:/,$str,2);
		$self->{'_summary'}->{'codonposition'}->[$pos-1]->{$base} = $freq;
	    }
	    $done = 1 if $pos == 3;
	} 
    }
    $done = 0;
    while( defined( $_ = $self->_readline) ) {
	if( /^Nei\s\&\sGojobori/ ) { $self->_pushback($_); last }
	last if ( $done );
	if( /^Codon frequencies under model, for use in evolver:/ ){
	    while( defined( $_ = $self->_readline) ) {
		last if( /^\s+$/ );
		s/^\s+//;
		s/\s+$//;
		push @{$self->{'_summary'}->{'codonfreqs'}},[split];
	    }
	    $done = 1;
	}
    }
}

sub _parse_aa_freqs {
    my ($self) = @_;
    my ($okay,$done,$header) = (0,0,0);
    my (@bases);
    my $numseqs = scalar @{$self->{'_summary'}->{'seqs'} || []};
    while( defined($_ = $self->_readline ) ) {
	if( /^TREE/o || /^AA distances/o ) { $self->_pushback($_); last }
	last if( $done);
	next if ( /^\s+$/o || /^\(Ambiguity/o );
	if( /^Frequencies\./o ) { 
	    $okay = 1;
	} elsif( ! $okay ) { # skip till we see 'Frequencies.
	    next;
	} elsif ( ! $header ) {
	    s/^\s+//;        # remove leading whitespace
	    @bases = split;  # get an array of the all the aa names
	    $header = 1;
	    $self->{'_summary'}->{'aafreqs'} = {}; # reset/clear values
	    next;
	} elsif( /^\#\s+constant\s+sites\:\s+
		 (\d+)\s+ # constant sites
		 \(\s*([\d\.]+)\s*\%\s*\)/ox){
	    $self->{'_summary'}->{'stats'}->{'constant_sites'} = $1;
	    $self->{'_summary'}->{'stats'}->{'constant_sites_percentage'} = $2;
	} elsif( /^ln\s+Lmax\s+\(unconstrained\)\s+\=\s+(\S+)/ox ) {
	    $self->{'_summary'}->{'stats'}->{'loglikelihood'} = $1;
	    $done = 1; # done for sure
	} else { 
	    my ($seqname,@freqs) = split;
	    my $basect = 0;
	    foreach my $f ( @freqs ) { 
		# this will also store 'Average'
		$self->{'_summary'}->{'aafreqs'}->{$seqname}->{$bases[$basect++]} = $f;
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
    my ($okay,$seen,$done) = (0,0,0);
    my (%matrix,@names,@values);
    my $numseqs = scalar @{$self->{'_summary'}->{'seqs'} || []};
    
    while( defined ($_ = $self->_readline ) ) {
	last if $done;
	if( /^TREE/ ) { $self->_pushback($_); last; }
	if( /^\s+$/ ) {
	    last if( $seen );
	    next;
	}
	$okay = 1, next if( /^(AA|ML) distances/ );
	if( $okay ) {
	    my ($seqname,@values) = split;
	    $seen = 1;
	    my $i = 0;
	    for my $s ( @names ) {
		$matrix{$seqname}->{$s} = $values[$i++];
	    }
	    push @names, $seqname;
	    $matrix{$seqname}->{$seqname} = 0;
	}
	$done = 1 if( scalar @names == $numseqs);
    }
    my %dist;
    my $i = 0;
    foreach my $lname ( @names ) {
	my @row;
	my $j = 0;
	foreach my $rname ( @names ) {
	    my $v = $matrix{$lname}->{$rname};
	    $v = $matrix{$rname}->{$lname} unless defined $v;
	    $matrix{$rname}->{$lname} = $matrix{$lname}->{$rname} = $v;
	    push @row, $v;
	    $dist{$lname}{$rname} = [$i,$j++];
	}
	$i++;
	push @values, \@row;
    }
    return new Bio::Matrix::PhylipDist(-program=> $self->{'_summary'}->{'seqtype'},
				       -matrix => \%dist,
				       -names  => \@names,
				       -values => \@values);

    1;
}

sub _parse_patterns { 
    my ($self) = @_;
    my ($patternct,@patterns,$ns,$ls);
    while( defined($_ = $self->_readline) ) {
	if( $patternct ) { 
#	    last unless ( @patterns == $patternct );
	    last if( /^\s+$/ );
	    s/^\s+//;
	    push @patterns, split;
	} elsif( /^ns\s+\=\s*(\d+)\s+ls\s+\=\s*(\d+)/ ) {
	    ($ns,$ls) = ($1,$2);
	} elsif( /^\# site patterns \=\s*(\d+)/ ) {
	    $patternct = $1;
	} else { 
#	    $self->debug("Unknown line: $_");
	}
    }
    $self->{'_summary'}->{'patterns'} = { -patterns => \@patterns,
					  -ns       => $ns,
					  -ls       => $ls};
}

sub _parse_seqs { 

    # this should in fact be packed into a Bio::SimpleAlign object instead of
    # an array but we'll stay with this for now 
    my ($self) = @_;
    my (@firstseq,@seqs);
    while( defined ($_ = $self->_readline) ) {
	last if( /^TREE/ );
	last if( /^\s+$/ && @seqs > 0 );
	next if ( /^\s+$/ );
	next if( /^\d+\s+$/ );

	my ($name,$seqstr) = split(/\s+/,$_,2);
	$seqstr =~ s/\s+//g; # remove whitespace 
	unless( @firstseq) {
	    @firstseq = split(//,$seqstr);
	    push @seqs, new Bio::PrimarySeq(-id  => $name,
					    -seq => $seqstr);
	} else { 

	    my $i = 0;
	    my $v;
	    while(($v = index($seqstr,'.',$i)) >= $i ) {
		# replace the '.' with the correct seq from the
		substr($seqstr,$v,1,$firstseq[$v]);
		$i = $v;
	    }
	    $self->debug( "adding seq $seqstr\n");
	    push @seqs, new Bio::PrimarySeq(-id  => $name,
					    -seq => $seqstr);
	}
    }
    $self->{'_summary'}->{'seqs'} = \@seqs;
    1;
}

sub _parse_codoncts { }

sub _parse_distmat { 
    my ($self) = @_;
    my @results;
    while( defined ($_ = $self->_readline) ) {
	next if/^\s+$/;
	last;
    }
    return unless (/^Nei\s*\&\s*Gojobori/);
    # skip the next 3 lines
    if( $self->{'_summary'}->{'seqtype'} eq 'CODONML' ) {
	$self->_readline;
	$self->_readline;
	$self->_readline;
    }
    my $seqct = 0;
    my @seqs;
    while( defined ($_ = $self->_readline ) ) {	
	last if( /^\s+$/ && exists $self->{'_summary'}->{'ngmatrix'} );
	next if( /^\s+$/ || /^NOTE:/i );
	chomp;
	my ($seq,$rest) = split(/\s+/,$_,2);
	$rest = '' unless defined $rest; # get rid of empty messages
	my $j = 0;
	if( $self->{'_summary'}->{'seqtype'} eq 'YN00') {
	    push @seqs, Bio::PrimarySeq->new(-display_id => $seq);
	}
	while ( $rest && $rest =~ 
	       /(\-?\d+(\.\d+)?)\s*\(\-?(\d+(\.\d+)?)\s+(\-?\d+(\.\d+)?)\)/g ) {
	    $self->{'_summary'}->{'ngmatrix'}->[$j++]->[$seqct] = 
	    { 'omega' => $1,
	      'dN'    => $3,
	      'dS'    => $5 };
	}
	$seqct++;
    }
    $self->{'_summary'}->{'seqs'} = \@seqs 
	if($self->{'_summary'}->{'seqtype'} eq 'YN00' && @seqs );
    1;
}


sub _parse_PairwiseCodon {
    my ($self) = @_;
    my @result;
    my ($a,$b,$log,$model);
    while( defined( $_ = $self->_readline) ) {
	if( /^pairwise comparison, codon frequencies\:\s*(\S+)\./) {
	    $model = $1;
	} elsif( /^(\d+)\s+\((\S+)\)\s+\.\.\.\s+(\d+)\s+\((\S+)\)/ ) {
	    ($a,$b) = ($1,$3);
	} elsif( /^lnL\s+\=\s*(\-?\d+(\.\d+)?)/ ) {
	    $log = $1;
	} elsif( m/^t\=\s*(\d+(\.\d+)?)\s+
		 S\=\s*(\d+(\.\d+)?)\s+
		 N\=\s*(\d+(\.\d+)?)\s+
		 dN\/dS\=\s*(\d+(\.\d+)?)\s+
		 dN\=\s*(\d+(\.\d+)?)\s+
		 dS\=\s*(\d+(\.\d+)?)/ox ) {
	    $result[$b-1]->[$a-1] = { 
		'lnL' => $log,
		't' => $1,
		'S' => $3,
		'N' => $5,
		'omega' => $7,
		'dN' => $9,
		'dS' => $11 };
	} elsif( /^\s+$/ ) { 
	    next; 
	} elsif( /^\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)/ ) {
	} else { 
	    $self->debug( "unknown line: $_");
	}
    }
    return ( -mlmatrix => \@result);
}

sub _parse_YN_Pairwise {
    my ($self) = @_;
    my @result;
    while( defined( $_ = $self->_readline) ) {
	last if( /^seq\.\s+seq\./);
    }
    while( defined( $_ = $self->_readline) ) {
	if( m/^\s+(\d+)\s+  # seq #
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
	    ) {
	    
	    $result[$2-1]->[$1-1] = { 
		'S' => $3,
		'N' => $5,
		't' => $7,
		'kappa' => $9,
		'omega' => $11,
		'dN' => $13,
		'dN_SE' => $15,
		'dS' => $17,
		'dS_SE' => $19,
	    };
	} elsif( /^\s+$/ ) { 
	    next; 
	}
    }
    return ( -mlmatrix => \@result);
}

sub _parse_Forestry {
    my ($self) = @_;
    my ($instancecount,$loglikelihood,$score,$done,$treelength) = (0,0,0,0,0);
    my @trees;
    my $okay = 0;
    while( defined ($_ = $self->_readline) ) {
	last if $done;	
	if( /^TREE/ ) {
	    ($score) = (/MP\s+score\:\s+(\S+)/ );
	} elsif( /^Node\s+\&/ || /^\s+N37/ || /^(CODONML|AAML|YN00|BASEML)/ ||
		 /^\*\*/ || /^Detailed output identifying parameters/) {
	    $self->_pushback($_);
	    $done = 1;
	    last;
	} elsif( /^tree\s+length\s+\=\s+(\S+)/ ) {
	    $treelength = $1; # not going to store this for now
            # as it is directly calculated from
	    # $tree->total_branch_length;
	} elsif( /^\s*lnL\(.+\)\:\s+(\S+)/ ) {
	    $loglikelihood = $1;
	} elsif( /^\(/) {
	    if( $okay > 0 ) {
		my $treestr = new IO::String($_);
		my $treeio = new Bio::TreeIO(-fh => $treestr,
					     -format => 'newick');
		my $tree = $treeio->next_tree;
		if( $tree ) {
		    $tree->score($loglikelihood);
		    push @trees, $tree;
		}
		last;
	    }
	    $okay++;
	} 
    }
    return @trees;
}

sub _parse_NSsitesBatch {
    my $self = shift;
    my %data;
    my ($okay,$done) =(0,0);
    while( defined($_ = $self->_readline) ) {
	last if $done;
	next if /^\s+$/;
	
	next unless( $okay || /^Model\s+\d+/ );
	if( /^Model\s+(\d+)/o ) {
	    if( $okay ) {
		# this only happens if $okay was already 1 and 
		# we hit a Model line
		$self->_pushback($_);
		$done = 1;
	    } else {
		chomp;
		$data{'-model_num'}        = $1;
		($data{'-model_description'}) = ( /\:\s+(.+)/ );
		$okay = 1;
	    }
	} elsif( /^Time used\:\s+(\S+)/o ) {
	    $data{'-time_used'} = $1;
	    $done = 1;
	} elsif( /^kappa\s+\(ts\/tv\)\s+\=\s+(\S+)/ ) { 	    
	    $data{'-kappa'} = $1;
	} elsif( /^TREE/ ) {
	    $data{'-trees'} = [$self->_parse_Forestry];
	    if( defined $data{'-trees'} && 
		scalar @{$data{'-trees'}} ) {
		$data{'-likelihood'}= $data{'-trees'}->[0]->score;
	    }
	} elsif( /^Positively selected sites/ ) {
	    $self->_pushback($_);
	    my @sites = $self->_parse_Pos_selected_sites;
	    $data{'-pos_sites'} = [@sites];
	} elsif( /^dN/oi ) {
	    if( /K\=(\d+)/ ) {
		$data{'-num_site_classes'} = $1;   
		my @p = split(/\s+/,$self->_readline);
		my @w = split(/\s+/,$self->_readline);
		shift @p;
		shift @w;
		$data{'-dnds_site_classes'} = { 'p' => \@p,
						'w' => \@w};
	    } elsif( /for each branch/ ) {
		my %branch_dnds = $self->_parse_branch_dnds;
		# These need to be added to the Node/branches
		$self->_pushback($_);
	    }
	} elsif( /^Parameters in beta:/ ) {
	    $_ = $self->_readline; # need the next line
	    if ( /p\=\s+(\S+)\s+q\=\s+(\S+)/ ) {
		$data{'-shape_params'} = { 
		    'shape' => 'beta',
		    'p'     => $1,
		    'q'     => $2 };
	    } else {
		$self->warn("unparseable beta paramters: $_");
	    }
	}  elsif( /^alpha\s+\(gamma\)\s+\=\s+(\S+)/ ) {
	    my $gamma = $1;
	    $_ = $self->_readline;
	    my (@r,@f);
	    if( s/^r\s+\(\s*\d+\)\:\s+// ) {
		@r = split;
	    }
	    $_ = $self->_readline;
	    if( s/^f\s*\:\s+// ) {
		@f = split;
	    }
	    $data{'-shape_params'} = { 
		'shape' => 'alpha',
		'gamma' => $gamma,
		'r'     => \@r,
		'f'     => \@f };
	}
    }
    return new Bio::Tools::Phylo::PAML::ModelResult(%data);
}

sub _parse_Pos_selected_sites {
    my $self = shift;
    my ($okay,$done) = (0,0);
    my @sites;
    while( defined($_ = $self->_readline) ) {

	next if ( /^\s+$/ );
	last if $done;
	next unless ( $okay || /^Positively selected sites Prob/);
	if( /^Positively selected sites Prob/ ) {
	    $okay = 1;
	} elsif( /^\s+(\d+)\s+(\S)\s+([\d\.\-\+]+)(\**)/ ) {
	    my $signif = $4; 
	    $signif = '' unless defined $signif;
	    push @sites, [$1,$2,$3,$signif];
	} else {
	    $self->_pushback($_);
	    last;
	}
    }
    return @sites;
}

sub _parse_branch_dnds { 
    my $self = shift;
    my ($okay) = (0);
    my %branch_dnds;
    while(defined($_ = $self->_readline ) ) {
	next if( /^\s+$/);
	next unless ( $okay || /^\s+branch\s+t/);
	if( /^\s+branch\s+t/ ) {
	    $okay = 1;
	} elsif( /^\s*(\d+)\.\.(\d+)/ ) {
	    my ($branch,$t,$S,$N,$dNdS,$dN,$dS,$S_dS,$N_dN) = split;
	    $branch_dnds{$branch} = [$branch,$t,$S,$N,$dNdS,$dN,
				     $dS,$S_dS,$N_dN];
	} else { 
	    $self->_pushback($_);
	    last;
	}
    }
    return %branch_dnds;
}
1;
