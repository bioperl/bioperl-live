
#
# BioPerl module for Bio::Search::Processor::Fasta
#
# Cared for by Aaron Mackey <amackey@virginia.edu>
#
# Copyright Aaron Mackey
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Search::Processor::Fasta - Processor of Fasta-generated data streams

=head1 SYNOPSIS

use Bio::Search::Processor

my $processor = new Bio::Search::Processor -file      => 'mysearchrun',
                                           -algorithm => 'Fasta';

=head1 DESCRIPTION

A Processor object is used to generate Bio::Search::Result::* objects, given
a source of Search data (a file or filehandle).  The Processor object works
very much like the SeqIO system: once initialized with the new() method, the
Processor object will continue to yield as many Result objects as are available
from the data source (for single "runs" this is often only one Result object).

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this
and other Bioperl modules. Send your comments and suggestions preferably
 to one of the Bioperl mailing lists.
Your participation is much appreciated.

  vsns-bcd-perl@lists.uni-bielefeld.de          - General discussion
  vsns-bcd-perl-guts@lists.uni-bielefeld.de     - Technically-oriented discussion
  http://bio.perl.org/MailList.html             - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
 the bugs and their resolution.
 Bug reports can be submitted via email or the web:

  bioperl-bugs@bio.perl.org
  http://bio.perl.org/bioperl-bugs/

=head1 AUTHOR - Aaron Mackey

Email amackey@virginia.edu

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Search::Processor::Fasta;
use vars qw($AUTOLOAD @ISA);
use strict;

# Object preamble - inherits from Bio::Root::Object

# "isa" modules
use Bio::Search::Processor::ProcessorI;

# "uses" modules
use Bio::Search::Result::Fasta;
use Bio::Search::Hit::Fasta; # will be contained in Result::Fasta
use Bio::SimpleAlign; # will be contained in Hit::Fasta

@ISA = qw(Bio::Search::Processor::ProcessorI Exporter);

my @AUTOLOAD_OK = qw(
		     _FILEHANDLE
		     _ALGORITHM
                     _SFSEP
		    );
my %AUTOLOAD_OK;

@AUTOLOAD_OK{@AUTOLOAD_OK} = 1;

# new() is inherited from Bio::Root::Object

# _initialize is where the heavy stuff will happen when new is called

sub _initialize {
    my($self,@args) = @_;

    my $make = $self->SUPER::_initialize();

    my ($fh, $file, $algorithm, $sfsep) =
	$self->_rearrange([qw(FH FILE ALGORITHM SFSEP)], @args);

    $fh = IO::File->new($file, "r") unless defined $fh;

    $self->_FILEHANDLE($fh);
    $self->_ALGORITHM($algorithm);
    $self->_SFSEP($sfsep || '|');

    return $make; # success - we hope!
}

=head2 next_result

 Title   : next_result
 Usage   : $result = $processor->next_result()
 Function: Used to obtain Result objects from a FASTA-generated data source
 Returns : a Bio::Search::Result::Fasta object
 Args    : <none>


=cut

sub next_result{

    my ($self,@args) = @_;

    my $fh = $self->_FILEHANDLE() or $self->throw("No data source!");

    my $result = new Bio::Search::Result::Fasta;

    return undef if eof $fh;

    # process data from $fh and set appropriate values in $result using
    # $result's publically available accessors, including creating and
    # storing Hit objects which themselves have created and contained
    # Alignment objects.  Sky's the limit here.

    local($/) = "\n";  # just to make sure ...
    my $interactive = 1;

    my $line = <$fh>; #" FASTA searches a protein or DNA sequence data bank"
    my ($algorithm_desc, $algorithm) = $line =~ m/^\s*((T?FAST\W).*?)\s*$/;

    $line = <$fh>; #" version 3.2t05  May 12, 1999"
    my ($version, $versiondate) = $line =~ m/^\s*version\s+(\S+)\s+(.*?)\s*$/;

    $line = <$fh>; # "Please cite:"
    $line = <$fh>; #" W.R. Pearson & D.J. Lipman PNAS (1988) 85:2444-2448"
    my ($citation) = $line =~ m/^\s*(.*)\s*$/;

    $line = <$fh>; # <blank line>
    $line = <$fh>; #" query.fa: 675 aa"
    my ($query_filename, $query_start, $query_end, $query_size, $query_type) =
	$line =~ m/^\s*(.*?):(?:(\d+)-(\d+):)?\s+(\d+)\s+(aa|nt)\s*$/;

    $line = <$fh>; #" >104K_THEPA |  12 34 | 104 KD MICRONEME-RHOPTRY ANTIGEN."
                   # ( bogus sfnum's btw )
    my $sfsep = $self->_SFSEP(); # default is '|', but could be anything!
    my ($query_id, $query_superfamilies, $query_desc) =
	$line =~ m/^\s*>(\S+)\s+         # id
                  (\Q$sfsep\E\s+         # opening sfsep, gotta metaquote it
                     [^\Q$sfsep\E]*\s*   # superfamily numbers
                   \Q$sfsep\E)?          # closing sfsep
                  \s*(.*?)\s*$/x;        # description
	$query_superfamilies = [ grep { length }
				 split(/\s+/, $query_superfamilies) ] ;

    $line = <$fh>; #" vs  /wrp_lib/spsinglecdspf1.fa library"
    my ($library_filename) = $line =~ m/^\s*vs\s+(.*)\s+library\s*$/;

    $line = <$fh>; #"searching /wrp_lib/spsinglecdspf1.fa library"
                   # no new information, so skip it.

    $line = <$fh>; # <blank line> or "..... ..... ..... Done!" if interactive
    if ( substr($line, 0, 1) eq '.' ) {
        $interactive = 1;
        until ($line =~ m/^[\.\s]*Done!\s*$/) {
	    $line = <$fh>;
        }
        $line = <$fh>; # now we get our <blank line>
    }

    $line = <$fh>; # histogram or "10821595 residues in 23981 sequences"
    my $histogram = '';
    if ($line =~ m/^\s+opt\s+E\(\)\s*$/) { # histogram, bleah
        $histogram .= $line;
        until ($line =~ m/^>/) {
	    $line = <$fh>;
	    $histogram .= $line;
        }
	$line = <$fh>; # now we get "10821595 residues in 23981 sequences"
    }
    my ($library_residues, $library_sequences) =
	$line =~ m/^\s*([\d,]+)\s+residues\s+in\s+(\d+)\s+sequences?\s*$/;

    $line = <$fh>; # First line of statistics, could be one of many:

    my $statistics = {};

 SWITCH: {

# <blank line>
	$line =~ m/^\s*$/ and do {
	    $statistics->{'description'} = 'no statistical estimates';
	  # We can't know $shuffled or not, FASTA output doesn't tell us.
	  # This is a bug with FASTA, not Bioperl
	    last SWITCH;
	};

#  unscaled statistics: mu= 46.7391  var=313.9553
# <blank line>
	$line =~ m/unscaled statistics/ and do {
	    $statistics->{'description'} = 'unscaled statistical estimates';
	    $statistics->{'shuffled'} = $line =~ m/^\s*(\(shuffled\))/; # scalar context here!
	    ($$statistics{'unscaled_mu'},
	     $$statistics{'unscaled_var'}) =
		 $line =~ m/mu=\s*
                          (-?\d*\.?\d*)\s*var=\s*
                          (-?\d*\.?\d*)/x;
	    if (length $histogram) {
		$line = <$fh>;
		($$statistics{'unscaled_ks'},
		 $$statistics{'unscaled_ks_N'},
		 $$statistics{'unscaled_ks_at'}) =
		     $line =~ m/^\s*Kolmogorov-Smirnov\s*statistic:\s*
                               (-?\d*\.?\d*)\s*\(N=(\d+)\)\s*at\s*(\d+)\s*$/x;
	    }
	    $line = <$fh>;
	    last SWITCH;
	};

#  Expectation_n fit: rho(ln(x))= 7.4444+/-0.00103; mu= 2.1186+/- 0.061;
# mean_var=211.1380+/-42.955, 0's: 0 Z-trim: 121  B-trim: 88 in 1/55
# Kolmogorov-Smirnov  statistic: 0.0587 (N=29) at  50
# <blank line>
	$line =~ m/Expectation_n\s+fit/ and do {
	    $statistics->{'description'} = 'regression-scaled estimates';
	    $statistics->{'shuffled'} = $line =~ m/^\s*(\(shuffled\))/; # scalar context here!
	    ($$statistics{'reg_rho'},
	     $$statistics{'reg_rho_sd'},
	     $$statistics{'reg_mu'},
	     $$statistics{'reg_mu_sd'}) =
		 $line =~ m/rho\(ln\(x\)\)=\s*
                         (-?\d*\.?\d*)\s*\+\/\-\s*
                         (\d*\.?\d*);\s*
                         mu=\s*
                         (-?\d*\.?\d*)\s*\+\/\-\s*
                         (\d*\.?\d*);\s*$/x;
	    $line = <$fh>;
	    ($$statistics{'reg_mean_var'},
	     $$statistics{'reg_mean_var_sd'},
	     $$statistics{'reg_zeros'},
	     $$statistics{'reg_ztrim'},
	     $$statistics{'reg_btrimmed_seqs'},
	     $$statistics{'reg_btrimmed_bins'},
	     $$statistics{'reg_total_bins'}) = 
		 $line =~ m/^\s*mean_var=\s*
                          (-?\d*\.?\d*)\s*\+\/\-\s*
                          (\d*\.?\d*),\s*0's:\s*(\d+)\s*Z-trim:\s*
                          (\d+)\s*B-trim:\s*(\d+)\s*in\s*(\d+)\/(\d+)\s*$/x;
	    if (length $histogram) {
		$line = <$fh>;
		($$statistics{'reg_ks'},
		 $$statistics{'reg_ks_N'},
		 $$statistics{'reg_ks_at'}) =
		     $line =~ m/^\s*Kolmogorov-Smirnov\s*statistic:\s*
                               (-?\d*\.?\d*)\s*\(N=(\d+)\)\s*at\s*(\d+)\s*$/x;
	   }
	    $line = <$fh>;
	    last SWITCH;
	};

# Ln statistics: mu= 41.7673  var=222.6012
# <blank line>
	$line =~ m/Ln statistics/ and do {
	    $statistics->{'description'} = 'log-corrected estimates';
	   # We can't know $shuffled or not, FASTA output doesn't tell us.
	   # This is a bug with FASTA, not Bioperl
	    ($$statistics{'log_mu'},
	     $$statistics{'log_var'}) =
		 $line =~ m/^\s*Ln\s+statistics:\s*mu=\s*
                          (-?\d*\.?\d*)\s*var=\s*
                          (-?\d*\.?\d*)/x;
	    if (length $histogram) {
		$line = <$fh>;
		($$statistics{'log_ks'},
		 $$statistics{'log_ks_N'},
		 $$statistics{'log_ks_at'}) =
		     $line =~ m/^\s*Kolmogorov-Smirnov\s*statistic:\s*
                               (-?\d*\.?\d*)\s*\(N=(\d+)\)\s*at\s*(\d+)\s*$/x;
	   }
	    $line = <$fh>;
	    last SWITCH;
	};
# Altschul/Gish params: n0: 675 Lambda: 0.158 K: 0.019 H: 0.100
# <blank line>
	$line =~ m/Altschul\/Gish/ and do {
	    $statistics->{'description'} = 'Altschul/Gish estimates';
	  # We can't know $shuffled or not, FASTA output doesn't tell us.
	  # This is a bug with FASTA, not Bioperl
	    ($$statistics{'ag_n0'},
	     $$statistics{'ag_lambda'},
	     $$statistics{'ag_K'},
	     $$statistics{'ag_H'}) =
		 $line =~ m/n0:\s*(\d+)\s*Lambda:\s*
                          (-?\d*\.?\d*)\s*K:\s*
                          (-?\d*\.?\d*)\s*H:\s*
                          (-?\d*\.?\d*)\s*$/x;
	    if (length $histogram) {
		$line = <$fh>;
		($$statistics{'ag_ks'},
		 $$statistics{'ag_ks_N'},
		 $$statistics{'ag_ks_at'}) =
		     $line =~ m/^\s*Kolmogorov-Smirnov\s*statistic:\s*
                               (-?\d*\.?\d*)\s*\(N=(\d+)\)\s*at\s*(\d+)\s*$/x;
	    }
	    $line = <$fh>;
	    last SWITCH;
	};

#  Expectation_i fit: rho(ln(x))= 7.5056+/-0.00104; mu= 1.9654+/- 0.061;
# mean_var=235.3261+/-46.802 0's: 0 Z-trim: 60 N-it: 1
# <blank line>
	$line =~ m/Expectation_i fit/ and do {
	    $statistics->{'description'} = 'regression-scaled estimates, -z 4';
	    $statistics->{'shuffled'} = $line =~ m/^\s*(\(shuffled\))/; # scalar context here!
	    ($$statistics{'reg4_rho'},
	     $$statistics{'reg4_rho_sd'},
	     $$statistics{'reg4_mu'},
	     $$statistics{'reg4_mu_sd'}) =
		 $line =~ m/rho\(ln\(x\)\)=\s*
                         (-?\d*\.?\d*)\s*\+\/\-\s*
                         (\d*\.?\d*);\s*
                         mu=\s*
                         (-?\d*\.?\d*)\s*\+\/\-\s*
                         (\d*\.?\d*);\s*$/x;
	    $line = <$fh>;
	    ($$statistics{'reg4_mean_var'},
	     $$statistics{'reg4_mean_var_sd'},
	     $$statistics{'reg4_zeros'},
	     $$statistics{'reg4_ztrim'},
	     $$statistics{'reg4_Nit'}) = 
		 $line =~ m/^\s*mean_var=\s*
                          (-?\d*\.?\d*)\s*\+\/\-\s*
                          (\d*\.?\d*),\s*0's:\s*(\d+)\s*Z-trim:\s*
                          (\d+)\s*N-it:\s*(\d+)\s*$/x;
	    if (length $histogram) {
		$line = <$fh>;
		($$statistics{'reg4_ks'},
		 $$statistics{'reg4_ks_N'},
		 $$statistics{'reg4_ks_at'}) =
		     $line =~ m/^\s*Kolmogorov-Smirnov\s*statistic:\s*
                               (-?\d*\.?\d*)\s*\(N=(\d+)\)\s*at\s*(\d+)\s*$/x;
	    }
	    $line = <$fh>;
	    last SWITCH;
	};

#  Expectation_v fit: rho(ln(x))= 7.4444+/-0.00103; mu= 2.1186+/- 0.061;
# rho2=  8.01; mu2= -261.80, 0's: 0 Z-trim: 121  B-trim: 88 in 1/55
# <blank line>
	$line =~ m/Expectation_v fit/ and do {
	    $statistics->{'description'} = 'regression-scaled estimates, -z 5';
	    $statistics->{'shuffled'} = $line =~ m/^\s*(\(shuffled\))/; # scalar context here!
	    ($$statistics{'reg5_rho'},
	     $$statistics{'reg5_rho_sd'},
	     $$statistics{'reg5_mu'},
	     $$statistics{'reg5_mu_sd'}) =
		 $line =~ m/rho\(ln\(x\)\)=\s*
                         (-?\d*\.?\d*)\s*\+\/\-\s*
                         (\d*\.?\d*);\s*
                         mu=\s*
                         (-?\d*\.?\d*)\s*\+\/\-\s*
                         (\d*\.?\d*);\s*$/x;
	    $line = <$fh>;
	    ($$statistics{'reg5_rho2'},
	     $$statistics{'reg5_rho2_sd'},
	     $$statistics{'reg5_zeros'},
	     $$statistics{'reg5_ztrim'},
	     $$statistics{'reg5_btrimmed_seqs'},
	     $$statistics{'reg5_btrimmed_bins'},
	     $$statistics{'reg5_total_bins'}) =
		 $line =~ m/^\s*mean_var=\s*
                          (-?\d*\.?\d*)\s*\+\/\-\s*
                          (\d*\.?\d*),\s*0's:\s*(\d+)\s*Z-trim:\s*
                          (\d+)\s*B-trim:\s*(\d+)\s*in\s*(\d+)\/(\d+)\s*$/x;
	    if (length $histogram) {
		$line = <$fh>;
		($$statistics{'reg5_ks'},
		 $$statistics{'reg5_ks_N'},
		 $$statistics{'reg5_ks_at'}) =
		     $line =~ m/^\s*Kolmogorov-Smirnov\s*statistic:\s*
                               (-?\d*\.?\d*)\s*\(N=(\d+)\)\s*at\s*(\d+)\s*$/x;
	    }
	    $line = <$fh>;
	    last SWITCH;
	};
    }

    $line = <$fh>;
#FASTA (3.2 December, 1998) function [optimized, BL50 matrix (15:-5)] ktup: 2
    my ($optimized,
	$matrix_name,
	$matrix_offset,
	$matrix_high_score,
	$matrix_low_score,
	$ktup) =
	    $line =~ m/\[(optimized,\s*)?
                       ([^\/]+)\/?(-?\d+)?\s*matrix\s*
                       \((-?\d+):(-?\d+)\]\s*
                       ktup:\s*(\d+)/x;

    $line = <$fh>;
# join: 38, opt: 26, gap-pen: -12/ -2, width:  16 reg.-scaled
    my ($join, $opt, $gap_open, $gap_extend, $width) =
	$line =~ m/^\s*join:\s*(\d+),\s*
                   opt:\s*(\d+),\s*
                   gap-pen:\s*(-?\d+)\s*\/\s*(-?\d+),\s*
                   width:\s*(\d+)/x;

    $line = <$fh>;
# Scan time:  0.390
    my ($scantime) = $line =~ m/(\d*\.?\d*)/;

    $line = <$fh>;
#The best  scores are:                             initn init1 opt z-sc E(23874)
    my ($app_lib_size) = $line =~ m/E\(\s*(\d+)\s*\)/;

#104K_THEPA 104 KD MICRONEME-RHOPTRY ANTI  ( 924) 4590 4590 4590 3639.4 1.6e-196
#SSGP_VOLCA SULFATED SURFACE GLYCOPROTEIN  ( 485)  184  184  247  162.3  0.0074
#EXTN_MAIZE EXTENSIN PRECURSOR (PROLINE-R  ( 267)  179  179  217  161.6  0.0081
#PRPM_HUMAN SALIVARY PROLINE-RICH PROTEIN  ( 234)  175   96  213  159.2   0.011
#PRP4_HUMAN SALIVARY PROLINE-RICH PROTEIN  ( 247)  173   94  205  157.3   0.014

    my @hits = ();
    while (1) {
	$line = <$fh>;
	if ($line =~ m/^\s*$/) {
	    last;
	} else {
	    my ($desc, $size, $initn, $init1, $opt, $zsc, $e_val) =
		$line =~ m/^(.{40})\s*
                           \(\s*(\d+)\s*\)\s*
                           (\d+)\s*(\d+)\s*(\d+)\s*
                           (?:(\d*\.?\d*)\s*(\d*[\.e]?\d*e?-?\d*)\s*)?$/x;
	    $desc =~ m/^(\S+)\s*(.*)\s*$/;
	    my $id = $1; $desc = $2;
	    my $hit = Bio::Search::Hit::Fasta->new('id' => $id,
						   'desc' => $desc,
						   'size' => $size,
						   'initn' => $initn,
						   'init1' => $init1,
						   'opt' => $opt,
						   'zsc' => $zsc,
						   'e_val' => $e_val);
	    push @hits, $hit;
	}
    }

    return $result;
}

sub AUTOLOAD {
    my ($self, $val) = @_;

    if ($AUTOLOAD_OK{$AUTOLOAD}) {
	$self->{$AUTOLOAD} = $val if defined $val;
	return $self->{$AUTOLOAD};
    } else {
	$self->throw("Unallowed accessor: $AUTOLOAD");
    }
}

1;

__END__

