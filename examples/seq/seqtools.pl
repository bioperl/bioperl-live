BEGIN {
        use vars qw($INSTALL_PATH);

        #### NOTE: $INSTALL_PATH needs to be customized
        ###         your system. It should point
        ###         to the directory right above Bio/ in order
        ###         for perl to be able to locate the .pm files. 

	$INSTALL_PATH = "/home/steve/perl/lib";

        ###
        ####
    }

#---------------------------------------------------------------------------
# LIBRARY      : seqtools.pl
# PURPOSE      : Provides a set of standard functions & variables for 
#                working with sets of sequences using Bio::SeqIO modules
#                and Bio::Seq objects.
#            Sequence data can be input and output in a variety of formats 
#            (raw embl fasta gcg genbank pir scf swiss).
#            See the Bio/SeqIO directory for the full list.
# AUTHOR       : Steve Chervitz <sac@bioper.org>
# CREATED      : 10 Apr 1998
# REVISION     : $Id$
# INSTALLATION : Edit $INSTALL_PATH to point to the directory containing
#                your Bioperl modules (the Bio/ directory).
#                Some operations require the Bio::Seq.pm module.
# USAGE:
#
#  #!/usr/bin/perl -w
#
#  ## A minimal script that uses this package:
#  ## Adjust path according to your system:
#  require "/home/me/bioperl/examples/seq/seqtools.pl"; 
#
#  &init_seq();
#  &load_ids();
#  &get_seq_objs();
#  &print_seqs();
#  &wrap_up();
#
# EXAMPLES : See the files seqs.pl, seqs2.pl, and seqs3.pl
#            in this dir for some working examples.
#
#
# PUBLIC FUNCTIONS DEFINED IN seqtools.pl:
#
#    init_seq()
#    load_ids()
#    load_seqs()
#    print_seq()
#    print_seqs()
#    write_file()
#    write_files()
#    wrap_up_seq()
#
#
# MODIFIED:
#  sac --- Wed Feb 23 00:32:20 2000
#      * Converted to use Bio::SeqIO instead of Bio::Tools::Fasta.pm
#  sac --- Thu Feb  4 06:48:42 1999
#      * Added the -wait command-line option to prevent read() timeouts when
#        using seqtools.pl in conjunction with blast_seq.pl.
#  sac -- Tue Nov 24 10:55:46 1998
#      * Usage text added for -incl/-excl options.
#  0.21, sac,  7 Sep 1998: 
#      * Added the -tag command-line option.
#        Allowed -incl and -excl options to accept a list on command-line
#        instead of requiring that they be specified in a file.
#  0.2,  sac,  4 Aug 1998: 
#      * Added print_params() function.
#      * Exits with non-zero status if any of the checks in wrap_up_seq() fail.
#      * Added private, up-cased versions of @ids_incl, %ids_incl, and xxx_excl
#        Private versions have '_uc' attached to their names.
#      * The $fh output FileHandle object variable is now public and can be 
#        used by scripts that require seqtools.pl.
#  0.12, sac, 17 Jul 1998: Added -col option and write_files() function.
#  0.1,  sac, 16 Jun 1998: Added alternate $INSTALL_PATH for testing purposes.
#
#---------------------------------------------------------------------------

use lib $INSTALL_PATH;
use lib '.','..';  # fail-safe incase you forget to edit $INSTALL_PATH.

use Bio::SeqIO;
use Bio::Root::Global qw(:devel $TIMEOUT_SECS);
use Getopt::Long;
use Carp;

select(STDOUT); $|=1;

my @SUPPORTED_FORMATS = qw(raw embl fasta gcg genbank pir scf swiss);

$VERSION    = 0.3;
$DESC       = "Provides standard functions and variables for working with sets\n".
              "of sequences and sequence objects using Bio::SeqIO modules.";
%seqParams  = ();  
$seqcount   = 0;
@seqs       = ();
@ids        = ();
@descs      = ();
$_count_processed = 0;
@errs    = ();

# Command-line options:
#$opt_strict     = 0;        # Verify seq data based on standard alphabet. 
$opt_nucl       = 0;
$opt_prot       = 0;
$opt_col        = 0;         # Column in the -incl or -excl files with the seq IDs.
$opt_exact      = 1;         # Require exact match when filtering based on seq id.
$opt_fmt        = 'fasta';   # Default input format
$opt_outfmt     = 'fasta';   # Default output format
$opt_seq        = undef;
$opt_incl       = undef;
$opt_excl       = undef;
$opt_out        = undef;
$opt_err        = undef;
$opt_tag        = undef;  # Optional tag string to be prepended to all descrptions.
$opt_write_files= undef;  # Directory name when writing individual files 
                          # for each input sequence. 
$opt_wait = $TIMEOUT_SECS;  # Seconds to wait for input before timing out. 

# General parameters
$opt_h      = 0;
$opt_eg     = 0;
$opt_mon    = 1;
$opt_debug  = 0;
@argv       = @ARGV;   # copy ARGV before GetOptions() massacres it.

@ids_incl    = ();
@ids_excl    = ();
%ids_incl    = ();
%ids_excl    = ();

# Private:
my (%not_filtered, $alphabet, %_tested, @_tested, $_dir, $seqout);
my (@ids_incl_uc, @ids_excl_uc, %ids_incl_uc, %ids_excl_uc);


#--------------------
sub print_seq_usage {
#--------------------
    &seq_usage;
    print STDERR "<RETURN> to view parameters."; <STDIN>;
    &seq_params;
    &seq_general_params;
}

#---------------
sub seq_usage {
#---------------
# Basic usage for a script that uses this library.

    print STDERR "$0, $VERSION\n$DESC\n";
    print STDERR <<"QQ_USAGE_QQ";

Usage: $0 [ parameters ] seq/seqs.fasta > out
       $0 [ parameters ] seq/*.fasta   > out
       gzcat *fasta.gz | $0 [ parameters ] > out

 input : Multiple Fasta-formatted sequences are read via STDIN
         or from a file(s) specified on the command line.
         
 Use the -eg option to see some examples.

QQ_USAGE_QQ

}

#---------------
sub seq_params {
#---------------

    print STDERR <<"QQ_PARAMS_QQ";

 SEQTOOLS PARAMETERS:
 --------------------
 -fmt <format> : Sequence format for readin sequence data (default = $opt_fmt).
                 Supported: @SUPPORTED_FORMATS
 -outfmt <format> : Sequence format for writing sequence data (default = $opt_outfmt).
                   Supported: @SUPPORTED_FORMATS
 -incl <file|list> : Filename containing list of sequence IDs to include
                     -or- list of IDs separated by whitespace.
                     For lines containing multiple ids, only
                     the first one is used.
 -excl <file|list> : Filename containing list of sequence IDs to exclude
                     -or- list of IDs separated by whitespace.
 -seq  file    : Filename containing Fasta formated sequences (if not
                 using STDIN). ('-seq' is optional but should be used
	         with scripts that accept different types of input files).
 -nucl         : Nucleotide sequences  (optional).
 -prot         : Peptide sequences (optional).
 -out <file>   : Filename for saving output (default = STDOUT).
 -tag <string> : Prepend the string to the descriptions of all sequences.
 -col <int>    : Column in -incl or -excl file containing sequence IDs. 
                 Default = 0 (left-most column).
 -noexact      : When using the -incl or -excl parameters, using -noexact
                 requires only that an ID in the Fasta file contains
                 an ID in the supplied list, not an exact match.
	         (useful for screening a GenBank Fasta file with a list of GIs;
	         the Fasta file contains the GI plus other identifiers).
	         Default = case-insensitive comparison.
 -write_files <dir> 
               : Directory for writing individual files for input sequences.
 -wait <int>   : Amount of seconds to wait before timing out when reading in 
                 sequence data (default = $opt_wait seconds).

QQ_PARAMS_QQ
}

#-----------------------
sub seq_general_params {
#-----------------------

    print STDERR <<"QQ_GENERAL_QQ";

 SEQTOOLS GENERAL PARAMETERS:
 ----------------------------
 -nomon     : Monitor off (does not print additional progress data)
 -err <file>: Save errors to the indicated file (append mode).
 -debug     : Debugging on.
 -h         : Print this help/usage info.
 -eg        : Print examples.

QQ_GENERAL_QQ
}

#---------------
sub init_seq {
#---------------
    ($usage_fref, @opts) = @_;

    &GetOptions('mon!', 'h!', 'debug!', 'incl=s', 'excl=s', 'wait=s', 
		'seq=s', 'nucl!', 'prot!', 'strict!', 'out=s',
		'exact!', 'outfmt=s', 'err=s', 'eg!', 'write_files=s', 'tag=s',
		@opts);
   
    $MONITOR && print STDERR "$0, v$VERSION\n",'-'x50,"\n";

    $opt_h and do{
	(ref($usage_fref) =~ /CODE/) ? &$usage_fref : &print_seq_usage;
	exit 1;
    };

    $opt_eg and do{
	eval{ print STDERR "\n$0 EXAMPLES:\n",'-'x30,"\n\n", &examples(); };
	$@ and print STDERR "\nSorry. No examples available.\n";
	exit 1;
    };

    monitor($opt_mon);
    debug($opt_debug);
    
    if(not $opt_seq and scalar(@ARGV) == 1) {
	# Mainly doing this to record the sequence file name.
	$opt_seq = $ARGV[0];
    }

    if($opt_err) {
	system('touch', $opt_err);
	open (STDERR, ">>$opt_err") or croak "*** $0: Can't open err file $opt_err: $!\n\n";
    }

    $alphabet = 'protein' if $opt_prot;
    $alphabet = 'dna' if $opt_nucl;

    if($DEBUG) {
	$alphabet or print STDERR "\n*** Sequence alphabet not specified (can use -nucl or -prot on command line)\n\n";
    }
    
    if(! grep /$opt_fmt/i, @SUPPORTED_FORMATS) {
	print STDERR "\n*** Input format not supported: $opt_fmt\n";
	print STDERR "\tValid formats: @SUPPORTED_FORMATS\n\n";
    }
    if(! grep /$opt_outfmt/i, @SUPPORTED_FORMATS) {
	print STDERR "\n*** Output format not supported: $opt_outfmt\n";
	print STDERR "\tValid formats: @SUPPORTED_FORMATS\n\n";
    }

    local($^W) = 0;

    $MONITOR && print STDERR "Output format: $opt_outfmt\n";

    if($opt_out) {
	$MONITOR and print STDERR "\nWriting sequence data to $opt_out\n\n";
    } 
    if($opt_write_files) {
	$MONITOR and print STDERR "\nWriting sequence files to $opt_write_files\n\n";
    } 

}

#----------------
sub print_params {
#----------------
    print STDERR "\nSEQTOOLS PARAMETERS:\n";
    print STDERR "-----------------------\n";
    printf STDERR "%-20s %s\n", 'incl', $opt_incl || '';
    printf STDERR "%-20s %s\n", 'excl', $opt_excl || '';
    printf STDERR "%-20s %s\n", 'seq', $opt_seq;
    printf STDERR "%-20s %s\n", 'prot/nucl', $opt_prot ? 'prot' : ($opt_nucl ? 'nucl' : '');
    printf STDERR "%-20s %s\n", 'outfmt', $opt_outfmt || '';
#    printf STDERR "%-20s %s\n", 'strict', $opt_strict;
    printf STDERR "%-20s %s\n", 'exact', $opt_exact;
    printf STDERR "%-20s %s\n", 'out', $opt_out || '';
    printf STDERR "%-20s %s\n", 'write_files', $opt_write_files || '';
}

#---------------
sub load_ids {
#---------------
    if($opt_incl) {
	&_load_list('incl');
#	  print STDERR "SOME IDS in hash ids_incl_uc:\n";
#	  my $i = 0;
#	  foreach(keys %ids_incl_uc) {
#	      print STDERR "  $_\n";
#	      last if $i++ == 10;
#	  }
	
    } else {
#	$MONITOR && print STDERR "\nNo included IDs.\n";
    }

    if($opt_excl) {
	&_load_list('excl');
    } else {
#	$MONITOR && print STDERR "\nNo excluded IDs.\n";
    }
}

#--------------
sub _load_list {
#--------------
# Modifies the sequence identifiers:
#  - Replaces '|' with '_' to avoid regexp problems
#  - Uppercases.

    my $type = shift;

    local(*file) = "opt_$type";
    local(*hash) = "ids_$type"; # Not used.

    # The following test of the *hash glob fails. Not sure why.
#    $hash = "ids_$type";
#    $MONITOR && print STDERR "\nLOADING LIST:\n";
#    foreach(keys %$hash) { print "  $_\n"; }
    
    my $count = 0;
    if(not -e $file) {
	# IDs were supplied on command line, not in file.
	my @ids = split(/\s+/, $file);
	my @ids_uc = map uc($_), @ids;  # Sequence IDs are uppercased for faster filtering.
	if($type eq 'incl') {
	    map $ids_incl_uc{$_}++, @ids_uc;
	    map $ids_incl{$_}++, @ids;
	    @ids_incl_uc = @ids_uc;
	    @ids_incl    = @ids;
	} else {
	    map $ids_excl_uc{$_}++, @ids_uc;
	    map $ids_excl{$_}++, @ids;
	    @ids_excl_uc = @ids_uc;
	    @ids_excl = @ids;
	}
	$count = scalar @ids;
    } else {
	
	open( LIST, $file) || croak "\n\a*** $0: Can't open ID file: $file: $!\n\n";
	
	my @set = ();
	my ($id, $id_uc);
	while(<LIST>) {
	    chomp();
	    next unless m/^\w/;
	    @set = split(/\s+/, $_);
	    $id = $set[$opt_col];
	    $id =~ s/\|/_/g;   # avoid regexp trouble.
	    $id_uc = uc($id);     # Sequence IDs are uppercased for filtering.
	    # Would like to use the *hash glob.
	    if($type eq 'incl') {
		$ids_incl_uc{$id_uc}++;
		push @ids_incl_uc, $id_uc;
		$ids_incl{$id}++;
		push @ids_incl, $id;
	    } else {
		$ids_excl_uc{$id_uc}++;
		push @ids_excl_uc, $id_uc;
		$ids_excl{$id}++;
		push @ids_excl, $id;
	    }
	    $count++;
	}
	close LIST;
    }
    $MONITOR && print STDERR "\n$count IDs to ${type}ude loaded.\n";
}


#----------------
sub load_seqs {
#----------------
    my $func_ref = shift;
    my %params;
    # Initialize input stream
    if($opt_seq) { 
        %params = ('-format' => $opt_fmt, '-file' => $opt_seq); 
    } else { 
        %params = ('-format' => $opt_fmt);  
    }

    my $seqin  = Bio::SeqIO->new(%params);

    # Initialize output stream
    if( $opt_out ) {
        %params = ('-format' => $opt_outfmt, '-file' => ">$opt_out");
    } else {
        %params = ('-format' => $opt_outfmt, '-fh' => \*STDOUT);
    }

    $seqout = Bio::SeqIO->new(%params);
    
    if(defined $alphabet) { $seqin->alphabet($alphabet); }

    $SIG{ALRM} = sub { die "Timed out!"; };

    my ($seq, $keep);
    eval {
        alarm($opt_wait);
        while ( $seq = $seqin->next_seq() ) {
	    alarm(0);  # Deactivate the alarm as soon as we start reading.
            $keep = 1;
            ## Should the sequence be filtered out?
            if($opt_incl or $opt_excl) {
                if(seq_filter($seq->length, $seq->id, $seq->desc)) {
                    $keep = 0;
                }
            }
            if(ref $func_ref eq "CODE") {
                &$func_ref($seq) if $keep;
            } else {
                push @seqs, $seq if $keep;
            }
        }
    };
    if($@ =~ /Timed out!/) {
	 croak "$0: Timed out while waiting for input.", " Timeout period = $opt_wait seconds.\nFor a longer time out period, supply a -wait <seconds> parameter\n";
    } elsif($@ =~ /\S/) {
         my $err = $@;
	 croak "$0: Unexpected error during read: $err";
    }

}


#---------------
sub seq_filter {
#---------------
# Returns true if the supplied sequence is to be filtered out.
    my($len, $id, $desc) = @_;
    my $filter_it = 0;
    my $id_uc = uc($id);  # Sequence IDs are uppercased  for filtering.

#    print "\nFILTERING $id_uc AGAINST:\n @ids_incl_uc[1..10]"; <STDIN>;

    if($opt_exact) {
	# Filter out ids NOT in the list of desireables.
	if ($opt_incl and not exists $ids_incl_uc{$id_uc}) {
	    $filter_it = 1;
	}
	# Filter out ids IN the list of undesireables.
	if ($opt_excl and exists $ids_excl_uc{$id_uc}) {
	    $filter_it = 1;
	}
    } else {
	# The hashes are faster but sometimes we need the 'contains' 
	# type of matching.
	$id_uc =~ s/\|/_/g;
	if ($opt_incl and not (exists $ids_incl_uc{$id_uc} or grep $id_uc =~ /$_/, @ids_incl_uc)) {
	    $filter_it = 1;
	} 
	if ($opt_excl and (exists $ids_excl_uc{$id_uc} or grep $id_uc =~ /$_/, @ids_excl_uc)) {
	    $filter_it = 1;
	}
    }

    $_tested{$id_uc}++;
    push @_tested, $id_uc;
    $not_filtered{$id_uc}++ if not $filter_it;
    
    $filter_it;
}
	

#---------------------
sub print_seq {
#---------------------
    my $seq = shift;

    # Prepend an optional tag to the sequence.
    $opt_tag and $seq->desc($opt_tag.' '.$seq->desc);

    $seqout->write_seq($seq);
    $_count_processed++;

}


#---------------------
sub print_seqs {
#---------------------

    foreach(@seqs) { 
	&print_seq($_);  
    }
}    


#----------------
sub write_file {
#----------------
# Write each sequence out to a separate file.
# The name of the file will be SEQ_ID.FORMAT (e.g., P34352.gcg).
# Any '|' in the ID are replaced by '_' since '|' can confuse the shell.

    my $seq = shift;

    &_set_dir unless $_dir;
 
    ($file = $seq->id) =~ s/\|/_/g;
    $file = $_dir . $file . substr(".\L$opt_outfmt\E", 0, 4);

    eval {
        my $out = Bio::SeqIO->new('-format' => $opt_outfmt, 
                                  '-file' => ">$file");

        $out->write_seq($seq);
        $_count_processed++;
        if(@ids_incl_uc) {
            return $count_processed >= scalar(@ids_incl_uc) ? 0 : 1;
        }
    };
    if($@) {
        my $err = $@;
        push @errs, $err;
    }

    1;

#    $MONITOR && print STDERR " wrote file: $file\n";
}    

#----------------
sub _set_dir {
#----------------
# Creates a new directory if the one specified by $opt_write_files
# does not already exist.
    $_dir = $opt_write_files || './';
    $_dir !~ /\/$/ and $_dir .= '/';
    system('mkdir', $_dir) unless -d $_dir;
}
    
#----------------
sub write_files {
#----------------
    $MONITOR && print STDERR "Writing individual files for sequences...\n";
    my $count = 0;
    foreach(@seqs) { 
	print STDERR $count%100 ? '' : '.';
	print STDERR $count%5000 ? '' : "\n";
	&write_file($_);  
    }
}

#----------------
sub wrap_up_seq {
#----------------
    # Want to verify that all sequences to be included, were.
    # Similarly for sequences to be excluded.
    # Also check for errors.
    my ($id, $count);
    my $trouble = 0;
    if(@errs) {
        printf STDERR "WARNING: %d errors were encountered.\n", scalar(@errs);
        print "View details of errors? [y|n] (y): ";
        if( ($response = <STDIN>) =~ /^y|^$/i) {
            foreach(@errs) { print STDERR $_; }
        }
    }

    if(%not_filtered) {
	$MONITOR && print STDERR "\nChecking for discrepencies.\n";
	$count = 0;
	foreach $id(keys %ids_incl_uc) {
	    if (!exists $_tested{$id}) {
		if(!$opt_exact and grep $_ =~ /$id/i, @_tested) {
		    last;
		}
		$count++;
		$MONITOR && print STDERR "NOT TESTED $count: $id\n";
		$trouble++;
	    }
	}
	$MONITOR && print STDERR "All tested.\n" unless $count;
	$count = 0;
	foreach $id(keys %ids_incl_uc) {
	    if (!exists $not_filtered{$id}) {
		if(!$opt_exact and grep $id =~ /$_/i, @ids_incl_uc) {
		    last;
		}
		$count++;
		$MONITOR && print STDERR "NOT INCLUDED $count: $id\n" unless exists $ids_excl_uc{$id};
		$trouble++;
	    }
	}
	$MONITOR && print STDERR "All included (as necessary).\n" unless $count;

	if(%ids_excl_uc) {
	    $count = 0;
	    foreach $id(keys %ids_excl_uc) {
		if (exists $not_filtered{$id}) {
		    if(!$opt_exact and grep $id =~ /$_/i, @ids_excl_uc) {
			last;
		    }
		    $count++;
		    $MONITOR && print STDERR "NOT EXCLUDED $count: $id\n";
		    $trouble++;
		}
	    }
	    $MONITOR && print STDERR "All excluded (as necessary).\n" unless $count;
	}
    }

    $MONITOR && print STDERR "Wrote files to $_dir\n" if $opt_write_files;
    local($^W) = 0;
    close $fh;
    $MONITOR && print STDERR "\n",'-'x50,"\nDone. $ID @argv\n";
    $MONITOR && print STDERR "$_count_processed sequences processed.\n";
    exit $trouble;
}

$VERSION;
