BEGIN {
        use vars qw($INSTALL_PATH);

        #### NOTE: $INSTALL_PATH needs to be customized
        ###         your system. It should point
        ###         to the directory right above Bio/ in order
        ###         for perl to be able to locate the .pm files. 

	$INSTALL_PATH = "/home/users/sac/perl/lib";

        ###
        ####
    }

#---------------------------------------------------------------------------
# LIBRARY      : seqtools.pl
# PURPOSE      : To provide a set of standard functions & variables for 
#                working with sets of sequences using Bioperl modules.
#                Currently focuses on the Bio::Tools::Fasta.pm module.
#                Thus, only Fasta-formatted sequence may be used as input.
# AUTHOR       : Steve A. Chervitz ((sac@genome.stanford.edu)
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
#  require "/home/me/perl/seq/seqtools.pl"; 
#            
#  &init_seq('myscript.pl', 0.1);
#  &load_ids();
#  &get_seq_objs();
#  &print_seqs();
#  &wrap_up();
#
# EXAMPLES : See the files seqs.pl, seqs2.pl, seqs3.pl, seqs4.pl
#            in this dir for some working examples.
#
#
# PUBLIC FUNCTIONS DEFINED IN seqtools.pl:
#
#    init_seq()
#    load_ids()
#    get_seq_objs()
#    get_seq_data()
#    print_seq()
#    print_seqs()
#    write_file()
#    write_files()
#    wrap_up_seq()
#
# COMMENTS :
#            The sequences can be output in a variety of formats 
#            (Fasta, Raw, GenBank, PIR, GCG, GCG_SEQ)
#            but the input format must be Fasta. The variety of output formats
#            will increase as Bio::Seq.pm matures.
#            
# BUGS:
#  FileHandle.pm may produce a harmless warning with the -w flag such as:
# "Close on unopened file <GEN0> at /usr/local/bin/perl/lib/FileHandle.pm line 255"
#
# MODIFIED:
#  sac --- Tue Nov 24 10:55:46 1998
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

use Bio::Tools::Fasta qw(:obj);
use Bio::Root::Global qw(:devel);
use Getopt::Long;
require FileHandle;
use Carp;

select(STDOUT); $|=1;

# Global:
$ID         = 'seqtools.pl';
$VERSION    = 0.21;
$DESC       = "Provides standard functions and variables for working with sets\n".
              "of sequences and sequence objects using Bioperl modules";
%seqParams  = ();  
$seqcount   = 0;
@seqs       = ();
@ids        = ();
@descs      = ();
$_count_processed = 0;

# Command-line options:
$opt_strict     = 0;
$opt_nucl       = 0;
$opt_prot       = 0;
$opt_eid        = 0; 
$opt_eseq       = 0; 
$opt_col        = 0;         # Column in the -incl or -excl files with the seq IDs.
$opt_exact      = 1;         # Require exact match when filtering based on seq id.
#$opt_fmt       = 'Fasta';   # Only working with Fasta input (for now)
$opt_outfmt     = 'Fasta';
$opt_seq        = undef;
$opt_incl       = undef;
$opt_excl       = undef;
$opt_out        = undef;
$opt_err        = undef;
$opt_tag        = undef;  # Optional tag string to be prepended to all descrptions.
$opt_write_files = undef; # Hold directory name for where to write individual files 
                          # for each input sequence. 

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
$fh          = undef;  # FileHandle object for output (STDOUT or the -out specified file).

# Private:
my (%not_filtered, $type, %_tested, @_tested, $_dir);
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

    print STDERR "$ID, $VERSION\n$DESC.\n";
    print STDERR <<"QQ_USAGE_QQ";

Usage: $ID [ parameters ] seq/seqs.fasta > out
       $ID [ parameters ] seq/*.fasta   > out
       gzcat *fasta.gz | $ID [ parameters ] > out

 input : Multiple Fasta-formatted sequences are read via STDIN
         or from a file(s) specified on the command line.
         
QQ_USAGE_QQ

}

#---------------
sub seq_params {
#---------------

    print STDERR <<"QQ_PARAMS_QQ";

 SEQTOOLS PARAMETERS:
 --------------------
 -incl <file|list> : Filename containing list of sequence IDs to include
                     -or- list of IDs separated by whitespace.
                     For lines containing multiple ids, only
                     the first one is used, unless the -last option is used.
 -excl <file|list> : Filename containing list of sequence IDs to exclude
                     -or- list of IDs separated by whitespace.
 -seq  file    : Filename containing Fasta formated sequences (if not
                 using STDIN). ('-seq' is optional but should be used
	         with scripts that accept different types of input files).
 -nucl         : Nucleotide sequences  (optional).
 -prot         : Peptide sequences (optional).
 -out <file>   : Filename for saving output (default = STDOUT).
 -outfmt <fmt> : Sequence format for saving sequence data (default = Fasta).
                   Also supported: raw, gcg, gcg_seq, gcg_ref.
 -eid          : Edit sequence IDs. 
                 (Uppercases and cleans up complex sequence identifiers:
	         gi|2980872|gnl|PID|e1283615 homeobox protein SHOTb
 	                    is converted to:
	         GI_2980872 homeobox protein SHOTb [ GNL PID e1283615 ]
                 Recommended when working with GenBank Fasta files.
	         Information will not be lost, just rearranged.
	         (Default editing: IDs are uppercased only.)
 -eseq         : Edit sequence (removes ambiguous characters at termini).
	         (Default editing: sequences will be uppercased and have all 
                  non-alphabetic, non-gap characters and white space removed.
		  Allowed gap characters = '.' and '-')
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
 -strict       : Create sequence objects in strict mode (default = nostrict).
 
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

    &GetOptions('last!', 'mon!', 'h!', 'debug!', 'incl=s', 'excl=s', 
		'seq=s', 'nucl!', 'prot!', 'strict!', 'out=s', 'eid!', 'eseq!', 
		'exact!', 'outfmt=s', 'err=s', 'eg!', 'write_files=s', 'tag=s',
		@opts);
   
    $MONITOR && print STDERR "$ID, v$VERSION\n",'-'x50,"\n";

    $opt_h and do{
	(ref($usage_fref) =~ /CODE/) ? &$usage_fref : &print_seq_usage;
	exit 1;
    };

    $opt_eg and do{
	eval{ print STDERR "\n$ID EXAMPLES:\n",'-'x30,"\n\n", &examples(); };
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

    $type = 'Amino' if $opt_prot;
    $type = 'Dna' if $opt_nucl;

    if($MONITOR) {
	$type or print STDERR "\n*** Sequence type not specified (can use -nucl or -prot on command line)\n\n";
    }
    
    if($opt_outfmt !~ /fasta|gcg|gcg_seq|raw/i ) {
	print STDERR "\n*** Output format not supported: $opt_outfmt\n";
	print STDERR "\tValid formats: Raw, Fasta, GCG, GCG_SEQ, GCG_REF\n\n";
    }

    local($^W) = 0;

    $MONITOR && print STDERR "Output format: $opt_outfmt\n";

    if($opt_out) {
	$MONITOR and print STDERR "\nWriting sequence data to $opt_out\n\n";
    } 
    if($opt_write_files) {
	$MONITOR and print STDERR "\nWriting sequence files in $opt_write_files\n\n";
    } 

    $fh = new FileHandle;

    if($opt_out) {
	open ($fh, ">$opt_out") || croak "\n\a*** $0: Can't open output file $opt_out: $!\n\n";
    } else {
	$fh = \*STDOUT;
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
    printf STDERR "%-20s %s\n", 'strict', $opt_strict;
    printf STDERR "%-20s %s\n", 'exact', $opt_exact;
    printf STDERR "%-20s %s\n", 'eid', $opt_eid;
    printf STDERR "%-20s %s\n", 'eseq', $opt_eseq || '';
    printf STDERR "%-20s %s\n", 'out', $opt_out || '';
    printf STDERR "%-20s %s\n", 'last', $opt_last || '';
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


#-----------------
sub get_seq_objs {
#-----------------
# get_seq_objs() may be called with an optional function ref as the sole argument.
# The supplied function should expect a single argument that is a
# object reference for Bio::Seq.pm object.
# If no function_ref is supplied, the sequence object will be loaded into
# the @seqs array.

    my $func_ref = shift;

    # The module now prints this.
#    $MONITOR && print STDERR "\nParsing Fasta sequence objects (100/dot, 5000/line).\n";

    # Use the static Fasta object provided by Fasta.pm.

    # As the sequence data is parsed, the Fasta object returns a set of
    # sequence objects, one for each sequence.
    # (Parameters tags can be lowercase if desired.)

    # The -SEQS parameter here serves only to indicate that we want 
    # to parse a Fasta sequence file (not a analysis report).
    # (Compare with get_seq_data() below).

    %params = (
	       -TYPE       => $type,
	    #  -PARSE      => 1,  # not needed since were calling parse() directly
	       -SEQS       => 1,
	       -EDIT_ID    => $opt_eid,
	       -EDIT_SEQ   => $opt_eseq,
	       -STRICT     => $opt_strict,
	       -FILT_FUNC  => ($opt_incl or $opt_excl) ? \&seq_filter : undef,
	       -EXEC_FUNC  => $func_ref || undef,
	       -SAVE_ARRAY => \@seqs,
	       );
    
    _parse_seqs();
}

#-----------------
sub get_seq_data {
#-----------------
# get_seq_data() may be called with an optional function ref:
#    &get_seqs(\&funct_ref)
# The supplied function should expect three string arguments:
#   id, description, and sequence  (in that order).

    my $func_ref = shift;

    $MONITOR && print STDERR "\nParsing Fasta sequence data (100/dot, 5000/line).\n";

    # Note that the array refs will not be used if a $func_ref is supplied.
    # Defining these parameters are neccessary, however, to indicate that 
    # we want data and not sequence objects.

    # The -SEQS parameter serves double duty here since it serves to indicate
    # that we want to parse a Fasta sequence file (not a analysis report)
    # and it provides a place to store the raw sequence data if a $func_ref
    # is not supplied.

    # Fasta object returns a set of sequence objects, one for each sequence.
    %params = (
	       -TYPE       => $type,
	    #  -PARSE      => 1,  # not needed since were calling parse() directly
	       -SEQS       => \@seqs,
	       -IDS        => \@ids,
	       -DESCS      => \@descs,
	       -EDIT_ID    => $opt_eid,
	       -EDIT_SEQ   => $opt_eseq,
	       -STRICT     => $opt_strict,
	       -FILT_FUNC  => ($opt_incl or $opt_excl) ? \&seq_filter : undef,
	       -EXEC_FUNC  => $func_ref || undef,
	       );
    
    _parse_seqs();
}

#------------------
sub _parse_seqs {
#------------------

    eval {
	if(@ARGV) {
	    foreach(@ARGV) {
		next unless -s;
		$MONITOR && print STDERR "\nParsing file $_\n";
		$params{-FILE} = $_;
		$seqcount += $Fasta->parse(%params);
		
	    }
	} elsif($opt_seq) {
	    $params{-FILE} = $opt_seq;
	    $seqcount = $Fasta->parse(%params);

	} else {
	    $seqcount = $Fasta->parse(%params);
	}
    };
    if($@) {
	print STDERR "\n*** TROUBLE:\n$@\n";
	exit 1;
    }

    $MONITOR && print STDERR "\n$seqcount sequence(s) loaded/processed.\n";
    exit 1 if !$seqcount;
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

    print $fh $seq->layout($opt_outfmt);
    print "\n" if $opt_outfmt =~ /raw/i; # raw format does not include new line.
    $_count_processed++;

    if(@ids_incl_uc) {
	return $count_processed >= scalar(@ids_incl_uc) ? 0 : 1;
    }
    1;
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
    open(OUT, ">$file") || croak "\n*** $0: Can't open file $file: $!\n\n";
    print OUT $seq->layout($opt_outfmt);
    close OUT;

    $_count_processed++;
    if(@ids_incl_uc) {
	return $count_processed >= scalar(@ids_incl_uc) ? 0 : 1;
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
    # There seems to be some discrepency.
    my ($id, $count);
    my $trouble = 0;
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
	$MONITOR && print STDERR "All included.\n" unless $count;

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
	    $MONITOR && print STDERR "All excluded.\n" unless $count;
	}
    }

    $MONITOR && print STDERR "Wrote files to $_dir\n" if $opt_write_files;
    local($^W) = 0;
    close $fh;
    $MONITOR && print STDERR "\n",'-'x50,"\nDone. $ID @argv\n";
    exit $trouble;
}

$VERSION;
