#
# BioPerl module for Bio::Tools::Run::StandAloneBlast
#
# Copyright Peter Schattner
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Tools::Run::StandAloneNCBIBlast - Object for the local execution 
of the NCBI BLAST program suite (blastall, blastpgp, bl2seq). With
experimental support for NCBI rpsblast.

=head1 SYNOPSIS

 # Do not use directly; see Bio::Tools::Run::StandAloneBlast

=head1 DESCRIPTION

See Bio::Tools::Run::StandAloneBlast

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists.  Your participation is much appreciated.

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
the bugs and their resolution.  Bug reports can be submitted via 
the web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Peter Schattner

Email schattner at alum.mit.edu

=head1 MAINTAINER - Torsten Seemann

Email torsten at infotech.monash.edu.au

=head1 CONTRIBUTORS

Sendu Bala  bix@sendu.me.uk (reimplementation)

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

package Bio::Tools::Run::StandAloneNCBIBlast;

use strict;
use warnings;

use base qw(Bio::Tools::Run::StandAloneBlast);

our $AUTOLOAD;
our $DEFAULTREADMETHOD = 'BLAST';

# If local BLAST databases are not stored in the standard
# /data directory, the variable BLASTDATADIR will need to be 
# set explicitly 
our $DATADIR = $Bio::Tools::Run::StandAloneBlast::DATADIR;

our %GENERAL_PARAMS  = (i => 'input',
                        o => 'outfile',
                        p => 'program',
                        d => 'database');
our @BLASTALL_PARAMS = qw(A B C D E F G K L M O P Q R S W X Y Z a b e f l m q r t v w y z n);
our @BLASTALL_SWITCH = qw(I g J T U n V s);
our @BLASTPGP_PARAMS = qw(A B C E F G H I J K L M N O P Q R S T U W X Y Z a b c e f h j k l m q s t u v y z);
our @RPSBLAST_PARAMS = qw(F I J L N O P T U V X Y Z a b e l m v y z);
our @BL2SEQ_PARAMS   = qw(A D E F G I J M S T U V W X Y a e g j m q r t);

our @OTHER_PARAMS = qw(_READMETHOD);


=head2 new

 Title   : new
 Usage   : my $obj = Bio::Tools::Run::StandAloneBlast->new();
 Function: Builds a newBio::Tools::Run::StandAloneBlast object 
 Returns : Bio::Tools::Run::StandAloneBlast
 Args    : -quiet => boolean # make program execution quiet
           -_READMETHOD => 'BLAST' (default, synonym 'SearchIO') || 'blast_pull'
                           # the parsing method, case insensitive

Essentially all BLAST parameters can be set via StandAloneBlast.pm.
Some of the most commonly used parameters are listed below. All
parameters have defaults and are optional except for -p in those programs that
have it. For a complete listing of settable parameters, run the relevant
executable BLAST program with the option "-" as in blastall -
Note that the input parameters (-i, -j, -input) should not be set directly by
you: this module sets them when you call one of the executable methods.

Blastall

  -p  Program Name [String]
        Input should be one of "blastp", "blastn", "blastx", 
        "tblastn", or "tblastx".
  -d  Database [String] default = nr
        The database specified must first be formatted with formatdb.
        Multiple database names (bracketed by quotations) will be accepted.
        An example would be -d "nr est"
  -e  Expectation value (E) [Real] default = 10.0
  -o  BLAST report Output File [File Out]  Optional,
	    default = ./blastreport.out ; set by StandAloneBlast.pm		
  -S  Query strands to search against database (for blast[nx], and tblastx). 3 is both, 1 is top, 2 is bottom [Integer]
	    default = 3

Blastpgp (including Psiblast)

  -j  is the maximum number of rounds (default 1; i.e., regular BLAST)
  -h  is the e-value threshold for including sequences in the
	    score matrix model (default 0.001)
  -c  is the "constant" used in the pseudocount formula specified in the paper (default 10)
  -B  Multiple alignment file for PSI-BLAST "jump start mode"  Optional
  -Q  Output File for PSI-BLAST Matrix in ASCII [File Out]  Optional

rpsblast

  -d  Database [String] default = (none - you must specify a database)
        The database specified must first be formatted with formatdb.
        Multiple database names (bracketed by quotations) will be accepted.
        An example would be -d "Cog Smart"
  -e  Expectation value (E) [Real] default = 10.0
  -o  BLAST report Output File [File Out]  Optional,
	    default = ./blastreport.out ; set by StandAloneBlast.pm		

Bl2seq

  -p  Program name: blastp, blastn, blastx. For blastx 1st argument should be nucleotide [String]
    default = blastp
  -o  alignment output file [File Out] default = stdout
  -e  Expectation value (E) [Real]  default = 10.0
  -S  Query strands to search against database (blastn only).  3 is both, 1 is top, 2 is bottom [Integer]
    default = 3

=cut

sub new {
    my ($caller, @args) = @_;
    my $self = $caller->SUPER::new(@args);
    
    # StandAloneBlast is special in that "one can modify the name of
    # the (ncbi) BLAST parameters as desired as long as the initial letter (and
    # case) of the parameter are preserved". We handle this by truncating input
    # args to their first char
    my %args = @args;
    @args = ();
    while (my ($attr, $value) = each %args) {
        $attr =~ s/^-//;
        $attr = substr($attr, 0, 1) unless $attr =~ /^_/;
        push(@args, $attr, $value);
    }
    
    $self->_set_from_args(\@args, -methods => {(map { $_ => $GENERAL_PARAMS{$_} } keys %GENERAL_PARAMS),
                                               (map { $_ => $_ } (@OTHER_PARAMS,
                                                                  @BLASTALL_PARAMS,
                                                                  @BLASTALL_SWITCH,
                                                                  @BLASTPGP_PARAMS,
                                                                  @RPSBLAST_PARAMS,
                                                                  @BL2SEQ_PARAMS))},
                                  -code => { map { $_ => 'my $self = shift;
                                                          if (@_) {
                                                              my $value = shift;
                                                              if ($value && $value ne \'F\') {
                                                                  $value = \'T\';
                                                              }
                                                              else {
                                                                  $value = \'F\';
                                                              }
                                                              $self->{\'_\'.$method} = $value;
                                                          }
                                                          return $self->{\'_\'.$method} || return;' } @BLASTALL_SWITCH },  # these methods can take boolean or 'T' and 'F'
                                  -create => 1,
                                  -force => 1,
                                  -case_sensitive => 1);
    
    my ($tfh, $tempfile) = $self->io->tempfile();
    my $outfile = $self->o || $self->outfile || $tempfile;
    $self->o($outfile);
    close($tfh);
    
    $self->_READMETHOD($DEFAULTREADMETHOD) unless $self->_READMETHOD;
    
    return $self;
}

# StandAloneBlast is special in that "one can modify the name of
# the (ncbi) BLAST parameters as desired as long as the initial letter (and
# case) of the parameter are preserved". We handle this with AUTOLOAD
# redirecting to the automatically created methods from _set_from_args() !
sub AUTOLOAD {
    my $self = shift;
    my $attr = $AUTOLOAD;
    $attr =~ s/.*:://;
    
    my $orig = $attr;
    
    $attr = substr($attr, 0, 1);
    
    $self->can($attr) || $self->throw("Unallowed parameter: $orig !");
    
    return $self->$attr(@_);
}

=head2  blastall

 Title   : blastall
 Usage   :  $blast_report = $factory->blastall('t/testquery.fa');
	or
	       $input = Bio::Seq->new(-id=>"test query",
				      -seq=>"ACTACCCTTTAAATCAGTGGGGG");
	       $blast_report = $factory->blastall($input);
	or 
	      $seq_array_ref = \@seq_array;  
         # where @seq_array is an array of Bio::Seq objects
	      $blast_report = $factory->blastall($seq_array_ref);
 Returns : Reference to a Blast object containing the blast report.
 Args    : Name of a file or Bio::Seq object or an array of 
           Bio::Seq object containing the query sequence(s). 
           Throws an exception if argument is not either a string 
           (eg a filename) or a reference to a Bio::Seq object 
           (or to an array of Seq objects).  If argument is string, 
           throws exception if file corresponding to string name can 
           not be found.

=cut

sub blastall {
    my ($self, $input1) = @_;
    $self->io->_io_cleanup();
    my $executable = 'blastall';
    
    # Create input file pointer
    my $infilename1 = $self->_setinput($executable, $input1) || $self->throw("$input1 not Bio::Seq object or array of Bio::Seq objects or file name!");
    $self->i($infilename1);
    
    my $blast_report = $self->_generic_local_blast($executable);
}

=head2  blastpgp

 Title   : blastpgp
 Usage   :  $blast_report = $factory-> blastpgp('t/testquery.fa');
	or
	       $input = Bio::Seq->new(-id=>"test query",
				      -seq=>"ACTADDEEQQPPTCADEEQQQVVGG");
	       $blast_report = $factory->blastpgp ($input);
	or
	      $seq_array_ref = \@seq_array;  
         # where @seq_array is an array of Bio::Seq objects
	      $blast_report = $factory-> blastpgp(\@seq_array);
 Returns : Reference to a Bio::SearchIO object containing the blast report 
 Args    : Name of a file or Bio::Seq object. In psiblast jumpstart 
           mode two additional arguments are required: a SimpleAlign 
           object one of whose elements is the query and a "mask" to 
           determine how BLAST should select scoring matrices see 
           DESCRIPTION above for more details.

           Throws an exception if argument is not either a string 
           (eg a filename) or a reference to a Bio::Seq object 
           (or to an array of Seq objects).  If argument is string, 
           throws exception if file corresponding to string name can 
           not be found.
 Returns : Reference to Bio::SearchIO object containing the blast report.

=cut

sub blastpgp {
    my $self = shift;
    my $executable = 'blastpgp';
    my $input1 = shift;
    my $input2 = shift;
    # used by blastpgp's -B option to specify which 
    # residues are position aligned
    my $mask = shift;
    
    my ($infilename1, $infilename2 ) = $self->_setinput($executable, 
                                                        $input1, $input2, 
                                                        $mask);
    if (!$infilename1) {$self->throw("$input1 not Bio::Seq object or array of Bio::Seq objects or file name!");}
    $self->i($infilename1);	# set file name of sequence to be blasted to inputfilename1 (-i param of blastpgp)
    if ($input2) {
        unless ($infilename2) {$self->throw("$input2 not SimpleAlign Object in pre-aligned psiblast\n");}
        $self->B($infilename2);	# set file name of partial alignment to inputfilename2 (-B param of blastpgp)
    }
    
    my $blast_report = $self->_generic_local_blast($executable);
}

=head2  rpsblast

 Title   : rpsblast
 Usage   :  $blast_report = $factory->rpsblast('t/testquery.fa');
	or
	       $input = Bio::Seq->new(-id=>"test query",
				      -seq=>"MVVLCRADDEEQQPPTCADEEQQQVVGG");
	       $blast_report = $factory->rpsblast($input);
	or
	      $seq_array_ref = \@seq_array;  
         # where @seq_array is an array of Bio::Seq objects
	      $blast_report = $factory->rpsblast(\@seq_array);
 Args    : Name of a file or Bio::Seq object or an array of 
           Bio::Seq object containing the query sequence(s). 
           Throws an exception if argument is not either a string 
           (eg a filename) or a reference to a Bio::Seq object 
           (or to an array of Seq objects).  If argument is string, 
           throws exception if file corresponding to string name can 
           not be found.
 Returns : Reference to a Bio::SearchIO object containing the blast report 

=cut

sub rpsblast {
    my ($self, $input1) = @_;
    $self->io->_io_cleanup();
    my $executable = 'rpsblast';
    
    # Create input file pointer
    my $infilename1 = $self->_setinput($executable, $input1) || $self->throw("$input1 not Bio::Seq object or array of Bio::Seq objects or file name!");
    $self->i($infilename1);
    
    my $blast_report = $self->_generic_local_blast($executable);
}

=head2   bl2seq

 Title   : bl2seq
 Usage   : $factory-> bl2seq('t/seq1.fa', 't/seq2.fa');
	or
	  $input1 = Bio::Seq->new(-id=>"test query1",
				  -seq=>"ACTADDEEQQPPTCADEEQQQVVGG");
	  $input2 = Bio::Seq->new(-id=>"test query2",
				  -seq=>"ACTADDEMMMMMMMDEEQQQVVGG");
	  $blast_report = $factory->bl2seq ($input1,  $input2);
 Returns : Reference to a BPbl2seq object containing the blast report.
 Args    : Names of 2 files  or 2 Bio::Seq objects containing the 
           sequences to be aligned by bl2seq.

           Throws an exception if argument is not either a pair of 
           strings (eg filenames) or references to Bio::Seq objects.  
           If arguments are strings, throws exception if files 
           corresponding to string names can not be found.

=cut

sub bl2seq {
    my $self = shift;
    my $executable = 'bl2seq';
    my $input1 = shift;
    my $input2 = shift;
    
    # Create input file pointer
    my ($infilename1, $infilename2 ) = $self->_setinput($executable, 
							  $input1, $input2);
    if (!$infilename1){$self->throw(" $input1  not Seq Object or file name!");}
    if (!$infilename2){$self->throw("$input2  not Seq Object or file name!");}
    
    $self->i($infilename1);	# set file name of first sequence to 
                            # be aligned to inputfilename1 
                            # (-i param of bl2seq)
    $self->j($infilename2);	# set file name of first sequence to 
                            # be aligned to inputfilename2 
                            # (-j param of bl2seq)
    
    my $blast_report = $self->_generic_local_blast($executable);   
}

=head2  _generic_local_blast

 Title   : _generic_local_blast
 Usage   : internal function not called directly
 Returns : Bio::SearchIO 
 Args    : Reference to calling object and name of BLAST executable 

=cut

sub _generic_local_blast {
    my $self = shift;
    my $executable = shift;
    
    # Create parameter string to pass to Blast program
    my $param_string = $self->_setparams($executable);
    
    # run Blast
    my $blast_report = $self->_runblast($executable, $param_string);
}

=head2  _runblast

 Title   :  _runblast
 Usage   :  Internal function, not to be called directly	
 Function:   makes actual system call to Blast program
 Example :
 Returns : Report Bio::SearchIO object in the appropriate format 
 Args    : Reference to calling object, name of BLAST executable, 
           and parameter string for executable 

=cut

sub _runblast {
	my ($self, $executable, $param_string) = @_;
	my ($blast_obj, $exe);
	if (! ($exe = $self->executable($executable)) ) {
		$self->warn("cannot find path to $executable");
		return;
	}
    
    # Use double quotes if executable path have empty spaces
    if ($exe =~ m/ /) {
        $exe = "\"$exe\"";
    }
	my $commandstring = $exe.$param_string;
    
	$self->debug("$commandstring\n");
	system($commandstring) && $self->throw("$executable call crashed: $? | $! | $commandstring\n");
    
    # set significance cutoff to set expectation value or default value
	# (may want to make this value vary for different executables)
	my $signif = $self->e() || 1e-5; 
    
    # get outputfilename
	my $outfile = $self->o();
    
    # this should allow any blast SearchIO parser (not just 'blast_pull' or 'blast',
    # but 'blastxml' and 'blasttable').  Fall back to 'blast' if not stipulated.
    my $method = $self->_READMETHOD;
	if ($method =~ /^(?:blast|SearchIO)/i )  {
        $method = 'blast' if $method =~ m{SearchIO}i;
		$blast_obj = Bio::SearchIO->new(-file => $outfile,
                                        -format => $method);
	}
    # should these be here?  They have been deprecated...
    elsif ($method =~ /BPlite/i ) {
		if ($executable =~ /bl2seq/i)  {
			# Added program info so BPbl2seq can compute strand info
			$self->throw("Use of Bio::Tools::BPbl2seq is deprecated; use Bio::SearchIO modules instead");
		}
        elsif ($executable =~ /blastpgp/i && defined $self->j() && $self->j() > 1) {
			$self->throw("Use of Bio::Tools::BPpsilite is deprecated; use Bio::SearchIO modules instead");
		}
        elsif ($executable =~ /blastall|rpsblast/i) { 
			$self->throw("Use of Bio::Tools::BPlite is deprecated; use Bio::SearchIO modules instead");
		}
        else { 
			$self->warn("Unrecognized executable $executable");
		}
	}
    else {
		$self->warn("Unrecognized readmethod $method");
	}
    
	return $blast_obj;
}

=head2  _setparams

 Title   : _setparams
 Usage   : Internal function, not to be called directly	
 Function: Create parameter inputs for Blast program
 Example :
 Returns : parameter string to be passed to Blast 
 Args    : Reference to calling object and name of BLAST executable

=cut

sub _setparams {
    my ($self, $executable) = @_;
    my ($attr, $value, @execparams);
    
    if    ($executable eq 'blastall') { @execparams = (@BLASTALL_PARAMS,
                                                       @BLASTALL_SWITCH); }
    elsif ($executable eq 'blastpgp') { @execparams =  @BLASTPGP_PARAMS;  }
    elsif ($executable eq 'rpsblast') { @execparams =  @RPSBLAST_PARAMS;  }
    elsif ($executable eq 'bl2seq'  ) { @execparams =  @BL2SEQ_PARAMS;    }
    
    # we also have all the general params
    push(@execparams, keys %GENERAL_PARAMS);
    
    my $database = $self->d;
    if ($database && $executable ne 'bl2seq') {
        # Need to prepend datadirectory to database name
        my @dbs = split(/ /, $database);
        for my $i (0..$#dbs) {
            # (works with multiple databases)
            if (! (-e $dbs[$i].".nin" || -e $dbs[$i].".pin") &&
                ! (-e $dbs[$i].".nal" || -e $dbs[$i].".pal") ) {
                $dbs[$i] = File::Spec->catdir($DATADIR, $dbs[$i]);
            }
        }
        $self->d('"'.join(" ", @dbs).'"');
    }
    
    # workaround for problems with shell metacharacters [bug 2707]
    # simply quoting does not always work!
    my $tmp = $self->o;
    $self->o(quotemeta($tmp)) if ($tmp && $^O !~ /^MSWin/);
    
    my $param_string = $self->SUPER::_setparams(-params => [@execparams],
                                                -dash => 1);
    
    $self->o($tmp) if ($tmp && $^O !~ /^MSWin/);

    $self->d($database) if $database;
    
    if ($self->quiet()) { 
        $param_string .= ' 2> '.File::Spec->devnull;
    }
    
    return $param_string;
}

1;
