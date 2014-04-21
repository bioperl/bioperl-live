#
# BioPerl module for Bio::Tools::Run::StandAloneBlast
#
# Copyright Peter Schattner
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Tools::Run::StandAloneBlast - Object for the local execution 
of the NCBI BLAST program suite (blastall, blastpgp, bl2seq). 
There is experimental support for WU-Blast and NCBI rpsblast.

=head1 SYNOPSIS

 # Local-blast "factory object" creation and blast-parameter
 # initialization:
 @params = (-database => 'swissprot', -outfile => 'blast1.out');
 $factory = Bio::Tools::Run::StandAloneBlast->new(@params);

 # Blast a sequence against a database:
 $str = Bio::SeqIO->new(-file=>'t/amino.fa', -format => 'Fasta');
 $input = $str->next_seq();
 $input2 = $str->next_seq();
 $blast_report = $factory->blastall($input);

 # Run an iterated Blast (psiblast) of a sequence against a database:
 $factory->j(3);    # 'j' is blast parameter for # of iterations
 $factory->outfile('psiblast1.out');
 $factory = Bio::Tools::Run::StandAloneBlast->new(@params);
 $blast_report = $factory->blastpgp($input);

 # Use blast to align 2 sequences against each other:
 $factory = Bio::Tools::Run::StandAloneBlast->new(-outfile => 'bl2seq.out');
 $factory->bl2seq($input, $input2);

 # Experimental support for WU-Blast 2.0
 my $factory = Bio::Tools::Run::StandAloneBlast->new(-program =>"wublastp",
                                                     -database =>"swissprot",
                                                     -e => 1e-20); 
 my $blast_report = $factory->wublast($seq);

 # Experimental support for NCBI rpsblast
 my $factory = Bio::Tools::Run::StandAloneBlast->new(-db => 'CDD/Cog', 
                                                     -expect => 0.001);
 $factory->F('T'); # turn on SEG filtering of query sequence
 my $blast_report = $factory->rpsblast($seq);

 # Use the experimental fast Blast parser, 'blast_pull'
 my $factory = Bio::Tools::Run::StandAloneBlast->new(-_READMETHOD =>'blast_pull',
                                                     @other_params);

 # Various additional options and input formats are available,
 # see the DESCRIPTION section for details.

=head1 DESCRIPTION

This DESCRIPTION only documents Bio::Tools::Run::StandAloneBlast, a
Bioperl object for running the NCBI standAlone BLAST package. Blast
itself is a large & complex program - for more information regarding
BLAST, please see the BLAST documentation which accompanies the BLAST
distribution. BLAST is available from ftp://ncbi.nlm.nih.gov/blast/.

A source of confusion in documenting a BLAST interface is that the
term "program" is used in - at least - three different ways in the
BLAST documentation. In this DESCRIPTION, "program" will refer to the
BLAST routine set by the BLAST C<-p> parameter that can be set to blastn,
blastp, tblastx etc. We will use the term Blast "executable" to refer
to the various different executable files that may be called - ie.
blastall, blastpgp or bl2seq. In addition, there are several BLAST
capabilities, which are also referred to as "programs", and are
implemented by using specific combinations of BLAST executables,
programs and parameters. They will be referred by their specific
names - eg PSIBLAST and PHIBLAST.

Before running StandAloneBlast it is necessary: to install BLAST 
on your system, to edit set the environmental variable $BLASTDIR 
or your $PATH variable to point to the BLAST directory, and to 
ensure that users have execute privileges for the BLAST program.  

If the databases which will be searched by BLAST are located in the 
data subdirectory of the blast program directory (the default 
installation location), StandAloneBlast will find them; however, 
if the database files are located in any other location, environmental 
variable $BLASTDATADIR will need to be set to point to that directory.

The use of the StandAloneBlast module is as follows: Initially, a
local blast "factory object" is created. The constructor may be passed
an optional array of (non-default) parameters to be used by the
factory, eg:

 @params = (-program => 'blastn', -database => 'ecoli.nt');
 $factory = Bio::Tools::Run::StandAloneBlast->new(@params);

Any parameters not explicitly set will remain as the defaults of the
BLAST executable. Note each BLAST executable has somewhat different
parameters and options. See the BLAST Documentation for a description
or run the BLAST executable from the command line followed solely with
a "-" to see a list of options and default values for that executable;
eg E<gt>blastall -.

BLAST parameters can be changed and/or examined at any time after the
factory has been created. The program checks that any
parameter/switch being set/read is valid. Except where specifically
noted, StandAloneBlast uses the same single-letter, case-sensitive
parameter names as the actual blast program. Currently no checks are
included to verify that parameters are of the proper type (e.g. string
or numeric) or that their values are within the proper range.

As an example, to change the value of the Blast parameter 'e' ('e' is
the parameter for expectation-value cutoff) 

  $expectvalue = 0.01;
  $factory->e($expectvalue);

Note that for improved script readibility one can modify the name of
the (ncbi) BLAST parameters as desired as long as the initial letter (and
case) of the parameter are preserved, e.g.:

  $factory->expectvalue($expectvalue);

Unfortunately, some of the BLAST parameters are not the single 
letter one might expect (eg "iteration round" in blastpgp is 'j'). 
Again one can check by using, for example:

  > blastpgp -

Wublast parameters need to be complete (ie. don't truncate them to their
first letter), but are case-insensitive.

Once the factory has been created and the appropriate parameters set,
one can call one of the supported blast executables. The input
sequence(s) to these executables may be fasta file(s) as described in
the BLAST documentation.

  $inputfilename = 't/testquery.fa';
  $blast_report = $factory->blastall($inputfilename);

In addition, sequence input may be in the form of either a Bio::Seq
object or (a reference to) an array of Bio::Seq objects, e.g.:

  $input = Bio::Seq->new(-id => "test query",
                         -seq => "ACTACCCTTTAAATCAGTGGGGG");
  $blast_report = $factory->blastall($input);

NOTE: Use of the BPlite method has been deprecated and is no longer supported.

For blastall and non-psiblast blastpgp runs, report object is a L<Bio::SearchIO>
object, selected by the user with the parameter _READMETHOD. The leading
underscore is needed to distinguish this option from options which are passed to
the BLAST executable. The default parser is Bio::SearchIO::blast. In any case,
the "raw" blast report is also available. The filename is set by the 'outfile'
parameter and has the default value of "blastreport.out".

For psiblast execution in the BLAST "jumpstart" mode, the program must
be passed (in addition to the query sequence itself) an alignment
containing the query sequence (in the form of a SimpleAlign object) as
well as a "mask" specifying at what residues position-specific scoring
matrices (PSSMs) are to used and at what residues default scoring
matrices (eg BLOSUM) are to be used. See psiblast documentation for
more details. The mask itself is a string of 0's and 1's which is the
same length as each sequence in the alignment and has a "1" at
locations where (PSSMs) are to be used and a "0" at all other
locations. So for example:

  $str = Bio::AlignIO->new(-file => "cysprot.msf", 
                           -format => 'msf');
  $aln = $str->next_aln();
  $len = $aln->length_aln();
  $mask = '1' x $len;
  # simple case where PSSM's to be used at all residues
  $report = $factory->blastpgp("cysprot1.fa", $aln, $mask);

For bl2seq execution, StandAloneBlast.pm can be combined with
AlignIO.pm to directly produce a SimpleAlign object from the alignment
of the two sequences produced by bl2seq as in:

  # Get 2 sequences
  $str = Bio::SeqIO->new(-file=>'t/amino.fa' , -format => 'Fasta');
  my $seq3 = $str->next_seq();
  my $seq4 = $str->next_seq();

  # Run bl2seq on them
  $factory = Bio::Tools::Run::StandAloneBlast->new(-program => 'blastp',
                                                   -outfile => 'bl2seq.out');
  my $bl2seq_report = $factory->bl2seq($seq3, $seq4);

  # Use AlignIO.pm to create a SimpleAlign object from the bl2seq report
  $str = Bio::AlignIO->new(-file=> 'bl2seq.out',-format => 'bl2seq');
  $aln = $str->next_aln();

For more examples of syntax and use of StandAloneBlast.pm, the user is
encouraged to run the scripts standaloneblast.pl in the bioperl
examples/tools directory and StandAloneBlast.t in the bioperl t/ 
directory.

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

package Bio::Tools::Run::StandAloneBlast;

use strict;
use warnings;

use Bio::Root::IO;
use Bio::Seq;
use Bio::SeqIO;
use Bio::SearchIO;
use File::Spec;

use base qw(Bio::Tools::Run::WrapperBase Bio::Factory::ApplicationFactoryI);

our $AUTOLOAD;
our $DEFAULTBLASTTYPE = 'NCBI';
our $DEFAULTREADMETHOD = 'BLAST';

# If local BLAST databases are not stored in the standard
# /data directory, the variable BLASTDATADIR will need to be 
# set explicitly 
our $DATADIR = $ENV{'BLASTDATADIR'} || $ENV{'BLASTDB'};
if (! defined $DATADIR && defined $ENV{'BLASTDIR'}) {
    my $dir = Bio::Root::IO->catfile($ENV{'BLASTDIR'}, 'data');
    if (-d $dir) {
        $DATADIR = $dir;
    }
    elsif ($ENV{'BLASTDIR'} =~ /bin/) {
        $dir = $ENV{'BLASTDIR'};
        $dir =~ s/bin/data/;
        $DATADIR = $dir if -d $dir;
    }
}

=head2 new

 Title   : new
 Usage   : my $obj = Bio::Tools::Run::StandAloneBlast->new();
 Function: Builds a newBio::Tools::Run::StandAloneBlast object 
 Returns : Bio::Tools::Run::StandAloneNCBIBlast or StandAloneWUBlast
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

WU-Blast

  -p Program Name [String] 
        Input should be one of "wublastp", "wublastn", "wublastx", 
        "wutblastn", or "wutblastx".
  -d  Database [String] default = nr
        The database specified must first be formatted with xdformat.
  -E  Expectation value (E) [Real] default = 10.0
  -o  BLAST report Output File [File Out]  Optional,
	    default = ./blastreport.out ; set by StandAloneBlast.pm		

=cut

sub new {
    my ($caller, @args) = @_;
    my $class = ref($caller) || $caller;
    
    # Because of case-sensitivity issues, ncbi and wublast methods are
    # mutually exclusive. We can't load ncbi methods if we start with wublast
    # (and vice versa) since wublast e() and E() should be the same thing,
    # whilst they must be different things in ncbi blast.
    #
    # Solution: split StandAloneBlast out into two more modules for NCBI and WU
    
    if ($class =~ /NCBI|WU/) {
        return $class->SUPER::new(@args);
    }
    
    my %args = @args;
    my $blasttype = $DEFAULTBLASTTYPE;
    while (my ($attr, $value) = each %args) {
        if ($attr =~/^-?\s*program\s*$|^-?p$/) {
            if ($value =~ /^wu*/) {
                $blasttype = 'WU';
            }
        }
    }
    
    my $module = "Bio::Tools::Run::StandAlone${blasttype}Blast";
    Bio::Root::Root->_load_module($module);
    return $module->new(@args);
}

=head2 executable

 Title   : executable
 Usage   : my $exe = $blastfactory->executable('blastall');
 Function: Finds the full path to the executable
 Returns : string representing the full path to the exe
 Args    : [optional] name of executable to set path to 
           [optional] boolean flag whether or not warn when exe is not found

=cut

sub executable {
    my ($self, $exename, $exe, $warn) = @_;
    $exename = 'blastall' unless (defined $exename || $self =~ /WUBlast/);
    $self->program_name($exename);
    
    if( defined $exe && -x $exe ) {
        $self->{'_pathtoexe'}->{$exename} = $exe;
    }
    unless( defined $self->{'_pathtoexe'}->{$exename} ) {
        my $f = $self->program_path($exename);	    
        $exe = $self->{'_pathtoexe'}->{$exename} = $f if(-e $f && -x $f );
        
        # This is how I meant to split up these conditionals --jason
        # if exe is null we will execute this (handle the case where
        # PROGRAMDIR pointed to something invalid)
        unless( $exe )  {  # we didn't find it in that last conditional
            if( ($exe = $self->io->exists_exe($exename)) && -x $exe ) {
                $self->{'_pathtoexe'}->{$exename} = $exe;
            }
            else { 
                $self->warn("Cannot find executable for $exename") if $warn;
                $self->{'_pathtoexe'}->{$exename} = undef;
            }
        }
   }
   return $self->{'_pathtoexe'}->{$exename};
}

=head2 program_dir

 Title   : program_dir
 Usage   : my $dir = $factory->program_dir();
 Function: Abstract get method for dir of program. 
 Returns : string representing program directory 
 Args    : none 

=cut

sub program_dir {
    my $self = shift;
    $self =~ /NCBIBlast/? $ENV{'BLASTDIR'}: $ENV{'WUBLASTDIR'};
}

sub program_name {
    my $self = shift;
    if (@_) { $self->{program_name} = shift }
    return $self->{program_name} || '';
}

sub program {
    my $self = shift;
    if( wantarray ) {
	return ($self->executable, $self->p());
    } else {
	return $self->executable(@_);
    }
}

=head2  _setinput

 Title   :  _setinput
 Usage   :  Internal function, not to be called directly	
 Function:   Create input file(s) for Blast executable
 Example :
 Returns : name of file containing Blast data input
 Args    : Seq object reference or input file name

=cut

sub _setinput {
	my ($self, $executable, $input1, $input2) = @_;
	my ($seq, $temp, $infilename1, $infilename2,$fh ) ;
	#  If $input1 is not a reference it better be the name of a file with
	#  the sequence/ alignment data...
	$self->io->_io_cleanup();

    SWITCH: {
        unless (ref $input1) {
			$infilename1 = (-e $input1) ? $input1 : 0 ;
			last SWITCH; 
        }
        
		# $input may be an array of BioSeq objects...
        if (ref($input1) =~ /ARRAY/i ) {
			($fh,$infilename1) = $self->io->tempfile();
			$temp =  Bio::SeqIO->new(-fh=> $fh, -format => 'fasta');
			foreach $seq (@$input1) {
				unless ($seq->isa("Bio::PrimarySeqI")) {return 0;}
				$seq->display_id($seq->display_id);
				$temp->write_seq($seq);
			}
			close $fh;
			$fh = undef;
			last SWITCH;
        }
        
		#  $input may be a single BioSeq object...
        elsif ($input1->isa("Bio::PrimarySeqI")) {
			($fh,$infilename1) = $self->io->tempfile();
            
			# just in case $input1 is taken from an alignment and has spaces (ie
			# deletions) indicated within it, we have to remove them - otherwise
			# the BLAST programs will be unhappy
			my $seq_string =  $input1->seq();
			$seq_string =~ s/\W+//g; # get rid of spaces in sequence
			$input1->seq($seq_string);
			$temp =  Bio::SeqIO->new(-fh=> $fh, '-format' => 'fasta');
			$temp->write_seq($input1);
			close $fh;
			undef $fh;
			last SWITCH;
        }
        
        $infilename1 = 0; # Set error flag if you get here
	}
	
    unless ($input2) { return $infilename1; }
    
    SWITCH2: {
        unless (ref $input2) {
			$infilename2 =   (-e $input2) ? $input2 : 0 ;
			last SWITCH2; 
        }
        if ($input2->isa("Bio::PrimarySeqI")  && $executable  eq 'bl2seq' ) {
			($fh,$infilename2) = $self->io->tempfile();
            
			$temp =  Bio::SeqIO->new(-fh=> $fh, '-format' => 'Fasta');
			$temp->write_seq($input2);
			close $fh;
			undef $fh;
			last SWITCH2;
        }
        
		# Option for using psiblast's pre-alignment "jumpstart" feature
        elsif ($input2->isa("Bio::SimpleAlign") && $executable eq 'blastpgp' ) {
			# a bit of a lie since it won't be a fasta file
            ($fh,$infilename2) = $self->io->tempfile(); 
            
            # first we retrieve the "mask" that determines which residues should
            # by scored according to their position and which should be scored
            # using the non-position-specific matrices
            my @mask = split("", shift );	#  get mask
            
            # then we have to convert all the residues in every sequence to upper
            # case at the positions that we want psiblast to use position specific
            # scoring
            foreach $seq ( $input2->each_seq() ) {
                my @seqstringlist = split("",$seq->seq());
                for (my $i = 0; $i < scalar(@mask); $i++) {
                    unless ( $seqstringlist[$i] =~ /[a-zA-Z]/ ) {next}
                    $seqstringlist[$i] = $mask[$i] ? uc $seqstringlist[$i]: lc $seqstringlist[$i] ;
                }
                my $newseqstring = join("", @seqstringlist);
                $seq->seq($newseqstring);
            }
            
            # Now we need to write out the alignment to a file 
            # in the "psi format" which psiblast is expecting
            $input2->map_chars('\.','-');
            $temp =  Bio::AlignIO->new(-fh=> $fh, '-format' => 'psi');
            $temp->write_aln($input2);
            close $fh;
            undef $fh;
            last SWITCH2;
        }
        
        $infilename2 = 0; # Set error flag if you get here
	}
    
	return ($infilename1, $infilename2);
}

=head1 Bio::Tools::Run::WrapperBase methods

=cut

=head2 no_param_checks

 Title   : no_param_checks
 Usage   : $obj->no_param_checks($newval)
 Function: Boolean flag as to whether or not we should
           trust the sanity checks for parameter values  
 Returns : value of no_param_checks
 Args    : newvalue (optional)

=cut

=head2 save_tempfiles

 Title   : save_tempfiles
 Usage   : $obj->save_tempfiles($newval)
 Function: 
 Returns : value of save_tempfiles
 Args    : newvalue (optional)

=cut

=head2 outfile_name

 Title   : outfile_name
 Usage   : my $outfile = $tcoffee->outfile_name();
 Function: Get/Set the name of the output file for this run
           (if you wanted to do something special)
 Returns : string
 Args    : [optional] string to set value to

=cut

=head2 tempdir

 Title   : tempdir
 Usage   : my $tmpdir = $self->tempdir();
 Function: Retrieve a temporary directory name (which is created)
 Returns : string which is the name of the temporary directory
 Args    : none

=cut

=head2 cleanup

 Title   : cleanup
 Usage   : $tcoffee->cleanup();
 Function: Will cleanup the tempdir directory after a PAML run
 Returns : none
 Args    : none

=cut

=head2 io

 Title   : io
 Usage   : $obj->io($newval)
 Function:  Gets a Bio::Root::IO object
 Returns : Bio::Root::IO
 Args    : none

=cut

1;
