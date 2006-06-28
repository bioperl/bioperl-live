# $Id$
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

 @params = (-database => 'swissprot',-outfile => 'blast1.out');
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

 # Various additional options and input formats are available,
 # see the DESCRIPTION section for details.

=head1 DESCRIPTION

This DESCRIPTION only documents Bio::Tools::Run::StandAloneBlast: - a
Bioperl object for running the NCBI standAlone BLAST package.  Blast,
itself, is a large & complex program - for more information regarding
BLAST, please see the BLAST documentation which accompanies the BLAST
distribution. BLAST is available from ftp://ncbi.nlm.nih.gov/blast/.

A source of confusion in documenting a BLAST interface is that the
term "program" is used in - at least - three different ways in the
BLAST documentation.  In this DESCRIPTION, "program" will refer to the
BLAST routine set by the BLAST C<-p> parameter that can be set to blastn,
blastp, tblastx etc.  We will use the term Blast "executable" to refer
to the various different executable files that may be called - ie
blastall, blastpgp or bl2seq.  In addition, there are several BLAST
capabilities, which are also referred to as "programs", and are
implemented by using specific combinations of BLAST executables,
programs and parameters.  They will be referred by their specific
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
BLAST executable.  Note each BLAST executable has somewhat different
parameters and options.  See the BLAST Documentation for a description
or run the BLAST executable from the command line followed solely with
a "-" to see a list of options and default values for that executable;
eg E<gt>blastall -.

BLAST parameters can be changed and/or examined at any time after the
factory has been created.  The program checks that any
parameter/switch being set/read is valid.  Except where specifically
noted, StandAloneBlast uses the same single-letter, case-sensitive
parameter names as the actual blast program.  Currently no checks are
included to verify that parameters are of the proper type (e.g. string
or numeric) or that their values are within the proper range.

As an example, to change the value of the Blast parameter 'e' ('e' is
the parameter for expectation-value cutoff) 

  $expectvalue = 0.01;
  $factory->e($expectvalue);

Note that for improved script readibility one can modify the name of
the BLAST parameters as desired as long as the initial letter (and
case) of the parameter are preserved, e.g.:

  $factory->expectvalue($expectvalue);

Unfortunately, some of the BLAST parameters are not the single 
letter one might expect (eg "iteration round" in blastpgp is 'j'). 
Again one can check by using, for example:

  > blastpgp - .

Once the factory has been created and the appropriate parameters set,
one can call one of the supported blast executables.  The input
sequence(s) to these executables may be fasta file(s) as described in
the BLAST documentation.

  $inputfilename = 't/testquery.fa';
  $blast_report = $factory->blastall($inputfilename);

In addition, sequence input may be in the form of either a Bio::Seq
object or or an array of Bio::Seq objects, e.g.:

  $input = Bio::Seq->new(-id => "test query",
                         -seq => "ACTACCCTTTAAATCAGTGGGGG");
  $blast_report = $factory->blastall($input);

For blastall and non-psiblast blastpgp runs, report object is either a
L<Bio::Tools::BPlite> or L<Bio::SearchIO> object, selected by the user 
with the parameter _READMETHOD.  The leading underscore is needed to
distinguish this option from options which are passed to the BLAST
executable. The default parser is Bio::SearchIO::blast.  If BPlite
method is selected, L<Bio::Tools::BPlite> objects will be returned for
standard blast and L<Bio::Tools::BPpsilite> for a multiple-iteration
blasts, and a L<Bio::Tools::BPbl2seq> for bl2seq.  In any case, the "raw"
blast report is also available. The filename is set by the in the
'outfile' parameter and has the default value of "blastreport.out".
The BPlite method is only provided to support legacy code since
the BPlite modules are no longer maintained - do not use BPlite
since these modules will be removed eventually.

For psiblast execution in the BLAST "jumpstart" mode, the program must
be passed (in addition to the query sequence itself) an alignment
containing the query sequence (in the form of a SimpleAlign object) as
well as a "mask" specifying at what residues position-specific scoring
matrices (PSSMs) are to used and at what residues default scoring
matrices (eg BLOSUM) are to be used. See psiblast documentation for
more details.  The mask itself is a string of 0's and 1's which is the
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

For more examples of syntax and use of Blast.pm, the user is
encouraged to run the scripts standaloneblast.pl in the bioperl
examples/tools directory and StandAloneBlast.t in the bioperl t/ 
directory.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists.  Your participation is much appreciated.

  bioperl-l@bioperl.org               - General discussion
  http://bio.perl.org/MailList.html   - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via 
the web:

  http://bugzilla.bioperl.org/

=head1 AUTHOR -  Peter Schattner

Email schattner at alum.mit.edu

=head1 MAINTAINER - Torsten Seemann

Email torsten at infotech.monash.edu.au

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

package Bio::Tools::Run::StandAloneBlast;

use vars qw($AUTOLOAD @ISA $PROGRAMDIR  $DATADIR $BLASTTYPE
	    @BLASTALL_PARAMS @BLASTPGP_PARAMS @WUBLAST_PARAMS 
		 @WUBLAST_SWITCH @RPSBLAST_PARAMS @BL2SEQ_PARAMS 
       @OTHER_PARAMS %OK_FIELD $DEFAULTREADMETHOD
	    );
		 
use strict;
use Bio::Root::Root;
use Bio::Root::IO;
use Bio::Seq;
use Bio::SeqIO;
use Bio::Tools::BPbl2seq;
use Bio::Tools::BPpsilite;
use Bio::SearchIO;
use Bio::Tools::Run::WrapperBase;
use Bio::Factory::ApplicationFactoryI;

BEGIN {
        @BLASTALL_PARAMS = qw(A B C D E F G I J K L M O P Q R S T U V W X Y Z a b d e f g i l m n o p q r s t v w y z);
        @BLASTPGP_PARAMS = qw(A B C E F G H I J K L M N O P Q R S T U W X Y Z a b c d e f h i j k l m o p q s t u v y z);
        @RPSBLAST_PARAMS = qw(F I J L N O P T U V X Y Z a b d e i l m o p v y z);
        @BL2SEQ_PARAMS = qw(A D E F G I J M S T U V W X Y a d e g i j m o p q r t);
	$DEFAULTREADMETHOD = 'BLAST';
	$BLASTTYPE = 'ncbi';
	@WUBLAST_PARAMS = 
	  qw( E S E2 S2 W T X M Y Z L K H V  B
     matrix Q R filter wordmask filter maskextra 
     hitdist wink ctxfactor gapE gapS gapE2 gapS2 gapW gapX olf golf 
     olmax golmax gapdecayrate topcomboN topcomboE sumstatsmethod 
     hspsepqmax hspsepsmax gapsepqmax gapsepsmax altscore hspmax gspmax 
     qoffset nwstart nwlen qrecmin qrecmax 
     dbrecmin dbrecmax vdbdescmax dbchunks sort_by_pvalue 
     cpus putenv getenv progress o database input);
	@WUBLAST_SWITCH = 
     qw(kap sump poissonp lcfilter lcmask echofilter stats nogap 
     gapall pingpong 
     nosegs postsw span2 span1 span prune consistency 
     links ucdb gi noseqs qtype qres sort_by_pvalue sort_by_count 
     sort_by_highscore sort_by_totalscore sort_by_subjectlength
     mmio nonnegok novalidctxok shortqueryok notes warnings 
     errors endputenv 
     getenv endgetenv abortonerror abortonfatal); 

	# Non BLAST parameters start with underscore to differentiate them
	# from BLAST parameters
	@OTHER_PARAMS = qw(_READMETHOD);

	# my @other_switches = qw(QUIET);

	# Authorize attribute fields
	foreach my $attr (@BLASTALL_PARAMS, @BLASTPGP_PARAMS, @RPSBLAST_PARAMS,
							@BL2SEQ_PARAMS, @OTHER_PARAMS ,@WUBLAST_PARAMS, 
                     @WUBLAST_SWITCH )
     { $OK_FIELD{$attr}++; }

        # You will need to enable Blast to find the Blast program.  
        # This can be done in at least two different ways:
        #  1. define an environmental variable blastDIR:
        #	export BLASTDIR=/home/peter/blast   or
        #  2. include a definition of an environmental variable 
        # BLASTDIR in every script that will use StandAloneBlast.pm.

	$PROGRAMDIR = $BLASTTYPE eq 'ncbi' ? $ENV{'BLASTDIR'}: $ENV{'WUBLASTDIR'};

	# If local BLAST databases are not stored in the standard
	# /data directory, the variable BLASTDATADIR will need to be 
	# set explicitly 
	$DATADIR =  $ENV{'BLASTDATADIR'} || $ENV{'BLASTDB'} || '';
}

@ISA = qw(Bio::Root::Root 
	  Bio::Tools::Run::WrapperBase 
	  Bio::Factory::ApplicationFactoryI);

=head1 BLAST parameters

Essentially all BLAST parameter can be set via StandAloneBlast.pm.
Some of the most commonly used parameters are listed below.  All
parameters have defaults and are optional (I think.)  For a complete
listing of settable parameters, run the relevant executable BLAST
program with the option "-" as in blastall -

=head2 Blastall

  -p  Program Name [String]
        Input should be one of "blastp", "blastn", "blastx", 
        "tblastn", or "tblastx".
  -d  Database [String] default = nr
        The database specified must first be formatted with formatdb.
        Multiple database names (bracketed by quotations) will be accepted.
        An example would be -d "nr est"
   -i  Query File [File In]   Set by StandAloneBlast.pm from script.
    default = stdin. The query should be in FASTA format.  If multiple FASTA entries are in the input
        file, all queries will be searched.
  -e  Expectation value (E) [Real] default = 10.0
  -o  BLAST report Output File [File Out]  Optional,
	default = ./blastreport.out ; set by StandAloneBlast.pm		
  -S  Query strands to search against database (for blast[nx], and tblastx).  3 is both, 1 is top, 2 is bottom [Integer]
	default = 3

=head2 Blastpgp (including Psiblast)

  -j   is the maximum number of rounds (default 1; i.e., regular BLAST)
  -h   is the e-value threshold for including sequences in the
	score matrix model (default 0.001)
  -c   is the "constant" used in the pseudocount formula specified in the paper (default 10)
  -B  Multiple alignment file for PSI-BLAST "jump start mode"  Optional
  -Q  Output File for PSI-BLAST Matrix in ASCII [File Out]  Optional

=head2 rpsblast

  -d  Database [String] default = (none - you must specify a database)
        The database specified must first be formatted with formatdb.
        Multiple database names (bracketed by quotations) will be accepted.
        An example would be -d "Cog Smart"
   -i  Query File [File In]   Set by StandAloneBlast.pm from script.
    default = stdin. The query should be in FASTA format.  If multiple FASTA entries are in the input
        file, all queries will be searched.
  -e  Expectation value (E) [Real] default = 10.0
  -o  BLAST report Output File [File Out]  Optional,
	default = ./blastreport.out ; set by StandAloneBlast.pm		

=head2 Bl2seq

  -i  First sequence [File In]
  -j  Second sequence [File In]
  -p  Program name: blastp, blastn, blastx. For blastx 1st argument should be nucleotide [String]
    default = blastp
  -o  alignment output file [File Out] default = stdout
  -e  Expectation value (E) [Real]  default = 10.0
  -S  Query strands to search against database (blastn only).  3 is both, 1 is top, 2 is bottom [Integer]
    default = 3

=head2 WU-Blast

  -p Program Name [String] 
        Input should be one of "wublastp", "wublastn", "wublastx", 
        "wutblastn", or "wutblastx".
  -d  Database [String] default = nr
        The database specified must first be formatted with xdformat.
  -i  Query File [File In]   Set by StandAloneBlast.pm from script.
    default = stdin. The query should be in FASTA format.  If multiple FASTA entries are in the input
        file, all queries will be searched.
  -E  Expectation value (E) [Real] default = 10.0
  -o  BLAST report Output File [File Out]  Optional,
	default = ./blastreport.out ; set by StandAloneBlast.pm		

=cut

sub new {
    my ($caller, @args) = @_;
    # chained new
    my $self = $caller->SUPER::new(@args);

    # to facilitiate tempfile cleanup
    my ($tfh,$tempfile) = $self->io->tempfile();
    close($tfh); # we don't want the filehandle, just a temporary name
    $self->o($tempfile) unless $self->o;
    $self->_READMETHOD($DEFAULTREADMETHOD);
    while (@args)  {
	my $attr =   shift @args;
    	my $value =  shift @args;
    	next if( $attr eq '-verbose');
    	# we allow both 'attr' and '-attr' options on the new() call
    	$attr =~ s/^-//;
    	# the workaround to deal with initializing
	if($attr =~/^\s*program\s*$|^p$/){
	    if($value =~/^wu*/){
		$BLASTTYPE="wublast";
	    }
	    $attr = 'p';
	}
	if($attr =~/outfile/){
	    $attr = 'o';
	}

	$self->$attr($value);
    }
    return $self;
}

=head2 quiet

 Title   : quiet
 Usage   : $obj->quiet($newval)
 Function: 
 Example : 
 Returns : value of quiet (a scalar)
 Args    : on set, new value (a scalar or undef, optional)


=cut

sub quiet{
    my $self = shift;
    return $self->{'_quiet'} = shift if @_;
    return $self->{'_quiet'};
}

sub AUTOLOAD {
    my $self = shift;
    my $attr = $AUTOLOAD;
    $attr =~ s/.*:://;    
    my $attr_letter = $BLASTTYPE eq 'ncbi' ? substr($attr, 0, 1) : $attr;

    # actual key is first letter of $attr unless first attribute
    # letter is underscore (as in _READMETHOD), the $attr is a BLAST
    # parameter and should be truncated to its first letter only
    $attr = ($attr_letter eq '_') ? $attr : $attr_letter;
    $self->throw("Unallowed parameter: $attr !") unless $OK_FIELD{$attr};
    $self->{$attr_letter} = shift if @_;
    return $self->{$attr_letter};
}

=head1 Methods

=head2 executable

 Title   : executable
 Usage   : my $exe = $blastfactory->executable('blastall');
 Function: Finds the full path to the 'codeml' executable
 Returns : string representing the full path to the exe
 Args    : [optional] name of executable to set path to 
           [optional] boolean flag whether or not warn when exe is not found


=cut

sub executable {
   my ($self, $exename, $exe,$warn) = @_;
   $exename = 'blastall' unless (defined $exename || $BLASTTYPE ne 'ncbi');

   if( defined $exe && -x $exe ) {
     $self->{'_pathtoexe'}->{$exename} = $exe;
   }
   unless( defined $self->{'_pathtoexe'}->{$exename} ) {
       my $f = $self->program_path($exename);	    
       $exe = $self->{'_pathtoexe'}->{$exename} = $f if(-e $f && -x $f );
        
       #  This is how I meant to split up these conditionals --jason
       # if exe is null we will execute this (handle the case where
       # PROGRAMDIR pointed to something invalid)
       unless( $exe )  {  # we didn't find it in that last conditional
	   if( ($exe = $self->io->exists_exe($exename)) && -x $exe ) {
	       $self->{'_pathtoexe'}->{$exename} = $exe;
	   } else { 
	       $self->warn("Cannot find executable for $exename") if $warn;
	       $self->{'_pathtoexe'}->{$exename} = undef;
	   }
       }
   }
   return $self->{'_pathtoexe'}->{$exename};
}


=head2 program_path

 Title   : program_path
 Usage   : my $path = $factory->program_path();
 Function: Builds path for executable 
 Returns : string representing the full path to the exe
 Args    : none

=cut

sub program_path {
    my ($self,$program_name) = @_;
    my @path;
    push @path, $self->program_dir if $self->program_dir;
    push @path, $program_name .($^O =~ /mswin/i ?'.exe':'');

    return Bio::Root::IO->catfile(@path);
}

=head2 program_dir

 Title   : program_dir
 Usage   : my $dir = $factory->program_dir();
 Function: Abstract get method for dir of program. 
 Returns : string representing program directory 
 Args    : none 

=cut

sub program_dir {
    $PROGRAMDIR;
}

sub program {
    my $self = shift;
    if( wantarray ) {
	return ($self->executable, $self->p());
    } else {
	return $self->executable(@_);
    }
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
	      $blast_report = $factory->blastall(\@seq_array);
 Returns : Reference to a Blast object or BPlite object 
           containing the blast report.
 Args    : Name of a file or Bio::Seq object or an array of 
           Bio::Seq object containing the query sequence(s). 
           Throws an exception if argument is not either a string 
           (eg a filename) or a reference to a Bio::Seq object 
           (or to an array of Seq objects).  If argument is string, 
           throws exception if file corresponding to string name can 
           not be found.

=cut

sub blastall {
    my ($self,$input1) = @_;
    $self->io->_io_cleanup();
    my $executable = 'blastall';
    my $input2;
# Create input file pointer
    my $infilename1 = $self->_setinput($executable, $input1);
    if (! $infilename1) {$self->throw(" $input1 ($infilename1) not Bio::Seq object or array of Bio::Seq objects or file name!");}

    $self->i($infilename1);	# set file name of sequence to be blasted to inputfilename1 (-i param of blastall)
    
    my $blast_report = &_generic_local_blast($self, $executable, 
					     $input1, $input2);
}

=head2  wublast

 Title   : wublast
 Usage   :  $blast_report = $factory->wublast('t/testquery.fa');
	or
	       $input = Bio::Seq->new(-id=>"test query",
				      -seq=>"ACTACCCTTTAAATCAGTGGGGG");
	       $blast_report = $factory->wublast($input);
	or 
	      $seq_array_ref = \@seq_array;  # where @seq_array is an array of Bio::Seq objects
	      $blast_report = $factory->wublast(\@seq_array);
 Returns :  Reference to a Blast object 
 Args    : Name of a file or Bio::Seq object or an array of 
           Bio::Seq object containing the query sequence(s). 
           Throws an exception if argument is not either a string 
           (eg a filename) or a reference to a Bio::Seq object 
           (or to an array of Seq objects).  If argument is string, 
           throws exception if file corresponding to string name can 
           not be found.

=cut

sub wublast {
  my ($self,$input1) = @_;
  $self->io->_io_cleanup();
  my $executable = 'wublast';
  my $infilename1 = $self->_setinput($executable, $input1);
  if (! $infilename1) {$self->throw(" $input1 ($infilename1) not Bio::Seq object or array of Bio::Seq objects or file name!");}
  $self->input($infilename1);	# set file name of sequence to be blasted to inputfilename1 (-i param of blastall)
  my $blast_report = &_generic_local_wublast($self, $executable, $input1);
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
 Returns : Reference to a Bio::SearchIO object or BPlite object 
           containing the blast report (BPlite only if you specify 
           _READMETHOD=> 'BPlite')
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
 Returns : Reference to Bio::SearchIO object 
           or Bio::Tools::BPpsilite if you specify 
           _READMETHOD => 'BPlite' object containing the blast report.

=cut

sub blastpgp {
    my $self = shift;
    my $executable = 'blastpgp';
    my $input1 = shift;
    my $input2 = shift;
	 # used by blastpgp's -B option to specify which 
	 # residues are position aligned
    my $mask = shift;

    my  ($infilename1, $infilename2 )  = $self->_setinput($executable, 
							  $input1, $input2, 
							  $mask);
    if (!$infilename1) {$self->throw(" $input1  not Bio::Seq object or array of Bio::Seq objects or file name!");}
    $self->i($infilename1);	# set file name of sequence to be blasted to inputfilename1 (-i param of blastpgp)
    if  ($input2) {
	unless ($infilename2) {$self->throw("$input2 not SimpleAlign Object in pre-aligned psiblast\n");}
	$self->B($infilename2);	# set file name of partial alignment to inputfilename2 (-B param of blastpgp)
    }
    my $blast_report = &_generic_local_blast($self, $executable, $input1, $input2);
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
 Returns : Reference to a Bio::SearchIO object or BPlite object 
           containing the blast report (BPlite only if you specify 
           _READMETHOD=> 'BPlite')

=cut

sub rpsblast {
    my ($self, $input1) = @_;
	
    my $executable = 'rpsblast';

    # Create input file pointer
    my $infilename1 = $self->_setinput($executable, $input1);
    if (! $infilename1) { 
	   $self->throw(" $input1 ($infilename1) not Bio::Seq object or array of Bio::Seq objects or file name!");
	 }
    $self->i($infilename1);	# set file name of sequence to be blasted to inputfilename1 (-i param of blastall)
    
	 # Run like a standard NCBI blast from this point
    my $blast_report = _generic_local_blast($self, $executable);
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
    my  ($infilename1, $infilename2 )  = $self->_setinput($executable, 
							  $input1, $input2);
    if (!$infilename1){$self->throw(" $input1  not Seq Object or file name!");}
    if (!$infilename2){$self->throw("$input2  not Seq Object or file name!");}

    $self->i($infilename1);	# set file name of first sequence to 
                                # be aligned to inputfilename1 
                                # (-i param of bl2seq)
    $self->j($infilename2);	# set file name of first sequence to 
                                # be aligned to inputfilename2 
                                # (-j param of bl2seq)
    my $blast_report = &_generic_local_blast($self, $executable);    
}
#################################################

=head2  _generic_local_blast

 Title   : _generic_local_blast
 Usage   : internal function not called directly
 Returns : Bio::SearchIO or Bio::Tools::BPlite object
 Args    : Reference to calling object and name of BLAST executable 

=cut

sub _generic_local_blast {
    my $self = shift;
    my $executable = shift;

    # Create parameter string to pass to Blast program
    my $param_string = $self->_setparams($executable);

    # run Blast
    my $blast_report = &_runblast($self, $executable, $param_string);
}


=head2  _generic_local_wublast

 Title   : _generic_local_wublast
 Usage   :  internal function not called directly
 Returns :  Blast object
 Args    :   Reference to calling object and name of BLAST executable 

=cut

sub _generic_local_wublast {
    my $self = shift;
    my $executable = shift;

    # Create parameter string to pass to Blast program
    my $param_string = $self->_setparams($executable);
    $param_string = " ".$self->database." ".$self->input." ".$param_string;

    # run Blast
    my $blast_report = &_runwublast($self, $executable, $param_string);
}

=head2  _runblast

 Title   :  _runblast
 Usage   :  Internal function, not to be called directly	
 Function:   makes actual system call to Blast program
 Example :
 Returns : Report object in the appropriate format (Bio::SearchIO)
           or if BPlite is requested: Bio::Tools::BPlite, 
           Bio::Tools::BPpsilite,or Bio::Tools::BPbl2seq)
 Args    : Reference to calling object, name of BLAST executable, 
           and parameter string for executable 

=cut

sub _runblast {
	my ($self,$executable,$param_string) = @_;
	my ($blast_obj,$exe);
	if( ! ($exe = $self->executable($executable)) ) {
		$self->warn("cannot find path to $executable");
		return;
	}
	my $commandstring = $exe. $param_string;

	# next line for debugging
	$self->debug( "$commandstring \n");

	my $status = system($commandstring);

	$self->throw("$executable call crashed: $? $commandstring\n")  
	  unless ($status==0) ;
	my $outfile = $self->o() ;	# get outputfilename
	my $signif = $self->e()  || 1e-5  ; 

	# set significance cutoff to set expectation value or default value
	# (may want to make this value vary for different executables)

	# If running bl2seq or psiblast (blastpgp with multiple iterations),
	# the specific parsers for these programs must be used (ie BPbl2seq or
	# BPpsilite).  Otherwise either the Blast parser or the BPlite
	# parsers can be selected.

	if ($self->_READMETHOD =~ /^(Blast|SearchIO)/i )  {
		$blast_obj = Bio::SearchIO->new(-file=>$outfile,
												  -format => 'blast' )  ;
	} elsif( $self->_READMETHOD =~ /BPlite/i ) {
		if ($executable =~ /bl2seq/i)  {
			# Added program info so BPbl2seq can compute strand info
			$blast_obj = Bio::Tools::BPbl2seq->new(-file => $outfile,
																-REPORT_TYPE => $self->p );
		} elsif ($executable =~ /blastpgp/i && defined $self->j() && 
					$self->j() > 1)  {
			$self->debug( "using psilite parser\n");
			$blast_obj = Bio::Tools::BPpsilite->new(-file => $outfile);
		} elsif( $executable =~ /blastall|rpsblast/i) { 
			$blast_obj = Bio::Tools::BPlite->new(-file=>$outfile);
		} else { 
			$self->warn("Unrecognized executable $executable");
		}
	} else {
		$self->warn("Unrecognized readmethod ".$self->_READMETHOD);
	}

	return $blast_obj;
}

=head2  _runwublast

 Title   :  _runwublast
 Usage   :  Internal function, not to be called directly	
 Function:   makes actual system call to WU-Blast program
 Example :
 Returns : Report Blast object
 Args    : Reference to calling object, name of BLAST executable, 
           and parameter string for executable 

=cut

sub _runwublast {
	my ($self,$executable,$param_string) = @_;
	my ($blast_obj,$exe);
	if( ! ($exe = $self->executable($self->p))){
            $self->warn("cannot find path to $executable");
            return;
	}
	my $commandstring = $exe.  " ".$param_string;

	# next line for debugging
	$self->debug( "$commandstring \n");

	my $status = system($commandstring);

	$self->throw("$executable call crashed: $? $commandstring\n")  
	  unless ($status==0) ;
	my $outfile = $self->o() ;	# get outputfilename
	$blast_obj = Bio::SearchIO->new(-file=>$outfile,
		                          			-format => 'blast') ;
	return $blast_obj;
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

 SWITCH:  {
      unless (ref $input1) {
			$infilename1 = (-e $input1) ? $input1 : 0 ;
			last SWITCH; 
      }
		#  $input may be an array of BioSeq objects...
      if (ref($input1) =~ /ARRAY/i ) {
			($fh,$infilename1) = $self->io->tempfile();
			$temp =  Bio::SeqIO->new(-fh=> $fh, 
											 -format => 'fasta');
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
      $infilename1 = 0;		# Set error flag if you get here
	}				# End SWITCH
	unless ($input2) { return $infilename1; }
 SWITCH2:  {
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
      elsif ($input2->isa("Bio::SimpleAlign")  && 
				 $executable  eq 'blastpgp' ) {
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
	  #  Now we need to write out the alignment to a file 
	  # in the "psi format" which psiblast is expecting
	  $input2->map_chars('\.','-');
	  $temp =  Bio::AlignIO->new(-fh=> $fh, '-format' => 'psi');
	  $temp->write_aln($input2);
	  close $fh;
	  undef $fh;
	  last SWITCH2;
  }
      $infilename2 = 0;		# Set error flag if you get here
	}				# End SWITCH2
	return ($infilename1, $infilename2);
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
    my ($self,$executable) = @_;
    my ($attr, $value, @execparams);

    if    ($executable eq 'blastall') { @execparams = @BLASTALL_PARAMS; }
    elsif ($executable eq 'blastpgp') { @execparams = @BLASTPGP_PARAMS; }
    elsif ($executable eq 'rpsblast') { @execparams = @RPSBLAST_PARAMS; }
    elsif ($executable eq 'bl2seq'  ) { @execparams = @BL2SEQ_PARAMS;   }
    elsif ($executable eq 'wublast' ) { @execparams = @WUBLAST_PARAMS;  }

    my $param_string = "";
    for $attr ( @execparams ) {
        $value = $self->$attr();
        next unless (defined $value);
        # Need to prepend datadirectory to database name
        if ($executable eq 'wublast') {
          next if $attr =~ /database|^d$/;
          next if $attr =~ /input|^i$/;
          $attr = 'o' if ($attr =~/outfile/);
        }

	if ($attr  eq 'd' && ($executable ne 'bl2seq')) { 
	    my @dbs = split(/ /, $value);
	    for (my $i = 0; $i < scalar(@dbs); $i++) {
		# moved the test for full path db to work with multiple databases
		if (! (-e $dbs[$i].".nin" || -e $dbs[$i].".pin") &&
		    ! (-e $dbs[$i].".nal" || -e $dbs[$i].".pal") ) {
		    $dbs[$i] = File::Spec->catdir($DATADIR, $dbs[$i]);
		}
		$value = '"'.join(" ", @dbs).'"';
	    }
	}
# put params in format expected by Blast
	$attr  = '-'. $attr ;       
	$param_string .= " $attr  $value ";
    }

    if ($self->quiet()) { 
      $param_string .= '  2>/dev/null';
    }
    if ($executable eq 'wublast') {
	foreach my $attr (@WUBLAST_SWITCH) {
	    my $value = $self->$attr();
	    next unless (defined $value);
	    my $attr_key = ' -'.(lc $attr);
	    $param_string .=$attr_key;
	}
    }
    return $param_string;
}


=head1 Bio::Tools::Run::Wrapper methods

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

__END__
