# $Id$
#
# BioPerl module for Bio::Tools::Run::Alignment::Clustalw
#
# Cared for by
#
# Copyright Peter Schattner
#
# You may distribute this module under the same terms as perl itself
# POD documentation - main docs before the code

=head1 NAME

Bio::Tools::Run::Alignment::Clustalw - Object for the calculation of a
multiple sequence alignment from a set of unaligned sequences or
alignments using the Clustalw program

=head1 SYNOPSIS

  #  Build a clustalw alignment factory
  @params = ('ktuple' => 2, 'matrix' => 'BLOSUM');
  $factory = Bio::Tools::Run::Alignment::Clustalw->new(@params);

  #  Pass the factory a list of sequences to be aligned.	
  $inputfilename = 't/cysprot.fa';
  $aln = $factory->align($inputfilename); # $aln is a SimpleAlign object.
  # or
  $seq_array_ref = \@seq_array;
  # where @seq_array is an array of Bio::Seq objects
  $aln = $factory->align($seq_array_ref);

  # Or one can pass the factory a pair of (sub)alignments
  #to be aligned against each other, e.g.:
  $aln = $factory->profile_align($aln1,$aln2);
  # where $aln1 and $aln2 are Bio::SimpleAlign objects.

  # Or one can pass the factory an alignment and one or more unaligned
  # sequences to be added to the alignment. For example: 	
  $aln = $factory->profile_align($aln1,$seq); # $seq is a Bio::Seq object.

There are various additional options and input formats available.  See
the DESCRIPTION section that follows for additional details.

=head1 DESCRIPTION

Note: this DESCRIPTION only documents the (Bio)perl interface to
Clustalw.  Clustalw, itself, is a large & complex program - for more
information regarding clustalw, please see the clustalw documentation
which accompanies the clustalw distribution. Clustalw is available
from (among others) ftp://ftp.ebi.ac.uk/pub/software/. Clustalw.pm has
been tested so far only under Linux. I expect that it should also work
under other Unix systems.  However, since the module is currently
implemented using (unix) system calls, extensive modification may be
necessary before Clustalw.pm would work under non-Unix operating
systems (eg Windows, MacOS).  Clustalw.pm has only been tested using
version 1.8 of clustalw.  Compatibility with earlier versions of the
clustalw program is currently unknown. Before running Clustalw.pm
successfully it will be necessary: to install clustalw on your system,
to edit the variable $clustdir in Clustalw.pm to point to the clustalw
program, and to ensure that users have execute privilieges for the
clustalw program.

Bio::Tools::Run::Alignment::Clustalw.pm: is an object for performing a
multiple sequence alignment from a set of unaligned sequences and/or
sub-alignments by means of the clustalw program.

Initially, a clustalw "factory object" is created. Optionally, the
factory may be passed most of the parameters or switches of the
clustalw program, e.g.:

	@params = ('ktuple' => 2, 'matrix' => 'BLOSUM');
	$factory = Bio::Tools::Run::Alignment::Clustalw->new(@params);

Any parameters not explicitly set will remain as the defaults of the
clustalw program.  Additional parameters and switches (not available
in clustalw) may also be set.  Currently, the only such parameter is
"quiet", which when set to a non-zero value, suppresses clustalw
terminal output. Not all clustalw parameters are supported at this
stage.

By default, Clustalw.pm output is returned solely in a the form of a
BioPerl Bio::SimpleAlign object which can then be printed and/or saved
in multiple formats using the AlignIO.pm module. Optionally the raw
clustalw output file can be saved if the calling script specifies an
output file (with the clustalw parameter OUTFILE).  Currently only the
GCG-MSF output file formats is supported.

Other parameters and features (such as those corresponding to tree
production) have not been implemented yet in Perl format.

Alignment parameters can be changed and/or examined at any time after
the factory has been created.  The program checks that any
parameter/switch being set/read is valid.  However, currently no
additional checks are included to check that parameters are of the
proper type (eg string or numeric) or that their values are within the
proper range.  As an example, to change the value of the clustalw
parameter ktuple to 3 and subsequently to check its value one would
write:

	$ktuple = 3;
	$factory->ktuple($ktuple);
 	$get_ktuple = $factory->ktuple();


Once the factory has been created and the appropriate parameters set,
one can call the method align() to align a set of unaligned sequences,
or call profile_align() to add one or more sequences or a second
alignment to an initial alignment.

Input to align() may consist of a set of unaligned sequences in the
form of the name of file containing the sequences. For example,
$inputfilename = 't/cysprot.fa'; $aln =
$factory-E<gt>align($inputfilename);

Alternately one can create an array of Bio::Seq objects somehow

	$str = Bio::SeqIO->new(-file=> 't/cysprot.fa', '-format' => 'Fasta');
	@seq_array =();
	while ( my $seq = $str->next_seq() ) {push (@seq_array, $seq) ;}

and pass the factory a reference to that array

	$seq_array_ref = \@seq_array;
	$aln = $factory->align($seq_array_ref);

In either case, align() returns a reference to a SimpleAlign object
which can then be displayed, stored, or converted to a UnivAlign
object for further manipulation.

Once an initial alignment exists, one can pass the factory additional
sequence(s) to be added (ie aligned) to the original alignment.  The
alignment can be passed as either an alignment file or a
Bio:SimpleAlign object.  The unaligned sequence(s) can be passed as a
filename or as an array of BioPerl sequence objects or as a single
BioPerl Seq object.  For example (to add a single sequence to an
alignment),

	$str = Bio::AlignIO->new(-file=> 't/cysprot1a.msf');
	$aln = $str->next_aln();
	$str1 = Bio::SeqIO->new(-file=> 't/cysprot1b.fa');
	$seq = $str1->next_seq();
	$aln = $factory->profile_align($aln,$seq);

In either case, profile_align() returns a reference to a SimpleAlign
object containing a new SimpleAlign object of the alignment with the
additional sequence(s) added in.

Finally one can pass the factory a pair of (sub)alignments to be
aligned against each other.  The alignments can be passed in the form
of either a pair of alignment files or a pair of Bio:SimpleAlign
objects. For example,

	$profile1 = 't/cysprot1a.msf';
	$profile2 = 't/cysprot1b.msf';
	$aln = $factory->profile_align($profile1,$profile2);
or
	$str1 = Bio::AlignIO->new(-file=> 't/cysprot1a.msf');
	$aln1 = $str1->next_aln();
	$str2 = Bio::AlignIO->new(-file=> 't/cysprot1b.msf');
	$aln2 = $str2->next_aln();
	$aln = $factory->profile_align($aln1,$aln2);

In either case, profile_align() returns a reference to a SimpleAlign
object containing an (super)alignment of the two input alignments.

For more examples of syntax and use of Clustalw.pm, the user is
encouraged to run the script Clustalw.t is the bioperl/t directory.

Note: Clustalw.pm is still under development. Various features of the
clustalw program have not yet been implemented.  If you would like
that a specific clustalw feature be added to this perl interface, let
me know.

These can be specified as paramters when instantiating a new TCoffee
object, or through get/set methods of the same name (lowercase).

=head1 PARAMETER FOR ALIGNMENT COMPUTATION

=head2 KTUPLE

 Title       : KTUPLE
 Description : (optional) set the word size to be used in the alignment
               This is the size of exactly matching fragment that is used.
               INCREASE for speed (max= 2 for proteins; 4 for DNA),
               DECREASE for sensitivity.
               For longer sequences (e.g. >1000 residues) you may
               need to increase the default

=head2 TOPDIAGS

 Title       : TOPDIAGS
 Description : (optional) number of best diagonals to use
               The number of k-tuple matches on each diagonal
               (in an imaginary dot-matrix plot) is calculated.
               Only the best ones (with most matches) are used in
               the alignment.  This parameter specifies how many.
               Decrease for speed; increase for sensitivity.

=head2 WINDOW

 Title       : WINDOW
 Description : (optional) window size
               This is the number of diagonals around each of the 'best'
               diagonals that will be used.  Decrease for speed;
               increase for sensitivity.

=head2 PAIRGAP

 Title       : PAIRGAP
 Description : (optional) gap penalty for pairwise alignments
               This is a penalty for each gap in the fast alignments.
               It has little affect on the speed or sensitivity except
               for extreme values.

=head2 FIXEDGAP

 Title       : FIXEDGAP
 Description : (optional) fixed length gap penalty

=head2 FLOATGAP

 Title       : FLOATGAP
 Description : (optional) variable length gap penalty

=head2 MATRIX

 Title       : MATRIX
 Default     : PAM100 for DNA - PAM250 for protein alignment
 Description : (optional) substitution matrix used in the multiple
               alignments. Depends on the version of clustalw as to
               what default matrix will be used

               PROTEIN WEIGHT MATRIX leads to a new menu where you are
               offered a choice of weight matrices. The default for
               proteins in version 1.8 is the PAM series derived by
               Gonnet and colleagues. Note, a series is used! The
               actual matrix that is used depends on how similar the
               sequences to be aligned at this alignment step
               are. Different matrices work differently at each
               evolutionary distance.

               DNA WEIGHT MATRIX leads to a new menu where a single
               matrix (not a series) can be selected. The default is
               the matrix used by BESTFIT for comparison of nucleic
               acid sequences.

=head2 TYPE

 Title       : TYPE
 Description : (optional) sequence type: protein or DNA. This allows
	       you to explicitly overide the programs attempt at
	       guessing the type of the sequence.  It is only useful
	       if you are using sequences with a VERY strange
	       composition.

=head2 OUTPUT

 Title       : OUTPUT
 Description : (optional) clustalw supports GCG or PHYLIP or PIR or
                Clustal format.  See the Bio::AlignIO modules for
                which formats are supported by bioperl.

=head2 OUTFILE

 Title       : OUTFILE
 Description : (optional) Name of clustalw output file. If not set
	       module will erase output file.  In any case alignment will
	       be returned in the form of SimpleAlign objects

=head2 TRANSMIT

 Title       : TRANSMIT
 Description : (optional) transitions not weighted.  The default is to
	       weight transitions as more favourable than other
	       mismatches in DNA alignments.  This switch makes all
	       nucleotide mismatches equally weighted.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists.  Your participation is much appreciated.

  bioperl-l@bioperl.org          - General discussion
  http://bio.perl.org/MailList.html             - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
 the bugs and their resolution.  Bug reports can be submitted via
 email or the web:

  bioperl-bugs@bio.perl.org
  http://bio.perl.org/bioperl-bugs/

=head1 AUTHOR -  Peter Schattner

Email schattner@alum.mit.edu

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

#'

package Bio::Tools::Run::Alignment::Clustalw;

use vars qw($AUTOLOAD @ISA $DEBUG $PROGRAM $PROGRAMDIR
	    $TMPDIR $TMPOUTFILE @CLUSTALW_SWITCHES @CLUSTALW_PARAMS
	    @OTHER_SWITCHES %OK_FIELD);
use strict;
use Bio::Seq;
use Bio::SeqIO;
use Bio::SimpleAlign;
use Bio::AlignIO;
use Bio::Root::RootI;
use Bio::Root::IO;

@ISA = qw(Bio::Root::RootI Bio::Root::IO);

BEGIN {

# You will need to enable Clustalw to find the clustalw program. This
# can be done in (at least) three ways:

# 1. Modify your $PATH variable to include your clustalw directory as
# in (for Linux):
# export PATH=$PATH:/home/peter/clustalw1.8
#
# 2. define an environmental variable CLUSTALDIR:
# export CLUSTALDIR=/home/peter/clustalw1.8
#
# 3. include a definition of an environmental variable CLUSTALDIR in
# every script that will use Clustal.pm.
# $ENV{CLUSTALDIR} = '/home/peter/clustalw1.8/';

    $PROGRAMDIR = $ENV{CLUSTALDIR} || '';
    $PROGRAM = Bio::Root::IO->catfile($PROGRAMDIR,'clustalw');

    @CLUSTALW_PARAMS = qw(KTUPLE TOPDIAGS WINDOW PAIRGAP FIXEDGAP
                   FLOATGAP MATRIX TYPE	TRANSIT DNAMATRIX OUTFILE
                   GAPOPEN GAPEXT MAXDIV GAPDIST HGAPRESIDUES PWMATRIX
                   PWDNAMATRIX PWGAPOPEN PWGAPEXT SCORE TRANSWEIGHT
                   SEED HELIXGAP OUTORDER STRANDGAP LOOPGAP TERMINALGAP
                   HELIXENDIN HELIXENDOUT STRANDENDIN STRANDENDOUT);

    @CLUSTALW_SWITCHES = qw(HELP CHECK OPTIONS NEGATIVE NOWEIGHTS ENDGAPS
                        NOPGAP NOHGAP NOVGAP KIMURA TOSSGAPS);

    @OTHER_SWITCHES = qw(QUIET);
    # Authorize attribute fields
    foreach my $attr ( @CLUSTALW_PARAMS, @CLUSTALW_SWITCHES,
		       @OTHER_SWITCHES ) { $OK_FIELD{$attr}++; }


}

sub new {
    my ($class,@args) = @_;
    my $self = $class->SUPER::new(@args);

    unless (&exists_clustal()) {
	warn "Clustalw program not found as $PROGRAM or not executable. \n  Clustalw can be obtained from eg- http://corba.ebi.ac.uk/Biocatalog/Alignment_Search_software.html/ \n";
    }

    my ($attr, $value);
    (undef,$TMPDIR) = $self->tempdir(CLEANUP=>1);
    (undef,$TMPOUTFILE) = $self->tempfile(-dir => $TMPDIR);
    while (@args)  {
	$attr =   shift @args;
	$value =  shift @args;
	next if( $attr =~ /^-/ ); # don't want named parameters
	$self->$attr($value);	
    }
    return $self;
}

sub AUTOLOAD {
    my $self = shift;
    my $attr = $AUTOLOAD;
    $attr =~ s/.*:://;
    $attr = uc $attr;
    $self->throw("Unallowed parameter: $attr !") unless $OK_FIELD{$attr};
    $self->{$attr} = shift if @_;
    return $self->{$attr};
}


=head2  exists_clustal()

 Title   : exists_clustal
 Usage   : $clustalfound = Bio::Tools::Run::Alignment::Clustalw->exists_clustal()
 Function: Determine whether clustalw program can be found on current host
 Example :
 Returns : 1 if clustalw program found at expected location, 0 otherwise.
 Args    :  none

=cut


sub exists_clustal {
    my $returnvalue = (-e $PROGRAM) ;
}


=head2  align

 Title   : align
 Usage   :
	$inputfilename = 't/cysprot.fa';
	$aln = $factory->align($inputfilename);
or
	$seq_array_ref = \@seq_array; @seq_array is array of Seq objs
	$aln = $factory->align($seq_array_ref);
 Function: Perform a multiple sequence alignment
 Example :
 Returns : Reference to a SimpleAlign object containing the
           sequence alignment.
 Args    : Name of a file containing a set of unaligned fasta sequences
           or else an array of references to Bio::Seq objects.

 Throws an exception if argument is not either a string (eg a
 filename) or a reference to an array of Bio::Seq objects.  If
 argument is string, throws exception if file corresponding to string
 name can not be found. If argument is Bio::Seq array, throws
 exception if less than two sequence objects are in array.

=cut

sub align {

    my ($self,$input) = @_;
    my ($temp,$infilename, $seq);
    my ($attr, $value, $switch);

# Create input file pointer
    $infilename = $self->_setinput($input);
    if (!$infilename) {$self->throw("Bad input data or less than 2 sequences in $input !");}

# Create parameter string to pass to clustalw program
    my $param_string = $self->_setparams();

# run clustalw
    my $aln = $self->_run('align', $infilename,$param_string);
}
#################################################

=head2  profile_align

 Title   : profile_align
 Usage   :
 Function: Perform an alignment of 2 (sub)alignments
 Example :
 Returns : Reference to a SimpleAlign object containing the (super)alignment.
 Args    : Names of 2 files containing the subalignments
         or references to 2 Bio::SimpleAlign objects.

Throws an exception if arguments are not either strings (eg filenames)
or references to SimpleAlign objects.


=cut

sub profile_align {

    my ($self,$input1,$input2) = @_;
    my ($temp,$infilename1,$infilename2,$input,$seq);


# Create input file pointer
    $infilename1 = $self->_setinput($input1,1);
    $infilename2 = $self->_setinput($input2,2);
    if (!$infilename1 || !$infilename2) {$self->throw("Bad input data: $input1 or $input2 !");}


# Create parameter string to pass to clustalw program
    my $param_string = $self->_setparams();

# run clustalw
    my $aln = $self->_run('profile-aln', $infilename1,
			  $infilename2, $param_string);

}
#################################################

=head2  _run

 Title   :  _run
 Usage   :  Internal function, not to be called directly	
 Function:   makes actual system call to clustalw program
 Example :
 Returns : nothing; clustalw output is written to a
           temporary file $TMPOUTFILE
 Args    : Name of a file containing a set of unaligned fasta sequences
           and hash of parameters to be passed to clustalw


=cut

sub _run {
    my ($self,$command,$infile1,$infile2,$param_string) = @_;
    my $instring;
    if ($command =~ /align/) {
	$instring =  "-infile=$infile1";
	$param_string = $infile2;
    }
    if ($command =~ /profile/) {
	$instring =  "-profile1=$infile1  -profile2=$infile2";
	$command = '-profile';
    }
    my $commandstring = $PROGRAM." $command"." $instring".
	" -output=gcg". " $param_string";

# next line is for debugging purposes
    if( $DEBUG ) {
	print "clustal command = $commandstring \n";
    }
    my $status = system($commandstring);
    $self->throw( "Clustalw call crashed: $? \n") unless $status==0;

    my $outfile = $self->outfile() || $TMPOUTFILE ;
# retrieve alignment (Note: MSF format for AlignIO = GCG format of clustalw)
    my $in  = Bio::AlignIO->new(-file => $outfile, '-format' => 'MSF');
    my $aln = $in->next_aln();

    # Clean up the temporary files created along the way...
    # Replace file suffix with dnd to find name of dendrogram file(s) to delete
    foreach my $f ( $infile1, $infile2 ) {
	$f =~ s/\.[^\.]*$// ;
	unlink $f .'.dnd' if( $f ne '' );
    }
    return $aln;
}


=head2  _setinput()

 Title   :  _setinput
 Usage   :  Internal function, not to be called directly	
 Function:   Create input file for clustalw program
 Example :
 Returns : name of file containing clustalw data input
 Args    : Seq or Align object reference or input file name


=cut

sub _setinput {
    my ($self, $input, $suffix) = @_;
    my ($infilename, $seq, $temp, $tfh);

    # suffix is used to distinguish alignment files If $input is not a
    # reference it better be the name of a file with the sequence/

    #  alignment data...

    unless (ref $input) {
	# check that file exists or throw
	$infilename = $input;
	unless (-e $input) {return 0;}
	return $infilename;
    }

    #  $input may be an array of BioSeq objects...
    if (ref($input) eq "ARRAY") {
        #  Open temporary file for both reading & writing of BioSeq array
	($tfh,$infilename) = $self->tempfile(-dir=>$TMPDIR);
	$temp =  Bio::SeqIO->new('-fh'=>$tfh,
				 '-format' =>'Fasta');

	# Need at least 2 seqs for alignment
	unless (scalar(@$input) > 1) {return 0;}

	foreach $seq (@$input) {
	    unless (ref($seq) eq "Bio::Seq")
	    {return 0;}
	    $temp->write_seq($seq);
	}
	$temp->close();
	return $infilename;
    }
#  $input may be a SimpleAlign object.
   elsif (ref($input) eq "Bio::SimpleAlign") {
	#  Open temporary file for both reading & writing of SimpleAlign object
	if ($suffix ==1 || $suffix== 2 ) {
	    ($tfh,$infilename) = $self->tempfile(-dir=>$TMPDIR);
	}
	$temp =  Bio::AlignIO->new('-fh'=> $tfh,
				   '-format' => 'Fasta');
	$temp->write_aln($input);
	return $infilename;
    }

#  or $input may be a single BioSeq object (to be added to a previous alignment)
    elsif (ref($input) && $input->isa("Bio::PrimarySeqI") && $suffix==2) {
        #  Open temporary file for both reading & writing of BioSeq object
	($tfh,$infilename) = $self->tempfile();
	$temp =  Bio::SeqIO->new(-fh=> $tfh, '-format' =>'Fasta');
	$temp->write_seq($input);
	return $infilename;
    }
    return 0;
}

=head2  _setparams()

 Title   :  _setparams
 Usage   :  Internal function, not to be called directly	
 Function:   Create parameter inputs for clustalw program
 Example :
 Returns : parameter string to be passed to clustalw
           during align or profile_align
 Args    : name of calling object

=cut

sub _setparams {
    my ($attr, $value, $self);

    $self = shift;

    my $param_string = "";
    for  $attr ( @CLUSTALW_PARAMS ) {
	$value = $self->$attr();
	next unless (defined $value);
	my $attr_key = lc $attr; #put params in format expected by clustalw
	$attr_key = ' -'.$attr_key;
	$param_string .= $attr_key.'='.$value;
    }

    for  $attr ( @CLUSTALW_SWITCHES) {
	$value = $self->$attr();
	next unless ($value);
	my $attr_key = lc $attr; #put switches in format expected by clustalw
	$attr_key = ' -'.$attr_key;
	$param_string .= $attr_key ;
#	$attr_key = '-'.$attr_key;
#	$param_string .= '"'.$attr_key.'",';
    }

# Set default output file if no explicit output file selected
    unless ($param_string =~ /outfile/) {
	$param_string .= " -outfile=$TMPOUTFILE" ;
    }

    if ($self->quiet()) { $param_string .= '  >/dev/null';}

    return $param_string;
}

1; # Needed to keep compiler happy
