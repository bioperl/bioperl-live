# $Id$
#
# BioPerl module for Bio::Tools::Run::Alignment::TCoffee
#
# Cared for by Jason Stajich
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Tools::Run::Alignment::TCoffee - Object for the calculation of a
multiple sequence alignment from a set of unaligned sequences or
alignments using the TCoffee program

=head1 SYNOPSIS

#  Build a tcoffee alignment factory
	@params = ('ktuple' => 2, 'matrix' => 'BLOSUM');
	$factory = Bio::Tools::Run::Alignment::TCoffee->new(@params);

#  Pass the factory a list of sequences to be aligned.	
	$inputfilename = 't/cysprot.fa';
	$aln = $factory->align($inputfilename); # $aln is a SimpleAlign object.
or
	$seq_array_ref = \@seq_array;  # where @seq_array is an array of Bio::Seq objects
	$aln = $factory->align($seq_array_ref);

#  Or one can pass the factory a pair of (sub)alignments to be aligned against each other, e.g.:	
	$aln = $factory->profile_align($aln1,$aln2); # where $aln1 and $aln2 are Bio::SimpleAlign objects.

#  Or one can pass the factory an alignment and one or more unaligned
#  sequences to be added to the alignment. For example: 	
	$aln = $factory->profile_align($aln1,$seq); # $seq is a Bio::Seq object.

There are various additional options and input formats available.  See
the DESCRIPTION section that follows for additional details.

=head1 DESCRIPTION

Note: this DESCRIPTION only documents the (Bio)perl interface to
TCoffee.  

=head1 DEVELOPERS NOTES


=head1 STILL TO WRITE


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

=head1 AUTHOR -  Jason Stajich

Email schattner@alum.mit.edu

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

package Bio::Tools::Run::Alignment::TCoffee;

use vars qw($AUTOLOAD @ISA $TMPOUTFILE $DEBUG $PROGRAM $PROGRAMDIR);
use strict;
use Bio::Seq;
use Bio::SeqIO;
use Bio::SimpleAlign;
use Bio::AlignIO;
use Bio::Root::RootI;

@ISA = qw(Bio::Root::RootI);

# You will need to enable TCoffee to find the tcoffee program. This can be done
# in (at least) three ways:
#  1. Modify your $PATH variable to include your tcoffee directory as in (for Linux):
#	export PATH=$PATH:/home/progs/tcoffee  or
#  2. define an environmental variable TCOFFEE:
#	export TCOFEEDIR=/home/progs/tcoffee   or
#  3. include a definition of an environmental variable TCOFFEEDIR in every script that will
#     use Bio::Tools::Run::Alignment::TCoffee.pm.
#	BEGIN {$ENV{TCOFFEEDIR} = '/home/progs/tcoffee'; }
$PROGRAMDIR = $ENV{TCOFFEEDIR} || '';
$PROGRAMDIR .= '/' if( substr($PROGRAMDIR, -1) ne '/' ); 
$PROGRAM =   $PROGRAMDIR.'t_coffee' ;
$DEBUG = 0;
unless (exists_tcoffee()) {
	warn "TCoffee program not found as $PROGRAM or not executable. \n  TCoffee can be obtained from eg- http://igs-server.cnrs-mrs.fr/~cnotred/Projects_home_page/t_coffee_home_page.html \n";
}
# Object preamble - inherits from Bio::Root::RootI


#***Pairwise alignments:***
#	KTUPLE      	#:(=n) word size
#	TOPDIAGS  	#:(=n) number of best diagonals
#	WINDOW   	#:(=n) window around best diagonals
#	PAIRGAP   	#:(=n)gap penalty

#**Multiple alignments:***
#	FIXEDGAP  	#:(=n)fixed length gap pen.
#	FLOATGAP  	#:(=n)variable length gap pen.
#	MATRIX     	#:= PAM100 or ID or file name. The default weight matrix
			#	for proteins is PAM 250.	
#	TYPE	 	#:(=p or d)type is protein or DNA.   This allows you to 	
			#	explicitly overide the programs attempt at guessing
			#	the type of the sequence.  It is only useful if you
			#	are using sequences with a VERY strange composition.
#	OUTPUT     	#:= tcoffee supports GCG or PHYLIP or PIR or Clustal format.
			# However currently only GCG format is supported by TCoffee.pm
			# (because AlignIO input modules haven't been written yet for the
                        # other formats)
#	OUTFILE     	#: Name of tcoffee's output file. If not set
			# module will erase output file.  In any case alignment will
			# be returned in the form of SimpleAlign objects
#	TRANSIT     	#:transitions not weighted.  The default is to weight
			#	transitions as more favourable than other mismatches
			#	in DNA alignments.  This switch makes all nucleotide
			#	mismatches equally weighted.

my @tcoffee_params = qw(IN TYPE PARAMETERS DO_NORMALIZE EXTEND
			DP_MODE KTUPLE NDIAGS DIAG_MODE SIM_MATRIX 
			MATRIX GAPOPEN GAPEXT COSMETIC_PENALTY TG_MODE
			WEIGHT SEQ_TO_ALIGN NEWTREE USETREE TREEMODE 
			OUTFILE OUTPUT CASE CPU OUT_LIB OUTORDER SEQNOS
			RUN_NAME CONVERT);

my @tcoffee_switches = qw(QUICKTREE);

my @other_switches = qw(QUIET ALIGN);
my %ok_field;


# Authorize attribute fields
foreach my $attr ( @tcoffee_params, @tcoffee_switches, @other_switches ) { 
    $ok_field{$attr}++; }

# new comes from RootI

sub _initialize {
    my($self,@args) = @_;
    my ($attr, $value);
    my $make = $self->SUPER::_initialize(@args);
    (undef,$TMPOUTFILE) = $self->tempfile();
    while (@args)  {
	$attr =  shift @args;
	$value =  shift @args;
	$self->$attr($value);
    }
    return $make;		# success - we hope!
}


sub AUTOLOAD {
    my $self = shift;
    my $attr = $AUTOLOAD;
    $attr =~ s/.*:://;
    $attr = uc $attr;
    $self->throw("Unallowed parameter: $attr !") unless $ok_field{uc $attr};

    $self->{uc $attr} = shift if @_;
    return $self->{uc $attr};
}


=head2  exists_tcoffee()

 Title   : exists_tcoffee
 Usage   : $coffeefound = Bio::Tools::Run::Alignment::TCoffee->exists_tcoffee()
 Function: Determine whether tcoffee program can be found on current host
 Example :
 Returns : 1 if tcoffee program found at expected location, 0 otherwise.
 Args    :  none

=cut


sub exists_tcoffee {
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

    my $self = shift;
    my $input = shift;
    my ($temp,$infilename, $seq);
    my ($attr, $value, $switch);

# Create input file pointer
    $infilename = $self->_setinput($input);
    if (!$infilename) {$self->throw("Bad input data or less than 2 sequences in $input !");}

# Create parameter string to pass to tcoffee program
    $self->{_in} = [];    

    my $param_string = $self->_setparams();

# run tcoffee
    my $aln = &_run($self, 'align', $infilename, $param_string);
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

    my $self = shift;
    my $input1 = shift;
    my $input2 = shift;
    my ($temp,$infilename1,$infilename2,$input,$seq);

# Create input file pointers
    $infilename1 = $self->_setinput($input1,1);
    $infilename2 = $self->_setinput($input2,2);
    if (!$infilename1 || !$infilename2) {$self->throw("Bad input data: $input1 or $input2 !");}

    # Create parameter string to pass to tcoffee program
    $self->{_in} = [];    
    my $param_string = $self->_setparams();

# run tcoffee
    my $aln = $self->_run('profile-aln', $infilename1, 
			  $infilename2, $param_string);

}
#################################################

=head2  _run

 Title   :  _run
 Usage   :  Internal function, not to be called directly	
 Function:  makes actual system call to tcoffee program
 Example :
 Returns : nothing; tcoffee output is written to a 
           temporary file $TMPOUTFILE
 Args    : Name of a file containing a set of unaligned fasta sequences
           and hash of parameters to be passed to tcoffee


=cut
sub _run {
    my ($infilename, $infile1,$infile2) = ('','','');
    my $self = shift;
    my $command = shift;
    if ($command =~ /align/) {
        $infilename = shift ;
	push @{$self->{_in}}, "$infilename";
    }
    if ($command =~ /profile/) {
	$infile1 = shift ;
        $infile2 = shift ;
	push @{$self->{_in}}, "$infile1", "$infile2";	
    }
    my $param_string = shift;
    my $instring = "-in=".join(",", @{$self->{_in}});
    my $commandstring = $PROGRAM." $instring".
	" -output=gcg". " $param_string";    
    # next line is for debugging purposes
    if( $DEBUG ) {
	print "tcoffee command = $commandstring \n";
    }
    
    my $status = system($commandstring);
    $self->throw( "TCoffee call crashed: $? \n") if(-z $TMPOUTFILE);

    my $outfile = $self->outfile() || $TMPOUTFILE;

    # retrieve alignment (Note: MSF format for AlignIO = GCG format of
    # tcoffee)

    my $in  = Bio::AlignIO->new(-file => $outfile, '-format' => 'MSF');
    my $aln = $in->next_aln();
   
    # Replace file suffix with dnd to find name of dendrogram file(s) to delete
    $infilename =~ s/\.[^\.]*$// ;
    $infile1 =~ s/\.[^\.]*$// ;
    $infile2 =~ s/\.[^\.]*$// ;
    unlink ( "$infilename.dnd", "$infile1.dnd", "$infile2.dnd" );
    return $aln;
}


=head2  _setinput()

 Title   :  _setinput
 Usage   :  Internal function, not to be called directly	
 Function:  Create input file for tcoffee program
 Example :
 Returns : name of file containing tcoffee data input
 Args    : Seq or Align object reference or input file name


=cut

sub _setinput {
    my ($self,$input, $suffix) = @_;
    my ($infilename, $seq, $temp, $tfh);    
# suffix used to distinguish alignment files
#  If $input is not a reference it better be the name of a file with the sequence/
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
	($tfh,$infilename) = $self->tempfile();
	$temp =  Bio::SeqIO->new(-fh => $tfh, 
				'-format' => 'Fasta');
	unless (scalar(@$input) > 1) {return 0;} # Need at least 2 seqs for alignment
	foreach $seq (@$input) {
	    unless (ref($seq) eq "Bio::Seq")
	    {return 0;}
	    $temp->write_seq($seq);
	}
	return $infilename;
    }
#  $input may be a SimpleAlign object.
    if (ref($input) eq "Bio::SimpleAlign") {
	#  Open temporary file for both reading & writing of SimpleAlign object
	($tfh, $infilename) = $self->tempfile() if ($suffix ==1 || $suffix
== 2 );
	$temp =  Bio::AlignIO->new(-fh=>$tfh, 
				   '-format' => 'Fasta');
	$temp->write_aln($input);
	return $infilename;
    }

#  or $input may be a single BioSeq object (to be added to a previous alignment)
    if (ref($input) eq "Bio::Seq" && $suffix==2) {
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
 Function:  Create parameter inputs for tcoffee program
 Example :
 Returns : parameter string to be passed to tcoffee 
           during align or profile_align
 Args    : name of calling object

=cut

sub _setparams {
    my ($self) = @_;
    my ($attr, $value,$param_string);

    my $laststr;
    for  $attr ( @tcoffee_params ) {
	$value = $self->$attr();
	next unless (defined $value);	
	my $attr_key = lc $attr;
	if( $attr_key =~ /matrix/ ) {
	    $self->{_in} = [ "X".lc($value) ];
	} else {
	    $attr_key = ' -'.$attr_key;
	    $param_string .= $attr_key .'='.$value; 
	}
   }
    for  $attr ( @tcoffee_switches) {
	$value = $self->$attr();
	next unless ($value);
	my $attr_key = lc $attr; #put switches in format expected by tcoffee
	$attr_key = ' -'.$attr_key;
	$param_string .= $attr_key ;
    }

# Set default output file if no explicit output file selected
    unless ($param_string =~ /outfile/) {
	$param_string .= " -outfile=$TMPOUTFILE" ;
    }

    if ($self->quiet()) { $param_string .= ' -quiet';}
    return $param_string;
}

1; # Needed to keep compiler happy
