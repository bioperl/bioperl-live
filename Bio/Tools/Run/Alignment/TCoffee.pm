# $Id$
#
# BioPerl module for Bio::Tools::Run::Alignment::TCoffee
#
# Cared for by Jason Stajich, Peter Schattner
#
# Copyright Jason Stajich, Peter Schattner
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Tools::Run::Alignment::TCoffee - Object for the calculation of a
multiple sequence alignment from a set of unaligned sequences or
alignments using the TCoffee program

=head1 SYNOPSIS

# Build a tcoffee alignment factory
    @params = ('ktuple' => 2, 'matrix' => 'BLOSUM');
    $factory = new Bio::Tools::Run::Alignment::TCoffee (@params);

# Pass the factory a list of sequences to be aligned.
    $inputfilename = 't/cysprot.fa';
# $aln is a SimpleAlign object.
    $aln = $factory->align($inputfilename); 

# or where @seq_array is an array of Bio::Seq objects
    $seq_array_ref = \@seq_array;  
    $aln = $factory->align($seq_array_ref);

# Or one can pass the factory a pair of (sub)alignments to be aligned against each other, e.g.:

# where $aln1 and $aln2 are Bio::SimpleAlign objects.
    $aln = $factory->profile_align($aln1,$aln2); 

# Or one can pass the factory an alignment and one or more unaligned sequences to be added to the alignment. For example:

# $seq is a Bio::Seq object.
    $aln = $factory->profile_align($aln1,$seq); 

There are various additional options and input formats available.  See
the DESCRIPTION section that follows for additional details.

=head1 DESCRIPTION

Note: this DESCRIPTION only documents the (Bio)perl interface to
TCoffee.  

There are a number of possible parameters one can pass in TCoffee.
One should really read the online manual for the best explaination of
all the features.  See
http://igs-server.cnrs-mrs.fr/~cnotred/Documentation/t_coffee/t_coffee_doc.html

These can be specified as paramters when instantiating a new TCoffee
object, or through get/set methods of the same name (lowercase).

=head1 PARAMETERS FOR ALIGNMENT COMPUTATION

=head2 IN

 Title       : IN
 Description : (optional) input filename, this is specified when
               align so should not use this directly unless one
               understand TCoffee program very well.

=head2 TYPE

 Title       : TYPE
 Args        : [string] DNA, PROTEIN
 Description : (optional) set the sequence type, guessed automatically
               so should not use this directly

=head2 PARAMETERS

 Title       : PARAMETERS
 Description : (optional) Indicates a file containing extra parameters

=head2 EXTEND

 Title       : EXTEND
 Args        : 0, 1, or positive value
 Default     : 1
 Description : Flag indicating that library extension should be
               carried out when performing multiple alignments, if set
               to 0 then extension is not made, if set to 1 extension
               is made on all pairs in the library.  If extension is
               set to another positive value, the extension is only
               carried out on pairs having a weigth value superior to
               the specified limit.

=head2 DP_NORMALISE

 Title       : DP_NORMALISE
 Args        : 0 or positive value
 Default     : 1000
 Description : When using a value different from 0, this flag sets the
               score of the highest scoring pair to 1000.

=head2 DP_MODE

 Title       : DP_MODE
 Args        : [string] gotoh_pair_wise, myers_miller_pair_wise,
               fasta_pair_wise cfasta_pair_wise
 Default     : cfast_fair_wise 
 Description : Indicates the type of dynamic programming used by 
               the program

    gotoh_pair_wise : implementation of the gotoh algorithm 
    (quadratic in memory and time)                

    myers_miller_pair_wise : implementation of the Myers and Miller
    dynamic programming algorithm ( quadratic in time and linear in
    space). This algorithm is recommended for very long sequences. It
    is about 2 time slower than gotoh. It only accepts tg_mode=1.

    fasta_pair_wise: implementation of the fasta algorithm. The
    sequence is hashed, looking for ktuples words. Dynamic programming
    is only carried out on the ndiag best scoring diagonals. This is
    much faster but less accurate than the two previous.

    cfasta_pair_wise : c stands for checked. It is the same
    algorithm. The dynamic programming is made on the ndiag best
    diagonals, and then on the 2*ndiags, and so on until the scores
    converge. Complexity will depend on the level of divergence of the
    sequences, but will usually be L*log(L), with an accuracy
    comparable to the two first mode ( this was checked on BaliBase).

=head2 KTUPLE

 Title       : KTUPLE 
 Args        : numeric value
 Default     : 1 or 2 (1 for protein, 2 for DNA ) 

 Description : Indicates the ktuple size for cfasta_pair_wise dp_mode
               and fasta_pair_wise. It is set to 1 for proteins, and 2
               for DNA. The alphabet used for protein is not the 20
               letter code, but a mildly degenerated version, where
               some residues are grouped under one letter, based on
               physicochemical properties:
               rk, de, qh, vilm, fy (the other residues are
               not degenerated). 

=head2 NDIAGS

 Title       : NDIAGS
 Args        : numeric value
 Default     : 0
 Description : Indicates the number of diagonals used by the
               fasta_pair_wise algorithm. When set to 0, 
               n_diag=Log (length of the smallest sequence)

=head2 DIAG_MODE

 Title       : DIAG_MODE
 Args        : numeric value
 Default     : 0


 Description : Indicates the manner in which diagonals are scored
              during the fasta hashing.  

              0 indicates that the score of a diagonal is equal to the
              sum of the scores of the exact matches it contains.  
              
            
              1 indicates that this score is set equal to the score of
              the best uninterrupted segment

              1 can be useful when dealing with fragments of sequences. 

=head2 SIM_MATRIX

 Title       : SIM_MATRIX
 Args        : string 
 Default     : vasiliky
 Description : Indicates the manner in which the amino acid is being
               degenerated when hashing. All the substitution matrix
               are acceptable. Categories will be defined as sub-group
               of residues all having a positive substitution score
               (they can overlap).

               If you wish to keep the non degenerated amino acid
               alphabet, use 'idmat'

=head2 MATRIX

 Title       : MATRIX
 Args        : 
 Default     :
 Description : This flag is provided for compatibility with
               ClustalW. Setting matrix = 'blosum' is equivalent to
               -in=Xblosum62mt , -matrix=pam is equivalent to
               in=Xpam250mt . Apart from this, the rules are similar
               to those applying when declaring a matrix with the
               -in=X fl

=head2 GAPOPEN

 Title       : GAPOPEN
 Args        : numeric
 Default     : 0
 Description : Indicates the penalty applied for opening a gap. The
               penalty must be negative. If you provide a positive
               value, it will automatically be turned into a negative
               number. We recommend a value of 10 with pam matrices,
               and a value of 0 when a library is used.

=head2 GAPEXT

 Title       : GAPEXT
 Args        : numeric
 Default     : 0
 Description : Indicates the penalty applied for extending a gap. 


=head2 COSMETIC_PENALTY

 Title       : COSMETIC_PENALTY
 Args        : numeric
 Default     : 100
 Description : Indicates the penalty applied for opening a gap. This
               penalty is set to a very low value. It will only have
               an influence on the portions of the alignment that are
               unalignable. It will not make them more correct, but
               only more pleasing to the eye ( i.e. Avoid stretches of
               lonely residues).

               The cosmetic penalty is automatically turned off if a
               substitution matrix is used rather than a library.

=head2 TG_MODE

 Title       : TG_MODE
 Args        : 0,1,2
 Default     : 1
 Description : (Terminal Gaps) 
               0: indicates that terminal gaps must be panelized with
                  a gapopen and a gapext penalty.
               1: indicates that terminal gaps must be penalized only
                  with a gapext penalty
               2: indicates that terminal gaps must not be penalized.

=head2 WEIGHT

 Title       : WEIGHT
 Args        : sim or sim_<matrix_name or matrix_file> or integer value 
 Default     : sim
 

 Description : Weight defines the way alignments are weighted when
               turned into a library.

               sim indicates that the weight equals the average
                   identity within the match residues.

               sim_matrix_name indicates the average identity with two
                   residues regarded as identical when their
                   substitution value is positive. The valid matrices
                   names are in matrices.h (pam250mt) . Matrices not
                   found in this header are considered to be
                   filenames. See the format section for matrices. For
                   instance, -weight=sim_pam250mt indicates that the
                   grouping used for similarity will be the set of
                   classes with positive substitutions. Other groups
                   include

                       sim_clustalw_col ( categories of clustalw
                       marked with :)

                       sim_clustalw_dot ( categories of clustalw
                       marked with .)


               Value indicates that all the pairs found in the
               alignments must be given the same weight equal to
               value. This is useful when the alignment one wishes to
               turn into a library must be given a pre-specified score
               (for instance if they come from a structure
               super-imposition program). Value is an integer:

                       -weight=1000 

  Note       : Weight only affects methods that return an alignment to
               T-Coffee, such as ClustalW. On the contrary, the
               version of Lalign we use here returns a library where
               weights have already been applied and are therefore
               insensitive to the -weight flag.

=head2 SEQ_TO_ALIGN

 Title       : SEQ_TO_ALIGN
 Args        : filename
 Default     : no file - align all the sequences

 Description : You may not wish to align all the sequences brought in
               by the -in flag. Supplying the seq_to_align flag allows
               for this, the file is simply a list of names in Fasta
               format.  

               However, note that library extension will be carried out 
               on all the sequences.

=head2 PARAMETERS FOR TREE COMPUTATION AND OUTPUT

=head2 NEWTREE

 Title       : NEWTREE
 Args        : treefile
 Default     : no file
 Description : Indicates the name of the new tree to compute. The
               default will be <sequence_name>.dnd, or <run_name.dnd>.
               Format is Phylips tree format

=head2 USETREE

 Title       : USETREE
 Args        : treefile 
 Default     : no file specified
 Description : This flag indicates that rather than computing a new
               dendrogram, t_coffee can use a pre-computed one. The
               tree files are in phylips format and compatible with
               ClustalW. In most cases, using a pre-computed tree will
               halve the computation time required by t_coffee. It is
               also possible to use trees output by ClustalW or
               Phylips. Format is Phylips tree format

=head2 TREE_MODE

 Title       : TREE_MODE
 Args        : slow, fast, very_fast
 Default     : very_fast
 Description : This flag indicates the method used for computing the 
               dendrogram. 
               slow : the chosen dp_mode using the extended library, 
               fast : The fasta dp_mode using the extended library.
               very_fast: The fasta dp_mode using pam250mt. 

=head2 QUICKTREE

 Title       : QUICKTREE
 Args        : 
 Default     : 
 Description : This flag is kept for compatibility with ClustalW. 
               It indicates that:  -tree_mode=very_fast

=head2 PARAMETERS FOR ALIGNMENT OUTPUT

=head2 OUTFILE

 Title       : OUTFILE
 Args        : out_aln file, default, no
 Default     : default ( yourseqfile.aln)
 Description : indicates name of output alignment file

=head2 OUTPUT

 Title       : OUTPUT
 Args        : format1, format2
 Default     : clustalw
 Description : Indicated format for outputting outputfile
               Supported formats are: 
               
               clustalw_aln, clustalw: ClustalW format. 
               gcg, msf_aln : Msf alignment. 
               pir_aln : pir alignment. 
               fasta_aln : fasta alignment. 
               phylip : Phylip format. 
               pir_seq : pir sequences (no gap). 
               fasta_seq : fasta sequences (no gap). 
    As well as: 
                score_html : causes the output to be a reliability 
                             plot in HTML 
                score_pdf : idem in PDF. 
                score_ps : idem in postscript. 

    More than one format can be indicated: 
                -output=clustalw,gcg, score_html 

=head2 CASE

 Title       : CASE
 Args        : upper, lower
 Default     : upper
 Description : triggers choice of the case for output

=head2 CPU

 Title       : CPU
 Args        : value
 Default     : 0
 Description : Indicates the cpu time (micro seconds) that must be
               added to the t_coffee computation time.

=head2 OUT_LIB

 Title       : OUT_LIB
 Args        : name of library, default, no 
 Default     : default
 Description : Sets the name of the library output. Default implies
               <run_name>.tc_lib

=head2 OUTORDER

 Title       : OUTORDER
 Args        : input or aligned
 Default     : input
 Description : Sets the name of the library output. Default implies
               <run_name>.tc_lib

=head2 SEQNOS

 Title       : SEQNOS
 Args        : on or off
 Default     : off
 Description : Causes the output alignment to contain residue numbers
               at the end of each line:

=head2 PARAMETERS FOR GENERIC OUTPUT

=head2 RUN_NAME

 Title       : RUN_NAME
 Args        : your run name
 Default     : 
 Description : This flag causes the prefix <your sequences> to be
               replaced by <your run name> when renaming the default
               files.

=head2 ALIGN

 Title       : ALIGN
 Args        : 
 Default     :
 Description : Indicates that the program must produce the
               alignment. This flag is here for compatibility with
               ClustalW

=head2 QUIET

 Title       : QUIET
 Args        : stderr, stdout, or filename, or nothing 
 Default     : stderr
 Description : Redirects the standard output to either a file. 
              -quiet on its own redirect the output to /dev/null. 

=head2 CONVERT

 Title       : CONVERT
 Args        : 
 Default     :
 Description : Indicates that the program must not compute the
               alignment but simply convert all the sequences,
               alignments and libraries into the format indicated with
               -output. This flag can also be used if you simply want
               to compute a library ( i.e. You have an alignment and
               you want to turn it into a library).

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

=head1 AUTHOR -  Jason Stajich, Peter Schattner

Email jason@chg.mc.duke.edu, schattner@alum.mit.edu

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

package Bio::Tools::Run::Alignment::TCoffee;

use vars qw($AUTOLOAD @ISA $TMPOUTFILE $PROGRAM $PROGRAMDIR
            @TCOFFEE_PARAMS @TCOFFEE_SWITCHES @OTHER_SWITCHES %OK_FIELD
            );
use strict;
use Bio::Seq;
use Bio::SeqIO;
use Bio::SimpleAlign;
use Bio::AlignIO;
use Bio::Root::RootI;
use Bio::Root::IO;

@ISA = qw(Bio::Root::RootI Bio::Root::IO);


BEGIN { 
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
    $PROGRAM = Bio::Root::IO->catfile($PROGRAMDIR,'t_coffee');
    
    @TCOFFEE_PARAMS = qw(IN TYPE PARAMETERS DO_NORMALISE EXTEND
			 DP_MODE KTUPLE NDIAGS DIAG_MODE SIM_MATRIX 
			 MATRIX GAPOPEN GAPEXT COSMETIC_PENALTY TG_MODE
			 WEIGHT SEQ_TO_ALIGN NEWTREE USETREE TREE_MODE 
			 OUTFILE OUTPUT CASE CPU OUT_LIB OUTORDER SEQNOS
			 RUN_NAME CONVERT);
    
    @TCOFFEE_SWITCHES = qw(QUICKTREE);
    
    @OTHER_SWITCHES = qw(QUIET ALIGN KEEPDND);
    
# Authorize attribute fields
    foreach my $attr ( @TCOFFEE_PARAMS, @TCOFFEE_SWITCHES, @OTHER_SWITCHES ) { 
	$OK_FIELD{$attr}++; }
}

sub new {
    my ($class,@args) = @_;
    my $self = $class->SUPER::new(@args);
    my ($attr, $value);
    (undef,$TMPOUTFILE) = $self->tempfile();
    while (@args)  {
	$attr =   shift @args;
	$value =  shift @args;
	next if( $attr =~ /^-/); # don't want named parameters 
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


=head2  exists_tcoffee

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
    my ($self,$input) = @_;
# Create input file pointer
    my $infilename = $self->_setinput($input);
    if (!$infilename) {$self->throw("Bad input data or less than 2 sequences in $input !");}

# Create parameter string to pass to tcoffee program
    $self->{'_in'} = [];    

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
    $self->{'_in'} = [];    
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
	push @{$self->{'_in'}}, "$infilename";
    }
    if ($command =~ /profile/) {
	$infile1 = shift ;
        $infile2 = shift ;
	push @{$self->{'_in'}}, "$infile1", "$infile2";	
    }
    my $param_string = shift;
    my $instring = "-in=".join(",", @{$self->{'_in'}});
    my $commandstring = $PROGRAM." $instring".
	" -output=gcg". " $param_string";    
    # next line is for debugging purposes
    if( $self->verbose > 0 ) {
	print "tcoffee command = $commandstring \n";
    }
    
    my $status = system($commandstring);
    $self->throw( "TCoffee call crashed: $? \n") if( -z $TMPOUTFILE );

    my $outfile = $self->outfile() || $TMPOUTFILE;

    # retrieve alignment (Note: MSF format for AlignIO = GCG format of
    # tcoffee)

    my $in  = Bio::AlignIO->new(-file => $outfile, '-format' => 'MSF');
    my $aln = $in->next_aln();
   
    # Replace file suffix with dnd to find name of dendrogram file(s) to delete
    if( ! $self->keepdnd ) {
	foreach my $f ( $infilename, $infile1, $infile2 ) {
	    next if( !defined $f || $f eq '');
	    $f =~ s/\.[^\.]*$// ;   
	    # because TCoffee writes these files to the CWD
	    if( $Bio::Root::IO::PATHSEP ) {
		my @line = split(/$Bio::Root::IO::PATHSEP/, $f);
		$f = pop @line;	    
	    } else { 		
		(undef, undef, $f) = File::Spec->splitpath($f);    
	    }
	    unlink $f .'.dnd' if( $f ne '' );
	}
    }
    return $aln;
}


=head2  _setinput

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
#  If $input is not a reference it better be the name of a 
# file with the sequence/ alignment data...
    if (! ref $input) {
	# check that file exists or throw
	$infilename = $input;
	unless (-e $input) {return 0;}
	return $infilename;
    }
#  $input may be an array of BioSeq objects...
    elsif (ref($input) =~ /ARRAY/i ) {
        #  Open temporary file for both reading & writing of BioSeq array
	($tfh,$infilename) = $self->tempfile();
	$temp =  Bio::SeqIO->new('-fh' => $tfh, 
				'-format' => 'Fasta');
	unless (scalar(@$input) > 1) {return 0;} # Need at least 2 seqs for alignment
	foreach $seq (@$input) {
	    return 0 if( !ref($seq) || ! $seq->isa("Bio::PrimarySeqI"));
	    $temp->write_seq($seq);
	}
	return $infilename;
    }
#  $input may be a SimpleAlign object.
    elsif ( $input->isa("Bio::SimpleAlign") ) {
	#  Open temporary file for both reading & writing of SimpleAlign object
	($tfh, $infilename) = $self->tempfile() if ($suffix ==1 
						    || $suffix == 2 );
	$temp =  Bio::AlignIO->new(-fh=>$tfh, 
				   '-format' => 'Fasta');
	$temp->write_aln($input);
	return $infilename;
    }

#  or $input may be a single BioSeq object (to be added to 
# a previous alignment)
    elsif ( $input->isa("Bio::PrimarySeqI") && $suffix==2) {
        #  Open temporary file for both reading & writing of BioSeq object
	($tfh,$infilename) = $self->tempfile();
	$temp =  Bio::SeqIO->new(-fh=> $tfh, '-format' =>'Fasta');
	$temp->write_seq($input);
	return $infilename;
    }
    return 0;
}


=head2  _setparams

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
    for  $attr ( @TCOFFEE_PARAMS ) {
	$value = $self->$attr();
	next unless (defined $value);	
	my $attr_key = lc $attr;
	if( $attr_key =~ /matrix/ ) {
	    $self->{'_in'} = [ "X".lc($value) ];
	} else {
	    $attr_key = ' -'.$attr_key;
	    $param_string .= $attr_key .'='.$value; 
	}
   }
    for  $attr ( @TCOFFEE_SWITCHES) {
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
