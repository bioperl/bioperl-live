# BioPerl module for Bio::Tools::Lagan
#
# Cared for by Stephen Montgomery <smontgom@bcgsc.bc.ca>
#
# Copyright Stephen Montgomery
#
# Special thanks to Peter Schattner.
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Tools::Lagan - Object for the local execution of the LAGAN suite of tools (including MLAGAN for multiple sequence alignments)

=head1 SYNOPSIS

To run mlagan/lagan, the executables "mlagan" and "lagan.pl" must be in your path or you must have an environment variable that points to the executable directory "LAGANDIR=/opt/lagan_executables/"

MLAGAN / LAGAN execution and alignment object creation.
	
	use Bio::Tools::Lagan;
	
	@params = (	'chaos' => "The contents of this string will be passed as args to chaos",
			#Read you chaos README file for more info/This functionality has not been tested and will be
			#integrated in future versions.			

			'order' => "-gs -7 -gc -2 -mt 2 -ms -1",
			#Where gap start penalty of- 7, gap continue of -2, match of 2, and mismatch of -1.
			
			'recurf1' => "(12,25),(7,25),(4,30)",
			#A list of (wordlength,score cutoff) pairs to be used in the recursive anchoring
			
			'tree' => "(sample1 (sample2 sample3))",
			#Used by mlagan / tree can also be passed when calling mlagan directly
			
			#SCORING PARAMETERS FOR MLAGAN
			'match' => 12,
			'mismatch' => -8,
			'gapstart' => -50,
			'gapend' => -50,
			'gapcont' => -2,
	);

	All lagan and mlagan parameters listed in their Readmes can be set except for the mfa flag which has been turned on by default to prevent parsing of the alignment format.

TO USE LAGAN:

	my $lagan = new Bio::Tools::Lagan(@params);
	my $report_out = $lagan->lagan($seq1, $seq2);

	A SimpleAlign object is returned	

TO USE MLAGAN:

	my $lagan = new Bio::Tools::Lagan();
	my $tree = "(($seqname1 $seqname2) $seqname3)";
	my @sequence_objs; 	#an array of bioperl Seq objects
	
	##If you use an unblessed seq array
	my $seq_ref = \@sequence_objs;
	bless $seq_ref, "ARRAY";

	my $report_out = $lagan->mlagan($seq_ref, $tree);

	A SimpleAlign object is returned	

Only basic mlagan/lagan functionality has been implemented due to the iterative development of their project.  Future maintenance upgrades will include enhanced features and scoring.

=head1 DESCRIPTION

A parser for Lagan output

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

=head1 AUTHOR - Stephen Montgomery

Email smontgom@bcgsc.bc.ca

Genome Sciences Centre in beautiful Vancouver, British Columbia CANADA

=head1 CONTRIBUTORS

MLagan/Lagan is the hard work of Michael Brudno et al.

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

package Bio::Tools::Lagan;

use vars qw(@ISA $PROGRAM_DIR @LAGAN_PARAMS @MLAGAN_PARAMS %OK_FIELD $AUTOLOAD);

use strict;
use Bio::Root::Root;
use Bio::Root::IO;
use Bio::Seq;
use Bio::SeqIO;
use Bio::AlignIO;
use Bio::AlignIO::fasta;
use Bio::SimpleAlign;
use Bio::Tools::Run::WrapperBase;

@ISA = qw(	Bio::Root::Root
		Bio::Tools::Run::WrapperBase);

BEGIN {
	@LAGAN_PARAMS = qw(chaos order recurse mfa out lazy maskedonly usebounds rc
			translate draft info fastreject);
	@MLAGAN_PARAMS = qw(nested postir lazy verbose tree match mismatch gapstart gapend gapcont
			out version);	
	#Not all of these parameters are useful in this context, care should be used in setting only standard ones

	#Authorize Attribute fields
	foreach my $attr (@LAGAN_PARAMS, @MLAGAN_PARAMS)
     		{ $OK_FIELD{$attr}++; }

	#The LAGANDIR environment variable should be set if the lagan executables aren't in your path.
	$PROGRAM_DIR = $ENV{'LAGANDIR'} || '';
}

sub new {
  	my($class, @args) = @_;
	
  	my $self = $class->SUPER::new(@args);
	my (undef, $tempfile) = $self->io->tempfile();
	$self->out($tempfile);
	while (@args) {
		my $attr = shift @args;
		my $value = shift @args;
		$self->$attr($value);
	}
  	return $self;
}

sub AUTOLOAD {
    my $self = shift;
    my $attr = $AUTOLOAD;
    $attr =~ s/.*:://;

    $self->throw("Unallowed parameter: $attr !") unless $OK_FIELD{$attr};
    $self->{$attr} = shift if @_;
    return $self->{$attr};
}

=head2 lagan

	Runs the Lagan pairwise alignment algorithm
	Inputs should be two PrimarySeq objects.
	Returns an SimpleAlign object / preloaded with the tmp file of the Lagan multifasta output.

=cut

sub lagan {
	my ($self, $input1, $input2) = @_;
	$self->io->_io_cleanup();
	my $executable = 'lagan.pl';
		
	#my (undef, $tempfile) = $self->io->tempfile();
        #$self->out($tempfile);

	my ($infile1, $infile2) = $self->_setinput($executable, $input1, $input2);
	my $lagan_report = &_generic_lagan(	$self,
						$executable,
						$infile1,
						$infile2 );
}

=head2 mlagan

        Runs the Mlagan multiple sequence alignment algorithm
        Inputs should be an Array of Primary Seq objects and a Phylogenetic Tree in String format
        Returns an SimpleAlign object / preloaded with the tmp file of the Mlagan multifasta output.

=cut

sub mlagan {
	my ($self, $input1, $tree) = @_;
	$self->io->_io_cleanup();
	my $executable = 'mlagan';
	my ($infiles, $tree) = $self->_setinput($executable, $input1, $tree);
	my $lagan_report = &_generic_lagan (	$self,
						$executable,
						$infiles,
						$tree );
}

=head2  _setinput

 Title   :  _setinput
 Usage   :  Internal function, not to be called directly
 Function:  Create input file(s) for Lagan executables
 Returns : name of files containing Lagan data input / or array of files and phylo tree for Mlagan data input

=cut


sub _setinput {
	my ($self, $executable, $input1, $input2) = @_;
	my ($fh, $infile1, $infile2, $temp1, $temp2, $seq1, $seq2);

	$self->io->_io_cleanup();
	
	SWITCH: {
		if ($input1->isa("Bio::PrimarySeqI")) {
			##INPUTS TO LAGAN
			($fh, $infile1) = $self->io->tempfile();

			#Want to make sure their are no white spaces in sequence.  Happens if input1 is taken
			#from an alignment.

			my $sequence = $input1->seq();
			$sequence =~ s/\W+//g;
			$input1->seq($sequence);
			$temp1 = Bio::SeqIO->new(	-fh => $fh,
							-format => 'Fasta' );
			$temp1->write_seq($input1);
			close $fh;
			undef $fh;
			last SWITCH;		
		}
		if (ref($input1) =~ /ARRAY/i) {
			##INPUTS TO MLAGAN / WILL hAVE TO BE CHANGED IF LAGAN EVER SUPPORTS MULTI-INPUT
			my @infilearr;
			foreach $seq1 (@$input1) {
				($fh, $infile1) = $self->io->tempfile();
				my $temp = Bio::SeqIO->new(	-fh => $fh,
								-format => 'Fasta' );
				unless ($seq1->isa("Bio::PrimarySeqI")) { return 0; }
				$temp->write_seq($seq1);
				close $fh;
			        undef $fh;
				push @infilearr, $infile1;
			}
			$infile1 = \@infilearr;
			last SWITCH;  
		}
	}
	SWITCH2: {
		if (ref($input2))
		{
			if ($input2->isa("Bio::PrimarySeqI")) {
                        	($fh, $infile2) = $self->io->tempfile();

                        	#Want to make sure their are no white spaces in sequence.  Happens if input2 is taken
                        	#from an alignment.

                        	my $sequence = $input2->seq();
                        	$sequence =~ s/\W+//g;
                        	$input2->seq($sequence);

                        	$temp2 = Bio::SeqIO->new(       -fh => $fh,
                                	                        -format => 'Fasta' );
                        	$temp2->write_seq($input2);
                        	close $fh;
                        	undef $fh;
                        	last SWITCH2;
                	}
		}
		else
		{
			$infile2 = $input2;
			##A tree as a scalar has been passed, pass it through
		}
        }
	return ($infile1, $infile2);
}

=head2  _generic_lagan

 Title   : _generic_lagan
 Usage   :  internal function not called directly
 Returns :  SimpleAlign object

=cut


sub _generic_lagan {
	my ($self, $executable, $input1, $input2) = @_;
	my $param_string = $self->_setparams($executable);
	my $lagan_report = &_runlagan($self, $executable, $param_string, $input1, $input2);	
}	

=head2  _setparams

 Title   : _setparams
 Usage   : Internal function, not to be called directly
 Function: Create parameter inputs for (m)Lagan program
 Returns : parameter string to be passed to Lagan
 Args    : Reference to calling object and name of (m)Lagan executable

=cut


sub _setparams {
	my ($self, $executable) = @_;
	my ($attr, $value, @execparams);

	if ($executable eq 'lagan.pl') { @execparams = @LAGAN_PARAMS; }
	if ($executable eq 'mlagan') { @execparams = @MLAGAN_PARAMS; }
	##EXPAND OTHER LAGAN SUITE PROGRAMS HERE

	my $param_string = "";
	for $attr (@execparams) {
		$value = $self->$attr();
		next unless (defined $value);
		$attr = '-' . $attr;
		$param_string .= " $attr $value ";
	}
	return $param_string . " -mfa ";
}	


=head2  _runlagan

 Title   :  _runlagan
 Usage   :  Internal function, not to be called directly
 Function:   makes actual system call to (m)Lagan program
 Example :
 Returns : Report object in the SimpleAlign object

=cut

sub _runlagan {
	my ($self, $executable, $param_string, $input1, $input2) = @_;
	my ($lagan_obj, $exe);
	if ( ! ($exe = $self->executable($executable)))  {
		$self->warn("cannot find path to $executable");
		return undef;
	}

	my $command_string;
	if ($executable eq 'lagan.pl')
	{
		$command_string = $exe . " " . $input1 . " " . $input2 . $param_string;
	}
	if ($executable eq 'mlagan')
	{
		$command_string = $exe;
		foreach my $tempfile (@$input1)
		{
			$command_string .= " " . $tempfile;
		}
		if (defined $input2)
		{
			$command_string .= " -tree " . "\"" . $input2 . "\"";
		}	
		$command_string .= " " . $param_string;
		print $command_string;
	}

	$self->debug("$command_string\n");
	my $status = system($command_string);
	my $outfile = $self->out();
	
	my $align = Bio::AlignIO->new(	'-file' => $outfile,
					'-format' => 'fasta' );
	my $aln = $align->next_aln();

	return $aln;
}   

=head2 executable

 Title   : executable
 Usage   : my $exe = $lagan->executable('mlagan');
 Function: Finds the full path to the 'lagan' executable
 Returns : string representing the full path to the exe
 Args    : [optional] name of executable to set path to
           [optional] boolean flag whether or not warn when exe is not found

 Thanks to Peter Schattner for providing the framework for this subroutine

=cut


sub executable {
   	my ($self, $exename, $exe, $warn) = @_;
   	$exename = 'lagan.pl' unless defined $exename;

   	if( defined $exe && -x $exe ) {
    	 	$self->{'_pathtoexe'}->{$exename} = $exe;
   	}
   	unless( defined $self->{'_pathtoexe'}->{$exename} ) {
       		my $f = $self->program_path($exename);
       		$exe = $self->{'_pathtoexe'}->{$exename} = $f if(-e $f && -x $f );

       		unless( $exe )  { 
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
 Usage   : my $path = $lagan->program_path();
 Function: Builds path for executable
 Returns : string representing the full path to the exe

 Thanks to Peter Schattner for providing the framework for this subroutine

=cut

sub program_path {
    my ($self,$program_name) = @_;
    my @path;
    push @path, $self->program_dir if $self->program_dir;
   	# push @path, $program_name .($^O =~ /mswin/i ?'':'');
	# Option for Windows variants / None so far

    return Bio::Root::IO->catfile(@path);
}

=head2 program_dir

 Title   : program_dir
 Usage   : my $dir = $lagan->program_dir();
 Function: Abstract get method for dir of program. To be implemented
           by wrapper.
 Returns : string representing program directory

 Thanks to Peter Schattner for providing the framework for this subroutine

=cut

sub program_dir {
    $PROGRAM_DIR;
}

