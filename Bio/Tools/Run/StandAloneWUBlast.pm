#
# BioPerl module for Bio::Tools::Run::StandAloneBlast
#
# Copyright Peter Schattner
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Tools::Run::StandAloneWUBlast - Object for the local execution 
of WU-Blast.

=head1 SYNOPSIS

 # Do not use directly; use Bio::Tools::Run::StandAloneBlast

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

package Bio::Tools::Run::StandAloneWUBlast;

use strict;

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
our @WUBLAST_PARAMS  = qw(e s e2 s2 w t x m y z l k h v b q r
    matrix filter wordmask filter maskextra  hitdist wink ctxfactor gape
    gaps gape2 gaps2 gapw gapx olf golf  olmax golmax gapdecayrate
    topcombon topcomboe sumstatsmethod hspsepqmax hspsepsmax gapsepqmax
    gapsepsmax altscore hspmax gspmax qoffset nwstart nwlen qrecmin qrecmax 
    dbrecmin dbrecmax vdbdescmax dbchunks sort_by_pvalue  cpus putenv
    getenv progress);
our @WUBLAST_SWITCH = qw(kap sump poissonp lcfilter lcmask echofilter
    stats nogap gapall pingpong nosegs postsw span2 span1 span prune
    consistency links ucdb gi noseqs qtype qres sort_by_pvalue
    sort_by_count sort_by_highscore sort_by_totalscore
    sort_by_subjectlength mmio nonnegok novalidctxok shortqueryok notes
    warnings errors endputenv getenv endgetenv abortonerror abortonfatal);

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
parameters have defaults and are optional except for -p.

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
    my $self = $caller->SUPER::new(@args);
    
    $self->_set_from_args(\@args, -methods => {(map { $_ => $GENERAL_PARAMS{$_} } keys %GENERAL_PARAMS),
                                               (map { $_ => $_ } (@OTHER_PARAMS,
                                                                  @WUBLAST_PARAMS,
                                                                  @WUBLAST_SWITCH))},
                                  -create => 1,
                                  -force => 1);
    
    my ($tfh, $tempfile) = $self->io->tempfile();
    my $outfile = $self->o || $self->outfile || $tempfile;
    $self->o($outfile);
    close($tfh);
    
    $self->_READMETHOD($DEFAULTREADMETHOD) unless $self->_READMETHOD;
    
    return $self;
}

# We let get/setter method names be case-insensitve
sub AUTOLOAD {
    my $self = shift;
    my $attr = $AUTOLOAD;
    $attr =~ s/.*:://;
    
    my $orig = $attr;
    
    $attr = lc($attr);
    
    $self->can($attr) || $self->throw("Unallowed parameter: $orig !");
    
    return $self->$attr(@_);
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
    my ($self, $input1) = @_;
    $self->io->_io_cleanup();
    my $executable = 'wublast';
    
    # Create input file pointer
    my $infilename1 = $self->_setinput($executable, $input1) || $self->throw("$input1 not Bio::Seq object or array of Bio::Seq objects or file name!");
    $self->i($infilename1);
    
    my $blast_report = $self->_generic_local_wublast($executable);
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
    my $blast_report = $self->_runwublast($executable, $param_string);
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
	my ($self, $executable, $param_string) = @_;
	my ($blast_obj, $exe);
	if (! ($exe = $self->executable($self->p))){
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
    
    # get outputfilename
	my $outfile = $self->o();	
	$blast_obj = Bio::SearchIO->new(-file => $outfile, -format => 'blast');
    
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
    
    @execparams = @WUBLAST_PARAMS;
    
    # of the general params, wublast only takes outfile at
    # this stage (we add in program, input and database manually elsewhere)
    push(@execparams, 'o');
    
    # workaround for problems with shell metacharacters [bug 2707]
    # simply quoting does not always work!
    # Fixed so Windows files are not quotemeta'd
    my $tmp = $self->o;
    $self->o(quotemeta($tmp)) if ($tmp && $^O !~ /^MSWin/);
    
    my $param_string = $self->SUPER::_setparams(-params => [@execparams],
                                                -switches => \@WUBLAST_SWITCH,
                                                -dash => 1);
    
    $self->o($tmp) if ($tmp && $^O !~ /^MSWin/);
    
    if ($self->quiet()) { 
        $param_string .= ' 2> '.File::Spec->devnull;
    }
    
    return $param_string;
}

1;
