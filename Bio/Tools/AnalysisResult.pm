
#
# BioPerl module for Bio::Tools::AnalysisResult
#
# Cared for by Hilmar Lapp <hlapp@gmx.net>
#
# Copyright Hilmar Lapp
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Tools::AnalysisResult - Base class for analysis result objects and parsers

=head1 SYNOPSIS

    # obtain a AnalysisResult derived object somehow
    print "Method ", $result->analysis_method(),
          ", version " $result->analysis_method_version(), 
          ", performed on ", $result->analysis_date(), "\n";
    # parse the result (utilizing SeqAnalysisParserI methods)
    $result->parse('-input' => 'myfile.out');
    # and annotate (e.g., a sequence)
    while($feat = $result->next_feature()) {
	$seq->add_SeqFeature($feat);
    }
    $result->close();
    # query object, e.g. a Bio::SeqI implementing object
    $queryseq = $result->analysis_query();
    # Subject of the analysis -- may be undefined. Refer to derived module
    # to find out what is returned.
    $subject = $result->analysis_subject();

=head1 DESCRIPTION

The AnalysisResult module is supposed to be the base class for modules
encapsulating parsers and interpreters for the result of a analysis that was
carried out with a query sequence.

The notion of an analysis represented by this base class is that of a unary or
binary operator, taking either one query or a query and a subject and producing
a result. The query is e.g. a sequence, and a subject is either a sequence,
too, or a database of sequences. 

This module also implements the Bio::SeqAnalysisParserI interface, and thus
can be used wherever such an object fits. See L<Bio::SeqAnalysisParserI>.
Developers will find a ready-to-use L<parse>() method, but need to implement 
next_feature() in an inheriting class. Support for initialization with input
file names and reading from streams is also ready to use.

Note that this module does not provide support for B<running> an analysis.
Rather, it is positioned in the subsequent parsing step (concerned with
turning raw results into BioPerl objects).

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists.  Your participation is much appreciated.

  bioperl-l@bioperl.org          - General discussion
  http://bio.perl.org/MailList.html             - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via email
or the web:

  bioperl-bugs@bio.perl.org
  http://bio.perl.org/bioperl-bugs/

=head1 AUTHOR - Hilmar Lapp

Email hlapp@gmx.net (or hilmar.lapp@pharma.novartis.com)

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Tools::AnalysisResult;
use vars qw(@ISA);
use strict;
use Symbol;

# Object preamble - inherits from Bio::Root::Object

use Bio::Root::Object;
use Bio::SeqAnalysisParserI;

@ISA = qw(Bio::Root::Object Bio::SeqAnalysisParserI);
# new() is inherited from Bio::Root::Object

# _initialize is where the heavy stuff will happen when new is called

sub _initialize {
  my($self,@args) = @_;

  my $make = $self->SUPER::_initialize(@args);

  $self->_initialize_state(@args);
  return $make; # success - we hope!
}

=head2 _initialize_state

 Title   : _initialize_state
 Usage   : n/a; usually called by _initialize()
 Function: This method is for BioPerl B<developers> only, as indicated by the
           leading underscore in its name.

           Performs initialization or reset of the state of this object. The
           difference to L<_initialize>() is that it may be called at any time,
           and repeatedly within the lifetime of this object. 

           This method is supposed to reset the state such that any 'history'
           is lost. State information that does not change during object
           lifetime is not considered as history, e.g. parent, name, etc shall
           not be reset. An inheriting object should only be concerned with
           state information it introduces itself, and for everything else
           call SUPER::_initialize_state(@args).

           An example is parsing an input file: a state reset implies
           discarding any unread input, and the actual input itself, followed
           by setting the new input.

           The argument syntax is the same as for L<new>() and L<_initialize>(),
           i.e., named parameters following the -name=>$value convention.
           The following parameters are dealt with by the implementation
           provided here:
              -INPUT, -FH, -FILE
           (tags are case-insensitive).
 Example :
 Returns :
 Args    :

=cut

sub _initialize_state {
    my ($self,@args) = @_;
    my ($input,$fh,$file) =
	$self->_rearrange([qw(INPUT
			      FH
			      FILE
			      )],
			  @args);

    $self->close();
    $self->{'readbuffer'} = "";
    $self->{'_analysis_sbjct'} = undef;
    $self->{'_analysis_query'} = undef;
    $self->{'_analysis_prog'} = undef;
    $self->{'_analysis_progVersion'} = undef;
    $self->{'_analysis_date'} = undef;

    # determine whether the input is a file(name) or a stream
    if($input) {
	if(ref(\$input) eq "SCALAR") {
	    # we assume that a scalar is a filename
	    if($file && ($file ne $input)) {
		$self->throw("input file given twice: $file and $input disagree");
	    }
	    $file = $input;
	} elsif(ref($input) eq "GLOB") {
	    # input is a stream
	    $fh = $input;
	} else {
	    # let's be strict for now
	    $self->throw("unable to determine type of input $input: ".
			 "not string and not GLOB");
	}
    }
    if( defined $fh && defined $file ) {
	$self->throw("You have defined both a filehandle and file to read from. Not good news!");
    }
    if((defined $file) && ($file ne '')) {
	$fh = Symbol::gensym();
	open ($fh,$file)
	    || $self->throw("Could not open $file for Fasta stream reading $!");
    }
    if((! defined($fh)) && ($file eq "")) {
	$fh = \*STDIN;
    }
    $self->_filehandle($fh) if defined $fh;
}

=head2 close

 Title   : close
 Usage   : $result->close()
 Function: Closes the file handle associated with this result file
 Example :
 Returns :
 Args    :

=cut

sub close {
   my ($self, @args) = @_;

   $self->{'_filehandle'} = undef;
}

=head2 parse

 Title   : parse
 Usage   : $obj->parse(-input=>$inputobj, [ -params=>[@params] ],
		       [ -method => $method ] )
 Function: Sets up parsing for feature retrieval from an analysis file, 
           or object

           This method is required for all classes implementing the
           SeqAnalysisParserI interface.

           The implementation provided here calls automatically
           _initialize_state() with passing on -input=>$inputobj and
           @params as final arguments.
 Example :
 Returns : void
 Args    : B<input>  - object/file where analysis are coming from
	   B<params> - parameter to use when parsing/running analysis
	   B<method> - method of analysis

=cut

sub parse {
    my ($self, @args) = @_;

    my ($input, $params, $method) = 
	$self->_rearrange([qw(INPUT
			      PARAMS
			      METHOD
			      )],
			  @args);

    # initialize with new input
    if($params) {
	$self->_initialize_state('-input' => $input, @$params);
    } else {
	$self->_initialize_state('-input' => $input);
    }
    $self->analysis_method($method) if $method;
}

=head2 analysis_query

 Usage     : $query_obj = $result->analysis_query();
 Purpose   : Set/Get the name of the query used to generate the result, that
             is, the entity on which the analysis was performed. Will mostly
             be a sequence object (Bio::PrimarySeq compatible).
 Argument  : 
 Returns   : The object set before. Mostly a Bio::PrimarySeq compatible object.

=cut

#--------
sub analysis_query {
    my ($self, $obj) = @_;
    if($obj) {
	$self->{'_analysis_query'} = $obj;
    }
    return $self->{'_analysis_query'};
}
#--------

=head2 analysis_subject

 Usage     : $result->analyis_subject();
 Purpose   : Set/Get the subject of the analysis against which it was
             performed. For similarity searches it will probably be a database,
             and for sequence feature predictions (exons, promoters, etc) it
             may be a collection of models or homologous sequences that were
             used, or undefined.
 Returns   : The object that was set before, or undef.
 Argument  : 

=cut

#---------------
sub analysis_subject { 
#---------------
    my ($self, $sbjct_obj) = @_; 
    if($sbjct_obj) {
	$self->{'_analysis_sbjct'} = $sbjct_obj;
    }
    return $self->{'_analysis_sbjct'};
}


=head2 analysis_date

 Usage     : $result->analysis_date();
 Purpose   : Set/Get the date on which the analysis was performed.
 Returns   : String
 Argument  : 
 Comments  : 

=cut

#----------
sub analysis_date {
    my ($self, $date) = @_;
    if($date) {
	$self->{'_analysis_date'} = $date;
    }
    return $self->{'_analysis_date'};
}
#----------

=head2 analysis_method

 Usage     : $result->analysis_method();
 Purpose   : Set/Get the name of the sequence analysis method that was used
             to produce this result (BLASTP, FASTA, etc.). May also be the
             actual name of a program.
 Returns   : String
 Argument  : n/a

=cut

#-------------
sub analysis_method { 
#-------------
    my ($self, $method) = @_;  
    if($method) {
	$self->{'_analysis_prog'} = $method;
    }
    return $self->{'_analysis_prog'}; 
}

=head2 analysis_method_version

 Usage     : $result->analysis_method_version();
 Purpose   : Set/Get the version string of the analysis program.
           : (e.g., 1.4.9MP, 2.0a19MP-WashU).
 Returns   : String
 Argument  : n/a

=cut

#---------------------
sub analysis_method_version {
#---------------------
    my ($self, $version) = @_; 
    if($version) {
	$self->{'_analysis_progVersion'} = $version;
    }
    return $self->{'_analysis_progVersion'}; 
}

=head2 _pushback

 Title   : _pushback
 Usage   : $obj->_pushback($newvalue)
 Function: puts a line previously read with _readline back into a buffer

           Copied from Bio/SeqIO.pm.
 Example :
 Returns :
 Args    : newvalue

=cut

sub _pushback {
  my ($obj, $value) = @_;
  $obj->{'readbuffer'} .= $value;
}


=head2 _filehandle

 Title   : _filehandle
 Usage   : $obj->_filehandle($newval)
 Function:
 Example :
 Returns : value of _filehandle
 Args    : newvalue (optional)


=cut

sub _filehandle {
    my ($obj, $value) = @_;
    if(defined $value) {
	$obj->{'_filehandle'} = $value;
    }
    return $obj->{'_filehandle'};
}


=head2 _readline

 Title   : _readline
 Usage   : $obj->_readline
 Function: Reads a line of input and returns it.

           Note that this method implicitely uses the value of $/ that is
           in effect when called.

           Note also that the current implementation does not handle pushed
           back input correctly unless the pushed back input ends with the
           value of $/.

           Copied from Bio/SeqIO.pm.
 Example :
 Returns : 

=cut

sub _readline {
  my $self = shift;
  my $fh = $self->_filehandle();
  my $line;

  # if the buffer been filled by _pushback then return the buffer
  # contents, rather than read from the filehandle
  if ( defined $self->{'readbuffer'} ) {
      $line = $self->{'readbuffer'};
      undef $self->{'readbuffer'};
  } else {
      $line = defined($fh) ? <$fh> : <>;
  }
  return $line;
}

sub DESTROY {
    my $self = shift;

    $self->close();
}

1;
