#
# BioPerl module for Bio::Tools::AnalysisResult
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Hilmar Lapp <hlapp-at-gmx.net>
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
          ", version ", $result->analysis_method_version(),
          ", performed on ", $result->analysis_date(), "\n";

    # annotate a sequence utilizing SeqAnalysisParserI methods
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
encapsulating parsers and interpreters for the result of a analysis
that was carried out with a query sequence.

The notion of an analysis represented by this base class is that of a
unary or binary operator, taking either one query or a query and a
subject and producing a result. The query is e.g. a sequence, and a
subject is either a sequence, too, or a database of sequences.

This module also implements the Bio::SeqAnalysisParserI interface, and
thus can be used wherever such an object fits.  See
L<Bio::SeqAnalysisParserI>.  Developers will
find a ready-to-use B<parse()> method, but need to implement
B<next_feature()> in an inheriting class. Support for initialization
with input file names and reading from streams is also ready to use.

Note that this module does not provide support for B<running> an
analysis.  Rather, it is positioned in the subsequent parsing step
(concerned with turning raw results into BioPerl objects).

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
the bugs and their resolution.  Bug reports can be submitted via the
web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Hilmar Lapp

Email hlapp-at-gmx.net

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Tools::AnalysisResult;
use strict;

use base qw(Bio::Root::Root Bio::SeqAnalysisParserI Bio::AnalysisResultI Bio::Root::IO);

sub new {
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(@args);
    $self->_initialize(@args);
    return $self;
}

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
           difference to _initialize() is that it may be called at any time,
           and repeatedly within the lifetime of this object. B<Note>, however,
           that this is potentially dangerous in a multi-threading
           environment. In general, calling this method twice is discouraged
           for this reason.

           This method is supposed to reset the state such that any 'history'
           is lost. State information that does not change during object
           lifetime is not considered as history, e.g. parent, name, etc shall
           not be reset. An inheriting object should only be concerned with
           state information it introduces itself, and for everything else
           call SUPER::_initialize_state(@args).

           An example is parsing an input file: a state reset implies
           discarding any unread input, and the actual input itself, followed
           by setting the new input.

           The argument syntax is the same as for L<new()|new> and L<_initialize()|_initialize>,
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

    $self->close();
    $self->_initialize_io(@args);

    $self->{'_analysis_sbjct'} = undef;
    $self->{'_analysis_query'} = undef;
    $self->{'_analysis_prog'} = undef;
    $self->{'_analysis_progVersion'} = undef;
    $self->{'_analysis_date'} = undef;

    return 1;
}

#  =head2 parse
#
#   Title   : parse
#   Usage   : $obj->parse(-input=>$inputobj, [ -params=>[@params] ],
#  		       [ -method => $method ] )
#   Function: Sets up parsing for feature retrieval from an analysis file,
#             or object.
#
#             This method was originally required by SeqAnalysisParserI, but
#             is now discouraged due to potential problems in a multi-
#             threading environment (CORBA!). If called only once, it doesn't
#             add any functionality to calling new() with the same
#             parameters.
#
#             The implementation provided here calls automatically
#             _initialize_state() and passes on -input=>$inputobj and
#             @params as final arguments.
#   Example :
#   Returns : void
#   Args    : B<input>  - object/file where analysis are coming from
#  	   B<params> - parameter to use when parsing/running analysis
#  	   B<method> - method of analysis
#
#  =cut

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


1;
