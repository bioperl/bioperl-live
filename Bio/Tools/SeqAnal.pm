#-------------------------------------------------------------------------------
# PACKAGE : Bio::Tools::SeqAnal.pm
# PURPOSE : To provide a base class for different sequence analysis tools.
# AUTHOR  : Steve A. Chervitz (sac@genome.stanford.edu)
# CREATED : 27 Mar 1998
# REVISION: $Id$
# STATUS  : Alpha
#
# For documentation, run this module through pod2html 
# (preferably from Perl v5.004 or better).
#-------------------------------------------------------------------------------

package Bio::Tools::SeqAnal;

use Bio::Root::Object ();
use Bio::Root::Global qw(:std);

use strict;
use vars qw($ID $VERSION @ISA);

@ISA        = qw( Bio::Root::Object );
$ID = 'Bio::Tools::SeqAnal';
$VERSION  = 0.011;  


## POD Documentation:

=head1 NAME

Bio::Tools::SeqAnal.pm - Bioperl sequence analysis base class.

=head1 SYNOPSIS

=head2 Object Creation

This module is an abstract base class. Perl will let you instantiate it, 
but it provides little functionality on its own. This module 
should be used via a specialized subclass. See L<_initialize>()
for a description of constructor parameters.

    require Bio::Tools::SeqAnal;

To run and parse a new report:

    $hit = new Bio::Tools::SeqAnal ( -run   => \%runParams,
				     -parse => 1);

To parse an existing report:

    $hit = new Bio::Tools::SeqAnal ( -file  => 'filename.data',
				     -parse => 1);

To run a report without parsing:

    $hit = new Bio::Tools::SeqAnal ( -run   => \%runParams
				     );

To read an existing report without parsing:

    $hit = new Bio::Tools::SeqAnal ( -file  => 'filename.data',
				     -read  => 1);


=head1 INSTALLATION

This module is included with the central Bioperl distribution:

   http://bio.perl.org/Core/Latest
   ftp://bio.perl.org/pub/DIST

Follow the installation instructions included in the README file.


=head1 DESCRIPTION

Bio::Tools::SeqAnal.pm is a base class for specialized
sequence analysis modules such as B<Bio::Tools::Blast.pm> and B<Bio::Tools::Fasta.pm>.
It provides some basic data and functionalities that are not unique to
a specialized module such as:

=over 4

=item * reading raw data into memory.

=item * storing name and version of the program.

=item * storing name of the query sequence.

=item * storing name and version of the database.

=item * storing & determining the date on which the analysis was performed.

=item * basic file manipulations (compress, uncompress, delete).

=back

Some of these functionalities (reading, file maipulation) are inherited from 
B<Bio::Root::Object.pm>, from which Bio::Tools::SeqAnal.pm derives.



=head1 RUN, PARSE, and READ

A SeqAnal.pm object can be created using one of three modes: run, parse, or read.

  MODE      DESCRIPTION
  -----     -----------
  run       Run a new sequence analysis report. New results can then
            be parsed or saved for analysis later.
 
  parse     Parse the data from a sequence analysis report loading it
            into the SeqAnal.pm object.
 
  read      Read in data from an existing raw analysis report without
            parsing it. In the future, this may also permit persistent 
            SeqAnal.pm objects. This mode is considered experimental.

The mode is set by supplying switches to the constructor, see L<_initialize>().



A key feature of SeqAnal.pm is the ability to access raw data in a generic
fashion. Regardless of what sequence analysis method is used, the raw data
always need to be read into memory.  The SeqAnal.pm class utilizes the L<read>() 
method inherited from B<Bio::Root::Object.pm> to permit the following:

=over 4

=item * read from a file or STDIN.

=item * read a single record or a stream containing multiple records.

=item * specify a record separator.

=item * store all input data in memory or process the data stream as it is being read.

=back

By permitting the parsing of data as it is being read, each record can be
analyzed as it is being read and saved or discarded as necessary.
This can be useful when cruching through thousands of reports.
For examples of this, see the L<parse>() methods defined in B<Bio::Tools::Blast.pm> and
B<Bio::Tools::Fasta.pm>.


=head2 Parsing & Running

Parsing and running of sequence analysis reports must be implemented for each 
specific subclass of SeqAnal.pm. No-op stubs ("virtual methods") are provided here for
the L<parse>() and L<run>() methods. See B<Bio::Tools::Blast.pm> and B<Bio::Tools::Fasta.pm>
for examples.


=head1 DEPENDENCIES

Bio::Tools::SeqAnal.pm is a concrete class that inherits from B<Bio::Root::Object.pm>.
This module also makes use of a number of functionalities inherited from
B<Bio::Root::Object.pm> (file manipulations such as reading, compressing, decompressing,
deleting, and obtaining date.


=head1 FEEDBACK

=head2 Mailing Lists 

User feedback is an integral part of the evolution of this and other Bioperl modules.
Send your comments and suggestions preferably to one of the Bioperl mailing lists.
Your participation is much appreciated.

    vsns-bcd-perl@lists.uni-bielefeld.de          - General discussion
    vsns-bcd-perl-guts@lists.uni-bielefeld.de     - Technically-oriented discussion
    http://bio.perl.org/MailList.html             - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track the bugs and 
their resolution. Bug reports can be submitted via email or the web:

    bioperl-bugs@bio.perl.org                   
    http://bio.perl.org/bioperl-bugs/           

=head1 AUTHOR

Steve A. Chervitz, sac@genome.stanford.edu

See the L<FEEDBACK> section for where to send bug reports and comments.

=head1 VERSION

Bio::Tools::SeqAnal.pm, 0.011

=head1 COPYRIGHT

Copyright (c) 1998 Steve A. Chervitz. All Rights Reserved.
This module is free software; you can redistribute it and/or 
modify it under the same terms as Perl itself.


=head1 SEE ALSO

 http://bio.perl.org/Projects/modules.html  - Online module documentation
 http://bio.perl.org/Projects/Blast/        - Bioperl Blast Project     
 http://bio.perl.org/                       - Bioperl Project Homepage


=cut



#
##
###
#### END of main POD documentation.
###
##
#

=head1 APPENDIX

Methods beginning with a leading underscore are considered private
and are intended for internal use by this module. They are
B<not> considered part of the public interface and are described here
for documentation purposes only.

=cut

##############################################################################
##                          CONSTRUCTOR                                     ##
##############################################################################


=head2 _initialize

 Usage     : n/a; automatically called by Bio::Root::Object::new()
 Purpose   : Calls private methods to extract the raw report data,
           : Calls superclass constructor first (Bio::Root::Object.pm).
 Returns   : string containing the make parameter value.
 Argument  : Named parameters (TAGS CAN BE ALL UPPER OR ALL LOWER CASE).
           : The SeqAnal.pm constructor only processes the following
           : parameters passed from new() 
           :     -RUN     => hash reference for named parameters to be used 
           :                 for running a sequence analysis program.
           :                 These are dereferenced and passed to the run() method.
	   :     -PARSE   => boolean, 
	   :     -READ    => boolean,
           :
           : If -RUN is HASH ref, the run() method will be called with the
           :   dereferenced hash.
           : If -PARSE is true, all parameters passed from new() are passed
           :   to the parse() method. This occurs after the run method call
           :   to enable combined running + parsing.
           : If -READ is true, all parameters passed from new() are passed
           :   to the read() method.
           : Either -PARSE or -READ should be true, not both.
 Comments  : Does not calls _rearrange() to handle parameters since only
           : a few are required and there may be potentially many.

See Also   : B<Bio::Root::Object::new()>, B<Bio::Root::Object::_rearrange()>

=cut
#-----------------
sub _initialize {
#-----------------
    my( $self, %param ) = @_;
    
    my $make = $self->SUPER::_initialize(%param);
    
    my($read, $parse, $runparam) = (
	($param{-READ}||$param{'-read'}), ($param{-PARSE}||$param{'-parse'}), 
	($param{-RUN}||$param{'-run'})
				    );

#	$self->_rearrange([qw(READ PARSE RUN)], @param);
	
    # Issue: How to keep all the arguments for running the analysis
    # separate from other arguments needed for parsing the results, etc?
    # Solution: place all the run arguments in a separate hash.

    $self->run(%$runparam) if ref $runparam eq 'HASH';

    if($parse) { $self->parse(%param); }
    elsif($read) { $self->read(%param) }
    
    $make;
}

#--------------
sub destroy {
#--------------
    my $self=shift; 
    $DEBUG==2 && print STDERR "DESTROYING $self ${\$self->name}";
    undef $self->{'_rawData'};
    $self->SUPER::destroy;
}


###############################################################################
#                                 ACCESSORS                                 
###############################################################################

# The mode of the SeqAnal object is no longer explicitly set.
# This simplifies the interface somewhat. 

##----------------------------------------------------------------------
#=head2 mode()

# Usage     : $object->mode();
#	    :
# Purpose   : Set/Get the mode for the sequence analysis object.
#	    :
# Returns   : String
#	    :
# Argument  : n/a
#	    :
#	    :
# Comments  : The mode specifies how much detail to extract from the
#	    : sequence analysis report. There are three modes:
#	    :
#	    :    'parse' -- Parse the sequence analysis output data.
#	    :
#	    :     'read' -- Reads in the raw report but does not 
#	    :               attempt to parse it. Useful when you just
#	    :               want to work with the output as-is 
#	    :               (e.g., create HTML-formatted output).
#	    :
#	    :     'run'  -- Generates a new report. 
#	    :
#	    : Allowable modes are defined by the exported package global array
#	    : @SeqAnal_modes. 
#
#See Also   : _set_mode()
#=cut
##----------------------------------------------------------------------
#sub mode { 
#    my $self = shift; 
#    if(@_) { $self->{'_mode'} = lc(shift); }
#    $self->{'_mode'}; 
#}
#
    

=head2 best

 Usage     : $object->best();
 Purpose   : Set/Get the indicator for processing only the best match.
 Returns   : Boolean (1 | 0)
 Argument  : n/a

=cut

#----------
sub best { 
#----------
    my $self = shift; 
    if(@_) { $self->{'_best'} = shift; }
    $self->{'_best'}; 
}



=head2 _set_db_stats

 Usage     : $object->_set_db_stats(<named parameters>);
 Purpose   : Set stats about the database searched.
 Returns   : String
 Argument  : named parameters:
           :   -LETTERS => <int>  (number of letters in db)
           :   -SEQS    => <int>  (number of sequences in db)

=cut

#-------------------
sub _set_db_stats { 
#-------------------
    my ($self, %param) = @_;

    $self->{'_db'}        ||= $param{-NAME}    || '';
    $self->{'_dbRelease'}   = $param{-RELEASE} || '';
    ($self->{'_dbLetters'}  = $param{-LETTERS} || 0)  =~ s/,//g;
    ($self->{'_dbSeqs'}     = $param{-SEQS}    || 0) =~ s/,//g;
    
}


=head2 database

 Usage     : $object->database();
 Purpose   : Set/Get the name of the database searched.
 Returns   : String
 Argument  : n/a

=cut

#---------------
sub database { 
#---------------
    my $self = shift; 
    if(@_) { $self->{'_db'} = shift; }
    $self->{'_db'};
}



=head2 database_release

 Usage     : $object->database_release();
 Purpose   : Set/Get the release date of the queried database.
 Returns   : String
 Argument  : n/a

=cut

#-----------------------
sub database_release { 
#-----------------------
    my $self = shift; 
    if(@_) { $self->{'_dbRelease'} = shift; }
    $self->{'_dbRelease'};
}


=head2 database_letters

 Usage     : $object->database_letters();
 Purpose   : Set/Get the number of letters in the queried database.
 Returns   : Integer
 Argument  : n/a

=cut

#----------------------
sub database_letters { 
#----------------------
    my $self = shift; 
    if(@_) { $self->{'_dbLetters'} = shift; }
    $self->{'_dbLetters'};
}



=head2 database_seqs

 Usage     : $object->database_seqs();
 Purpose   : Set/Get the number of sequences in the queried database.
 Returns   : Integer
 Argument  : n/a

=cut

#------------------
sub database_seqs { 
#------------------
    my $self = shift; 
    if(@_) { $self->{'_dbSeqs'} = shift; }
    $self->{'_dbSeqs'};
}



=head2 set_date

 Usage     : $object->set_date([<string>]);
 Purpose   : Set the name of the date on which the analysis was performed.
 Argument  : The optional string argument ca be the date or the 
           : string 'file' in which case the date will be obtained from
           : the report file
 Returns   : String
 Throws    : Exception if no date is supplied and no file exists.
 Comments  : This method attempts to set the date in either of two ways:
           :   1) using data passed in as an argument,
           :   2) using the Bio::Root::Utilities.pm file_date() method
           :      on the output file. 
           : Another way is to extract the date from the contents of the
           : raw output data. Such parsing will have to be specialized
           : for different seq analysis reports. Override this method
           : to create such custom parsing code if desired.

See Also   : L<date>(), B<Bio::Root::Object::file_date()>

=cut

#---------------
sub set_date { 
#---------------
    my $self = shift; 
    my $date = shift;
    my ($file);

    if( !$date and ($file = $self->file)) { 
	# If no date is passed and a file exists, determine date from the file.
	# (provided by superclass Bio::Root::Object.pm)
	eval {
	    $date = $self->SUPER::file_date(-FMT => 'd m y');  
	};
	if($@) {
	    $date = 'UNKNOWN';
	    $self->warn("Can't set date of report.");
	}
    }
    $self->{'_date'} = $date;
}



=head2 date

 Usage     : $object->date();
 Purpose   : Get the name of the date on which the analysis was performed.
 Returns   : String
 Argument  : n/a
 Comments  : This method is not a combination set/get, it only gets.
           
See Also   : L<set_date>()

=cut

#----------
sub date {  my $self = shift;  $self->{'_date'}; }
#----------




=head2 length

 Usage     : $object->length();
 Purpose   : Set/Get the length of the query sequence (number of monomers).
 Returns   : Integer
 Argument  : n/a
 Comments  : Developer note: when using the built-in length function within
           : this module, call it as CORE::length().
=cut

#------------
sub length { 
#------------
    my $self = shift;  
    if(@_) { $self->{'_length'} = shift; }
    $self->{'_length'}; 
}



=head2 program

 Usage     : $object->program();
 Purpose   : Set/Get the name of the sequence analysis (BLASTP, FASTA, etc.)
 Returns   : String
 Argument  : n/a

=cut

#-------------
sub program { 
#-------------
    my $self = shift;  
    if(@_) { $self->{'_prog'} = shift; }
    $self->{'_prog'}; 
}



=head2 program_version

 Usage     : $object->program_version();
 Purpose   : Set/Get the version number of the sequence analysis program.
           : (e.g., 1.4.9MP, 2.0a19MP-WashU).
 Returns   : String
 Argument  : n/a

=cut

#---------------------
sub program_version {
#---------------------
    my $self = shift; 
    if(@_) { $self->{'_progVersion'} = shift; }
    $self->{'_progVersion'}; 
}


=head2 query

 Usage     : $name = $object->query();
 Purpose   : Get the name of the query sequence used to generate the report.
 Argument  : n/a
 Returns   : String
 Comments  : Equivalent to $object->name().

=cut

#--------
sub query { my $self = shift; $self->name; }
#--------


=head2 query_desc

 Usage     : $object->desc();
 Purpose   : Set/Get the description of the query sequence for the analysis.
 Returns   : String
 Argument  : n/a

=cut

#--------------
sub query_desc { 
#--------------
    my $self = shift;  
    if(@_) { $self->{'_qDesc'} = shift; }
    $self->{'_qDesc'};
}




=head2 display

 Usage     : $object->display(<named parameters>);
 Purpose   : Display information about Bio::Tools::SeqAnal.pm data members.
           : Overrides Bio::Root::Object::display().
 Example   : $object->display(-SHOW=>'stats');
 Argument  : Named parameters: -SHOW  => 'file' | 'stats'
           :                   -WHERE => filehandle (default = STDOUT)
 Returns   : n/a
 Status    : Experimental

See Also   : L<_display_stats>(), L<_display_file>(), B<Bio::Root::Object::display>

=cut

#---------------
sub display {
#---------------
    my( $self, %param ) = @_;
    
    $self->SUPER::display(%param);
    
    my $OUT = $self->fh();
    $self->show =~ /file/i and $self->_display_file($OUT);
    1;
}



=head2 _display_file

 Usage     : n/a; called automatically by display()
 Purpose   : Print the contents of the raw report file.
 Example   : n/a
 Argument  : one argument = filehandle object.
 Returns   : true (1)
 Status    : Experimental

See Also   : L<display>()  

=cut

#------------------
sub _display_file {
#------------------
    my( $self, $OUT) = @_;

    print $OUT scalar($self->read);
    1;
}



=head2 _display_stats

 Usage     : n/a; called automatically by display()
 Purpose   : Display information about Bio::Tools::SeqAnal.pm data members.
           : Prints the file name, program, program version, database name, 
           : database version, query name, query length, 
 Example   : n/a
 Argument  : one argument = filehandle object.
 Returns   : printf call.
 Status    : Experimental

See Also   : B<Bio::Root::Object.pm>::display()  

=cut

#--------------------
sub _display_stats {
#--------------------
    my( $self, $OUT ) = @_;
    
    printf( $OUT "\n%-15s: %s\n", "QUERY NAME", $self->query ||'UNKNOWN' );
    printf( $OUT "%-15s: %s\n", "QUERY DESC", $self->query_desc || 'UNKNOWN');
    printf( $OUT "%-15s: %s\n", "LENGTH", $self->length || 'UNKNOWN');
    printf( $OUT "%-15s: %s\n", "FILE", $self->file || 'STDIN');
    printf( $OUT "%-15s: %s\n", "DATE", $self->date || 'UNKNOWN');
    printf( $OUT "%-15s: %s\n", "PROGRAM", $self->program || 'UNKNOWN');
    printf( $OUT "%-15s: %s\n", "VERSION", $self->program_version || 'UNKNOWN');
    printf( $OUT "%-15s: %s\n", "DB-NAME", $self->database || 'UNKNOWN');
    printf( $OUT "%-15s: %s\n", "DB-RELEASE", ($self->database_release || 'UNKNOWN'));
    printf( $OUT "%-15s: %s\n", "DB-LETTERS", ($self->database_letters) ? $self->database_letters : 'UNKNOWN');
    printf( $OUT "%-15s: %s\n", "DB-SEQUENCES", ($self->database_seqs) ? $self->database_seqs : 'UNKNOWN');
}


#####################################################################################
##                                 VIRTUAL METHODS                                 ##
#####################################################################################

=head1 VIRTUAL METHODS

=head2 parse

 Usage     : $object->parse( %named_parameters )
 Purpose   : Parse a raw sequence analysis report.
 Returns   : Integer (number of sequence analysis reports parsed).
 Argument  : Named parameters.
 Throws    : Exception: virtual method not defined.
           : Propagates any exception thrown by read()
 Status    : Virtual
 Comments  : This is virtual method that should be overridden to
           : parse a specific type of data.

See Also   : B<Bio::Root::Object::read()>

=cut

#---------
sub parse {
#---------
    my ($self, @param) = @_;

    $self->throw("Virtual method parse() not defined ${ref($self)} objects.");

    # The first step in parsing is reading in the data:
    $self->read(@param);
}



=head2 run

 Usage     : $object->run( %named_parameters )
 Purpose   : Run a sequence analysis program on one or more sequences.
 Returns   : n/a
           : Run mode should be configurable to return a parsed object or
           : the raw results data.
 Argument  : Named parameters:
 Throws    : Exception: virtual method not defined.
 Status    : Virtual

=cut

#--------
sub run {
#--------
    my ($self, %param) = @_;
    $self->throw("Virtual method run() not defined ${ref($self)} objects.");
}


1;
__END__

#####################################################################################
#                                END OF CLASS                                       #
#####################################################################################

 
=head1 FOR DEVELOPERS ONLY

=head2 Data Members

Information about the various data members of this module is provided for those 
wishing to modify or understand the code. Two things to bear in mind: 

=over 4

=item 1 Do NOT rely on these in any code outside of this module. 

All data members are prefixed with an underscore to signify that they are private.
Always use accessor methods. If the accessor doesn't exist or is inadequate, 
create or modify an accessor (and let me know, too!). 

=item 2 This documentation may be incomplete and out of date.

It is easy for these data member descriptions to become obsolete as 
this module is still evolving. Always double check this info and search 
for members not described here.

=back

An instance of Bio::Tools::SeqAnal.pm is a blessed reference to a hash containing
all or some of the following fields:

 FIELD           VALUE
 --------------------------------------------------------------
  _file            Full path to file containing raw sequence analysis report.

  _mode            Affects how much detail to extract from the raw report.
 		   Future mode will also distinguish 'running' from 'parsing'

 
 THE FOLLOWING MAY BE EXTRACTABLE FROM THE RAW REPORT FILE:
 
  _prog            Name of the sequence analysis program.

  _progVersion     Version number of the program.
 		   
  _db              Database searched. 

  _dbRelease       Version or date of the database searched.
 		   
  _dbLetters       Total number of letters in the database.
 		   
  _dbSequences     Total number of sequences in the database.
 		   
  _query           Name of query sequence.

  _length          Length of the query sequence. 

  _date            Date on which the analysis was performed.

 
  INHERITED DATA MEMBERS 
 
  _name            From Bio::Root::Object.pm. String representing the name of the query sequence.
 		   Typically obtained from the report file.

  _parent          From Bio::Root::Object.pm. This member contains a reference to the
 		   object to which this seq anal report belongs. Optional & experimenta.
                   (E.g., a protein object could create and own a Blast object.)

=cut

1;
