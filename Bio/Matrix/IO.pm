#
# BioPerl module for Bio::Matrix::IO
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Jason Stajich <jason-at-bioperl-dot-org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Matrix::IO - A factory for Matrix parsing

=head1 SYNOPSIS

  use Bio::Matrix::IO;
  my $parser = Bio::Matrix::IO->new(-format => 'scoring',
                                    -file   => 'BLOSUMN50');

  my $matrix = $parser->next_matrix;

=head1 DESCRIPTION

This is a general factory framework for writing parsers for Matricies.
This includes parsing output from distance output like PHYLIP's
ProtDist.  Additionally it should be possible to fit parsers for PWM
and PSSMs once their Matrix objects are written.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

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
of the bugs and their resolution. Bug reports can be submitted via
the web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Jason Stajich

Email jason-at-bioperl-dot-org

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Matrix::IO;
use strict;


use base qw(Bio::Root::IO);

=head2 new

 Title   : new
 Usage   : my $obj = Bio::Matrix::IO->new();
 Function: Builds a new Bio::Matrix::IO object 
 Returns : an instance of Bio::Matrix::IO
 Args    :


=cut

sub new { 
  my($caller,@args) = @_;
  my $class = ref($caller) || $caller;
    
    # or do we want to call SUPER on an object if $caller is an
    # object?
    if( $class =~ /Bio::Matrix::IO::(\S+)/ ) {
	my ($self) = $class->SUPER::new(@args);	
	$self->_initialize(@args);
	return $self;
    } else { 

	my %param = @args;
	@param{ map { lc $_ } keys %param } = values %param; # lowercase keys
	my $format = $param{'-format'} || 
	    $class->_guess_format( $param{'-file'} || $ARGV[0] ) ||
		'scoring';
	$format = "\L$format";	# normalize capitalization to lower case

	# normalize capitalization
	return unless( $class->_load_format_module($format) );
	return "Bio::Matrix::IO::$format"->new(@args);
    }
}

=head2 newFh

 Title   : newFh
 Usage   : $fh = Bio::Matrix::IO->newFh(-file=>$filename,-format=>'Format')
 Function: does a new() followed by an fh()
 Example : $fh = Bio::Matrix::IO->newFh(-file=>$filename,-format=>'Format')
           $matrix = <$fh>;   # read a matrix object
           print $fh $matrix; # write a matrix object
 Returns : filehandle tied to the Bio::SeqIO::Fh class
 Args    :

=cut

sub newFh {
  my $class = shift;
  return unless my $self = $class->new(@_);
  return $self->fh;
}

=head2 fh

 Title   : fh
 Usage   : $obj->fh
 Function: Get a filehandle type access to the matrix parser
 Example : $fh = $obj->fh;      # make a tied filehandle
           $matrix = <$fh>;     # read a matrix object
           print $fh $matrix;   # write a matrix object
 Returns : filehandle tied to Bio::Matrix::IO class
 Args    : none

=cut


sub fh {
  my $self = shift;
  my $class = ref($self) || $self;
  my $s = Symbol::gensym;
  tie $$s,$class,$self;
  return $s;
}


=head2 format

 Title   : format
 Usage   : $format = $obj->format()
 Function: Get the matrix format
 Returns : matrix format
 Args    : none

=cut

# format() method inherited from Bio::Root::IO


=head2 next_matrix

 Title   : next_matrix
 Usage   : my $matrix = $matixio->next_matrix;
 Function: Parse the next matrix from the data stream
 Returns : L<Bio::Matrix::MatrixI> type object or undef when finished
 Args    : none


=cut

sub next_matrix{
   my ($self) = @_;
   $self->throw_not_implemented();
}

=head2 write_matrix

 Title   : write_matrix
 Usage   : $io->write_matrix($matrix)
 Function: Writes a matrix out to the data stream
 Returns : none
 Args    : Array of Bio::Matrix::MatrixI object
          - note that not all matricies can be converted to 
            each format, beware with mixing matrix types and output formats

=cut

sub write_matrix{
   my ($self) = @_;
   $self->throw_not_implemented();
}

sub _initialize {
    my ($self,@args) = @_;
    $self->_initialize_io(@args);
}

=head2 _load_format_module

 Title   : _load_format_module
 Usage   : *INTERNAL Matrix::IO stuff*
 Function: Loads up (like use) a module at run time on demand

=cut

sub _load_format_module {
  my ($self,$format) = @_;
  my $module = "Bio::Matrix::IO::" . $format;
  my $ok;
  
  eval {
      $ok = $self->_load_module($module);
  };
  if ( $@ ) {
    print STDERR <<END;
$self: $format cannot be found
Exception $@
For more information about the Matrix::IO system please see the
Matrix::IO docs.  This includes ways of checking for formats at
compile time, not run time
END
  ;
  }
  return $ok;
}


=head2 _guess_format

 Title   : _guess_format
 Usage   : $obj->_guess_format($filename)
 Returns : guessed format of filename (lower case)
 Args    : filename

=cut

sub _guess_format {
   my $class = shift;
   return unless $_ = shift;
   return 'scoring'   if /BLOSUM|PAM$/i;
   return 'phylip'   if /\.dist$/i;
}

sub DESTROY {
    my $self = shift;
    $self->close();
}

sub TIEHANDLE {
  my $class = shift;
  return bless {'matrixio' => shift},$class;
}

sub READLINE {
  my $self = shift;
  return $self->{'matrixio'}->next_tree() || undef unless wantarray;
  my (@list,$obj);
  push @list,$obj  while $obj = $self->{'treeio'}->next_tree();
  return @list;
}

sub PRINT {
  my $self = shift;
  $self->{'matrixio'}->write_tree(@_);
}

1;
