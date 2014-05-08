#
# BioPerl module for Bio::PopGen::IO
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Jason Stajich <jason-at-bioperl.org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::PopGen::IO - Input individual,marker,allele information

=head1 SYNOPSIS
 
  use Bio::PopGen::IO;
  my $io = Bio::PopGen::IO->new(-format => 'csv',
                                -file   => 'data.csv');

  # Some IO might support reading in a population at a time

  my @population;
  while( my $ind = $io->next_individual ) {
      push @population, $ind;
  }


=head1 DESCRIPTION

This is a generic interface to reading in population genetic data (of
which there really isn't too many standard formats).  This implementation
makes it easy to provide your own parser for the data.  You need to
only implement one function next_individual.  You can also implement 
next_population if your data has explicit information about population
memberhsip for the indidviduals.

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

Email jason-at-bioperl.org

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...
#TODO 
# Set the Individual creation as a factory rather than
# hardcoded

package Bio::PopGen::IO;
use strict;

# Object preamble - inherits from Bio::Root::Root

use Bio::Root::Root;

use base qw(Bio::Root::IO);

=head2 new

 Title   : new
 Usage   : my $obj = Bio::PopGen::IO->new();
 Function: Builds a new Bio::PopGen::IO object 
 Returns : an instance of Bio::PopGen::IO
 Args    :


=cut

sub new {
  my($class,@args) = @_;

  if( $class =~ /Bio::PopGen::IO::(\S+)/ ) {
    my ($self) = $class->SUPER::new(@args);	
    $self->_initialize(@args);
    return $self;
  } else { 
    my %param = @args;
    @param{ map { lc $_ } keys %param } = values %param; # lowercase keys
    my $format = $param{'-format'} ||
      $class->_guess_format( $param{'-file'} || $ARGV[0] ) || 'csv';

    # normalize capitalization to lower case
    $format = "\L$format";
    
    return unless( $class->_load_format_module($format) );
    return "Bio::PopGen::IO::${format}"->new(@args);
  }
}


=head2 format

 Title   : format
 Usage   : $format = $stream->format()
 Function: Get the PopGen format
 Returns : PopGen format
 Args    : none

=cut

# format() method inherited from Bio::Root::IO


# _initialize is chained for all PopGen::IO classes

sub _initialize {
    my($self, @args) = @_;
#    my ($indfact, $popfact) = $self->_rearrange([qw(INDIVIDUAL_FACTORY
#						    POPULATION_FACTORY)],
#						@args);
#    $indfact = Bio::PopGen::IndividualBuilder->new() unless $indfact;
#    $indfact = Bio::PopGen::PopulationBuilder->new() unless $indfact;

    # initialize the IO part
    $self->_initialize_io(@args);
    return 1;
}

=head2 next_individual

 Title   : next_individual
 Usage   : my $ind = $popgenio->next_individual;
 Function: Retrieve the next individual from a dataset
 Returns : L<Bio::PopGen::IndividualI> object
 Args    : none


=cut

sub next_individual{
    my ($self) = @_;
    $self->throw_not_implemented();
}


=head2 next_population

 Title   : next_population
 Usage   : my $pop = $popgenio->next_population;
 Function: Retrieve the next population from a dataset
 Returns : L<Bio::PopGen::PopulationI> object
 Args    : none
 Note    : Many implementation will not implement this

=cut

sub next_population{
    my ($self) = @_;
    $self->throw_not_implemented();
}

=head2 write_individual

 Title   : write_individual
 Usage   : $popgenio->write_individual($ind);
 Function: Write an individual out in the implementation format
 Returns : none
 Args    : L<Bio::PopGen::PopulationI> object(s)

=cut

sub write_individual{
    my ($self) = @_;
    $self->throw_not_implemented();
}



=head2 write_population

 Title   : write_population
 Usage   : $popgenio->write_population($pop);
 Function: Write a population out in the implementation format
 Returns : none
 Args    : L<Bio::PopGen::PopulationI> object(s)
 Note    : Many implementation will not implement this

=cut

sub write_population{
    my ($self) = @_;
    $self->throw_not_implemented();
}


=head2 newFh

 Title   : newFh
 Usage   : $fh = Bio::SeqIO->newFh(-file=>$filename,-format=>'Format')
 Function: does a new() followed by an fh()
 Example : $fh = Bio::SeqIO->newFh(-file=>$filename,-format=>'Format')
           $sequence = <$fh>;   # read a sequence object
           print $fh $sequence; # write a sequence object
 Returns : filehandle tied to the Bio::SeqIO::Fh class
 Args    :

See L<Bio::SeqIO::Fh>

=cut

sub newFh {
  my $class = shift;
  return unless my $self = $class->new(@_);
  return $self->fh;
}

=head2 fh

 Title   : fh
 Usage   : $obj->fh
 Function:
 Example : $fh = $obj->fh;      # make a tied filehandle
           $sequence = <$fh>;   # read a sequence object
           print $fh $sequence; # write a sequence object
 Returns : filehandle tied to Bio::SeqIO class
 Args    : none

=cut


sub fh {
  my $self = shift;
  my $class = ref($self) || $self;
  my $s = Symbol::gensym;
  tie $$s,$class,$self;
  return $s;
}

=head2 _load_format_module

 Title   : _load_format_module
 Usage   : *INTERNAL Bio::PopGen::IO stuff*
 Function: Loads up (like use) a module at run time on demand
 Example : 
 Returns : 
 Args    : 

=cut

sub _load_format_module {
  my ($self,$format) = @_;
  my $module = "Bio::PopGen::IO::" . $format;
  my $ok;
  
  eval {
      $ok = $self->_load_module($module);
  };
  if ( $@ ) {
      print STDERR <<END;
$self: $format cannot be found
Exception $@
For more information about the Bio::PopGen::IO system please see the 
Bio::PopGen::IO docs.  This includes ways of checking for formats at 
compile time, not run time
END
  ;
  }
  return $ok;
}


=head2 _guess_format

 Title   : _guess_format
 Usage   : $obj->_guess_format($filename)
 Function:
 Example :
 Returns : guessed format of filename (lower case)
 Args    :

=cut


sub _guess_format {
   my $class = shift;
   return unless $_ = shift;
   return 'csv'   if (/csv/i or /\.dat\w$/i);
}

sub close { 
    my $self = shift;
    $self->SUPER::close(@_);
}

sub DESTROY {
    my $self = shift;
    $self->close();
}

sub TIEHANDLE {
  my $class = shift;
  return bless {processor => shift}, $class;
}

sub READLINE {
  my $self = shift;
  return $self->{'processor'}->next_result() || undef unless wantarray;
  my (@list, $obj);
  push @list, $obj while $obj = $self->{'processor'}->next_result();
  return @list;
}

sub PRINT {
  my $self = shift;
  $self->{'processor'}->write_result(@_);
}

1;
