
#
# BioPerl module for Bio::Search::Processor::Fasta
#
# Cared for by Aaron Mackey <amackey@virginia.edu>
#
# Copyright Aaron Mackey
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Search::Processor::Fasta - Processor of Fasta-generated data streams

=head1 SYNOPSIS

use Bio::Search::Processor

my $processor = new Bio::Search::Processor -file      => 'mysearchrun',
                                           -algorithm => 'Fasta';

=head1 DESCRIPTION

A Processor object is used to generate Bio::Search::Result::* objects, given
a source of Search data (a file or filehandle).  The Processor object works
very much like the SeqIO system: once initialized with the new() method, the
Processor object will continue to yield as many Result objects as are available
from the data source (for single "runs" this is often only one Result object).

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this
and other Bioperl modules. Send your comments and suggestions preferably
 to one of the Bioperl mailing lists.
Your participation is much appreciated.

  vsns-bcd-perl@lists.uni-bielefeld.de          - General discussion
  vsns-bcd-perl-guts@lists.uni-bielefeld.de     - Technically-oriented discussion
  http://bio.perl.org/MailList.html             - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
 the bugs and their resolution.
 Bug reports can be submitted via email or the web:

  bioperl-bugs@bio.perl.org
  http://bio.perl.org/bioperl-bugs/

=head1 AUTHOR - Aaron Mackey

Email amackey@virginia.edu

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Search::Processor::Fasta;
use vars qw($AUTOLOAD @ISA);
use strict;

# Object preamble - inherits from Bio::Root::Object

# "isa" modules
use Bio::Search::Processor::ProcessorI;

# "uses" modules
use Bio::Search::Result::Fasta;
use Bio::Search::Hit::Fasta; # will be contained in Result::Fasta
use Bio::SimpleAlign; # will be contained in Hit::Fasta
                      # should be Bio::Alignment, my 2 cents -AJM

@ISA = qw(Bio::Search::Processor::ProcessorI Exporter);

# new() is inherited from Bio::Root::Object

# _initialize is where the heavy stuff will happen when new is called

sub _initialize {
  my($self,@args) = @_;

  # do any Fasta-specific stuff with @args,
  # then pass the rest back to ProcessorI !
  my $make = $self->SUPER::_initialize(@args);

 return $make; # success - we hope!

}

=head2 next_result

 Title   : next_result
 Usage   : $result = $processor->next_result()
 Function: Used to obtain Result objects from a Fasta-generated data source
 Returns : a Bio::Search::Result::Fasta object
 Args    : <none>


=cut

sub next_result{
   my ($self,@args) = @_;

   my $fh = $self->_filehandle() or $self->throw("No data source!");

   my $result = new Bio::Search::Result::Fasta;

   # process data from $fh and set appropriate values in $result using
   # $result's publically available accessors, including creating and
   # storing Hit objects which themselves have created and contained
   # Alignment objects.  Sky's the limit here.

   return $result;
}

1;

__END__
