
#
# BioPerl module for Bio::Search::Processor::ProcessorI
#
# Cared for by Aaron Mackey <amackey@virginia.edu>
#
# Copyright Aaron Mackey
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Search::Processor::ProcessorI - Abstract Interface for Processor Objects

=head1 SYNOPSIS

  use Bio::Search::Processor

  my $processor = new Bio::Search::Processor -file      => 'mysearchrun',
                                             -algorithm => 'Blast';

  while ($result = $processor->next_result()) {

    $id = $result->get_query_id();
    $lib = $result->get_library_name();
    $size = $result->get_library_size();

    foreach $hit ( $result->get_hits() ) {
        $matchid = $hit->get_id();
        $matchdesc = $hit->get_desc();
        # etc, etc, do stuff.
    }
  }

=head1 DESCRIPTION

A Processor object is used to generate Bio::Search::Result::* objects, given
a source of Search data (a file or filehandle).  The Processor object works
very much like the SeqIO system: once initialized with the new() method, the
Processor object will continue to yield as many Result objects as are available
from the data source (for single "runs" this is often only one Result object).

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists.  Your participation is much appreciated.

  bioperl-l@bioperl.org            - General discussion
  http://bio.perl.org/MailList.html - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via email
or the web:

  bioperl-bugs@bio.perl.org
  http://bio.perl.org/bioperl-bugs/

=head1 AUTHOR - Aaron Mackey

Email amackey@virginia.edu

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Search::Processor::ProcessorI;

use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::Root::Object

use Bio::Root::Root;

@ISA = qw(Bio::Root::Root);

=head2 next_result

 Title   : next_result
 Usage   : $result = $processor->next_result()
 Function: Returns the next Bio::Search::Result::* object available from
           the provided data stream or undef if no more are available.
 Returns : Bio::Search::Result object
 Args    : <none>


=cut

sub next_result{
   my ($self,@args) = @_;

   $self->throw("Abstract processor call of next_result() - your processor type has not implemented this method!");

}

__END__


