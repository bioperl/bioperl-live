
#
# BioPerl module for Bio::Search::Result::ResultI.pm
#
# Cared for by Aaron Mackey <amackey@virginia.edu>
#
# Copyright Aaron Mackey
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

    Bio::Search::Result::ResultI.pm - Abstract interface to Result objects

=head1 SYNOPSIS

    These objects are generated automatically by Bio::Search::Processor
objects, and wouldn't be used directly.

=head1 DESCRIPTION

Bio::Search::Result::* objects are data structures containing the results from
the execution of a search algorithm.  As such, it may contain various
algorithm specific information as well as details of the execution, but will
contain a few fundamental elements, including the ability to return
Bio::Search::Hit::* objects.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this
and other Bioperl modules. Send your comments and suggestions preferably
 to one of the Bioperl mailing lists.
Your participation is much appreciated.

  bioperl-l@bioperl.org              - General discussion
  http://bio.perl.org/MailList.html  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
 the bugs and their resolution.
 Bug reports can be submitted via email or the web:

  bioperl-bugs@bio.perl.org
  http://bio.perl.org/bioperl-bugs/

=head1 AUTHOR - Aaron Mackey

Email amackey@virginia.edu

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut

#'
# Let the code begin...


package Bio::Search::Result::ResultI;

use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::Root::Object

use Bio::Root::RootI;

@ISA = qw(Bio::Root::RootI);


=head2 get_query_id

 Title   : get_query_id
 Usage   : $id = $result->get_query_id();
 Function: Used to obtain the string identifier of the query used by the
           algorithm.
 Returns : a scalar string.
 Args    : <none>

=cut

sub get_query_id {
   my ($self,@args) = @_;

   return $self->throw("Abstract interface call.");
}

=head2 get_library_name

 Title   : get_library_name
 Usage   : $name = $result->get_library_name()
 Function: Used to obtain the name of the library that the query was searched
           against by the algorithm.
 Returns : a scalar string
 Args    : <none>

=cut

sub get_library_name {
   my ($self,@args) = @_;

   return $self->throw("Abstract interface call.");
}

=head2 get_library_size

 Title   : get_library_size
 Usage   : $size = $result->get_library_size()
 Function: Used to obtain the size of library that was searched against.
 Returns : a scalar integer (units specific to algorithm, but probably the
           total number of residues in the library, if available) or undef if
           the information was not available to the Processor object.
 Args    : <none>


=cut

sub get_library_size {
   my ($self,@args) = @_;

   return $self->throw("Abstract interface call.");
}

=head2 get_library_count

 Title   : get_library_count
 Usage   : $count = $result->get_library_count()
 Function: Used to obtain the number of entries contained in the library.
 Returns : a scalar integer representing the number of entities in the library
           or undef if the information was not available.
 Args    : <none>


=cut

sub get_library_count {
   my ($self,@args) = @_;

   return $self->throw("Abstract interface call.");
}

=head2 get_hits

 Title   : get_hits
 Usage   : @hits = $result->get_hits();
 Function: Used to obtain the array of hit objects, representing potential
           matches between the query and various entities from the library.
 Returns : an array of Bio::Search::Hit::* object (specific type depends on
           algorithm), or undef if there were none.
 Args    : <none>


=cut

sub get_hits {
   my ($self,@args) = @_;

   return $self->throw("Abstract interface call.");
}


1;


