# $Id$
#
# BioPerl module for Bio::Search::ReportI
#
# Cared for by Jason Stajich <jason@bioperl.org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Search::ReportI - DESCRIPTION of Interface

=head1 SYNOPSIS

Give standard usage here

=head1 DESCRIPTION

Describe the interface here

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
  http://bioperl.org/bioperl-bugs/

=head1 AUTHOR - Jason Stajich

Email jason@bioperl.org

Describe contact details here

=head1 CONTRIBUTORS

Additional contributors names and emails here

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Search::ReportI;
use vars qw(@ISA);
use strict;
use Carp;

sub _abstractDeath {
  my $self = shift;
  my $package = ref $self;
  my $caller = (caller)[1];
  confess "Abstract method '$caller' defined in interface Bio::Search::ReportI not implemented by pacakge $package. Not your fault - author of $package should be blamed!";
}

=head2 next_subject

 Title   : next_subject
 Usage   : my $subject = $report->next_subject;
 Function: Returns the next Subject from a search
 Returns : Bio::Search::SubjectI object
 Args    : none

=cut

sub next_subject{
   my ($self) = @_;
   $self->_abstractDeath;

}

=head2 database_name

 Title   : database_name
 Usage   : my $name = $report->database_name;
 Function: Returns the name of database searched 
 Returns : string
 Args    : none

=cut

sub database_name{
   my ($self,@args) = @_;
   $self->_abstractDeath;
}

=head2 database_size

 Title   : database_size
 Usage   : my $size = $report->database_size
 Function: Returns the size of the database searched
 Returns : integer
 Args    : none


=cut

sub database_size{
   my ($self,@args) = @_;
   $self->_abstractDeath;
}

=head2 query_name

 Title   : query_name
 Usage   : my $q_name = $report->query_name
 Function: Returns the name of the query sequence used to search the database
 Returns : string
 Args    : none


=cut

sub query_name{
   my ($self,@args) = @_;
   $self->_abstractDeath;
}

=head2 query_size

 Title   : query_size
 Usage   : my $qsize = $report->query_size;
 Function: Returns the size of the query sequence used to search the database
 Returns : integer
 Args    : none


=cut

sub query_size{
   my ($self,@args) = @_;
   $self->_abstractDeath;
}

=head2 program_name

 Title   : program_name
 Usage   : my $prog_name = $report->program_name
 Function: Returns the full name of the program that generated this report
 Returns : String
 Args    : none


=cut

sub program_name{
   my ($self,@args) = @_;
   $self->_abstractDeath;

}

=head2 program_version

 Title   : program_version
 Usage   : my $version = $report->program_version
 Function: Returns the version number of the program which generated 
           this report
 Returns : String
 Args    : none


=cut

sub program_version{
   my ($self,@args) = @_;
   $self->_abstractDeath;   

}

1;
