#-----------------------------------------------------------------
# $Id$
#
# BioPerl module Bio::Search::Result::ResultI
#
# Cared for by Steve Chervitz <sac@bioperl.org>
#
# Originally created by Aaron Mackey <amackey@virginia.edu>
#
# You may distribute this module under the same terms as perl itself
#-----------------------------------------------------------------

# POD documentation - main docs before the code

=head1 NAME

Bio::Search::Result::ResultI - Abstract interface to Search Result objects

=head1 SYNOPSIS

Bio::Search::Result::ResultI objects cannot be instantiated since this
module defines a pure interface.

Given an object that implements the Bio::Search::Result::ResultI  interface,
you can do the following things with it:

    while( $hit = $result->next_hits()) { ... }

    $id = $result->query_name();

    $id = $result->query_description();

    $name = $result->database_name()

    $size = $result->database_letters()

    $num_entries = $result->database_entries()

    my $gap_ext = $result->get_parameter('gapext')

    my @params = $result->available_parameters

    my $gap_ext = $result->get_statistic('kappa')

    my @statnames = $result->available_statistics

=head1 DESCRIPTION

Bio::Search::Result::ResultI objects are data structures containing
the results from the execution of a search algorithm.  As such, it may
contain various algorithm specific information as well as details of
the execution, but will contain a few fundamental elements, including
the ability to return Bio::Search::Hit::HitI objects.

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

=head1 AUTHOR 

Aaron Mackey E<lt>amackey@virginia.eduE<gt>  (original author)

Steve Chervitz E<lt>sac@bioperl.orgE<gt>

See L<the FEEDBACK section | FEEDBACK> for where to send bug reports and comments.

=head1 COPYRIGHT

Copyright (c) 1999-2001 Aaron Mackey, Steve Chervitz. All Rights Reserved.

=head1 DISCLAIMER

This software is provided "as is" without warranty of any kind.

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut

#'
# Let the code begin...


package Bio::Search::Result::ResultI;

use strict;
use vars qw(@ISA);

use Bio::AnalysisResultI;

@ISA = qw( Bio::AnalysisResultI );


=head2 next_hit

 Title   : next_hit
 Usage   : while( $hit = $result->next_hit()) { ... }
 Function: Returns the next available Hit object, representing potential
           matches between the query and various entities from the database.
 Returns : a Bio::Search::Hit::HitI object or undef if there are no more.
 Args    : none


=cut

sub next_hit {
    my ($self,@args) = @_;
    $self->throw_not_implemented;
}

=head2 query_name

 Title   : query_name
 Usage   : $id = $result->query_name();
 Function: Get the string identifier of the query used by the
           algorithm that performed the search.
 Returns : a string.
 Args    : none

=cut

sub query_name {
    my ($self,@args) = @_;
    $self->throw_not_implemented;
}

=head2 query_accession

 Title   : query_accession
 Usage   : $id = $result->query_accession();
 Function: Get the accession (if available) for the query sequence
 Returns : a string
 Args    : none

=cut

sub query_accession {
    my ($self,@args) = @_;
    $self->throw_not_implemented;
}


=head2 query_length

 Title   : query_length
 Usage   : $id = $result->query_length();
 Function: Get the length of the query sequence
           used in the search.
 Returns : a number
 Args    : none

=cut

sub query_length {
    my ($self,@args) = @_;
    $self->throw_not_implemented;
}

=head2 query_description

 Title   : query_description
 Usage   : $id = $result->query_description();
 Function: Get the description of the query sequence
           used in the search.
 Returns : a string
 Args    : none

=cut

sub query_description {
    my ($self,@args) = @_;
    $self->throw_not_implemented;
}


=head2 database_name

 Title   : database_name
 Usage   : $name = $result->database_name()
 Function: Used to obtain the name of the database that the query was searched
           against by the algorithm.
 Returns : a scalar string
 Args    : none

=cut

sub database_name {
    my ($self,@args) = @_;

    $self->throw_not_implemented;
}

=head2 database_letters

 Title   : database_letters
 Usage   : $size = $result->database_letters()
 Function: Used to obtain the size of database that was searched against.
 Returns : a scalar integer (units specific to algorithm, but probably the
           total number of residues in the database, if available) or undef if
           the information was not available to the Processor object.
 Args    : none


=cut

sub database_letters {
    my ($self,@args) = @_;
    $self->throw_not_implemented();
}

=head2 database_entries

 Title   : database_entries
 Usage   : $num_entries = $result->database_entries()
 Function: Used to obtain the number of entries contained in the database.
 Returns : a scalar integer representing the number of entities in the database
           or undef if the information was not available.
 Args    : none


=cut

sub database_entries {
    my ($self,@args) = @_;

    $self->throw_not_implemented();
}

=head2 get_parameter

 Title   : get_parameter
 Usage   : my $gap_ext = $result->get_parameter('gapext')
 Function: Returns the value for a specific parameter used
           when running this result
 Returns : string
 Args    : name of parameter (string)

=cut

sub get_parameter{
   my ($self,@args) = @_;
   $self->throw_not_implemented();
}

=head2 available_parameters

 Title   : available_parameters
 Usage   : my @params = $result->available_parameters
 Function: Returns the names of the available parameters
 Returns : Return list of available parameters used for this result
 Args    : none

=cut

sub available_parameters{
   my ($self) = @_;
   $self->throw_not_implemented();
}

=head2 get_statistic

 Title   : get_statistic
 Usage   : my $gap_ext = $result->get_statistic('kappa')
 Function: Returns the value for a specific statistic available 
           from this result
 Returns : string
 Args    : name of statistic (string)

=cut

sub get_statistic{
   my ($self,@args) = @_;
   $self->throw_not_implemented();
}

=head2 available_statistics

 Title   : available_statistics
 Usage   : my @statnames = $result->available_statistics
 Function: Returns the names of the available statistics
 Returns : Return list of available statistics used for this result
 Args    : none

=cut

sub available_statistics{
   my ($self) = @_;
   $self->throw_not_implemented();
}

1;


