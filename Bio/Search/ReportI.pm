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

Bio::Search::ReportI - Interface for a Search Report (BLAST, FASTA, HMMER)

=head1 SYNOPSIS
{
    # get a Search Report somehow
    # like from a SearchIO
    my $searchio = new Bio::SearchIO(-format => 'blastxml', -file => 'report.xml');
    my $report = $searchio->next_report;
     
    my @searchparams = $report->available_parameters;
    print "program was ", $report->program_name, " ", $report->program_version,
    "\n";

    foreach my $p ( @searchparams ) {
	print "$p = ", $report->get_parameter($p), "\n";
    }
    print "db was ", $report->database_name, 
    ($report->database_size) ? " and size was ".$report->database_size : '',
    "\n";
    print "query was ", $report->query_name, " query length was ", 
    $report->query_size, "\n";
    while( my $subject = $report->next_subject ) {
	# process a Bio::Search::SubjectI object
    }
    
=head1 DESCRIPTION

This is an interface describing the minimal information for a
describing a Search Report.

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

=head2 get_parameter

 Title   : get_parameter
 Usage   : my $gap_ext = $report->get_parameter('gapext')
 Function: Returns the value for a specific parameter used
           when running this report
 Returns : string
 Args    : name of parameter (string)

=cut

sub get_parameter{
   my ($self,@args) = @_;
   $self->_abstractDeath;
}

=head2 available_parameters

 Title   : available_parameters
 Usage   : my @params = $report->available_parameters
 Function: Returns the names of the available parameters
 Returns : Return list of available parameters used for this report
 Args    : none

=cut

sub available_parameters{
   my ($self) = @_;
   $self->_abstractDeath;
}

=head2 get_statistic

 Title   : get_statistic
 Usage   : my $gap_ext = $report->get_statistic('kappa')
 Function: Returns the value for a specific statistic available 
           from this report
 Returns : string
 Args    : name of statistic (string)

=cut

sub get_statistic{
   my ($self,@args) = @_;
   $self->_abstractDeath;
}

=head2 available_statistics

 Title   : available_statistics
 Usage   : my @statnames = $report->available_statistics
 Function: Returns the names of the available statistics
 Returns : Return list of available statistics used for this report
 Args    : none

=cut

sub available_statistics{
   my ($self) = @_;
   $self->_abstractDeath;
}

1;
