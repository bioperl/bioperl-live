# $Id$
#
# BioPerl module for Bio::SearchIO::SearchResultEventBuilder
#
# Cared for by Jason Stajich <jason@bioperl.org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::SearchIO::SearchResultEventBuilder - DESCRIPTION of Object

=head1 SYNOPSIS

Give standard usage here

=head1 DESCRIPTION

Describe the object here

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


package Bio::SearchIO::SearchResultEventBuilder;
use vars qw(@ISA %KNOWNEVENTS);
use strict;

use Bio::Root::RootI;
use Bio::SearchIO::EventHandlerI;
use Bio::Search::Report;
use Bio::Search::Subject;
use Bio::Search::HSP;


@ISA = qw(Bio::Root::RootI Bio::SearchIO::EventHandlerI);

=head2 new

 Title   : new
 Usage   : my $obj = new Bio::SearchIO::SearchResultEventBuilder();
 Function: Builds a new Bio::SearchIO::SearchResultEventBuilder object 
 Returns : Bio::SearchIO::SearchResultEventBuilder
 Args    :


=cut

sub new {
  my($class,@args) = @_;

  my $self = $class->SUPER::new(@args);
  return $self;
}

=head2 will_handle

 Title   : will_handle
 Usage   : if( $handler->will_handle($event_type) ) { ... }
 Function: Tests if this event builder knows how to process a specific event
 Returns : boolean
 Args    : event type name


=cut

sub will_handle{
   my ($self,$type) = @_;
   return ( $type eq 'hsp' || $type eq 'subject' || $type eq 'report' );
}

=head2 SAX methods

=head2 start_report

 Title   : start_report
 Usage   : $handler->start_report($reporttype)
 Function: Begins a report event cycle
 Returns : none 
 Args    : Type of Report

=cut

sub start_report {
   my ($self,$type) = @_;
   $self->{'_reporttype'} = $type;
   $self->{'_subjects'} = [];
   
   return;
}

=head2 end_report

 Title   : end_report
 Usage   : my @reports = $parser->end_report
 Function: Finishes a report handler cycle
 Returns : A Bio::Search::ReportI
 Args    : none

=cut

sub end_report {
    my ($self,$type,$data) = @_;
    my $report = new Bio::Search::Report
	('-db_name'     => $data->{'dbname'},
	 '-db_size'     => $data->{'dbsize'},
	 '-query_name'  => $data->{'queryname'},
	 '-query_size'  => $data->{'querylen'},
	 '-program_name'=> $data->{'programname'},
	 '-program_version'=> $data->{'programver'},
	 '-report_type' => $type,
	 '-subjects'    => $self->{'_subjects'} );
    $self->{'_subjects'} = [];
    return $report;
}

=head2 start_hsp

 Title   : start_hsp
 Usage   : $handler->start_hsp($name,$data)
 Function: Processes a HSP
 Returns : none
 Args    : type of element 
           associated data (hashref)

=cut

sub start_hsp {
    my ($self,@args) = @_;
   return;
}

=head2 end_hsp

 Title   : end_hsp
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub end_hsp{
   my ($self,$type,$data) = @_;

   my $hsp = new Bio::Search::HSP
       (
	'-report_type'   => $type,
	'-score'         => $data->{'score'},
	'-bits'          => $data->{'bits'},
	'-match'         => $data->{'match'},
	'-hsp_length'    => $data->{'hsplen'},
	'-positive'      => $data->{'positive'},
	'-gaps'          => $data->{'gaps'},
	'-evalue'        => $data->{'evalue'},
	'-query_begin'   => $data->{'querystart'},
	'-query_end'     => $data->{'queryend'},
	'-subject_begin' => $data->{'subjectstart'},
	'-subject_end'   => $data->{'subjectend'},
	'-query_seq'     => $data->{'queryseq'},
	'-subject_seq'   => $data->{'subjectseq'},
	'-homology_seq'  => $data->{'homolseq'},
	'-query_length'  => $data->{'querylen'},
	'-subject_length'=> $data->{'subjectlen'},
	'-query_name'    => $data->{'queryname'},
	'-subject_name'  => $data->{'subjectname'},
	'-query_frame'   => $data->{'queryframe'},
	'-subject_frame' => $data->{'subjectframe'},
	);
   push @{$self->{'_hsps'}}, $hsp;
   return;
}

=head2 start_subject

 Title   : start_subject
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub start_subject{
    my ($self,$type) = @_;
    $self->{'_hsps'} = [];    
    return;
}

=head2 end_subject

 Title   : end_subject
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub end_subject{
   my ($self,$type,$data) = @_;
    my $subject = new Bio::Search::Subject
	( '-report_type' => $type,
	  '-name'        => $data->{'subjectname'}, 
	  '-length'      => $data->{'subjectlen'},
	  '-accession'   => $data->{'subjectacc'},
	  '-desc'        => $data->{'subjectdesc'},
	  '-hsps'        => $self->{'_hsps'});
   push @{$self->{'_subjects'}}, $subject;
   $self->{'_hsps'} = [];
   return;
}

1;
