# $Id$
#
# BioPerl module for Bio::Search::Subject
#
# Cared for by Jason Stajich <jason@bioperl.org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Search::Subject - DESCRIPTION of Object

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


package Bio::Search::Subject;
use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::Root::Root

use Bio::Root::Root;
use Bio::Search::SubjectI;

@ISA = qw(Bio::Root::Root Bio::Search::SubjectI );

=head2 new

 Title   : new
 Usage   : my $obj = new Bio::Search::Subject(-name => $name, -type => $type);
 Function: Builds a new Bio::Search::Subject object 
 Returns : Bio::Search::Subject
 Args    : -name      => Subject name
           -type      => Report Type
           -length    => Subject length
           -accession => Subject Accession
           -desc      => Subject Description 
           -hsps => array ref of Bio::Search::HSPIs

=cut


sub new {
  my($class,@args) = @_;

  my $self = $class->SUPER::new(@args);
  my $hsps;
  ($self->{'_name'},  $self->{'_type'},
      $self->{'_length'}, $self->{'_accession'},
      $self->{'_description'},
      $hsps) = $self->_rearrange([qw(NAME
				     REPORT_TYPE
				     LENGTH
				     ACCESSION
				     DESC
				     HSPS)],
				 @args);
  $self->{'_hspindex'} = 0;
  $self->{'_hsps'} = [];
  if( defined $hsps  ) {
      $self->throw("Must define arrayref of HSPIs when initializing a $class\n") 
	  unless ref($hsps) =~ /array/i;
      
      foreach my $h ( @$hsps ) {
	  $self->add_hsp($h);
      }
  }
  return $self;
}

=head2 Bio::Search::SubjectI methods

=cut

=head2 name

 Title    : name
 Usage    : $name = $obj->name();
 Function : returns the name of the Sbjct 
 Example  : 
 Returns  : name of the Sbjct 
 Args     :

=cut

sub name { 
    my ($self) = @_;
    return $self->{'_name'}; 
}

=head2 report_type

 Title    : report_type
 Usage    : $type = $sbjct->report_type()
 Function : Returns the type of report from which this subject was obtained.
            This usually pertains only to BLAST and friends reports, for which
            the report type denotes what type of sequence was aligned against
            what (BLASTN: dna-dna, BLASTP prt-prt, BLASTX translated dna-prt, 
            TBLASTN prt-translated dna, TBLASTX translated dna-translated dna).
 Example  : 
 Returns  : A string (BLASTN, BLASTP, BLASTX, TBLASTN, TBLASTX, UNKNOWN)
 Args     : a string on set (you should know what you are doing)

=cut

sub report_type {
    my ($self) = @_;
    return $self->{'_type'}; 
}


=head2 length

 Title   : length
 Usage   : my $len = $subject->length
 Function: Returns the length of the Subject (accessible from HSP as well)
 Returns : integer
 Args    : none

=cut

sub length{
   my ($self) = @_;
   
   return $self->{'_length'} || 0;
}

=head2 accession

 Title   : accession
 Usage   : my $accession = $subject->accession
 Function: Retrieve the accession for the subject (if any)
 Returns : string
 Args    : none

=cut

sub accession{
   my ($self,@args) = @_;
   return $self->{'_accession'} || '';
}

=head2 description

 Title   : description
 Usage   : my $desc = $subject->description
 Function: Retrieve the description for the subject (if any)
 Returns : string
 Args    : none

=cut

sub description{
   my ($self,@args) = @_;

   return $self->{'_description'} || '';
}

=head2 nextFeaturePair

 Title    : nextFeaturePair
 Usage    : $name = $obj->nextFeaturePair();
 Function : same as the nextHSP function 
 Example  : 
 Returns  : next FeaturePair 
 Args     :

=cut

sub nextFeaturePair { shift->next_hsp } # just another name

=head2 next_hsp

 Title    : next_hsp
 Usage    : $hsp = $obj->next_hsp();
 Function : returns the next available High Scoring Pair
 Example  : 
 Returns  : Bio::Search::HSPI  or null if finished
 Args     :

=cut

sub next_hsp {
    my ($self) = @_;
    my $index = $self->_nexthspindex;
    return undef if ($index > @{$self->{'_hsps'}} ); 
    return $self->{'_hsps'}->[$index];
}

=head2 Bio::Search::Subject methods

=cut

=head2 add_hsp

 Title   : add_hsp
 Usage   : $subject->add_hsp($hsp)
 Function: Adds a HSPI to the stored list
 Returns : number of HSPs
 Args    : Bio::Search::HSPI

=cut

sub add_hsp{
    my ($self,$h) = @_;
    if( $h->isa('Bio::Search::HSPI') ) { 
	push @{$self->{'_hsps'}}, $h;
    } else {
	$self->warn("Passed in " .ref($h). " as a HSP which is not a Bio::Search::HSPI... skipping");
    }
    return scalar @{$self->{'_hsps'}};
}

=head2 rewind

 Title   : rewind
 Usage   : $hsp->rewind;
 Function: Allow one to reset the HSP iteration to the beginning
           Since this is an in-memory implementation
 Returns : none
 Args    : none

=cut

sub rewind{
   my ($self) = @_;
   $self->{'_hspindex'} = 0;
}


=head2 _nexthspindex

 Title   : _nexthspindex
 Usage   : private

=cut

sub _nexthspindex{
   my ($self,@args) = @_;
   return $self->{'_hspindex'}++;
}

1;
