# $Id$
#
# BioPerl module for Bio::Search::SubjectI
#
# Cared for by Jason Stajich <jason@bioperl.org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Search::SubjectI - A Subject (Hit) from a Sequence Database Search

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


package Bio::Search::SubjectI;
use vars qw(@ISA);
use strict;
use Carp;

sub _abstractDeath {
  my $self = shift;
  my $package = ref $self;
  my $caller = (caller)[1];
  confess "Abstract method '$caller' defined in interface Bio::Search::SubjectI not implemented by pacakge $package. Not your fault - author of $package should be blamed!";
}

=head2 name

 Title    : name
 Usage    : $name = $obj->name();
 Function : returns the name of the Sbjct 
 Example  : 
 Returns  : name of the Sbjct 
 Args     :

=cut

sub name { $_[0]->_abstractDeath; }

=head2 length

 Title   : length
 Usage   : my $len = $subject->length
 Function: Returns the length of the Subject (accessible from HSP as well)
 Returns : integer
 Args    : none

=cut

sub length{
   my ($self,@args) = @_;
   $self->_abstractDeath;
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
   $self->_abstractDeath;
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
   $self->_abstractDeath;
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
    $_[0]->_abstractDeath; 
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
    $_[0]->_abstractDeath;
}

1;
