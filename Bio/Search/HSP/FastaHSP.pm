# $Id$
#
# BioPerl module for Bio::Search::HSP::FastaHSP
#
# Cared for by Jason Stajich <jason@bioperl.org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Search::HSP::FastaHSP - HSP object for FASTA specific data

=head1 SYNOPSIS

# get a FastaHSP from a SearchIO stream
my $in = new Bio::SearchIO(-format => 'fasta', -file => 'filename.fasta');

while( my $r = $in->next_result) {
    while( my $hit = $r->next_result ) {
         while( my $hsp = $hit->next_hsp ) {
            print "smith-waterman score (if available): ", 
                  $hsp->sw_score(),"\n";
         }
    }
}

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


package Bio::Search::HSP::FastaHSP;
use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::Root::Root

use Bio::Search::HSP::GenericHSP;

@ISA = qw(Bio::Search::HSP::GenericHSP );

=head2 new

 Title   : new
 Usage   : my $obj = new Bio::Search::HSP::FastaHSP();
 Function: Builds a new Bio::Search::HSP::FastaHSP object 
 Returns : Bio::Search::HSP::FastaHSP
 Args    :


=cut

sub new {
  my($class,@args) = @_;

  my $self = $class->SUPER::new(@args);
  
  my ($swscore) = $self->_rearrange([qw(SWSCORE)], @args);

  defined $swscore && $self->sw_score($swscore);

  return $self;
}


=head2 sw_score

 Title   : sw_score
 Usage   : $obj->sw_score($newval)
 Function: Get/Set Smith-Waterman score
 Returns : value of sw_score
 Args    : newvalue (optional)


=cut

sub sw_score{
    my ($self,$value) = @_;
    if( defined $value || ! defined $self->{'_sw_score'} ) {
	$value = 0 unless defined $value; # default value
	$self->{'_sw_score'} = $value;
    }
    return $self->{'_sw_score'};
}


1;
