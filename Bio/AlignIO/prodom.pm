#
# BioPerl module for Bio::AlignIO::prodom

#	based on the Bio::SeqIO::prodom module
#       by Ewan Birney <birney@sanger.ac.uk>
#       and Lincoln Stein  <lstein@cshl.org>
#
#       and the SimpleAlign.pm module of Ewan Birney
#
# Copyright Peter Schattner
#
# You may distribute this module under the same terms as perl itself
# _history
# September 5, 2000
# POD documentation - main docs before the code

=head1 NAME

Bio::AlignIO::prodom - prodom sequence input/output stream

=head1 SYNOPSIS

Do not use this module directly.  Use it via the Bio::AlignIO class.

=head1 DESCRIPTION

This object can transform Bio::SimpleAlign objects to and from prodom flat
file databases.

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

=head1 AUTHORS - Peter Schattner

Email: schattner@alum.mit.edu


=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::AlignIO::prodom;
use vars qw(@ISA);
use strict;
# Object preamble - inherits from Bio::Root::Object

use Bio::AlignIO;

@ISA = qw(Bio::AlignIO);
# new() is inherited from Bio::Root::Object

# _initialize is where the heavy stuff will happen when new is called

sub _initialize {
  my($self,@args) = @_;
  return unless my $make = $self->SUPER::_initialize(@args);
}

=head2 next_aln

 Title   : next_aln
 Usage   : $aln = $stream->next_aln()
 Function: returns the next alignment in the stream.
 Returns : SimpleAlign object
 Args    : NONE

=cut

sub next_aln {
    my $self = shift;
    my $entry;
    my ($acc, $fake_id, $start, $end, $seq, $add, %names);



    my $aln =  Bio::SimpleAlign->new();


    while( $entry = $self->_readline) {

       if ($entry =~ /^AC\s+(\S+)\s*$/) {         #ps 9/12/00
	   $aln->id( $1 );
       }
       elsif ($entry =~ /^AL\s+(\S+)\|(\S+)\s+(\d+)\s+(\d+)\s+\S+\s+(\S+)\s*$/){    #ps 9/12/00
	   $acc=$1;
	   $fake_id=$2;  # Accessions have _species appended
	   $start=$3;
	   $end=$4;
	   $seq=$5;
	
	   $names{'fake_id'} = $fake_id;

	   $add = new Bio::LocatableSeq('-seq'=>$seq,
			       '-id'=>$acc,
			       '-start'=>$start,
			       '-end'=>$end,
			       );
	
	   $aln->addSeq($add);
       }
       elsif ($entry =~ /^CO/) {
	   # the consensus line marks the end of the alignment part of the entry
	   last;
       }
   }

#  If $end <= 0, we have either reached the end of
#  file in <> or we have encountered some other error
#
   if ($end <= 0) { undef $aln;}


   return $aln;
}



=head2 write_aln

 Title   : write_aln
 Usage   : $stream->write_aln(@aln)
 Function: writes the $aln object into the stream in prodom format  ###Not yet implemented!###
 Returns : 1 for success and 0 for error
 Args    : Bio::SimpleAlign object


=cut

sub write_aln {
    my ($self,@aln) = @_;

    $self->throw("Sorry: prodom-format output, not yet implemented! /n");
}

