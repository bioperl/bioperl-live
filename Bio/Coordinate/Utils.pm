# $Id$
#
# BioPerl module for Bio::Coordinate::Utils
#
# Cared for by Heikki Lehvaslaiho <heikki@ebi.ac.uk>
#
# Copyright Heikki Lehvaslaiho
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Coordinate::Utils - Additional methods to create Bio::Coordinate objects

=head1 SYNOPSIS

    use Bio::Coordinate::Utils;
    # get a Bio::Align::AlignI compliant object, $aln, somehow
    # it could be a Bio::SimpleAlign

    $mapper = Bio::Coordinate::Utils->from_align($aln, 1);

=head1 DESCRIPTION

This class is a holder of methods that work on or create
Bio::Coordinate::MapperI- compliant objects. . These methods are not
part of the Bio::Coordinate::MapperI interface and should in general
not be essential to the primary function of sequence objects. If you
are thinking of adding essential functions, it might be better to
create your own sequence class.  See L<Bio::PrimarySeqI>,
L<Bio::PrimarySeq>, and L<Bio::Seq> for more.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists.  Your participation is much appreciated.

  bioperl-l@bioperl.org                        - General discussion
  http://bio.perl.org/MailList.html             - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via email
or the web:

  bioperl-bugs@bio.perl.org
  http://bugzilla.bioperl.org/

=head1 AUTHOR - Heikki Lehvaslaiho

Email:  heikki@ebi.ac.uk
Address:

     EMBL Outstation, European Bioinformatics Institute
     Wellcome Trust Genome Campus, Hinxton
     Cambs. CB10 1SD, United Kingdom

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Coordinate::Utils;
use vars qw(@ISA);

use Bio::Location::Simple;
use Bio::Coordinate::Pair;
use Bio::Coordinate::Collection;

use strict;

@ISA = qw(Bio::Root::Root);
# new inherited from Root

=head2 from_align

 Title   : from_align
 Usage   : $mapper = Bio::Coordinate::Utils->from_align($aln, 1);
 Function:
           Create a mapper out of an alignment.
           The mapper will return a value only when both ends of
           the input range find a match.

           Note: This implementation works only on pairwise alignments
           and is not yet well tested!

 Returns : A Bio::Coordinate::MapperI
 Args    : Bio::Align::AlignI object
           Id for the reference sequence, optional

=cut

sub from_align {
   my ($self, $aln, $ref ) = @_;

   $aln->isa('Bio::Align::AlignI') ||
       $self->throw('Not a Bio::Align::AlignI object but ['. ref($self). ']');

   # default reference sequence to the first sequence
   $ref ||= 1;
   
   my $collection = Bio::Coordinate::Collection->new(-return_match=>1);

   # this works only for pairs, so split the MSA
   # take the ref
   #foreach remaining seq in aln, do:
   $aln->map_chars('\.','-');
   my $cs = $aln->gap_line;
   my $seq1 = $aln->get_seq_by_pos(1);
   my $seq2 = $aln->get_seq_by_pos(2);   
   while ( $cs =~ /([^\-]+)/g) {
       # alignment coordinates
       my $lenmatch = length($1);
       my $start = pos($cs) - $lenmatch +1;
       my $end   = $start + $lenmatch -1;
       my $match1 = Bio::Location::Simple->new
	   (-seq_id => $seq1->id,
	    -start  => $seq1->location_from_column($start)->start,
	    -end    => $seq1->location_from_column($end)->start,
	    -strand => $seq1->strand );

       my $match2 = Bio::Location::Simple->new
	   (-seq_id => $seq2->id,
	    -start  => $seq2->location_from_column($start)->start,
	    -end    => $seq2->location_from_column($end)->start,
	    -strand => $seq2->strand );
       
       my $pair = Bio::Coordinate::Pair->new
	   (-in  => $match1,
	    -out => $match2
	    );
       if(! $pair->test ) {
	   $self->warn("pair align did not pass test:\n",
		       "\tm1=",$match1->to_FTstring(), " len=",$match1->length, 
		       " m2=", $match2->to_FTstring()," len=",$match2->length,"\n");
       }
       $collection->add_mapper($pair);
   }
   return ($collection->each_mapper)[0] if $collection->mapper_count == 1;
   return $collection;

}



1;
