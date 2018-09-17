#
# BioPerl module for Bio::Search::Hit::GenericHit
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Jason Stajich <jason@bioperl.org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Search::Hit::BlastHit - Blast-specific subclass of Bio::Search::Hit::GenericHit

=head1 SYNOPSIS

    use Bio::Search::Hit::BlastHit;
    my $hit = Bio::Search::Hit::BlastHit->new(-algorithm => 'blastp');

# See Bio::Search::Hit::GenericHit for information about working with Hits.

# TODO: Describe how to configure a SearchIO stream so that it generates
#       GenericHit objects.

=head1 DESCRIPTION

This object is a subclass of Bio::Search::Hit::GenericHit
and provides some operations that facilitate working with BLAST
and PSI-BLAST Hits.

For general information about working with Hits, see 
Bio::Search::Hit::GenericHit.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

  bioperl-l@bioperl.org                  - General discussion
  http://bioperl.org/wiki/Mailing_lists  - About the mailing lists

=head2 Support 

Please direct usage questions or support issues to the mailing list:

I<bioperl-l@bioperl.org>

rather than to the module maintainer directly. Many experienced and 
reponsive experts will be able look at the problem and quickly 
address it. Please include a thorough description of the problem 
with code and data examples if at all possible.

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution. Bug reports can be submitted via the
web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Jason Stajich and Steve Chervitz

Email jason@bioperl.org
Email sac@bioperl.org

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Search::Hit::BlastHit;
use strict;

use Bio::Search::SearchUtils;

use base qw(Bio::Search::Hit::GenericHit);

=head2 new

 Title   : new
 Usage   : my $obj = Bio::Search::Hit::GenericHit->new();
 Function: Builds a new Bio::Search::Hit::GenericHit object 
 Returns : Bio::Search::Hit::GenericHit
 Args    : See Bio::Search::Hit::GenericHit() for other args.
           Here are the BLAST-specific args that can be used when
           creating BlastHit objects:
           -iteration    => integer for the PSI-Blast iteration number
           -found_again  => boolean, true if hit appears in a 
                            "previously found" section of a PSI-Blast report.

=cut

sub new {
  my($class,@args) = @_;

  my $self = $class->SUPER::new(@args);
  my ($iter,$found) = $self->_rearrange([qw(ITERATION 
                                            FOUND_AGAIN
                                           )], @args);

  defined $iter   && $self->iteration($iter);
  defined $found  && $self->found_again($found);

  return $self;
}

=head2 iteration

 Usage     : $hit->iteration( $iteration_num );
 Purpose   : Gets the iteration number in which the Hit was found.
 Example   : $iteration_num = $sbjct->iteration();
 Returns   : Integer greater than or equal to 1
             Non-PSI-BLAST reports will report iteration as 1, but this number
             is only meaningful for PSI-BLAST reports.
 Argument  : iteration_num (optional, used when setting only)
 Throws    : none

See Also   : L<found_again()|found_again>

=cut

sub iteration{
   my ($self,$value) = @_;
   if( defined $value) {
      $self->{'_psiblast_iteration'} = $value;
    }
    return $self->{'_psiblast_iteration'};
}

=head2 found_again

 Title     : found_again
 Usage     : $hit->found_again;
             $hit->found_again(1);
 Purpose   : Gets a boolean indicator whether or not the hit has
             been found in a previous iteration.
             This is only applicable to PSI-BLAST reports.

              This method indicates if the hit was reported in the 
              "Sequences used in model and found again" section of the
              PSI-BLAST report or if it was reported in the
              "Sequences not found previously or not previously below threshold"
              section of the PSI-BLAST report. Only for hits in iteration > 1.

 Example   : if( $hit->found_again()) { ... };
 Returns   : Boolean, true (1) if the hit has been found in a 
             previous PSI-BLAST iteration.
             Returns false (0 or undef) for hits that have not occurred in a
             previous PSI-BLAST iteration.
 Argument  : Boolean (1 or 0). Only used for setting.
 Throws    : none

See Also   : L<iteration()|iteration>

=cut

sub found_again {
   my $self = shift;
   return $self->{'_found_again'} = shift if @_;
   return $self->{'_found_again'};
}


sub expect { shift->significance(@_) }


1;
