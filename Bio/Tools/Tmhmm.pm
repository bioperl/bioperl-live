#
# BioPerl module for Bio::Tools::Tmhmm
#
# Original copyright Balamurugan Kumarasamy
# Re-written cleanly by Torsten Seemann, Sep 2006
#
# Copyright:
# You may distribute this module under the same terms as Perl itself

=head1 NAME

Bio::Tools::Tmhmm - parse TMHMM output (TransMembrane HMM)

=head1 SYNOPSIS

  use Bio::Tools::Tmhmm;
  my $parser = Bio::Tools::Tmhmm->new(-fh => $filehandle );
  while ( my $tmhmm_feat = $parser->next_result ) {
     # do something, e.g.
     push @tmhmm_feat, $tmhmm_feat;
  }

=head1 DESCRIPTION

TMHMM is software for the prediction of transmembrane helices in proteins.
See  L<http://www.cbs.dtu.dk/services/TMHMM/> for more details.

This module parses the "long output" format of TMHMM 2.0 and
creates a Bio:SeqFeature::Generic object for each C<TMHelix> feature found
from lines like this:

  my_sequence_id  TMHMM2.0  TMhelix     54    76


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

=head1 AUTHOR - Torsten Seemann

Email torsten.seemann AT infotech.monash.edu.au

=head1 CONTRIBUTOR - Bala

Email savikalpa@fugu-sg.org

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

package Bio::Tools::Tmhmm;

use strict;

use Bio::Tools::AnalysisResult;
use Bio::Root::Root;
use Bio::Root::IO;

use base qw(Bio::Root::Root Bio::Root::IO Bio::Tools::AnalysisResult);

use Bio::SeqFeature::Generic;


=head2 new

 Title   : new
 Usage   : my $obj = Bio::Tools::Tmhmm->new();
 Function: Builds a new Bio::Tools::Tmhmm object
 Returns : Bio::Tools::Tmhmm
 Args    : Either of the following as per L<Bio::Root::IO> interface
             -fh   => $filehandle 
             -file => $filename

=cut

sub new {
  my($class,@args) = @_;
  my $self = $class->SUPER::new(@args);
  $self->_initialize_io(@args);
  return $self;
}


=head2 next_result

 Title   : next_result
 Usage   : my $feat = $Tmhmm->next_result
 Function: Get the next result set from parser data
 Returns : Bio::SeqFeature::Generic
 Args    : none

=cut

sub next_result {
  my $self = shift;

  # # my_sequence_id Length: 178
  # my_sequence_id  TMHMM2.0  outside      1    53
  # my_sequence_id  TMHMM2.0  TMhelix     54    76
  # my_sequence_id  TMHMM2.0  inside      77   115

  while (my $line = $self->_readline) { 
    if ( $line =~ m/^(\S+)\s+(\S+)\s+(TMhelix)\s+(\d+)\s+(\d+)$/i ) {
       return Bio::SeqFeature::Generic->new(
	 -primary => 'transmembrane',
         -seq_id  => $1,
	 -source  => $2,
	 -start   => $4,
	 -end     => $5,
       );
    }
  }
}

1;


