# $Id: bioperl.lisp 15559 2009-02-23 12:11:20Z maj $
#
# BioPerl module for Bio::Search::Result::hmmer3Result
#
# Please direct questions and support issues to <bioperl-l@bioperl.org>
#
# Cared for by Thomas Sharpton <thomas.sharpton@gmail.com>
#
# Copyright Thomas Sharpton
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Search::Result::hmmer3Result - DESCRIPTION of Object

=head1 SYNOPSIS

Give standard usage here

=head1 DESCRIPTION

Describe the object here

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

  bioperl-l@bioperl.org                  - General discussion
  http://bioperl.org/wiki/Mailing_lists  - About the mailing lists

=head2 Support

Please direct usage questions or support issues to the mailing list:

L<bioperl-l@bioperl.org>

rather than to the module maintainer directly. Many experienced and
reponsive experts will be able look at the problem and quickly
address it. Please include a thorough description of the problem
with code and data examples if at all possible.

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution. Bug reports can be submitted via
the web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Thomas Sharpton

Email thomas.sharpton@gmail.com

Describe contact details here

=head1 CONTRIBUTORS

Additional contributors names and emails here

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Search::Result::hmmer3Result;
use strict;

use base qw(Bio::Search::Result::GenericResult);

=head2 new

 Title   : new
 Usage   : my $obj = new Bio::Search::Result::hmmer3Result.pm();
 Function: Builds a new Bio::Search::Result::hmmer3Result.pm object 
 Returns : an instance of Bio::Search::Result::hmmer3Result.pm
 Args    : -hmm_name => string, name of hmm file
           -sequence_file => name of the sequence file

=cut

sub new {
  my($class,@args) = @_;
  my $self = $class->SUPER::new(@args);

  my ($hmm,$seqfile) = $self->_rearrange([qw(HMM_NAME SEQUENCE_FILE)],
					 @args);
  defined( $seqfile ) && $self->sequence_file( $seqfile );
  defined( $hmm ) && $self->hmm_name( $hmm );

  return $self;
}

=head2 hmm_name

 Title   : hmm_name
 Usage   : $obj->hmm_name($newval)
 Function: Get/Set the value of hmm_name
 Returns : value of hmm_name
 Args    : newvalue (optional)


=cut

sub hmm_name{
   my ($self,$value) = @_;
   if( defined $value) {
      $self->{'_hmm_name'} = $value;
    }
    return $self->{'_hmm_name'};
}

=head2 sequence_file

 Title   : sequence_file
 Usage   : $obj->sequence_file($newval)
 Function: Get/Set the value of sequence_file
 Returns : value of sequence_file
 Args    : newvalue (optional)


=cut

sub sequence_file{
   my ($self,$value) = @_;
   if( defined $value) {
      $self->{'_sequence_file'} = $value;
    }
    return $self->{'_sequence_file'};

}

=head2 next_model

 Title   : next_model
 Usage   : my $domain = $result->next_model
 Function: Returns the next domain - this
           is an alias for next_hit
 Returns : L<Bio::Search::Hit::HitI> object
 Args    : none


=cut

sub next_model{ shift->next_hit }

=head2 models

 Title   : models
 Usage   : my @domains = $result->models;
 Function: Returns the list of HMM models seen - this
           is an alias for hits()
 Returns : Array of L<Bio::Search::Hit::HitI> objects
 Args    : none


=cut

sub models{ shift->hits } 

=head2 rewind

 Title   : rewind
 Usage   : $result->rewind;
 Function: Allow one to reset the Hit iteration to the beginning
           Since this is an in-memory implementation
 Returns : none
 Args    : none

=cut

sub rewind{
   my ($self) = @_;
   $self->{'_hitindex'} = 0;
}



1;
