#
# BioPerl module for Bio::SeqFeature::Primer
#
# Cared for by Chad Matsalla
#chad
# Copyright Chad Matsalla
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::SeqFeature::Primer - Primer Generic SeqFeature

=head1 SYNOPSIS

=head1 DESCRIPTION

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists.  Your participation is much appreciated.

  bioperl-l@bioperl.org          - General discussion
  http://bio.perl.org/MailList.html             - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via email
or the web:

  bioperl-bugs@bio.perl.org
  http://bio.perl.org/bioperl-bugs/

=head1 AUTHOR - Chad Matsalla

Chad Matsalla <lt>bioinformatics1@dieselwurks.comE<gt>

=head1 DEVELOPERS

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::SeqFeature::Primer;
use vars qw(@ISA);
use strict;

use Bio::Root::Root;
use Bio::SeqFeatureI;
use Bio::SeqFeature::Generic;
     # I am really not sure which of these I need so I included them at this time
use Bio::Annotation;
use Bio::Location::Simple;
use Bio::Seq;
use Dumpvalue qw(dumpValue);

my $dumper = new Dumpvalue();


@ISA = qw(Bio::Root::Root Bio::SeqFeature::Generic);



=head2 new()

 Title   : new()
 Usage   :
 Function:
 Example :
 Returns :
 Args    : 
Devel notes: I think that I want to accept a hash

=cut


sub new {
        my ( $class, @args) = @_;  
          my %arguments = @args;
        my $self = $class->SUPER::new(@args);
          # these are from generic.pm, with which i started
        $self->{'_parse_h'}       = {};
        $self->{'_gsf_tag_hash'}  = {}; 
               # things that belong with the primer
          my ($sequence, $primer_sequence_id, $id) = 
               $self->_rearrange([qw(SEQUENCE PRIMER_SEQUENCE_ID ID)],@args);
          $self->{'_gsf_seqname'} = $self->{primer_sequence_id};
          foreach my $argument (sort keys %arguments) {
                  if ($argument eq "-sequence") {
                         print("Processing sequence\n");
                          if (!%arguments->{-sequence}) {
                               # print("No sequence!!!!\n");
                          }
                          elsif (ref($sequence) eq "Bio::Seq") {
                              print("Seq isa reference\n");
                               $self->{seq} = %arguments->{-sequence};
                                    # i don't know what to do about the id() and display_id() thing
                               $self->{id} = $self->{seq}->display_id();
                          }
                          else {
                              print("Seq was not a reference. creating a new object.\n");
                               $self->{seq} = new Bio::Seq( -seq => %arguments->{-sequence} );
                          }
                    }
                    else {
                         $self->{$argument} = %arguments->{$argument};
                    }
          }
               # what about the ids?
          $self->{id} = $id ? $id : $primer_sequence_id;
          $self->{primer_sequence_id} = $primer_sequence_id ? $primer_sequence_id : $id;
          if (!$self->{seq}->display_id()) {
               $self->{seq}->display_id($self->{id});
          }
               # a bunch of things now need to be set for this SeqFeature
               # things like:
                    # TARGET=513,26
                    # PRIMER_FIRST_BASE_INDEX=1
                    # PRIMER_LEFT=484,20
          my $target = $self->{TARGET} ? $self->{TARGET} : $self->{target};
          if ($target ) {
                 $$target = m/(\d)\,(\d)/;
                  $self->start() = $1;
                  $self->end() = $1 + $2;
          }

          return $self;
}


=head2 seq()

 Title   : seq()
 Usage   : $seq = $primer->seq();
 Function: Return the _entire_ sequence associated with this Primer. 
 Returns : A Bio::Seq object
 Args    : None.

=cut

sub seq() {
     my $self = shift;
     return $self->{seq};
}


=head2 get_feature()

 Title   : get_feature()
 Usage   : $seq = $primer->get_feature("PENALTY");
 Function: When the Bio::SeqFeature::Primer object is created, the user can
     pass in an arbitrary hash containing key->value pairs. This is allowed
     because I didn't want to assume that the user was trying to model a
     primer3 construct. 
 Returns : A scalar returning the value of the provided key or null if that
     key was not defined in the hash used to construct this Bio::SeqFeature::
     Primer object.
 Args    : A string representing a hash key.

=cut

sub get_feature {
     my ($self,$key) = @_;
     return $self->{$key}; 



}




1;
