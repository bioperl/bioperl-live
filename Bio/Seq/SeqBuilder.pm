# $Id$
#
# BioPerl module for Bio::Seq::SeqBuilder
#
# Cared for by Jason Stajich <jason@bioperl.org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Seq::SeqBuilder - Instantiates a new Bio::PrimarySeqI (or derived class) through a factory

=head1 SYNOPSIS

    use Bio::Seq::SeqBuilder;
    my $factory = new Bio::Seq::SeqBuilder(-type => 'Bio::PrimarySeq');
    my $seq = $factory->new_sequence(-seq => 'WYRAVLC',
				     -id  => 'name');

=head1 DESCRIPTION

This object will build Bio::Seq objects generically.

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


package Bio::Seq::SeqBuilder;
use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::Root::Root

use Bio::Root::Root;
use Bio::Factory::SequenceBuilderI;

@ISA = qw(Bio::Root::Root Bio::Factory::SequenceBuilderI);

=head2 new

 Title   : new
 Usage   : my $obj = new Bio::Seq::SeqBuilder();
 Function: Builds a new Bio::Seq::SeqBuilder object 
 Returns : Bio::Seq::SeqBuilder
 Args    :


=cut

sub new {
  my($class,@args) = @_;
  my $self = $class->SUPER::new(@args);
  my ($type) = $self->_rearrange([qw(TYPE)], @args);
  if( ! defined $type ) { 
      $type = 'Bio::PrimarySeq';
  }
  eval "require $type";
  if( $@ ) { $self->throw("$@: Unrecognized Sequence type for SeqBuilder '$type'");}
  $self->type($type);
  return $self;
}


=head2 new_sequence

 Title   : new_sequence
 Usage   : my $seq = $seqbuilder->new_sequence(-seq => 'CAGT', -id => 'name');
 Function: Instantiates new Bio::SeqI (or one of its child classes)
           This object allows us to genericize the instantiation of sequence
           objects.
 Returns : Bio::SeqI
 Args    : initialization parameters specific to the type of sequence
           object we want.  Typically 
           -seq        => $str,
           -display_id => $name

=cut

sub new_sequence{
   my ($self,@args) = @_;
   return $self->type->new(-verbose => $self->verbose, @args);
}

=head2 type

 Title   : type
 Usage   : $obj->type($newval)
 Function: 
 Returns : value of type
 Args    : newvalue (optional)


=cut

sub type{
   my ($self,$value) = @_;
   if( defined $value) {
      $self->{'type'} = $value;
    }
    return $self->{'type'};
}

1;
