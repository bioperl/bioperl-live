#
# BioPerl module for Bio::Seq::SeqFactory
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

Bio::Seq::SeqFactory - Instantiation of generic Bio::PrimarySeqI (or derived) objects through a factory

=head1 SYNOPSIS

    use Bio::Seq::SeqFactory;
    my $factory = Bio::Seq::SeqFactory->new();
    my $primaryseq = $factory->create( -seq => 'WYRAVLC',
                                       -id  => 'name'     );

    # Create Bio::Seq instead of Bio::PrimarySeq objects:
    my $factory = Bio::Seq::SeqFactory->new( -type => 'Bio::Seq' );


=head1 DESCRIPTION

This object will build L<Bio::PrimarySeqI> and L<Bio::SeqI> objects
generically.

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

=head1 AUTHOR - Jason Stajich

Email jason@bioperl.org

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Seq::SeqFactory;
use strict;


use base qw(Bio::Root::Root Bio::Factory::SequenceFactoryI);

=head2 new

 Title   : new
 Usage   : my $obj = Bio::Seq::SeqFactory->new();
 Function: Builds a new Bio::Seq::SeqFactory object 
 Returns : Bio::Seq::SeqFactory
 Args    : -type => string, name of a PrimarySeqI derived class
                    This is optional. Default=Bio::PrimarySeq.

=cut

sub new {
  my($class,@args) = @_;
  my $self = $class->SUPER::new(@args);
  my ($type) = $self->_rearrange([qw(TYPE)], @args);
  if( ! defined $type ) { 
      $type = 'Bio::PrimarySeq';
  }
  $self->type($type);
  return $self;
}


=head2 create

 Title   : create
 Usage   : my $seq = $seqbuilder->create(-seq => 'CAGT', -id => 'name');
 Function: Instantiates new Bio::SeqI (or one of its child classes)
           This object allows us to genericize the instantiation of sequence
           objects.
 Returns : Bio::PrimarySeq object (default)
           The return type is configurable using new(-type =>"...").
 Args    : initialization parameters specific to the type of sequence
           object we want.  Typically 
           -seq        => $str,
           -display_id => $name

=cut

sub create {
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

sub type {
   my ($self, $value) = @_;
   if (defined $value) {
       eval "require $value";
       if( $@ ) { $self->throw("$@: Unrecognized Sequence type for SeqFactory '$value'");}
       
       my $a = bless {},$value;
       unless( $a->isa('Bio::PrimarySeqI') ||
               $a->isa('Bio::Seq::QualI' ) ) {
           $self->throw("Must provide a valid Bio::PrimarySeqI or Bio::Seq::QualI or child class to SeqFactory Not $value");
       }
       $self->{'type'} = $value;
    }
    return $self->{'type'};
}

1;
