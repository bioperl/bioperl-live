#
# BioPerl module for Bio::Annotation::Comment
#
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Annotation::Comment - A Comment on an Annotation

=head1 SYNOPSIS

    # comment objects attached to annotations
    foreach my $comment ( $seq->annotation->each_Comment() ) {
	# comment object currently pretty stupid. Just gives back 
	# text as a string
	$text = $comment->text();
    }

=head1 DESCRIPTION

A comment object is meant to represent one logical comment in a 
piece of annotation (common CC lines in EMBL files etc). At the moment
is a very simple object, but this will give us a placeholder for 
more advanced things, eg, able to provide some XML stuff and eventually
things like authorship tracking.

This object originally came from the Pfam annotation object

=head1 CONTACT

Ewan Birney <birney@ebi.ac.uk>

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::Annotation::Comment;
use vars qw(@ISA);
use strict;

# Object preamble - inheriets from Bio::Root::Object

use Bio::Root::Object;

@ISA = qw(Bio::Root::Object);
# new() is inherited from Bio::Root::Object

# _initialize is where the heavy stuff will happen when new is called

sub _initialize {
  my($self,@args) = @_;

  my $make = $self->SUPER::_initialize;

  $self->{'flat'} = [];
# set stuff in self from @args
 return $make; # success - we hope!
}

=head2 text

 Title   : text
 Usage   : $self->text($newval)
 Function: 
 Example : 
 Returns : value of text
 Args    : newvalue (optional)

=cut

sub text{
   my ($self,$value) = @_;
   if( defined $value) {
      $self->{'text'} = $value;
    }
    return $self->{'text'};

}

1;
