
#
# BioPerl module for Bio::Annotation::SimpleValue
#
# Cared for by bioperl <bioperl-l@bio.perl.org>
#
# Copyright bioperl
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Annotation::SimpleValue - A simple scalar 

=head1 SYNOPSIS

   use Bio::Annotation::SimpleValue;
   use Bio::Annotation::Collection;

   my $col = new Bio::Annotation::Collection;
   my $sv = new Bio::Annotation::SimpleValue(-value => 'someval');   
   $col->add_Annotation('tagname', $sv);

=head1 DESCRIPTION

Scalar value annotation object 

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists. Your participation is much appreciated.

  bioperl-l@bioperl.org              - General discussion
  http://bio.perl.org/MailList.html  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via email
or the web:

  bioperl-bugs@bioperl.org
  http://bio.perl.org/bioperl-bugs/

=head1 AUTHOR - bioperl

Email bioperl-l@bio.perl.org

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Annotation::SimpleValue;
use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::Root::Root

use Bio::AnnotationI;
use Bio::Root::Root;

@ISA = qw(Bio::Root::Root Bio::AnnotationI);

=head2 new

 Title   : new
 Usage   : my $sv = new Bio::Annotation::SimpleValue;
 Function: Instantiate a new SimpleValue object
 Returns : Bio::Annotation::SimpleValue object
 Args    : -value => $value to initialize the object data field [optional]

=cut

sub new{
   my ($class,@args) = @_;

   my $self = $class->SUPER::new(@args);

   my ($value) = $self->_rearrange([qw(VALUE)], @args);
   defined $value  && $self->value($value);
   return $self;
}


=head2 AnnotationI implementing functions

=cut

=head2 as_text

 Title   : as_text
 Usage   : my $text = $obj->as_text
 Function: return the string "Value: $v" where $v is the value 
 Returns : string
 Args    : none


=cut

sub as_text{
   my ($self) = @_;

   return "Value: ".$self->value;
}

=head2 hash_tree

 Title   : hash_tree
 Usage   : my $hashtree = $value->hash_tree
 Function: For supporting the AnnotationI interface just returns the value
           as a hashref with the key 'value' pointing to the value
 Returns : hashrf
 Args    : none


=cut

sub hash_tree{
   my ($self) = @_;
   
   my $h = {};
   $h->{'value'} = $self->value;
}

=head2 Specific accessors for SimpleValue

=cut

=head2 value

 Title   : value
 Usage   : $obj->value($newval)
 Function: Get/Set the value for simplevalue
 Returns : value of value
 Args    : newvalue (optional)


=cut

sub value{
   my ($self,$value) = @_;
   
   if( defined $value) {
      $self->{'value'} = $value;
    }
    return $self->{'value'};
}

1;
