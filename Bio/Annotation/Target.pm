#
# BioPerl module for Bio::Annotation::Target
#
# Cared for by Scott Cain <cain@cshl.org>
#
# Copyright Scott Cain
#
# Based on the Bio::Annotation::DBLink by Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Annotation::Target 

=head1 SYNOPSIS

   $target1 = new Bio::Annotation::Target(-target_id  => 'F321966.1',
                                          -tstart     => 1,
                                          -tend       => 200,
                                          -tstrand    => 1,    # or -1
					 );

   #or 

   $target2 = new Bio::Annotation::Target();
   $target2->target_id('Q75IM5');
   $target2->tstart(7);
   ... etc ...

   # Target is-a Bio::AnnotationI object, can be added to annotation
   # collections, e.g. the one on features or seqs
   $feat->annotation->add_Annotation('Target', $target2);


=head1 DESCRIPTION

Provides an object which represents a target (ie, a similarity hit) from
one object to something in another database without prescribing what is
in the other database

=head1 AUTHOR - Scott Cain

Scott Cain - cain@cshl.org

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...

package Bio::Annotation::Target;
use vars qw(@ISA);
use strict;
# use overload '""' => \&as_text; 

use Bio::Root::Root;
use Bio::AnnotationI;
use Bio::Location::Simple;

@ISA = qw(Bio::Root::Root Bio::AnnotationI);


sub new {
  my($class,@args) = @_;

  my $self = $class->SUPER::new(@args);

  my ($target_id, $tstart, $tend, $tstrand) =
      $self->_rearrange([qw(
                             TARGET_ID,
                             TSTART,
                             TEND,
                             TSTRAND,
			    )], @args);
  
  $target_id    && $self->target_id($target_id);
  $tstart       && $self->tstart($tstart);
  $tend         && $self->tend($tend);
  $tstrand      && $self->tstrand($tstrand);

  my $location  = Bio::Location::Simple->new(
                     -seq_id       => $target_id,
                     -start    => $tstart,
                     -end      => $tend,
                     -is_remote=> 1,
                  );
  $location->strand($tstrand) if $tstrand;
  $self->location($location);

  return $self;
}

=head1 AnnotationI implementing functions

=cut


=head2 as_text

 Title   : as_text
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub as_text{
   my ($self) = @_;

   return "Target=".$self->target_id." ".$self->tstart." ".
                    $self->tend." ".$self->strand;
}

=head2 tagname

 Title   : tagname
 Usage   : $obj->tagname($newval)
 Function: Get/set the tagname for this annotation value.

           Setting this is optional. If set, it obviates the need to
           provide a tag to Bio::AnnotationCollectionI when adding
           this object. When obtaining an AnnotationI object from the
           collection, the collection will set the value to the tag
           under which it was stored unless the object has a tag
           stored already.

 Example : 
 Returns : value of tagname (a scalar)
 Args    : new value (a scalar, optional)


=cut

sub tagname{
    my ($self,$value) = @_;
    if( defined $value) {
	$self->{'tagname'} = $value;
    }
    return $self->{'tagname'};
}

=head1 Specific accessors for Targets

=cut

=head2 target_id

=over

=item Usage

  $obj->target_id()        #get existing value
  $obj->target_id($newval) #set new value

=item Function

=item Returns

value of target_id (a scalar)

=item Arguments

new value of target_id (to set)

=back

=cut

sub target_id {
    my $self = shift;
    return $self->{'target_id'} = shift if defined(@_);
    return $self->{'target_id'};
}

=head2 tstart

=over

=item Usage

  $obj->tstart()        #get existing value
  $obj->tstart($newval) #set new value

=item Function

=item Returns

value of tstart (a scalar)

=item Arguments

new value of tstart (to set)

=back

=cut

sub tstart {
    my $self = shift;
    return $self->{'tstart'} = shift if defined(@_);
    return $self->{'tstart'};
}

=head2 tend

=over

=item Usage

  $obj->tend()        #get existing value
  $obj->tend($newval) #set new value

=item Function

=item Returns

value of tend (a scalar)

=item Arguments

new value of tend (to set)

=back

=cut

sub tend {
    my $self = shift;
    return $self->{'tend'} = shift if defined(@_);
    return $self->{'tend'};
}

=head2 tstrand

=over

=item Usage

  $obj->tstrand()        #get existing value
  $obj->tstrand($newval) #set new value

=item Function

=item Returns

value of tstrand (a scalar)

=item Arguments

new value of tstrand (to set)

=back

=cut

sub tstrand {
    my $self = shift;
    return $self->{'tstrand'} = shift if defined(@_);
    return $self->{'tstrand'};
}

=head2 location

=over

=item Usage

  $obj->location()        #get existing value
  $obj->location($newval) #set new value

=item Function

=item Returns

value of location (a scalar)

=item Arguments

new value of location (to set)

=back

=cut

sub location {
    my $self = shift;
    return $self->{'location'} = shift if defined(@_);
    return $self->{'location'};
}


1;
