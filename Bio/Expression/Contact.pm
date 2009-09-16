# $Id$
#
# BioPerl module for Bio::Expression::Contact
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Allen Day <allenday@ucla.edu>
#
# Copyright Allen Day
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Expression::Contact - DESCRIPTION of Object

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

I<bioperl-l@bioperl.org>

rather than to the module maintainer directly. Many experienced and 
reponsive experts will be able look at the problem and quickly 
address it. Please include a thorough description of the problem 
with code and data examples if at all possible.

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution. Bug reports can be submitted via
the web:

  http://bugzilla.open-bio.org/

=head1 AUTHOR - Allen Day

Email allenday@ucla.edu

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Expression::Contact;
use strict;
use base qw(Bio::Root::Root);

=head2 new()

 Usage   : my $obj = Bio::Expression::Contact->new();
 Function: Builds a new Bio::Expression::Contact object 
 Returns : an instance of Bio::Expression::Contact
 Args    :


=cut

sub new {
  my($class,@args) = @_;

  my $self = $class->SUPER::new(@args);
  $self->_initialize(@args);
  return $self;
}

=head2 _initialize()

 Usage   : $obj->_initialize(%arg);
 Function: Internal method to initialize a new Bio::Expression::Contact object
 Returns : true on success
 Args    : passed through to new()

=cut

sub _initialize {
  my($self,%arg) = @_;

  foreach my $arg (keys %arg){
    my $marg = $arg;
    $marg =~ s/^-//;
    $self->$marg($arg{$arg}) if $self->can($marg);
  }

  return 1;
}

=head2 source()

 Usage   : $obj->source($newval)
 Function: 
 Example : 
 Returns : value of source (a scalar)
 Args    : on set, new value (a scalar or undef, optional)


=cut

sub source {
  my($self,$val) = @_;
  $self->{'source'} = $val if defined($val);
  return $self->{'source'};
}

=head2 accession()

 Usage   : $obj->accession($newval)
 Function: 
 Example : 
 Returns : value of accession (a scalar)
 Args    : on set, new value (a scalar or undef, optional)


=cut

sub accession {
  my($self,$val) = @_;
  $self->{'accession'} = $val if defined($val);
  return $self->{'accession'};
}

=head2 name()

 Usage   : $obj->name($newval)
 Function: 
 Example : 
 Returns : value of name (a scalar)
 Args    : on set, new value (a scalar or undef, optional)


=cut

sub name {
  my($self,$val) = @_;
  $self->{'name'} = $val if defined($val);
  return $self->{'name'};
}

=head2 db()

 Usage   : $obj->db($newval)
 Function: 
 Example : 
 Returns : value of db (a scalar)
 Args    : on set, new value (a scalar or undef, optional)


=cut

sub db {
  my($self,$val) = @_;
  $self->{'db'} = $val if defined($val);
  return $self->{'db'};
}


1;
