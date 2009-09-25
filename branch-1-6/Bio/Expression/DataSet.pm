# $Id$
#
# BioPerl module for Bio::Expression::DataSet
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

Bio::Expression::DataSet - DESCRIPTION of Object

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


package Bio::Expression::DataSet;
use strict;
use base qw(Bio::Root::Root);

=head2 new()

 Usage   : my $obj = Bio::Expression::DataSet->new();
 Function: Builds a new Bio::Expression::DataSet object 
 Returns : an instance of Bio::Expression::DataSet
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
 Function: Internal method to initialize a new Bio::Expression::DataSet object
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

 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub name {
  my($self,$val) = @_;
  $self->_load();
  $self->{'db'} = $val if defined($val);
  return $self->{'db'};
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

=head2 pubmed_id()

 Usage   : $obj->pubmed_id($newval)
 Function: 
 Example : 
 Returns : value of pubmed_id (a scalar)
 Args    : on set, new value (a scalar or undef, optional)


=cut

sub pubmed_id {
  my($self,$val) = @_;
  $self->_load();
  $self->{'pubmed_id'} = $val if defined($val);
  return $self->{'pubmed_id'};
}

=head2 web_link()

 Usage   : $obj->web_link($newval)
 Function: 
 Example : 
 Returns : value of web_link (a scalar)
 Args    : on set, new value (a scalar or undef, optional)


=cut

sub web_link {
  my($self,$val) = @_;
  $self->_load();
  $self->{'web_link'} = $val if defined($val);
  return $self->{'web_link'};
}

=head2 contact()

 Usage   : $obj->contact($newval)
 Function: 
 Example : 
 Returns : value of contact (a scalar)
 Args    : on set, new value (a scalar or undef, optional)


=cut

sub contact {
  my($self,$val) = @_;
  $self->_load();
  $self->{'contact'} = $val if defined($val);
  return $self->{'contact'};
}

=head2 samples()

 Usage   : $obj->samples($newval)
 Function: 
 Example : 
 Returns : value of samples (a scalar)
 Args    : on set, new value (a scalar or undef, optional)


=cut

sub samples {
  my($self,$val) = @_;
  $self->_load();
  $self->{'samples'} = $val if defined($val);
  return $self->{'samples'};
}

=head2 description()

 Usage   : $obj->description($newval)
 Function: 
 Example : 
 Returns : value of description (a scalar)
 Args    : on set, new value (a scalar or undef, optional)


=cut

sub description {
  my($self,$val) = @_;
  $self->_load();
  $self->{'description'} = $val if defined($val);
  return $self->{'description'};
}

=head2 design()

 Usage   : $obj->design($newval)
 Function: 
 Example : 
 Returns : value of design (a scalar)
 Args    : on set, new value (a scalar or undef, optional)


=cut

sub design {
  my($self,$val) = @_;
  $self->_load();
  $self->{'design'} = $val if defined($val);
  return $self->{'design'};
}

=head2 design_description()

 Usage   : $obj->design_description($newval)
 Function: 
 Example : 
 Returns : value of design_description (a scalar)
 Args    : on set, new value (a scalar or undef, optional)


=cut

sub design_description {
  my($self,$val) = @_;
  $self->_load();
  $self->{'design_description'} = $val if defined($val);
  return $self->{'design_description'};
}













=head2 get_samples()

 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_samples {
  my ($self,@args) = @_;
  if ( $self->samples() ) {
    return @{ $self->samples() };
  }
  else {
    return ();
  }
}







sub _load {
  my $self = shift;
  if ( $self->{'_load'} ) {
    return 1;
  }
  $self->{'_load'}++;
  $self->db->fill_dataset( $self );
  return $self->{'_load'};
}

1;
