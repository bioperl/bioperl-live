# $Id$
#
# BioPerl module for Bio::Expression::Platform
#
# Cared for by Allen Day <allenday@ucla.edu>
#
# Copyright Allen Day
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Expression::Platform - DESCRIPTION of Object

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


package Bio::Expression::Platform;
use strict;
use base qw(Bio::Root::Root);
use Bio::DB::Taxonomy;

=head2 new()

 Usage   : my $obj = Bio::Expression::Platform->new();
 Function: Builds a new Bio::Expression::Platform object 
 Returns : an instance of Bio::Expression::Platform
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
 Function: Internal method to initialize a new Bio::Expression::Platform object
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

  $self->taxdb( Bio::DB::Taxonomy->new(-source => 'entrez') );
  return 1;
}



=head2 get_datasets()

 Usage   :
 Function:
 Example :
 Returns : 
 Args    :

=cut

sub get_datasets {
  my ($self,@args) = @_;
  my $db = $self->db();

  my @datasets = $db->get_datasets( $self );

  return @datasets;
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

=head2 taxon()

 Usage   : $obj->taxon()
 Function: 
 Example : 
 Returns : A Bio::Taxonomy::Node object
 Args    : none


=cut

sub taxon {
  my($self) = @_;
  if ( ! $self->{'taxon'} ) {
    $self->{'taxon'} = $self->taxdb->get_Taxonomy_Node( $self->_taxon_id() );
  }
  return $self->{'taxon'};
}

=head2 contact()

 Usage   : $obj->contact($newval)
 Function: 
 Example : 
 Returns : a Bio::Expression::Contact object
 Args    : on set, new value (a scalar or undef, optional)


=cut

sub contact {
  my($self,$val) = @_;
  $self->{'contact'} = $val if defined($val);
  return $self->{'contact'};
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

=head2 _taxon_id()

 Usage   : $obj->_taxon_id($newval)
 Function: 
 Example : 
 Returns : value of _taxon_id (a scalar)
 Args    : on set, new value (a scalar or undef, optional)


=cut

sub _taxon_id {
  my($self,$val) = @_;
  $self->{'_taxon_id'} = $val if defined($val);
  return $self->{'_taxon_id'};
}

=head2 taxdb()

 Usage   : $obj->taxdb($newval)
 Function: 
 Example : 
 Returns : a Bio::DB::Taxonomy object
 Args    : on set, new value (a scalar or undef, optional)


=cut

sub taxdb {
  my($self,$val) = @_;
  $self->{'taxdb'} = $val if defined($val);
  return $self->{'taxdb'};
}



1;
