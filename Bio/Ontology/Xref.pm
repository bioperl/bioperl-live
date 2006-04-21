# $Id$
#
# This GO module is maintained by Chris Mungall <cjm@fruitfly.org>
#
# see also - http://www.geneontology.org
#          - http://www.godatabase.org/dev
#
# You may distribute this module under the same terms as perl itself

package Bio::Ontology::Xref;

=head1 NAME

  GO::Model::Xref;

=head1 SYNOPSIS

  my $xrefs = $term->dbxref_list();
  foreach my $xref (@$xrefs) P
    printf "Term %s has an xref %s:%s\n",
            $term->name, $xref->xref_key, $xref->dbname;
  }

=head1 DESCRIPTION

represents a cross reference to an external database. an Xref is made
up of a key (ie the accession number, or whatever the value of the
unique field being keyed off of is) and a database name. this should
theorerically be enough to uniquely identify any databased entity.

=head1 NOTES

Like all the GO::Model::* classes, this uses accessor methods to get
or set the attributes. by using the accessor method without any
arguments gets the value of the attribute. if you pass in an argument,
then the attribuet will be set according to that argument.

for instance

  # this sets the value of the attribute
  $my_object->attribute_name("my value");

  # this gets the value of the attribute
  $my_value = $my_object->attribute_name();

=cut


use Carp qw(cluck confess);
#use Exporter;
use GO::Utils qw(rearrange);
use GO::Model::Root;
use strict;
use vars qw(@ISA);

#@ISA = qw(GO::Model::Root Exporter);
@ISA = qw( Bio::Root::Root );


sub new {

    my( $class,@args ) = @_;

    my $self = $class->SUPER::new( @args );
    my ( $accession,
         $dbname,
         $description,
         $type
       )
        = $self->_rearrange( [ qw( ACCESSION
                                   DBNAME
                                   DESCRIPTION
                                   TYPE
       ) ], @args );

    defined($accession)   && $self->accession( $accession );
    defined($dbname)         && $self->dbname( $dbname );
    defined($description)   && $self->description( $description );
    defined($type)   && $self->description( $type );

    return $self;
} # new


=head2 accession

  Usage   -
  Returns -
  Args    -

 accessor: gets/sets the key/id of the cross reference

=cut

sub accession {
    my $self = shift;
    $self->{accession} = shift if @_;
    return $self->{accession};
}


=head2 type

  Usage   -
  Returns -
  Args    -

 accessor: gets/sets the key/id type of the cross reference


=cut

sub type {
    my $self = shift;
    $self->{type} = shift if @_;
    return $self->{type};
}


=head2 formatted_xref

  Usage   -
  Returns -
  Args    -

=cut

sub formatted_xref_ {
    my $self=shift;
    return $self->dbname().":".$self->accession();
}


=head2 dbname

  Alias   - dbname
  Usage   -
  Returns -
  Args    -

 accessor: gets/sets the database name of the cross reference

must be a valid database name

=cut

sub dbname {
    my $self = shift;
    $self->{dbname} = shift if @_;
    return $self->{dbname};
}


=head2 description

  Alias   - name
  Usage   -
  Returns -
  Args    -

=cut

sub description {
    my $self = shift;
    $self->{description} = shift if @_;
    return $self->{description};
}


sub to_idl_struct {
    my $self = shift;
    return
      {
       dbname=>$self->dbname,
       keyvalue=>$self->accession
      };
}





1;