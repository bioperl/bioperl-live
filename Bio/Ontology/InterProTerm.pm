# $Id$
#
# BioPerl module for Bio::Ontology::InterProTerm
#
# Cared for by Peter Dimitrov <dimitrov@gnf.org>
#
# Copyright Peter Dimitrov
# (c) Peter Dimitrov, dimitrov@gnf.org, 2002.
# (c) GNF, Genomics Institute of the Novartis Research Foundation, 2002.
#
# You may distribute this module under the same terms as perl itself.
# Refer to the Perl Artistic License (see the license accompanying this
# software package, or see http://www.perl.com/language/misc/Artistic.html)
# for the terms under which you may use, modify, and redistribute this module.
#
# THIS PACKAGE IS PROVIDED "AS IS" AND WITHOUT ANY EXPRESS OR IMPLIED
# WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED WARRANTIES OF
# MERCHANTIBILITY AND FITNESS FOR A PARTICULAR PURPOSE.
#
# POD documentation - main docs before the code

=head1 NAME

InterProTerm - Implementation of InterProI term interface

=head1 SYNOPSIS

  my $term = Bio::Ontology::InterProTerm->new( -interpro_id => "IPR000001",
					       -name => "Kringle",
					       -definition => "Kringles are autonomous structural domains ...",
					       -category => "Domain"
					     );
  print $term->interpro_id(), "\n";
  print $term->name(), "\n";
  print $term->definition(), "\n";
  print $term->is_obsolete(), "\n";
  print $term->category(), "\n";

=head1 DESCRIPTION

This is a simple implementation for InterPro terms.

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

=head1 AUTHOR - Peter Dimitrov

Email dimitrov@gnf.org

=head1 CONTRIBUTORS

Additional contributors names and emails here

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Ontology::InterProTerm;
use vars qw(@ISA);
use strict;
use Carp;
use Bio::Ontology::Term;
use Bio::Annotation::Reference;

use constant INTERPRO_ID_DEFAULT => "IPR000000";

@ISA = qw( Bio::Ontology::Term );

=head2 new

 Title   : new
 Usage   : $term = Bio::Ontology::InterProTerm->new( -interpro_id => "IPR000002",
						     -name => "Cdc20/Fizzy",
						     -definition => "The Cdc20/Fizzy region is almost always ...",
						     -category => "Domain"
						     -comment => ""
						   );
 Function: Creates a new Bio::Ontology::InterProTerm.
 Example :
 Returns : A new Bio::Ontology::InterProTerm object.
 Args    :
  -interpro_id => the InterPro ID of the term. Has the form IPRdddddd, where dddddd is a zero-padded six digit number
  -name => the name of this InterPro term [scalar]
  -definition => the definition/abstract of this InterPro term [scalar]
  -category => type of InterPro term: "Domain", "Family", "Repeat", "PTM" (post-translational modification) [scalar]
  -comment => a comment [scalar]

=cut

sub new{
  my ($class, @args) = @_;
  my $self = $class->SUPER::new(@args);

  my ( $interpro_id,
       $name,
       $short_name,
       $definition,
       $category,
       $comment )
    = $self->_rearrange( [ qw( INTERPRO_ID
			       NAME
			       SHORT_NAME
			       DEFINITION
			       CATEGORY
			       COMMENT ) ], @args );

  $self->init();

  $interpro_id && $self->interpro_id( $interpro_id );
  $name && $self->name( $name );
  $short_name && $self->short_name( $name );
  $definition && $self->definition( $definition );
  $category && $self->category( $category );
  $comment && $self->comment( $comment  );
  $self->identifier( $interpro_id );

  return $self;
}

=head2 init

 Title   : init
 Usage   : $term->init();
 Function: Initializes this InterProTerm to all "" and empty lists.
 Example :
 Returns : 
 Args    :


=cut

sub init{
  my ($self) = @_;

  $self->interpro_id( INTERPRO_ID_DEFAULT );
  $self->name("");
  $self->short_name("");
  $self->definition( "" );
  $self->comment( "" );
  $self->category( "" );

}

=head2 _check_interpro_id

 Title   : _check_interpro_id
 Usage   :
 Function: Performs simple check in order to validate that its argument has the form IPRdddddd, where dddddd is a zero-padded six digit number.
 Example :
 Returns : Returns its argument if valid, otherwise throws exception.
 Args    : String


=cut

sub _check_interpro_id{
  my ($self, $value) = @_;

  $self->throw( "InterPro ID ".$value." is incorrect\n" )
    unless ( $value =~ /^IPR\d{6}$/ ||
	     $value eq INTERPRO_ID_DEFAULT );

  return $value;
}

=head2 interpro_id

 Title   : interpro_id
 Usage   : $obj->interpro_id($newval)
 Function: Set/get for the interpro_id of this InterProTerm
 Example : 
 Returns : value of interpro_id (a scalar)
 Args    : new value (a scalar, optional)


=cut

sub interpro_id{
  my ($self, $value) = @_;

  if( defined $value) {
    $value = $self->_check_interpro_id($value);
  }

  return $self->identifier($value);
}

=head2 short_name

 Title   : short_name
 Usage   : $obj->short_name($newval)
 Function: Set/get for the short name of this InterProTerm.
 Example : 
 Returns : value of short_name (a scalar)
 Args    : new value (a scalar, optional)


=cut

sub short_name{
  my ($self, $value) = @_;

  if( defined $value) {
    $self->{'short_name'} = $value ? $value : undef;
  }

  return $self->{'short_name'};
}

=head2 protein_count

 Title   : protein_count
 Usage   : $obj->protein_count($newval)
 Function: Set/get for the protein count of this InterProTerm.
 Example : 
 Returns : value of protein_count (a scalar)
 Args    : new value (a scalar, optional)


=cut

sub protein_count{
  my ($self,$value) = @_;

  if( defined $value) {
    $self->{'protein_count'} = $value ? $value : undef;
  }

  return $self->{'protein_count'};
}

=head2 references

 Title   : references
 Usage   : $obj->references($newval)
 Function: 
 Example : 
 Returns : reference to an array of Bio::Annotation::Reference objects
 Args    : reference to an array of Bio::Annotation::Reference objects


=cut

sub references{
  my ($self, $value) = @_;

  if( defined $value) {
    $self->{'references'} = $value;
  }

  return $self->{'references'};
}

=head2 member_list

 Title   : member_list
 Usage   : $obj->member_list($newval)
 Function: Set/get for member list element of the InterPro xml schema
 Example : 
 Returns : reference to an array of Bio::Annotation::DBLink objects
 Args    : reference to an array of Bio::Annotation::DBLink objects


=cut

sub member_list{
  my ($self, $value) = @_;

  if( defined $value) {
    $self->{'member_list'} = $value;
  }

  return $self->{'member_list'};
}

=head2 examples

 Title   : examples
 Usage   : $obj->examples($newval)
 Function: Set/get for example list element of the InterPro xml schema
 Example : 
 Returns : reference to an array of Bio::Annotation::DBLink objects
 Args    : reference to an array of Bio::Annotation::DBLink objects


=cut

sub examples{
  my ($self, $value) = @_;

  if( defined $value) {
    $self->{'examples'} = $value;
  }

  return $self->{'examples'};
}

=head2 external_doc_list

 Title   : external_doc_list
 Usage   : $obj->external_doc_list($newval)
 Function: Set/get for external doc list element of the InterPro xml schema
 Example : 
 Returns : reference to an array of Bio::Annotation::DBLink objects
 Args    : reference to an array of Bio::Annotation::DBLink objects


=cut

sub external_doc_list{
  my ($self, $value) = @_;

  if( defined $value) {
    $self->{'external_doc_list'} = $value;
  }

  return $self->{'external_doc_list'};
}

=head2 class_list

 Title   : class_list
 Usage   : $obj->class_list($newval)
 Function: Set/get for class list element of the InterPro xml schema
 Example : 
 Returns : reference to an array of Bio::Annotation::DBLink objects
 Args    : reference to an array of Bio::Annotation::DBLink objects


=cut

sub class_list{
  my ($self, $value) = @_;

  if( defined $value) {
    $self->{'class_list'} = $value;
  }

  return $self->{'class_list'};
}

=head2 secondary_ids

 Title   : secondary_ids
 Usage   : $obj->secondary_ids($newval)
 Function: Set/get for secondary accession numbers
 Example : 
 Returns : reference to an array of strings
 Args    : reference to an array of strings


=cut

sub secondary_ids{
  my ($self, $value) = @_;

  if( defined $value) {
    $self->{'secondary_ids'} = $value;
  }

  return $self->{'secondary_ids'};
}

=head2 to_string

 Title   : to_string()
 Usage   : print $term->to_string();
 Function: to_string method for InterPro terms.
 Returns : A string representation of this InterPro term.
 Args    :

=cut

sub to_string {
  my($self) = @_;
  my $s = "";

  $s .= "-- InterPro id:\n";
  $s .= $self->interpro_id()."\n";
  if (defined $self->name) {
    $s .= "-- Name:\n";
    $s .= $self->name()."\n";
    $s .= "-- Definition:\n";
    $s .= $self->definition()."\n";
    $s .= "-- Category:\n";
    if ( defined( $self->category() ) ) {
      $s .= $self->category()->name()."\n";
    } else {
      $s .= "\n";
    }
    $s .= "-- Version:\n";
    $s .= $self->version()."\n";
    $s .= "-- Is obsolete:\n";
    $s .= $self->is_obsolete()."\n";
    $s .= "-- Comment:\n";
    $s .= $self->comment()."\n"; 
    if (defined $self->references) {
      $s .= "-- References:\n";
      foreach my $ref ( @{$self->references} ) {
	$s .= $ref->authors."\n".$ref->title."\n".$ref->location."\n\n";
      };
      $s .= "\n";
    }
    if (defined $self->member_list) {
      $s .= "-- Member List:\n";
      foreach my $ref ( @{$self->member_list} ) {
	$s .= $ref->database."\t".$ref->primary_id."\n";
      };
      $s .= "\n";
    }
    if (defined $self->external_doc_list) {
      $s .= "-- External Document List:\n";
      foreach my $ref ( @{$self->external_doc_list} ) {
	$s .= $ref->database."\t".$ref->primary_id."\n";
      };
      $s .= "\n";
    }
    if (defined $self->examples) {
      $s .= "-- Examples:\n";
      foreach my $ref ( @{$self->examples} ) {
	$s .= $ref->database."\t".$ref->primary_id."\t".$ref->comment."\n";
      };
      $s .= "\n";
    }
    if (defined $self->class_list) {
      $s .= "-- Class List:\n";
      foreach my $ref ( @{$self->class_list} ) {
	$s .= $ref->primary_id."\n";
      };
      $s .= "\n";
    }
    if (defined $self->secondary_ids) {
      $s .= "-- Secondary IDs:\n";
      foreach my $ref ( @{$self->secondary_ids} ) {
	$s .= $ref."\n";
      };
      $s .= "\n";
    }
  }
  else {
    $s .= "InterPro term not fully instantiated\n";
  }
  return $s;
}

1;
