#
# BioPerl module for Bio::Ontology::InterProTerm
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
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

Bio::Ontology::InterProTerm - Implementation of InterProI term interface

=head1 SYNOPSIS

  my $term = Bio::Ontology::InterProTerm->new( 
        -interpro_id => "IPR000001",
		  -name => "Kringle",
		  -definition => "Kringles are autonomous structural domains ...",
		  -ontology => "Domain"
					     );
  print $term->interpro_id(), "\n";
  print $term->name(), "\n";
  print $term->definition(), "\n";
  print $term->is_obsolete(), "\n";
  print $term->ontology->name(), "\n";

=head1 DESCRIPTION

This is a simple extension of L<Bio::Ontology::Term> for InterPro terms.

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

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Peter Dimitrov

Email dimitrov@gnf.org

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Ontology::InterProTerm;
use strict;

use Bio::Annotation::Reference;

use constant INTERPRO_ID_DEFAULT => "IPR000000";

use base qw(Bio::Ontology::Term);

=head2 new

 Title   : new
 Usage   : $term = Bio::Ontology::InterProTerm->new( -interpro_id => "IPR000002",
						     -name => "Cdc20/Fizzy",
						     -definition => "The Cdc20/Fizzy region is almost always ...",
						     -ontology => "Domain"
						   );

 Function: Creates a new Bio::Ontology::InterProTerm.
 Example :
 Returns : A new Bio::Ontology::InterProTerm object.
 Args    :
  -interpro_id => the InterPro ID of the term. Has the form IPRdddddd, where dddddd is a zero-padded six digit number
  -name => the name of this InterPro term [scalar]
  -definition => the definition/abstract of this InterPro term [scalar]
  -ontology => ontology of InterPro terms [Bio::Ontology::OntologyI]
  -comment => a comment [scalar]

=cut

sub new{
  my ($class, @args) = @_;
  my $self = $class->SUPER::new(@args);

  my ( $interpro_id,
       $short_name)
      = $self->_rearrange( [qw( INTERPRO_ID
				SHORT_NAME
				)
			    ], @args );

  $interpro_id && $self->interpro_id( $interpro_id );
  $short_name && $self->short_name( $short_name );

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
    my $self = shift;

    # first call the inherited version to properly chain up the hierarchy
    $self->SUPER::init(@_);

    # then only initialize what we implement ourselves here
    $self->interpro_id( INTERPRO_ID_DEFAULT );
    $self->short_name("");

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
    return $self->identifier($value);
  }

  return $self->identifier();
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

=head2 get_references

 Title   : get_references
 Usage   :
 Function: Get the references for this InterPro term.
 Example :
 Returns : An array of L<Bio::Annotation::Reference> objects
 Args    :


=cut

# Defined in parent class

#sub get_references{
#    my $self = shift;
#
#    return @{$self->{"_references"}} if exists($self->{"_references"});
#    return ();
#}

=head2 add_reference

 Title   : add_reference
 Usage   :
 Function: Add one or more references to this InterPro term.
 Example :
 Returns : 
 Args    : One or more L<Bio::Annotation::Reference> objects.


=cut

# Defined in parent class

#sub add_reference{
#    my $self = shift;
#
#    $self->{"_references"} = [] unless exists($self->{"_references"});
#    push(@{$self->{"_references"}}, @_);
#}

=head2 remove_references

 Title   : remove_references
 Usage   :
 Function: Remove all references for this InterPro term.
 Example :
 Returns : The list of previous references as an array of
           L<Bio::Annotation::Reference> objects.
 Args    :


=cut

# Defined in parent class
#sub remove_references{
#    my $self = shift;
#
#    my @arr = $self->get_references();
#    $self->{"_references"} = [];
#    return @arr;
#}

=head2 get_members

 Title   : get_members
 Usage   : @arr = get_members()
 Function: Get the list of member(s) for this object.
 Example :
 Returns : An array of Bio::Annotation::DBLink objects
 Args    :


=cut

sub get_members{
    my $self = shift;
    return $self->get_dbxrefs('member_list');
}

=head2 add_member

 Title   : add_member
 Usage   :
 Function: Add one or more member(s) to this object.
 Example :
 Returns : 
 Args    : One or more Bio::Annotation::DBLink objects.


=cut

sub add_member{
    my $self = shift;
    $self->add_dbxref(-dbxrefs => \@_, -context => 'member_list');
}

=head2 remove_members

 Title   : remove_members
 Usage   :
 Function: Remove all members for this class.
 Example :
 Returns : The list of previous members as an array of
           Bio::Annotation::DBLink objects.
 Args    :


=cut

sub remove_members{
    my $self = shift;
    return $self->remove_dbxrefs('member_list');
}

=head2 get_examples

 Title   : get_examples
 Usage   : @arr = get_examples()
 Function: Get the list of example(s) for this object.

           This is an element of the InterPro xml schema.

 Example :
 Returns : An array of Bio::Annotation::DBLink objects
 Args    :


=cut

sub get_examples{
    my $self = shift;
    return $self->get_dbxrefs('example_list');
}

=head2 add_example

 Title   : add_example
 Usage   :
 Function: Add one or more example(s) to this object.

           This is an element of the InterPro xml schema.

 Example :
 Returns : 
 Args    : One or more Bio::Annotation::DBLink objects.


=cut

sub add_example{
    my $self = shift;
    return $self->add_dbxref(-dbxrefs => \@_, -context => 'example_list');
}

=head2 remove_examples

 Title   : remove_examples
 Usage   :
 Function: Remove all examples for this class.

           This is an element of the InterPro xml schema.

 Example :
 Returns : The list of previous examples as an array of
           Bio::Annotation::DBLink objects.
 Args    :


=cut

sub remove_examples{
    my $self = shift;
    return $self->remove_dbxrefs('example_list');
}

=head2 get_external_documents

 Title   : get_external_documents
 Usage   : @arr = get_external_documents()
 Function: Get the list of external_document(s) for this object.

           This is an element of the InterPro xml schema.

 Example :
 Returns : An array of Bio::Annotation::DBLink objects
 Args    :


=cut

sub get_external_documents{
    my $self = shift;
    return $self->get_dbxrefs('external_doc_list');
}

=head2 add_external_document

 Title   : add_external_document
 Usage   :
 Function: Add one or more external_document(s) to this object.

           This is an element of the InterPro xml schema.

 Example :
 Returns : 
 Args    : One or more Bio::Annotation::DBLink objects.


=cut

sub add_external_document{
    my $self = shift;
    return $self->add_dbxref(-dbxrefs => \@_, -context => 'external_doc_list');
}

=head2 remove_external_documents

 Title   : remove_external_documents
 Usage   :
 Function: Remove all external_documents for this class.

           This is an element of the InterPro xml schema.

 Example :
 Returns : The list of previous external_documents as an array of
           Bio::Annotation::DBLink objects.
 Args    :


=cut

sub remove_external_documents{
    my $self = shift;
    return $self->remove_dbxrefs('external_doc_list');
}

=head2 class_list

 Title   : class_list
 Usage   : $obj->class_list($newval)
 Function: Set/get for class list element of the InterPro xml schema
 Example : 
 Returns : reference to an array of Bio::Annotation::DBLink objects
 Args    : reference to an array of Bio::Annotation::DBLink objects

=cut

# this is inconsistent with the above, but we work around it and hope nothing
# breaks

sub class_list{
    my ($self, $value) = @_;
    if( defined $value && ref $value eq 'ARRAY') {
        if (!@$value) {
            # passing an empty array ref is essentially same as remove_dbxrefs,
            # so do that
            $self->remove_dbxrefs('class_list');
        } else {
            $self->add_dbxref(-dbxrefs => $value, -context => 'class_list');
        }
    }
    return [$self->get_dbxrefs('class_list')];
}

=head2 to_string

 Title   : to_string()
 Usage   : print $term->to_string();
 Function: to_string method for InterPro terms.
 Returns : A string representation of this InterPro term.
 Args    :

=cut

sub to_string {
    my ($self) = @_;
    my $s = "";

    $s .= "-- InterPro id:\n";
    $s .= $self->interpro_id() . "\n";
    if ( defined $self->name ) {
        $s .= "-- Name:\n";
        $s .= $self->name() . "\n";
        $s .= "-- Definition:\n";
        $s .= ( $self->definition() || '' ) . "\n";
        $s .= "-- Category:\n";
        if ( defined( $self->ontology() ) ) {
            $s .= $self->ontology()->name() . "\n";
        }
        else {
            $s .= "\n";
        }
        $s .= "-- Version:\n";
        $s .= ( $self->version() || '' ) . "\n";
        $s .= "-- Is obsolete:\n";
        $s .= $self->is_obsolete() . "\n";
        $s .= "-- Comment:\n";
        $s .= ( $self->comment() || '' ) . "\n";
        if ( defined $self->get_references ) {
            $s .= "-- References:\n";
            foreach my $ref ( $self->get_references ) {
                $s .=
                    $ref->authors . "\n"
                  . $ref->title . "\n"
                  . $ref->location . "\n\n";
            }
            $s .= "\n";
        }
        if ( defined $self->get_members ) {
            $s .= "-- Member List:\n";
            foreach my $ref ( $self->get_members ) {
                $s .= $ref->database . "\t" . $ref->primary_id . "\n";
            }
            $s .= "\n";
        }
        if ( defined $self->get_external_documents ) {
            $s .= "-- External Document List:\n";
            foreach my $ref ( $self->get_external_documents ) {
                $s .= $ref->database . "\t" . $ref->primary_id . "\n";
            }
            $s .= "\n";
        }
        if ( defined $self->get_examples ) {
            $s .= "-- Examples:\n";
            foreach my $ref ( $self->get_examples ) {
                $s .= join( "\t",
                    map { $ref->$_ || '' } qw(database primary_id comment) )
                  . "\n";
            }
            $s .= "\n";
        }
        if ( defined $self->class_list ) {
            $s .= "-- Class List:\n";
            foreach my $ref ( @{ $self->class_list } ) {
                $s .= $ref->primary_id . "\n";
            }
            $s .= "\n";
        }
        if ( $self->get_secondary_ids ) {
            $s .= "-- Secondary IDs:\n";
            foreach my $ref ( $self->get_secondary_ids() ) {
                # TODO: getting undef here in some cases, needs to be checked
                next unless defined ($ref);
                $s .= $ref . "\n";
            }
            $s .= "\n";
        }
    }
    else {
        $s .= "InterPro term not fully instantiated\n";
    }
    return $s;
}

=head1 Deprecated methods

These are here for backwards compatibility.

=cut

=head2 secondary_ids

 Title   : secondary_ids
 Usage   : $obj->secondary_ids($newval)
 Function: This is deprecated. Use get_secondary_ids() or 
           add_secondary_id() instead.
 Example : 
 Returns : reference to an array of strings
 Args    : reference to an array of strings


=cut

sub secondary_ids{
    my $self = shift;
    my @ids;

    $self->warn("secondary_ids is deprecated. Use ".
		"get_secondary_ids/add_secondary_id instead.");

    # set mode?
    if(@_) {
	my $sids = shift;
	if($sids) {
	    $self->add_secondary_id(@$sids);
	    @ids = @$sids; 
	} else {
	    # we interpret setting to undef as removing the array
	    $self->remove_secondary_ids();
	}
    } else {
	# no; get mode
	@ids = $self->get_secondary_ids();
    }
    return \@ids;
}

1;
