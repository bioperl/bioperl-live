#
# BioPerl module for Bio::Ontology::RelationshipType  
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Christian M. Zmasek <czmasek-at-burnham.org> or <cmzmasek@yahoo.com>
#
# (c) Christian M. Zmasek, czmasek-at-burnham.org, 2002.
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
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Ontology::RelationshipType  - a relationship type for an ontology

=head1 SYNOPSIS

  #

=head1 DESCRIPTION

This class can be used to model various types of relationships
(such as "IS_A", "PART_OF", "CONTAINS", "FOUND_IN", "RELATED_TO").

This class extends L<Bio::Ontology::Term>, so it essentially is-a
L<Bio::Ontology::TermI>. In addition, all methods are overridden such
as to make the object immutable.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to the 
Bioperl mailing lists  Your participation is much appreciated.

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
the bugs and their resolution.  Bug reports can be submitted via
the web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR

Christian M. Zmasek

Email: czmasek-at-burnham.org  or  cmzmasek@yahoo.com

WWW:   http://monochrome-effect.net/

Address:

  Genomics Institute of the Novartis Research Foundation
  10675 John Jay Hopkins Drive
  San Diego, CA 92121

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...

package Bio::Ontology::RelationshipType;
use strict;


use constant PART_OF    => "PART_OF";
use constant RELATED_TO => "RELATED_TO";
use constant IS_A       => "IS_A";
use constant CONTAINS   => "CONTAINS";
use constant FOUND_IN   => "FOUND_IN";
use constant REGULATES   => "REGULATES";
use constant POSITIVELY_REGULATES   => "POSITIVELY_REGULATES";
use constant NEGATIVELY_REGULATES   => "NEGATIVELY_REGULATES";


use base qw(Bio::Ontology::Term);


#
# cache for terms
#
my %term_name_map = ();


=head2 get_instance

 Title   : get_instance
 Usage   : $IS_A       = Bio::Ontology::RelationshipType->get_instance( "IS_A" );
           $PART_OF    = Bio::Ontology::RelationshipType->get_instance( "PART_OF" );
           $RELATED_TO = Bio::Ontology::RelationshipType->get_instance( "RELATED_TO" );
           $CONTAINS   = Bio::Ontology::RelationshipType->get_instance( "CONTAINS" );
           $FOUND_IN   = Bio::Ontology::RelationshipType->get_instance( "FOUND_IN" );
 Function: Factory method to create instances of RelationshipType
 Returns : [Bio::Ontology::RelationshipType]
 Args    : "IS_A" or "PART_OF" or "CONTAINS" or "FOUND_IN" or 
           "RELATED_TO" [scalar]
           the ontology [Bio::Ontology::OntologyI] (optional)

=cut

sub get_instance {
    my ( $class, $name, $ont ) = @_;

    $class->throw("must provide predicate name") unless $name;

    # is one in the cache?
    my $reltype = $term_name_map{$name};

    if($reltype &&
       # check whether ontologies match
       (($ont && $reltype->ontology() &&
	 ($ont->name() eq $reltype->ontology->name())) ||
	(! ($reltype->ontology() || $ont)))) {
	# we're done, return cached type
	return $reltype;
    }
    # valid relationship type?

#
#see the cell ontology.  this code is too strict, even for dag-edit files. -allen
#
#    if ( ! (($name eq IS_A) || ($name eq PART_OF) ||
#	    ($name eq CONTAINS) || ( $name eq FOUND_IN ))) {
#        my $msg = "Found unknown type of relationship: [" . $name . "]\n";
#        $msg .= "Known types are: [" . IS_A . "], [" . PART_OF . "], [" . CONTAINS . "], [" . FOUND_IN . "]";
#        $class->throw( $msg );
#    }
    # if we get here we need to create the rel.type
    $reltype = $class->new(-name     => $name,
			   -ontology => $ont);
    # cache it (FIXME possibly overrides one from another ontology)
    $term_name_map{$name} = $reltype;
    return $reltype;
} # get_instance


=head2 init

 Title   : init()
 Usage   : $type->init();
 Function: Initializes this to all undef and empty lists.
 Returns :
 Args    :

=cut

sub init {
    my $self = shift;

    $self->SUPER::init();

    # at this point we don't really need to do anything special for us
} # init


=head2 equals

 Title   : equals
 Usage   : if ( $type->equals( $other_type ) ) { ...
 Function: Compares this type to another one, based on string "eq" of
           the "identifier" field, if at least one of the two types has
           the identifier set, or string eq of the name otherwise.
 Returns : true or false
 Args    : [Bio::Ontology::RelationshipType]

=cut

sub equals {
    my( $self, $type ) = @_;

    $self->_check_class( $type, "Bio::Ontology::RelationshipType" );

    if ( $self->identifier() xor $type->identifier() ) {
        $self->warn("comparing relationship types when only ".
		    "one has an identifier will always return false" );
    }

    return
 	($self->identifier() || $type->identifier()) ?
	$self->identifier() eq $type->identifier() :
	$self->name() eq $type->name();
	
} # equals


=head2 identifier

 Title   : identifier
 Usage   : $term->identifier( "IS_A" );
           or
           print $term->identifier();
 Function: Set/get for the immutable identifier of this Type.
 Returns : The identifier [scalar].
 Args    : The identifier [scalar] (optional).

=cut

sub identifier {
    my $self = shift;
    my $ret = $self->SUPER::identifier();
    if(@_) {
	$self->throw($self->veto_change("identifier",$ret,$_[0]))
	    if $ret && ($ret ne $_[0]);
	$ret = $self->SUPER::identifier(@_);
    }
    return $ret;
} # identifier


=head2 name

 Title   : name
 Usage   : $term->name( "is a type" );
           or
           print $term->name();
 Function: Set/get for the immutable name of this Type.
 Returns : The name [scalar].
 Args    : The name [scalar] (optional).

=cut

sub name {
    my $self = shift;
    my $ret = $self->SUPER::name();
    if(@_) {
	$self->throw($self->veto_change("name",$ret,$_[0]))
	    if $ret && ($ret ne $_[0]);
	$ret = $self->SUPER::name(@_);
    }
    return $ret;
} # name





=head2 definition

 Title   : definition
 Usage   : $term->definition( "" );
           or
           print $term->definition();
 Function: Set/get for the immutable definition of this Type.
 Returns : The definition [scalar].
 Args    : The definition [scalar] (optional).

=cut

sub definition {
    my $self = shift;
    my $ret = $self->SUPER::definition();
    if(@_) {
	$self->veto_change("definition",$ret,$_[0]) 
	    if $ret && ($ret ne $_[0]);
	$ret = $self->SUPER::definition(@_);
    }
    # let's be nice and return something readable here
    return $ret if $ret;
    return $self->name()." relationship predicate (type)" if $self->name();
} # definition



=head2 ontology

 Title   : ontology
 Usage   : $term->ontology( $top );
           or
           $top = $term->ontology();
 Function: Set/get for the ontology this relationship type lives in.
 Returns : The ontology [Bio::Ontology::OntologyI].
 Args    : On set, the ontology [Bio::Ontology::OntologyI] (optional).

=cut

sub ontology {
    my $self = shift;
    my $ret = $self->SUPER::ontology();
    if(@_) {
	my $ont = shift;
	if($ret) {
	    $self->throw($self->veto_change("ontology",$ret->name,
					    $ont ? $ont->name : $ont))
		unless $ont && ($ont->name() eq $ret->name());
	}
	$ret = $self->SUPER::ontology($ont,@_);
    }
    return $ret;
} # category



=head2 version

 Title   : version
 Usage   : $term->version( "1.00" );
           or
           print $term->version();
 Function: Set/get for immutable version information.
 Returns : The version [scalar].
 Args    : The version [scalar] (optional).

=cut

sub version {
    my $self = shift;
    my $ret = $self->SUPER::version();
    if(@_) {
	$self->throw($self->veto_change("version",$ret,$_[0]))
	    if $ret && ($ret ne $_[0]);
	$ret = $self->SUPER::version(@_);
    }
    return $ret;
} # version



=head2 is_obsolete

 Title   : is_obsolete
 Usage   : $term->is_obsolete( 1 );
           or
           if ( $term->is_obsolete() )
 Function: Set/get for the immutable obsoleteness of this Type.
 Returns : the obsoleteness [0 or 1].
 Args    : the obsoleteness [0 or 1] (optional).

=cut

sub is_obsolete {
    my $self = shift;
    my $ret = $self->SUPER::is_obsolete();
    if(@_) {
	$self->throw($self->veto_change("is_obsolete",$ret,$_[0]))
	    if $ret && ($ret != $_[0]);
	$ret = $self->SUPER::is_obsolete(@_);
    }
    return $ret;
} # is_obsolete


=head2 comment

 Title   : comment
 Usage   : $term->comment( "..." );
           or
           print $term->comment();
 Function: Set/get for an arbitrary immutable comment about this Type.
 Returns : A comment.
 Args    : A comment (optional).

=cut

sub comment {
    my $self = shift;
    my $ret = $self->SUPER::comment();
    if(@_) {
	$self->throw($self->veto_change("comment",$ret,$_[0]))
	    if $ret && ($ret ne $_[0]);
	$ret = $self->SUPER::comment(@_);
    }
    return $ret;
} # comment

=head1 Private methods 

May be overridden in a derived class, but should never be called from
outside.

=cut

sub _check_class {
    my ( $self, $value, $expected_class ) = @_;

    if ( ! defined( $value ) ) {
        $self->throw( "Found [undef] where [$expected_class] expected" );
    }
    elsif ( ! ref( $value ) ) {
        $self->throw( "Found [scalar] where [$expected_class] expected" );
    }
    elsif ( ! $value->isa( $expected_class ) ) {
        $self->throw( "Found [" . ref( $value ) . "] where [$expected_class] expected" );
    }

} # _check_type

=head2 veto_change

 Title   : veto_change
 Usage   :
 Function: Called if an attribute is changed. Setting an attribute is
           considered a change if it had a value before and the attempt
           to set it would change the value.

           This method returns the message to be printed in the exception.

 Example :
 Returns : A string
 Args    : The name of the attribute that was attempted to change.
           Optionally, the old value and the new value for reporting
           purposes only.

=cut

sub veto_change{
    my ($self,$attr,$old,$new) = @_;

    my $changetype = $old ? ($new ? "change" : "unset") : "change";
    my $msg = "attempt to $changetype attribute $attr in ".ref($self).
    ", which is immutable";
    $msg .= " (\"$old\" to \"$new\")" if $old && $new;
    return $msg;
}

1;
