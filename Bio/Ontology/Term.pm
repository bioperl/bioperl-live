#
# BioPerl module for Bio::Ontology::Term
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

Bio::Ontology::Term - implementation of the interface for ontology terms

=head1 SYNOPSIS

#get Bio::Ontology::TermI somehow.

  print $term->identifier(), "\n";
  print $term->name(), "\n";
  print $term->definition(), "\n";
  print $term->is_obsolete(), "\n";
  print $term->comment(), "\n";

  foreach my $synonym ( $term->each_synonym() ) {
      print $synonym, "\n";
  }

=head1 DESCRIPTION

This is a simple implementation for ontology terms providing basic
methods (it provides no functionality related to graphs). It
implements the L<Bio::Ontology::TermI> interface.

This class also implements L<Bio::IdentifiableI> and
L<Bio::DescribableI>.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists.  Your participation is much appreciated.

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
the bugs and their resolution.  Bug reports can be submitted via the web:

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
methods.

=cut


# Let the code begin...

package Bio::Ontology::Term;
use strict;
use Bio::Ontology::Ontology;
use Bio::Ontology::OntologyStore;
use Bio::Annotation::DBLink;
use Data::Dumper;

use constant TRUE    => 1;
use constant FALSE   => 0;

use base qw(Bio::Root::Root Bio::Ontology::TermI Bio::IdentifiableI Bio::DescribableI);

=head2 new

 Title   : new
 Usage   : $term = Bio::Ontology::Term->new(
                -identifier  => "16847",
                -name        => "1-aminocyclopropane-1-carboxylate synthase",
                -definition  => "Catalysis of ...",
                -is_obsolete => 0,
                -comment     => "" );
 Function: Creates a new Bio::Ontology::Term.
 Returns : A new Bio::Ontology::Term object.
 Args    : -identifier            => the identifier of this term [scalar]
           -name                  => the name of this term [scalar]
           -definition            => the definition of this term [scalar]
           -ontology              => the ontology this term lives in
                                     (a Bio::Ontology::OntologyI object)
           -version               => version information [scalar]
           -is_obsolete           => the obsoleteness of this term [0 or 1]
           -comment               => a comment [scalar]
           -dblinks               => Bio::Annotation::DBLink objects
                                     [reference to array]
           -references            => Bio::Annotation::Reference objects
                                     [reference to array]

See L<Bio::Ontology::OntologyI>, L<Bio::Annotation::Reference>,
L<Bio::Annotation::DBLink>.

=cut

sub new {

    my( $class,@args ) = @_;

    my $self = $class->SUPER::new( @args );
    my ( $identifier, $name, $definition, $category, $ont, $version,
        $is_obsolete, $comment, $dblinks, $dbxrefs, $references)
        = $self->_rearrange( [
        qw(IDENTIFIER NAME DEFINITION CATEGORY ONTOLOGY VERSION IS_OBSOLETE
        COMMENT DBLINKS DBXREFS REFERENCES) ], @args );

    $self->init();

    defined($identifier)   && $self->identifier( $identifier );
    defined($name)         && $self->name( $name );
    defined($definition)   && $self->definition( $definition );
    defined($category)     && $self->category( $category );
    defined($ont)          && $self->ontology( $ont );
    defined($version)      && $self->version( $version );
    defined($is_obsolete)  && $self->is_obsolete( $is_obsolete );
    defined($comment)      && $self->comment( $comment  );
    defined($dbxrefs)      && $self->add_dbxref(-dbxrefs => $dbxrefs);
    # deprecated methods, allow to pass on to get the dep. notification
    ref($dblinks)          && $self->add_dblink(@$dblinks);
    ref($references)       && $self->add_reference(@$references);

    return $self;
} # new



sub init {

    my $self = shift;

    $self->identifier(undef);
    $self->name(undef);
    $self->comment(undef);
    $self->definition(undef);
    $self->ontology(undef);
    $self->is_obsolete(0);
    $self->remove_synonyms();
    $self->remove_dbxrefs();
    $self->remove_references;
    $self->remove_secondary_ids();

} # init



=head2 identifier

 Title   : identifier
 Usage   : $term->identifier( "GO:0003947" );
           or
           print $term->identifier();
 Function: Set/get for the identifier of this Term.
 Returns : The identifier [scalar].
 Args    : The identifier [scalar] (optional).

=cut

sub identifier {
    my $self = shift;

    return $self->{'identifier'} = shift if @_;
    return $self->{'identifier'};
} # identifier


=head2 name

 Title   : name
 Usage   : $term->name( "N-acetylgalactosaminyltransferase" );
           or
           print $term->name();
 Function: Set/get for the name of this Term.
 Returns : The name [scalar].
 Args    : The name [scalar] (optional).

=cut

sub name {
    my $self = shift;

    return $self->{'name'} = shift if @_;
    return $self->{'name'};
} # name


=head2 definition

 Title   : definition
 Usage   : $term->definition( "Catalysis of ..." );
           or
           print $term->definition();
 Function: Set/get for the definition of this Term.
 Returns : The definition [scalar].
 Args    : The definition [scalar] (optional).

=cut

sub definition {
    my $self = shift;

    return $self->{'definition'} = shift if @_;
    return $self->{'definition'};
} # definition


=head2 ontology

 Title   : ontology
 Usage   : $ont = $term->ontology();
           or
           $term->ontology( $ont );
 Function: Get the ontology this term is in.

           Note that with the ontology in hand you can query for all
           related terms etc.

 Returns : The ontology of this Term as a Bio::Ontology::OntologyI
           implementing object.
 Args    : On set, the  ontology of this Term as a Bio::Ontology::OntologyI
           implementing object or a string representing its name.

See L<Bio::Ontology::OntologyI>.

=cut

sub ontology {
    my $self = shift;
    my $ont;

    if(@_) {
        $ont = shift;
        if($ont) {
            $ont = Bio::Ontology::Ontology->new(-name => $ont) if ! ref($ont);
            if(! $ont->isa("Bio::Ontology::OntologyI")) {
                $self->throw(ref($ont)." does not implement ".
                             "Bio::Ontology::OntologyI. Bummer.");
            }
        }
        return $self->{"_ontology"} = $ont;
    }
    return $self->{"_ontology"};
} # ontology

=head2 version

 Title   : version
 Usage   : $term->version( "1.00" );
           or
           print $term->version();
 Function: Set/get for version information.
 Returns : The version [scalar].
 Args    : The version [scalar] (optional).

=cut

sub version {
    my $self = shift;

    return $self->{'version'} = shift if @_;
    return $self->{'version'};
} # version

=head2 is_obsolete

 Title   : is_obsolete
 Usage   : $term->is_obsolete( 1 );
           or
           if ( $term->is_obsolete() )
 Function: Set/get for the obsoleteness of this Term.
 Returns : the obsoleteness [0 or 1].
 Args    : the obsoleteness [0 or 1] (optional).

=cut

sub is_obsolete{
    my $self = shift;

    return $self->{'is_obsolete'} = shift if @_;
    return $self->{'is_obsolete'};
} # is_obsolete


=head2 comment

 Title   : comment
 Usage   : $term->comment( "Consider the term ..." );
           or
           print $term->comment();
 Function: Set/get for an arbitrary comment about this Term.
 Returns : A comment.
 Args    : A comment (optional).

=cut

sub comment{
    my $self = shift;

    return $self->{'comment'} = shift if @_;
    return $self->{'comment'};
} # comment

=head2 get_synonyms

 Title   : get_synonyms
 Usage   : @aliases = $term->get_synonyms;
 Function: Returns a list of aliases of this Term.
 Returns : A list of aliases [array of [scalar]].
 Args    :

=cut

sub get_synonyms {
    my $self = shift;

    return @{ $self->{ "_synonyms" } } if exists($self->{ "_synonyms" });
    return ();
} # get_synonyms


=head2 add_synonym

 Title   : add_synonym
 Usage   : $term->add_synonym( @asynonyms );
           or
           $term->add_synonym( $synonym );
 Function: Pushes one or more synonyms into the list of synonyms.
 Returns :
 Args    : One synonym [scalar] or a list of synonyms [array of [scalar]].

=cut

sub add_synonym {
    my ( $self, @values ) = @_;

    return unless( @values );

    # avoid duplicates
    foreach my $syn (@values) {
        next if grep { $_ eq $syn; } @{$self->{ "_synonyms" }};
        push( @{ $self->{ "_synonyms" } }, $syn );
    }

} # add_synonym


=head2 remove_synonyms

 Title   : remove_synonyms()
 Usage   : $term->remove_synonyms();
 Function: Deletes (and returns) the synonyms of this Term.
 Returns : A list of synonyms [array of [scalar]].
 Args    :

=cut

sub remove_synonyms {
    my ( $self ) = @_;

    my @a = $self->get_synonyms();
    $self->{ "_synonyms" } = [];
    return @a;

} # remove_synonyms

=head2 get_dblinks

 Title   : get_dblinks()
 Usage   : @ds = $term->get_dblinks();
 Function: Returns a list of each dblinks of this GO term.
 Returns : A list of dblinks [array of [scalars]].
 Args    : A scalar indicating the context (optional).
           If omitted, all dblinks will be returned.
 Note    : deprecated method due to past use of mixed data types; use
           get_dbxrefs() instead, which handles both strings and DBLink
           instances

=cut

sub get_dblinks {
    my ($self, $context) = @_;
    $self->deprecated("Use of get_dblinks is deprecated.  Note that prior use\n".
                      "of this method could return either simple scalar values\n".
                      "or Bio::Annotation::DBLink instances; only \n".
                      "Bio::Annotation::DBLink is now supported.\n ".
                      "Use get_dbxrefs() instead");
    $self->get_dbxrefs($context);
} # get_dblinks

=head2 get_dbxrefs

 Title   : get_dbxrefs()
 Usage   : @ds = $term->get_dbxrefs();
 Function: Returns a list of each link for this term.

           If an implementor of this interface permits modification of
           this array property, the class should define at least
           methods add_dbxref() and remove_dbxrefs(), with obvious
           functionality.

 Returns : A list of L<Bio::Annotation::DBLink> instances
 Args    : [optional] string which specifies context (default : returns all dbxrefs)

=cut

sub get_dbxrefs {
    my ($self, $context) = shift;
    my @dbxrefs;
    if (defined($context)) {
        if (exists($self->{_dblinks}->{$context})) {
            @dbxrefs =  @{$self->{_dblinks}->{$context}};
        }
    } else {
        @dbxrefs = map { @$_ } values %{$self->{_dblinks}} ;
    }
    return @dbxrefs;
} # get_dbxrefs

=head2 get_dblink_context

  Title   : get_dblink_context
  Usage   : @context = $term->get_dblink_context;
  Function: Return all context existing in Term
  Returns : a list of scalar
  Args    : [none]
  Note    : deprecated method due to past use of mixed data types; use
            get_dbxref_context() instead

=cut

sub get_dblink_context {
    my $self=shift;
    $self->deprecated("Use of get_dblink_context() is deprecated; use get_dbxref_context() instead");
    return $self->get_dbxref_context(@_);
}

=head2 get_dbxref_context

  Title   : get_dbxref_context
  Usage   : @context = $term->get_dbxref_context;
  Function: Return all context strings existing in Term
  Returns : a list of scalars
  Args    : [none]

=cut

sub get_dbxref_context {
    my $self=shift;
    return keys %{$self->{_dblinks}};
}

=head2 add_dblink

 Title   : add_dblink
 Usage   : $term->add_dblink( @dbls );
           or
           $term->add_dblink( $dbl );
 Function: Pushes one or more dblinks onto the list of dblinks.
 Returns :
 Args    : One or more L<Bio::Annotation::DBLink> instances
 Note    : deprecated method due to past use of mixed data types; use
           add_dbxref() instead, which handles both strings and
           DBLink instances

=cut

sub add_dblink {
    my $self = shift;
    $self->deprecated("Use of simple strings and add_dblink() is deprecated; use\n".
                      "Bio::Annotation::DBLink instances and add_dbxref() instead");
    # here we're assuming the data is in a simple DB:ID format
    my @dbxrefs;
    for my $string (@_) {
        my ($db, $id) = split(':',$string);
        push @dbxrefs, Bio::Annotation::DBLink->new(-database => $db, -primary_id => $id);
    }
    return $self->add_dbxref(-dbxrefs => \@dbxrefs, -context => '_default');
} # add_dblink

=head2 add_dbxref

 Title   : add_dbxref
 Usage   : $term->add_dbxref( @dbls );
           or
           $term->add_dbxref( $dbl );
 Function: Pushes one or more dblinks onto the list of dblinks.
 Returns :
 Args    : -dbxrefs : array ref of Bio::Annotation::DBLink instances
           -context : string designating the context for the DBLink
                       (default : '_default' - contextless)

=cut

sub add_dbxref {
    my $self = shift;
    my ($links, $context) = $self->_rearrange([qw(DBXREFS CONTEXT)],@_);
    return unless defined $links;
    $context ||= '_default';
    $self->throw("DBLinks must be passed as an array reference") if ref $links ne 'ARRAY';
    foreach my $dbxref (@{$links}) {
        $self->throw("$dbxref is not a DBLink") unless ref $dbxref &&
            $dbxref->isa('Bio::Annotation::DBLink');
        $self->throw("'all' is a reserved word for context.") if $context eq 'all';
        if (! exists($self->{_dblinks}->{$context})) {
            $self->{_dblinks}->{$context} = [];
        }
        my $linktext = ref $dbxref ? $dbxref->display_text : $dbxref;
        if (grep {$_->display_text eq $linktext}
            @{$self->{_dblinks}->{$context}})
        {
            $self->warn("DBLink exists in the dblink of $context");
        }
        push @{$self->{_dblinks}->{$context}}, $dbxref;
    }
} # add_dbxref

# alias, for consistency
*add_dbxrefs = \&add_dbxref;

=head2 has_dblink

  Title   : has_dblink
  Usage   : $term->has_dblink($dblink);
  Function: Checks if a DBXref is already existing in the OBOterm object
  Return  : TRUE/FALSE
  Args    : [arg1] A DBxref identifier
  Note    : deprecated method due to past use of mixed data types; use
            has_dbxref() instead, which handles both strings and
            DBLink instances

=cut

sub has_dblink {
    my ( $self, $value ) = @_;
    $self->deprecated("use of has_dblink() is deprecated; use has_dbxref() instead");
    return $self->has_dbxref($value);
}

=head2 has_dbxref

  Title   : has_dbxref
  Usage   : $term->has_dbxref($dbxref);
  Function: Checks if a dbxref string is already existing in the OBOterm object
  Return  : TRUE/FALSE
  Args    : [arg1] A DBxref identifier (string).
            Bio::Annotation::DBLink::display_text() is used for comparison
            against the string.

=cut

sub has_dbxref {
    my ( $self, $value ) = @_;
    return unless defined $value;
    my $context = "_default";
    $self->throw("'all' is a reserved word for context.") if $context eq 'all';
    $context ||= '_default';
    if ( ( $self->{_dblinks}->{$context} ) &&
        grep { $_->display_text eq $value } 
        @{ $self->{_dblinks}->{$context} } )
    {
        return TRUE;
    }
    else {
        return FALSE;
    }
}

=head2 add_dblink_context

  Title   : add_dblink_context
  Usage   : $term->add_dblink_context($db, $context);
  Function: add a dblink with its context
  Return  : [none]
  Args    : [arg1] a Bio::Annotation::DBLink instance
            [arg2] a string for context; if omitted, the
                   default/context-less one will be used.
  Note    : deprecated method due to past use of mixed data types; use
            add_dbxref() instead

=cut

sub add_dblink_context {
    my ($self, $value, $context) = @_;
    $self->deprecated("Use of simple strings and add_dblink_context() is deprecated; use\n
                      Bio::Annotation::DBLink instances and add_dbxref() instead");
    return $self->add_dbxref([$value],$context);
}

=head2 remove_dblinks

 Title   : remove_dblinks()
 Usage   : $term->remove_dblinks();
 Function: Deletes (and returns) the definition references of this GO term.
 Returns : A list of definition references [array of [scalars]].
 Args    : Context. If omitted or equal to 'all', all dblinks
           will be removed.
 Note    : deprecated method due to past use of mixed data types; use
           remove_dblinks() instead, which handles both strings and
           DBLink instances

=cut

sub remove_dblinks {
    my ($self, $context) = @_;
    $self->deprecated("use of remove_dblinks() is deprecated; use remove_dbxrefs() instead");
    return $self->remove_dbxrefs(@_);
} # remove_dblinks

=head2 remove_dbxrefs

 Title   : remove_dbxrefs()
 Usage   : $term->remove_dbxrefs();
 Function: Deletes (and returns) the definition references of this GO term.
 Returns : A list of definition references [array of [scalars]].
 Args    : Context. If omitted or equal to 'all', all dblinks
           will be removed.

=cut

sub remove_dbxrefs {
    my ($self, $context) = @_;
    $context = undef if $context && ($context eq "all");
    my @old = $self->get_dbxrefs($context);
    if (defined($context)) {
        $self->{_dblinks}->{$context}=[];
    } else {
        $self->{_dblinks} = {};
    }
    return @old;
} # remove_dbxrefs

=head2 get_references

  Title   : get_references
  Usage   : @references = $self->get_references
  Fuctnion: Returns a list of references
  Return  : A list of objects
  Args    : [none]

=cut

sub get_references {
    my $self=shift;
    return @{$self->{_references}} if exists $self->{_references};
    return ();
}

=head2 add_reference

  Title   : add_reference
  Usage   : $self->add_reference($reference);
            $self->add_reference($reference1, $reference2);
  Fuctnion: Add one or more references
  Returns : [none]

=cut

sub add_reference {
    my ($self, @values) =@_;
    return unless @values;
    # avoid duplicates and undefs
    foreach my $reference (@values){
        $self->throw("Passed data not an Bio::Annotation::Reference") unless ref $reference &&
            $reference->isa('Bio::AnnotationI');
        next unless defined $reference;
        next if grep{$_ eq $reference} @{$self->{_references}};
        push @{$self->{_references}}, $reference;
    }
}

=head2 remove_references

  Title   : remove_references
  Usage   : $self->remove_references;
  Function: Deletes (and returns) all references
  Returns : A list of references
  Args    : [none]

=cut

sub remove_references {
    my $self=shift;
    my @references=$self->get_references;
    $self->{_references}=[];
    return @references;
}

=head2 get_secondary_ids

 Title   : get_secondary_ids
 Usage   : @ids = $term->get_secondary_ids();
 Function: Returns a list of secondary identifiers of this Term.

           Secondary identifiers mostly originate from merging terms,
           or possibly also from splitting terms.

 Returns : A list of secondary identifiers [array of [scalar]]
 Args    :

=cut

sub get_secondary_ids {
    my $self = shift;

    return @{$self->{"_secondary_ids"}} if exists($self->{"_secondary_ids"});
    return ();
} # get_secondary_ids


=head2 add_secondary_id

 Title   : add_secondary_id
 Usage   : $term->add_secondary_id( @ids );
           or
           $term->add_secondary_id( $id );
 Function: Adds one or more secondary identifiers to this term.
 Returns :
 Args    : One or more secondary identifiers [scalars]

=cut

sub add_secondary_id {
    my $self = shift;

    return unless @_;

    # avoid duplicates
    foreach my $id (@_) {
        next if grep { !$_ or $_ eq $id; } @{$self->{ "_secondary_ids" }};
        push( @{ $self->{ "_secondary_ids" } }, $id );
    }

} # add_secondary_id


=head2 remove_secondary_ids

 Title   : remove_secondary_ids
 Usage   : $term->remove_secondary_ids();
 Function: Deletes (and returns) the secondary identifiers of this Term.
 Returns : The previous list of secondary identifiers [array of [scalars]]
 Args    :

=cut

sub remove_secondary_ids {
    my $self = shift;

    my @a = $self->get_secondary_ids();
    $self->{ "_secondary_ids" } = [];
    return @a;

} # remove_secondary_ids


# Title   :_is_true_or_false
# Function: Checks whether the argument is TRUE or FALSE.
# Returns :
# Args    : The value to be checked.
sub _is_true_or_false {
    my ( $self, $value ) = @_;
    unless ( $value !~ /\D/ && ( $value == TRUE || $value == FALSE ) ) {
        $self->throw( "Found [" . $value
        . "] where " . TRUE . " or " . FALSE . " expected" );
    }
} # _is_true_or_false

=head1 Methods implementing L<Bio::IdentifiableI> and L<Bio::DescribableI>

=cut

=head2 object_id

 Title   : object_id
 Usage   : $string    = $obj->object_id()
 Function: a string which represents the stable primary identifier
           in this namespace of this object.

           This is a synonym for identifier().

 Returns : A scalar

=cut

sub object_id {
    return shift->identifier(@_);
}

=head2 authority

 Title   : authority
 Usage   : $authority    = $obj->authority()
 Function: a string which represents the organisation which
           granted the namespace, written as the DNS name for
           organisation (eg, wormbase.org)

           This forwards to ontology()->authority(). Note that you
           cannot set the authority before having set the ontology or
           the namespace (which will set the ontology).

 Returns : A scalar
 Args    : on set, the new value (a scalar)

=cut

sub authority {
    my $self = shift;
    my $ont = $self->ontology();

    return $ont->authority(@_) if $ont;
    $self->throw("cannot manipulate authority prior to ".
                 "setting the namespace or ontology") if @_;
    return;
}


=head2 namespace

 Title   : namespace
 Usage   : $string    = $obj->namespace()
 Function: A string representing the name space this identifier
           is valid in, often the database name or the name
           describing the collection.

           This forwards to ontology() (set mode) and
           ontology()->name() (get mode). I.e., setting the namespace
           will set the ontology to one matching that name in the
           ontology store, or to one newly created.

 Returns : A scalar
 Args    : on set, the new value (a scalar)

=cut

sub namespace {
    my $self = shift;

    $self->ontology(@_) if(@_);
    my $ont = $self->ontology();
    return defined($ont) ? $ont->name() : undef;
}

=head2 display_name

 Title   : display_name
 Usage   : $string    = $obj->display_name()
 Function: A string which is what should be displayed to the user.

           The definition in Bio::DescribableI states that the
           string should not contain spaces. As this is not very
           sensible for ontology terms, we relax this here. The
           implementation just forwards to name().

 Returns : A scalar
 Args    : on set, the new value (a scalar)

=cut

sub display_name {
    return shift->name(@_);
}


=head2 description

 Title   : description
 Usage   : $string    = $obj->description()
 Function: A text string suitable for displaying to the user a
           description. This string is likely to have spaces, but
           should not have any newlines or formatting - just plain
           text.

           This forwards to definition(). The caveat is that the text
           will often be longer for ontology term definitions than the
           255 characters stated in the definition in
           Bio::DescribableI.

 Returns : A scalar
 Args    : on set, the new value (a scalar)

=cut

sub description {
    return shift->definition(@_);
}

#################################################################
# aliases or forwards to maintain backward compatibility
#################################################################

=head1  Deprecated methods

Used for looking up the methods that supercedes them.

=cut

sub each_dblink {shift->throw("use of each_dblink() is deprecated; use get_dbxrefs() instead")} 
sub add_dblinks {shift->throw("use of add_dblinks() is deprecated; use add_dbxref() instead")}
*each_synonym = \&get_synonyms;
*add_synonyms = \&add_synonym;

1;
