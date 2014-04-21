#
# BioPerl module for Bio::Phenotype::MeSH::Term
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Heikki Lehvaslaiho, heikki-at-bioperl-dot-org
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Phenotype::MeSH::Term - A MeSH term

=head1 SYNOPSIS

  use Bio::Phenotype::MeSH::Term;

  # create a term object
  my $term = Bio::Phenotype::MeSH::Term->new
      (-id => 'D000001',
       -name => 'Dietary Fats',
       -description => 'dietary fats are...'
      );

  # get a Bio::Phenotype::MeSH::Twig somehow...
  $term->add_twig($twig1);


=head1 DESCRIPTION

This class keeps information about MeSH terms. MeSH stands for Medical
Subject Headings and is one of the ways for annotaing biomedical
literature.  The terminology is maintained by National Library of
Medicine of USA . See http://www.nlm.nih.gov/mesh/meshhome.html.

In addition to id, name and description a term can know about its
surrounding terms (Bio::Phenotype::MeSH::Twig) in the term hierarchy.

This class is mainly used from Bio::DB::MeSH which retrieves terms
over the Web.

=head1 SEE ALSO

L<Bio::DB::MeSH>, 
L<Bio::Phenotype::MeSH::Twig>

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

report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via the
web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR

Heikki Lehvaslaiho, heikki-at-bioperl-dot-org

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Phenotype::MeSH::Term;
use strict;


use base qw(Bio::Root::Root);

sub new {

    my( $class,@args ) = @_;
    my $self = $class->SUPER::new( @args );

    my ( $id, $name, $description, $comment ) = $self->_rearrange
        ( [ qw( ID
                NAME
                DESCRIPTION
                SPECIES
                COMMENT
              ) ],
          @args );

    $self->{"_twigs"} = [];

    $id            && $self->id( $id );
    $name          && $self->name( $name );
    $description   && $self->description( $description );

    return $self;
}


=head2 id

 Title   : id
 Usage   : $obj->id( "r1" );
           or
           print $obj->id();
 Function: Set/get for the id.
 Returns : A id [scalar].
 Args    : A id [scalar] (optional).

=cut

sub id {
    my ( $self, $value ) = @_;
    $self->{ "_id" } = $value if defined $value;
    return $self->{ "_id" };
}

=head2 name

 Title   : name
 Usage   : $obj->name( "r1" );
           or
           print $obj->name();
 Function: Set/get for the name.
 Returns : A name [scalar].
 Args    : A name [scalar] (optional).

=cut

sub name {
    my ( $self, $value ) = @_;
    $self->{ "_name" } = $value if defined $value;
    return $self->{ "_name" };
}

=head2 description

 Title   : description
 Usage   : $obj->description( "r1" );
           or
           print $obj->description();
 Function: Set/get for the description.
 Returns : A description [scalar].
 Args    : A description [scalar] (optional).

=cut

sub description {
    my ( $self, $value ) = @_;
    $self->{ "_description" } = $value if defined $value;
    return $self->{ "_description" };
}


=head2 add_synonym

 Title   : add_synonym
 Usage   : $obj->add_synonym( @synonyms );
           or
           $obj->add_synonym( $synonym );
 Function: Pushes one or more synonyms for the term  term
           into the list of synonyms.
 Returns : 
 Args    : scalar(s).

=cut

sub add_synonym {
    my ( $self, @values ) = @_;
    push( @{ $self->{ "_synonyms" } }, @values );
}

=head2 each_synonym

 Title   : each_synonym()
 Usage   : @gs = $obj->each_synonym();
 Function: Returns a list of gene symbols [scalars, most likely Strings]
           associated with this phenotype.
 Returns : A list of scalars.
 Args    :

=cut

sub each_synonym {
    my ( $self ) = shift;
    return @{ $self->{ "_synonyms" } };
}

=head2 purge_synonyms

 Usage   : $obj->purge_synonym();
 Function: Deletes  the list of synonyms to this term.
 Returns : A list of scalars.
 Args    :

=cut

sub purge_synonyms {
    my ( $self ) = @_;
    $self->{ "_synonyms" } = [];
}


=head2 Twig management

Each MeSH term belongs to a complex tree like hierarchy of terms where
each term can appear multiple times. The immediately surrounding nodes
of the tree are modelled in twigs.

See: L<Bio::Phenotype::MeSH::Twig>.

=cut

=head2 add_twig

 Title   : add_twig
 Usage   : $obj->add_twig( @twigs );
           or
           $obj->add_twig( $twig );
 Function: Pushes one or more twig term names [scalars, most likely Strings]
           into the list of twigs.
 Returns : 
 Args    : scalar(s).

=cut

sub add_twig {
    my ( $self, @values ) = @_;
    foreach my $twig (@values) {
        $self->warn ("Not a MeSH twig [$twig]")
            unless $twig->isa('Bio::Phenotype::MeSH::Twig');
        $twig->term($self);
        push( @{ $self->{ "_twigs" } }, $twig );
    }
    1;
}

=head2 each_twig

 Title   : each_twig()
 Usage   : @gs = $obj->each_twig();
 Function: Returns a list of gene symbols [scalars, most likely Strings]
           associated with this phenotype.
 Returns : A list of scalars.
 Args    :

=cut

sub each_twig {
    my ( $self ) = shift;
    return @{ $self->{ "_twigs" } };
}

=head2 purge_twigs

 Usage   : $obj->purge_twig();
 Function: Deletes  the list of twigs associated with this term.
 Returns : A list of scalars.
 Args    :

=cut

sub purge_twigs {
    my ( $self ) = @_;
    $self->{ "_twigs" } = [];
}


=head2 each_parent

 Title   : each_parent()
 Usage   : @gs = $obj->each_parent();
 Function: Returns a list of names of parents for this term
 Returns : A list of scalars.
 Args    :

=cut

sub each_parent {
    my ( $self ) = shift;
    return map {$_->parent()} @{ $self->{ "_twigs" } };
}

1;
