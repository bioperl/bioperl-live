#
# BioPerl module for Bio::Phenotype::MeSH::Twig
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Heikki Lehvaslaiho, heikki-at-bioperl-dot-org
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Phenotype::MeSH::Twig - Context for a MeSH term

=head1 SYNOPSIS

  use Bio::Phenotype::MeSH::Twig
  # create a twig object
  my $twig = Bio::Phenotype::MeSH::Twig->new();

  # the term has only one parent in any twig
  $twig->parent('Fats');


  # a twig makeas sense only in the context of a term
  # which is a  Bio::Phenotype::MeSH::Term object

  # a term can have many twigs i.e. it can appear in many places in
  # the hierarchy
  #
  $ term->add_twig($twig);

  # adding the twig into a term adds a link into into it 
  $twig->term eq $term;

  # a twig can know about other terms under the parant node
  $twig->add_sister('Bread', 'Candy', 'Cereals');
  print join ( ', ', $twig->each_sister()), "\n";

  # a twig can know about other terms under this term
  $twig->add_child('Butter', 'Margarine');
  print join ( ', ', $twig->each_child()), "\n";



=head1 DESCRIPTION

This class represents the immediate surrounding of a MeSH term. It
keeps track on nodes names above the current node ('parent') other
nodes at the same level ('sisters') and nodes under it ('children').
Note that these are name strings, not objects.

Each twig can be associated with only one term, but term can have
multiple twigs. (Twigs can be though to be roles for a term.)

=head1 SEE ALSO

L<Bio::Phenotype::MeSH::Term>

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

package Bio::Phenotype::MeSH::Twig;
use strict;


use base qw(Bio::Root::Root);


sub new {

    my( $class,@args ) = @_;
    my $self = $class->SUPER::new( @args );

    my ($term, $parent ) = $self->_rearrange
        ( [ qw(
               TERM
               PARENT
              ) ],
          @args );

    $self->{"_children"} = [];
    $self->{"_sisters"} = [];

    $term && $self->term($term );
    $parent  && $self->parent($parent );
    return $self;
}


=head2 parent

 Title   : parent
 Usage   : $obj->parent( "r1" );
           or
           print $obj->parent();
 Function: Set/get for the parent.
 Returns : A parent [scalar].
 Args    : A parent [scalar] (optional).

=cut

sub parent {
    my ( $self, $value ) = @_;
    $self->{ "_parent" } = $value if defined $value;
    return $self->{ "_parent" };
}

=head2 term

 Title   : term
 Usage   : $obj->term( "r1" );
           or
           print $obj->term();
 Function: Set/get for the term.
 Returns : A term [scalar].
 Args    : A term [scalar] (optional).

=cut

sub term {
    my ( $self, $value ) = @_;
    if (defined $value) {
        $self->throw ("Not a MeSH term [$value]")
            unless $value->isa('Bio::Phenotype::MeSH::Term');
        $self->{ "_term" } = $value
    }
    return $self->{ "_term" };
}


=head2 add_child

 Title   : add_child
 Usage   : $obj->add_child( @children );
           or
           $obj->add_child( $child );
 Function: Pushes one or more child term names [scalars, most likely Strings]
           into the list of children.
 Returns : 
 Args    : scalar(s).

=cut

sub add_child {
    my ( $self, @values ) = @_;
    push( @{ $self->{ "_children" } }, @values );
    return scalar @values;
}

=head2 each_child

 Title   : each_child()
 Usage   : @gs = $obj->each_child();
 Function: Returns a list of gene symbols [scalars, most likely Strings]
           associated with this phenotype.
 Returns : A list of scalars.
 Args    :

=cut

sub each_child {
    my ( $self ) = shift;
    return @{ $self->{ "_children" } };
}

=head2 purge_children

 Usage   : $obj->purge_child();
 Function: Deletes  the list of children associated with this term.
 Returns : A list of scalars.
 Args    :

=cut

sub purge_children {
    my ( $self ) = @_;
    $self->{ "_children" } = [];
}


=head2 add_sister

 Title   : add_sister
 Usage   : $obj->add_sister( @sisters );
           or
           $obj->add_sister( $sister );
 Function: Pushes one or more sister term names [scalars, most likely Strings]
           into the list of sisters.
 Returns : 
 Args    : scalar(s).

=cut

sub add_sister {
    my ( $self, @values ) = @_;
    push( @{ $self->{ "_sisters" } }, @values );
    return scalar @values;
}

=head2 each_sister

 Title   : each_sister()
 Usage   : @gs = $obj->each_sister();
 Function: Returns a list of gene symbols [scalars, most likely Strings]
           associated with this phenotype.
 Returns : A list of scalars.
 Args    :

=cut

sub each_sister {
    my ( $self ) = shift;
    return @{ $self->{ "_sisters" } };
}

=head2 purge_sisters

 Usage   : $obj->purge_sister();
 Function: Deletes  the list of sisters associated with this term.
 Returns : A list of scalars.
 Args    :

=cut

sub purge_sisters {
    my ( $self ) = @_;
    $self->{'_sisters'} = [];
}

1;
