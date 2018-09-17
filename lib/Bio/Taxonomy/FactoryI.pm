#
#
# BioPerl interface of Bio::Taxnomoy::FactoryI
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Juguang Xiao
#
# You may distribute this module under the same terms as Perl itself
#
# POD documentation - main does before the code

=head1 NAME

Bio::Taxonomy::FactoryI - interface to define how to access NCBI Taxonoy

=head1 DESCRIPTION

NB: This module has been deprecated.

$factory-E<gt>fetch is a general method to fetch Taxonomy by either NCBI
taxid or any types of names.

$factory-E<gt>fetch_parent($taxonomy), returns a Taxonomy that is
one-step higher rank of the taxonomy specified as argument.

$factory-E<gt>fetch_children($taxonomy), reports an array of Taxonomy
those are one-step lower rank of the taxonomy specified as the
argument.

=head1 AUTHOR - Juguang Xiao

juguang@tll.org.sg

=head1 CONTRIBUTORS

Additional contributors' names and emails here

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

package Bio::Taxonomy::FactoryI;
use strict;


use base qw(Bio::Root::Root);

=head2 fetch

  Title:    fetch
  Usage:    my $taxonomy = $factory->fetch(-taxon_id => 9605);
            my $taxonomy = $factory->fetch(-common_name => 'mammals');
  Fuctnion: Fetch taxonomy by taxon_id, common name or scientific name.
  Returns:  an instance of Bio::Taxonomy
  Args:     -taxon_id => NCBI taxonomy ID
            -common_name => comon name, such as 'human', 'mammals'
            -scientifc_name => specitic name, such as 'sapiens', 'Mammalia'

=cut

sub fetch {
    shift->throw_not_implemented;
}

=head2 fuzzy_fetch

  Title:    fuzzy_fetch
  Usage:    my @taxonomy = $factory->fuzzy_fetch(-name => 'mouse');
  Function: Fuzzy fetch by name, or any text information found in DB
  Returns:  an array reference of Bio::Taxonomy objects
  Args:     -name => any name, such as common name, variant, scientific name
            -description, or -desc => any text information

=cut

sub fuzzy_fetch {
    shift->throw_not_implemented;
}

=head2 fetch_parent

  Title:    fetch_parent
  Usage:    my $parent_taxonomy = $factory->fetch_parent($taxonomy);
  Function: Fetch the parent that is one-rank higher than the argument.
  Returns:  an instance of Bio::Taxonomy, or undef if the arg is the top one.
  Args:     a Bio::Taxonomy object.

=cut

sub fetch_parent {
    shift->throw_not_implemented;
}

=head2 fetch_children

  Title:    fetch_children
  Usage:    my @children_taxonomy = $factory->fetch_children($taxonomy);
  Function: Fetch all children those are one-rank lower than the argument.
  Returns:  an array reference of Bio::Taxonomy objects
  Args:     a Bio::Taxonomy object.

=cut

sub fetch_children {
    shift->throw_not_implemented;
}

1;
