# $Id: NodePhyloXML.pm 11508 2007-06-23 01:38:32Z jason $
#
# BioPerl module for Bio::Tree::NodePhyloXML
#
# Cared for by Mira Han <mirhan@indiana.edu>
#
# Copyright Mira Han
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Tree::NodePhyloXML - A Simple Tree Node with support for PhyloXML tags

=head1 SYNOPSIS

    use Bio::Tree::NodePhyloXML;
    my $nodeA = Bio::Tree::NodePhyloXML->new();
    my $nodeL = Bio::Tree::NodePhyloXML->new();
    my $nodeR = Bio::Tree::NodePhyloXML->new();

    my $node = Bio::Tree::NodePhyloXML->new();
    $node->add_Descendents($nodeL);
    $node->add_Descendents($nodeR);

    print "node is not a leaf \n" if( $node->is_leaf);

=head1 DESCRIPTION

Makes a Tree Node with PhyloXML tags, suitable for building a Tree.  See
L<Bio::Tree::Node> for a full list of functionality.

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

=head1 AUTHOR - Mira Han

Email mirhan@indiana.edu

=head1 CONTRIBUTORS

The PhyloXML format was created by Chris Zmasek,
and is described at:

  http://www.phyloxml.org/

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...

package Bio::Tree::NodePhyloXML;
use strict;


use base qw(Bio::Tree::Node);

=head2 new

 Title   : new
 Usage   : my $obj = Bio::Tree::NodePhyloXML->new();
 Function: Builds a new Bio::Tree::NodePhyloXML object
 Returns : Bio::Tree::NodePhyloXML
 Args    : -left          => pointer to Left descendent (optional)
           -right         => pointer to Right descenent (optional)
	         -branch_length => branch length [integer] (optional)
           -bootstrap     => bootstrap value (string)
           -description   => description of node
           -id            => unique id for node
           -user_tag      => hashref of PhyloXML tags and values

=cut

sub new {
  my ($class,@args) = @_;
  my $self = $class->SUPER::new(@args);
  $self->debug("new NodePhyloXML\n");
  my ($user_tag) = $self->_rearrange([qw(PhyloXML)], @args);
  $self->_tag($user_tag);
  return $self;
}

sub DESTROY {
    my ($self) = @_;
    # try to insure that everything is cleaned up
    $self->SUPER::DESTROY();
    if( defined $self->{'_desc'} &&
	ref($self->{'_desc'}) =~ /ARRAY/i ) {
	while( my ($nodeid,$node) = each %{ $self->{'_desc'} } ) {
	    $node->{'_ancestor'} = undef; # insure no circular references
	    $node->DESTROY();
	    $node = undef;
	}
	$self->{'_desc'} = {};
    }
}

sub to_string{
   my ($self) = @_;
   my @tags = $self->get_all_tags;
   my $tagstr = '';
   if( @tags ) {
       $tagstr = '[' . join(":", "&&PhyloXML", 
			    map { "$_=" .join(',',
					      $self->get_tag_values($_))}
			    @tags ) . ']';
   }
   return sprintf("%s%s%s",
		  defined $self->id ? $self->id : '',
		  defined $self->branch_length ? ':' . 
		  $self->branch_length : ' ',
		  $tagstr);
}

=head2 _tag

 Title   : _tag
 Usage   : my $tag = $nodephyloXML->_tag(%tags);
 Function: Set tag-value pairs for PhyloXML nodes
 Returns : none
 Args    : hashref to update the tags/value pairs
           OR 
           with a scalar value update the bootstrap value by default


=cut

sub _tag 
{
  my ($self, $tags) = @_;
  if (defined $tags && (ref($tags) =~ /HASH/i)) 
  {
    while( my ($tag,$val) = each %$tags ) 
    {
      if( ref($val) =~ /ARRAY/i ) 
      {
        for my $v ( @$val ) 
        {
          $self->add_tag_value($tag,$v);
        }
      } 
      else {
        $self->add_tag_value($tag,$val);
      }
    }
    if (exists $tags->{'B'}) 
    {
      $self->bootstrap($tags->{'B'});
    }
  } 
  elsif (defined $tags and ! ref ($tags)) 
  {
# bootstrap by default
    $self->bootstrap($tags);
  }
}

1;
