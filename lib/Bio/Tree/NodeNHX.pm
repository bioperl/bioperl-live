#
# BioPerl module for Bio::Tree::NodeNHX
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Aaron Mackey <amackey@virginia.edu>
#
# Copyright Aaron Mackey
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Tree::NodeNHX - A Simple Tree Node with support for NHX tags

=head1 SYNOPSIS

    use Bio::Tree::NodeNHX;
    my $nodeA = Bio::Tree::NodeNHX->new();
    my $nodeL = Bio::Tree::NodeNHX->new();
    my $nodeR = Bio::Tree::NodeNHX->new();

    my $node = Bio::Tree::NodeNHX->new();
    $node->add_Descendents($nodeL);
    $node->add_Descendents($nodeR);

    print "node is not a leaf \n" if( $node->is_leaf);

=head1 DESCRIPTION

Makes a Tree Node with NHX tags, suitable for building a Tree.  See
L<Bio::Tree::Node> for a full list of functionality.

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

=head1 AUTHOR - Aaron Mackey

Email amackey@virginia.edu

=head1 CONTRIBUTORS

The NHX (New Hampshire eXtended) format was created by Chris Zmasek,
and is described at:

  http://sourceforge.net/projects/forester-atv/

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...

package Bio::Tree::NodeNHX;
use strict;


use base qw(Bio::Tree::Node);

=head2 new

 Title   : new
 Usage   : my $obj = Bio::Tree::NodeNHX->new();
 Function: Builds a new Bio::Tree::NodeNHX object
 Returns : Bio::Tree::NodeNHX
 Args    : -left          => pointer to Left descendent (optional)
           -right         => pointer to Right descenent (optional)
	   -branch_length => branch length [integer] (optional)
           -bootstrap     => bootstrap value (string)
           -description   => description of node
           -id            => unique id for node
           -nhx           => hashref of NHX tags and values

=cut

sub new {
  my($class,@args) = @_;

  my $self = $class->SUPER::new(@args);
  my ($nhx) = $self->_rearrange([qw(NHX)], @args);
  $self->nhx_tag($nhx);
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
   if( scalar(@tags) > 0 ) {
       $tagstr = '[' . join(":", "&&NHX", 
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

=head2 nhx_tag

 Title   : nhx_tag
 Usage   : my $tag = $nodenhx->nhx_tag(%tags);
 Function: Set tag-value pairs for NHX nodes
 Returns : none
 Args    : hashref to update the tags/value pairs
           OR 
           with a scalar value update the bootstrap value by default


=cut

sub nhx_tag {
    my ($self, $tags) = @_;
    if (defined $tags && (ref($tags) =~ /HASH/i)) {
	while( my ($tag,$val) = each %$tags ) {
	    if( ref($val) =~ /ARRAY/i ) {
		for my $v ( @$val ) {
		    $self->add_tag_value($tag,$v);
		}
	    } else {
		$self->add_tag_value($tag,$val);
	    }
	}
	if (exists $tags->{'B'}) {
	    $self->bootstrap($tags->{'B'});
	}
    } elsif (defined $tags and ! ref ($tags)) {
	$self->debug( "here with $tags\n");
        # bootstrap by default
	$self->bootstrap($tags);
    }
}

1;
