# $Id$
#
# BioPerl module for Bio::Tree::NodeNHX
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
    my $nodeA = new Bio::Tree::NodeNHX();
    my $nodeL = new Bio::Tree::NodeNHX();
    my $nodeR = new Bio::Tree::NodeNHX();

    my $node = new Bio::Tree::NodeNHX();
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

  bioperl-l@bioperl.org              - General discussion
  http://bioperl.org/MailList.shtml  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution. Bug reports can be submitted via
email or the web:

  bioperl-bugs@bioperl.org
  http://bioperl.org/bioperl-bugs/

=head1 AUTHOR - Aaron Mackey

Email amackey@virginia.edu

=head1 CONTRIBUTORS

The NHX (New Hampshire eXtended) format was created by Chris Zmasek,
and is described at:

  http://www.genetics.wustl.edu/eddy/forester/NHX.html

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...

package Bio::Tree::NodeNHX;
use vars qw(@ISA);
use strict;

use Bio::Tree::Node;

@ISA = qw(Bio::Tree::Node);

=head2 new

 Title   : new
 Usage   : my $obj = new Bio::Tree::NodeNHX();
 Function: Builds a new Bio::Tree::NodeNHX object
 Returns : Bio::Tree::NodeNHX
 Args    : -left   => pointer to Left descendent (optional)
           -right  => pointer to Right descenent (optional)
	   -branch_length => branch length [integer] (optional)
           -bootstrap => value   bootstrap value (string)
           -desc      => description of node
           -id        => unique id for node
           -nhx       => hashref of NHX tags and values
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
   return sprintf("%s%s%s",
		  defined $self->id ? $self->id : '',
		  defined $self->branch_length ? ':' . $self->branch_length : ' ',
		  '[' . join(":", "&&NHX", map { "$_=" . $self->nhx_tag($_) } keys %{$self->nhx_tag || {}}) . ']'
		 );
}

sub nhx_tag {

    my ($self, $tags) = @_;

    if (defined $self->bootstrap) {
	$self->{_nhx_tag}->{B} = $self->bootstrap;
    }

    if (defined $tags && (ref $tags eq 'HASH')) {
	$self->{_nhx_tag} = $tags;
	if (exists $self->{_nhx_tag}->{B}) {
	    if (defined $self->bootstrap &&
		($self->bootstrap != $self->{_nhx_tag}) ) {
		$self->warn("bootstrap value (" . $self->bootstrap . ") being overwritten by NHX B: value ($self->{_nhx_tag}->{B})!");
	    }
	    $self->bootstrap($self->{_nhx_tag}->{B});
	}
    } elsif (defined $tags and ! ref $tags) {
	return $self->{_nhx_tag}->{$tags};
    }
        return $self->{_nhx_tag};
}

1;

