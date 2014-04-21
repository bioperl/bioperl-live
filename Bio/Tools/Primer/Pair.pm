# BioPerl module for Bio::Tools::Primer::Pair
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Ewan Birney <birney@ebi.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Tools::Primer::Pair - two primers on left and right side

=head1 SYNOPSIS

    use Bio::Tools::Primer::Pair;

    my $pair = Bio::Tools::Primer::Pair->new( -left => $leftp , -right => $rightp);

    # helper functions

    print "GC percentage different",$pf->gc_difference(),"\n";
    print "product length is ",$pf->product_length,"\n";



=head1 DESCRIPTION

Primer Pairs represents one primer in a primer pair. This object is mainly for
designing primers, and probably principly used in the primer design system

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
of the bugs and their resolution. Bug reports can be submitted via the
web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Ewan Birney

Email birney-at-ebi.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...



package Bio::Tools::Primer::Pair;

use base qw(Bio::Root::Root);

sub new {
    my ( $caller, @args) = @_;   
    my ($self) = $caller->SUPER::new(@args); 

    my ($left,$right) = $self->_rearrange([qw(LEFT RIGHT)],@args);

    if( !defined $left || !defined $right ) {
	$self->throw("Pair must be initialised with left and right primers");
    }

    $self->left($left);
    $self->right($right);

    # done - we hope
    return $self;
}

sub left {
    my $self = shift;
    my $left = shift;

    if( defined $left ) {
	if( !ref $left || !$left->isa("Bio::Tools::Primer::Feature") ) {
	    $self->throw("left primer must be a Bio::Tools::Primer::Feature, not $left");
	}
	$self->{'left'} = $left;
    }

    return $self->{'left'};
}


sub right {
    my $self = shift;
    my $right = shift;

    if( defined $right ) {
	if( !ref $right || !$right->isa("Bio::Tools::Primer::Feature") ) {
	    $self->throw("right primer must be a Bio::Tools::Primer::Feature, not $right");
	}
	$self->{'right'} = $right;
    }

    return $self->{'right'};
}

sub gc_difference {
    my $self = shift;

    return abs ( $self->left->gc_percent - $self->right->gc_percent );
}

sub product_length {
    my $self = shift;

    return $self->right->end - $self->left->start +1;
}


1;
