
#
# BioPerl module for Bio::Search::Hit::HitI
#
# Cared for by Aaron Mackey <amackey@virginia.edu>
#
# Copyright Aaron Mackey
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Search::Hit::HitI - Abstract interface and base class for Hit objects

=head1 SYNOPSIS

This is a base class and should not be used directly.

=head1 DESCRIPTION

    Bio::Search::Hit::* objects are data structures that contain information
about specific hits obtained during a library search.  Some information will
be algorithm-specific, but others will be generally defined.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this
and other Bioperl modules. Send your comments and suggestions preferably
 to one of the Bioperl mailing lists.
Your participation is much appreciated.

  vsns-bcd-perl@lists.uni-bielefeld.de          - General discussion
  vsns-bcd-perl-guts@lists.uni-bielefeld.de     - Technically-oriented discussion
  http://bio.perl.org/MailList.html             - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
 the bugs and their resolution.
 Bug reports can be submitted via email or the web:

  bioperl-bugs@bio.perl.org
  http://bio.perl.org/bioperl-bugs/

=head1 AUTHOR - Aaron Mackey

Email amackey@virginia.edu

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::Search::Hit::HitI;

use vars qw($AUTOLOAD @ISA);
use strict;

# Object preamble - inherits from Bio::Root::Object

use Bio::Root::Object;

@ISA = qw(Bio::Root::Object Exporter);

my @AUTOLOAD_OK = qw( _HIT_ID
		      _HIT_DESC
		      _ALGORITHM
		      _EVALUE
		    );
my %AUTOLOAD_OK = ();
@AUTOLOAD_OK{@AUTOLOAD_OK} = 1;

# new() is inherited from Bio::Root::Object

# _initialize is where the heavy stuff will happen when new is called

sub _initialize {
    my($self,@args) = @_;

    my $make = $self->SUPER::_initialize;

    # set stuff in self from @args
    while (@args > 1) {
	my ($arg, $val) = splice(@args, 0, 2);
	$arg =~ s/^-?(.*)/_\U$1\E/; # we require hash parameters, no arrays!
        $self->$arg($val);
    }

    $self->throw("Ill-defined arguments to Hit's new() function!") if @args;

    return $make; # success - we hope!
}

=head2 get_hit_id

 Title   : get_hit_id
 Usage   : $id = $hit->get_hit_id();
 Function: Used to obtain the id of the matched entity.
 Returns : a scalar string
 Args    : <none>

=cut

sub get_hit_id{
    my ($self,@args) = @_;

    return $self->_HIT_ID();
}

=head2 get_hit_desc

 Title   : get_hit_desc
 Usage   : $desc = $hit->get_hit_desc();
 Function: Used to obtain the description of a matched entity
 Returns : a scalar string
 Args    : <none>

=cut

sub get_hit_desc{
    my ($self,@args) = @_;

    return $self->_HIT_DESC();
}

=head2 get_algorithm

 Title   : get_algorithm
 Usage   : $which = $hit->get_algorithm();
 Function: Used to obtain the algorithm specification used to obtain the hit
 Returns : a scalar string 
 Args    : <none>

=cut

sub get_algorithm{
    my ($self,@args) = @_;

    return $self->_ALGORITHM();
}

=head2 get_expectation_value

 Title   : get_expectation_value
 Usage   : $evalue = $hit->get_expectation_value();
 Function: Used to obtain the E() value of a hit, i.e. the probability that
           this particular hit was obtained purely by random chance.
 Returns : a scalar value
 Args    : <none>

=cut

sub get_expectation_value{
    my ($self,@args) = @_;

    return $self->_EVALUE();
}

=head2 get_evalue

 Title   : get_evalue
 Usage   : $evalue = $hit->get_evalue();
 Function: same as get_expectation_value (just an alias)
 Returns : see above
 Args    : <none>


=cut

sub get_evalue{
    my ($self,@args) = @_;

    return $self->get_expectation_value(@args);
}

=head2 get_alignments

 Title   : get_alignments
 Usage   : @aligns = $hit->get_alignments()
 Function: Used to obtain an array of Bio::Alignment objects corresponding to
           the Hit.
 Returns : an array of Alignment objects, if more than one Alignments exist and
           called in wantarray context.  In scalar context will return a
           reference to an array of Alignment objects, if more than one
           Alignment objects exist, otherwise will return the unique
           Alignment object.
 Args    :


=cut

sub get_alignments{
    my ($self,@args) = @_;

    return unless defined wantarray;

    my @align = @{$self->{'_ALIGNMENTS'}};
    return wantarray ? @align : ( @align > 1 ? \@align : $align[0] );
} 

=head2 _add_alignment

 Title   : _add_alignment
 Usage   : $hit->_add_alignment($alignobject);
 Function: Internal method, used by Processor objects to load Alignment objects
           into Hit objects.
 Returns : reference to new array of alignment objects
 Args    : a Bio::Alignment object

=cut

sub _add_alignment{
    my ($self, $alignment) = @_;

    my @align = @{$self->{'_ALIGNMENTS'}} || ();
    push @align, $alignment;
    return $self->{'_ALIGNMENTS'} = \@align;
}

sub AUTOLOAD {
    my ($self, $val) = @_;

    if ( $AUTOLOAD_OK{$AUTOLOAD} ) {
        $self->{$AUTOLOAD} = $val if defined $val;
        return $self->{$AUTOLOAD};
    } else {
        my $check = $self->SUPER::AUTOLOAD($val);
        $self->throw("Failed accessor: $AUTOLOAD !") unless defined $check;
        return $check;
    }
}

1;
