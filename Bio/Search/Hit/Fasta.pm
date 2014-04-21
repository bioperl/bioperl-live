#
# BioPerl module for Bio::Search::Hit::Fasta
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

Bio::Search::Hit::Fasta - Hit object specific for Fasta-generated hits

=head1 SYNOPSIS

  # You wouldn't normally create these manually; 
  # instead they would be produced by Bio::SearchIO::fasta

  use Bio::Search::Hit::Fasta;
  my $hit = Bio::Search::Hit::Fasta->new(id=>'LBL_6321', desc=>'lipoprotein', e_val=>0.01);

=head1 DESCRIPTION

L<Bio::Search::Hit::HitI> objects are data structures that contain information
about specific hits obtained during a library search.  Some information will
be algorithm-specific, but others will be generally defined, such as the
ability to obtain alignment objects corresponding to each hit.

=head1 SEE ALSO

L<Bio::Search::Hit::HitI>,
L<Bio::Search::Hit::GenericHit>,
L<Bio::SearchIO::fasta>.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this
and other Bioperl modules. Send your comments and suggestions preferably
 to one of the Bioperl mailing lists.
Your participation is much appreciated.

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
the bugs and their resolution.  Bug reports can be submitted via the
web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Aaron Mackey

Email amackey-at-virginia.edu

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

#'
# Let the code begin...

package Bio::Search::Hit::Fasta;

use vars qw($AUTOLOAD);
use strict;

use base qw(Bio::Search::Hit::HitI);

my @AUTOLOAD_OK = qw(_ID _DESC _SIZE _INITN _INIT1 _OPT _ZSC _E_VAL);

my %AUTOLOAD_OK = ();
@AUTOLOAD_OK{@AUTOLOAD_OK} = (1) x @AUTOLOAD_OK;

=head2 _initialize

 Function: where the heavy stuff will happen when new is called

=cut

sub _initialize {
    my($self, %args) = @_;

    my $make = $self->SUPER::_initialize(%args);

    while (my ($key, $val) = each %args) {
	$key = '_' . uc($key);
	$self->$key($val);
    }

    return $make; # success - we hope!
}

=head2 AUTOLOAD

 Function: Provide getter/setters for ID,DESC,SIZE,INITN,INIT1,OPT,ZSC,E_VAL

=cut

sub AUTOLOAD {
    my ($self, $val) = @_;

    $AUTOLOAD =~ s/.*:://;

    if ( $AUTOLOAD_OK{$AUTOLOAD} ) {
        $self->{$AUTOLOAD} = $val if defined $val;
        return $self->{$AUTOLOAD};
    } else {
        $self->throw("Unallowed accessor: $AUTOLOAD !");
    }
}

1;

