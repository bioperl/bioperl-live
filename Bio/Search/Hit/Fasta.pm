
#
# BioPerl module for Bio::Search::Hit::Fasta
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

    These objects are generated automatically by Bio::Search::Processor::Fasta,
and shouldn't be used directly.


=head1 DESCRIPTION

    Bio::Search::Hit::* objects are data structures that contain information
about specific hits obtained during a library search.  Some information will
be algorithm-specific, but others will be generally defined, such as the
ability to obtain alignment objects corresponding to each hit.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this
and other Bioperl modules. Send your comments and suggestions preferably
 to one of the Bioperl mailing lists.
Your participation is much appreciated.

  bioperl-l@bioperl.org             - General discussion
  http://bio.perl.org/MailList.html - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via email
or the web:

  bioperl-bugs@bio.perl.org
  http://bugzilla.bioperl.org/

=head1 AUTHOR - Aaron Mackey

Email amackey@virginia.edu

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

#'
# Let the code begin...

package Bio::Search::Hit::Fasta;

use vars qw($AUTOLOAD @ISA);
use strict;

# Object preamble - inherits from Bio::Root::Object

use Bio::Search::Hit::HitI;

@ISA = qw(Bio::Search::Hit::HitI);

my @AUTOLOAD_OK = qw(        _ID
                             _DESC
                             _SIZE
                             _INITN
                             _INIT1
                             _OPT
                             _ZSC
                             _E_VAL
                    );

my %AUTOLOAD_OK = ();
@AUTOLOAD_OK{@AUTOLOAD_OK} = (1) x @AUTOLOAD_OK;

# new() is inherited from Bio::Root::Object

# _initialize is where the heavy stuff will happen when new is called

sub _initialize {
    my($self, %args) = @_;

    my $make = $self->SUPER::_initialize(%args);

    while (my ($key, $val) = each %args) {
	$key = '_' . uc($key);
	$self->$key($val);
    }

    return $make; # success - we hope!
}

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

__END__
