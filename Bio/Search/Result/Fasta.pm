
#
# BioPerl module for Bio::Search::Result::Fasta
#
# Cared for by Aaron Mackey <amackey@virginia.edu>
#
# Copyright Aaron Mackey
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Search::Result::Fasta - Result object for FASTA-generated data sources

=head1 SYNOPSIS

    These objects are generated automatically by Bio::Search::Processor::Fasta,
and shouldn't be used directly.

=head1 DESCRIPTION

Bio::Search::Result::* objects are data structures containing the results from
the execution of a search algorithm.  As such, it may contain various
algorithm specific information as well as details of the execution, but will
contain a few fundamental elements, including the ability to return
Bio::Search::Hit::* objects.

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


package Bio::Search::Result::Fasta;

use vars qw($AUTOLOAD @ISA);
use strict;

# Object preamble - inherits from Bio::Root::Object

use Bio::Search::Result::ResultI;

@ISA = qw(Bio::Search::Result::ResultI);

my @AUTOLOAD_OK = qw(
                         _INTERACTIVE
                         _ALGORITHM_DESC
                         _ALGORITHM
                         _VERSION
                         _VERSION_DATE
                         _CITATION
                         _QUERY_FILENAME
                         _QUERY_START
                         _QUERY_END
                         _QUERY_SIZE
                         _QUERY_TYPE
                         _QUERY_SUPERFAMILIES
                         _QUERY_DESC
                         _LIBRARY_FILENAME
                         _LIBRARY_RESIDUES
                         _LIBRARY_SEQUENCES
                         _STATISTICS
                         _OPTIMIZED
                         _MATRIX_NAME
                         _MATRIX_OFFSET
                         _MATRIX_HIGH_SCORE
                         _MATRIX_LOW_SCORE
                         _KTUP
                         _JOIN
                         _OPT
                         _GAP_OPEN
                         _GAP_EXTEND
                         _WIDTH
                         _SCANTIME
                         _DISPLAY_TIME
                         _START_TIME
                         _END_TIME
                         _APP_LIB_SIZE
                         _HITS
                         _COMPLIB
                         _COMPLIB_VERSION
                         _COMPLIB_VERSION_DATE
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
        $self->{$AUTOLOAD} = $val if defined $val &&
	    caller eq 'Bio::Search::Result::Fasta';
        return $self->{$AUTOLOAD};
    } else {
        $self->throw("Unallowed accessor: $AUTOLOAD !");
    }
}

1;

__END__
