#
# $Id$
#
# BioPerl module for Bio::DB::Flat::BDB
#
# Cared for by Lincoln Stein <lstein@cshl.org>
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::DB::Flat::BDB::fasta - fasta adaptor for Open-bio standard BDB-indexed flat file

=head1 SYNOPSIS

See Bio::DB::Flat.

=head1 DESCRIPTION

This module allows fasta files to be stored in Berkeley DB flat files
using the Open-Bio standard BDB-indexed flat file scheme.  You should
not be using this directly, but instead use it via Bio::DB::Flat.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists.  Your participation is much appreciated.

  bioperl-l@bioperl.org             - General discussion
  http://bioperl.org/MailList.shtml - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via
email or the web:

  bioperl-bugs@bio.perl.org
  http://bio.perl.org/bioperl-bugs/

=head1 SEE ALSO

L<Bio::DB::Flat>,

=head1 AUTHOR - Lincoln Stein

Email - lstein@cshl.org

=cut

package Bio::DB::Flat::BDB::fasta;

use strict;
use Bio::DB::Flat::BDB;
use vars '@ISA';

@ISA = qw(Bio::DB::Flat::BDB);

sub parse_one_record {
  my $self  = shift;
  my $fh    = shift;

  undef $self->{fasta_stored_id} if exists $self->{fasta_stored_fh}
    && $fh ne $self->{fasta_stored_fh} ;
  $self->{fasta_stored_fh} = $fh;

  while (<$fh>) {		# don't try this at home
    if (/^>(\S+)/) {
      my $id = $self->{fasta_stored_id};
      $self->{fasta_stored_id} = $1;
      next unless defined $id;
      return ($id,-length($_));
    }
  }
  # we get here at the end of the file
  return $self->{fasta_stored_id};
}

sub default_file_format { "fasta" }

1;
