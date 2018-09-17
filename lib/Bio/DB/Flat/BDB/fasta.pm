#
#
# BioPerl module for Bio::DB::Flat::BDB
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
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

=head1 SEE ALSO

L<Bio::DB::Flat>,

=head1 AUTHOR - Lincoln Stein

Email - lstein@cshl.org

=cut

package Bio::DB::Flat::BDB::fasta;

use strict;

use base qw(Bio::DB::Flat::BDB);

sub default_file_format { "fasta" }

sub seq_to_ids {
  my $self = shift;
  my $seq  = shift;
  my %ids;
  $ids{$self->primary_namespace} = $seq->primary_id;
  \%ids;
}

sub parse_one_record {
    my $self  = shift;
    my $fh    = shift;
    
    # fasta parses by changing $/ to '\n>', need to adjust accordingly
    my $adj = -1;
    my $parser =
      $self->{cached_parsers}{fileno($fh)}
        ||= Bio::SeqIO->new(-fh=>$fh,-format=>$self->default_file_format);
    my $seq = $parser->next_seq or return;
    $self->{flat_alphabet} ||= $seq->alphabet;
    my $ids = $self->seq_to_ids($seq);
    return ($ids, $adj);
}

1;
