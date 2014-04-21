#
# BioPerl module for Bio::SearchIO::hmmer
#
# Please direct questions and support issues to <bioperl-l@bioperl.org>
#
# Cared for by Kai Blin <kai.blin@biotech.uni-tuebingen.de>
#
# Copyright Kai Blin
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=head1 NAME

Bio::SearchIO::hmmer - A parser for HMMER2 and HMMER3 output (hmmscan, hmmsearch, hmmpfam)

=head1 SYNOPSIS

    # do not use this class directly it is available through Bio::SearchIO
    use Bio::SearchIO;
    my $in = Bio::SearchIO->new(-format => 'hmmer',
                               -file   => 't/data/L77119.hmmer');
    while( my $result = $in->next_result ) {
        # this is a Bio::Search::Result::HMMERResult object
        print $result->query_name(), " for HMM ", $result->hmm_name(), "\n";
        while( my $hit = $result->next_hit ) {
            print $hit->name(), "\n";
            while( my $hsp = $hit->next_hsp ) {
                print "length is ", $hsp->length(), "\n";
            }
        }
    }

=head1 DESCRIPTION

This object implements a parser for HMMER output. It works with both HMMER2 and HMMER3

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

=head1 AUTHOR - Kai Blin

Email kai.blin-at-biotech.uni-tuebingen.de

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::SearchIO::hmmer;

use strict;

use Bio::Factory::ObjectFactory;

use base qw(Bio::SearchIO);

sub new {
  my ( $caller, @args ) = @_;
  my $class = ref($caller) || $caller;

  my $self = $class->SUPER::new(@args);
  $self->_initialize(@args);

  # Try to guess the hmmer format version if it's not specified.
  my $version;
  my %param = @args;

  @param{ map { lc $_ } keys %param } = values %param; # lowercase keys

  # If the caller specified a version, go for that
  if (defined($param{"-version"})) {
    $version = $param{"-version"};
  } else {

    # read second line of the file
    my $first_line = $self->_readline;
    $_ = $self->_readline;

    if ( m/HMMER\s3/ ) {
      $version = "3";
    } else {
      $version = "2";
    }

    $self->_pushback($_);
    $self->_pushback($first_line);
  }

  my $format = "hmmer$version";
  return unless( $class->_load_format_module($format) );

  bless($self, "Bio::SearchIO::$format");

  return $self;
}

sub _initialize {
    my ( $self, @args ) = @_;
    $self->SUPER::_initialize(@args);
    my $handler = $self->_eventHandler;
    $handler->register_factory(
        'result',
        Bio::Factory::ObjectFactory->new(
            -type      => 'Bio::Search::Result::HMMERResult',
            -interface => 'Bio::Search::Result::ResultI'
        )
    );

    $handler->register_factory(
        'hit',
        Bio::Factory::ObjectFactory->new(
            -type      => 'Bio::Search::Hit::HMMERHit',
            -interface => 'Bio::Search::Hit::HitI'
        )
    );

    $handler->register_factory(
        'hsp',
        Bio::Factory::ObjectFactory->new(
            -type      => 'Bio::Search::HSP::HMMERHSP',
            -interface => 'Bio::Search::HSP::HSPI'
        )
    );
}

1;
