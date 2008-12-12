# $Id$
#
# BioPerl module for Bio::DB::EUtilities::Query
#
# Cared for by Chris Fields
#
# Copyright Chris Fields
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code
# 
# Part of the EUtilities BioPerl package

=head1 NAME

Bio::DB::EUtilities::Query

=head1 SYNOPSIS

  #### should not create instance directly; Bio::Tools::EUtilities does this ####

  # can also use '-response' (for HTTP::Response objects) or '-fh' (for filehandles)
  my $info = Bio::Tools::EUtilities->new(-eutil => 'esearch',
                                         -file => 'esearch.xml');

  # esearch

  # esearch with history

  # egquery

  # espell (just for completeness, really)

=head1 DESCRIPTION

...

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the
evolution of this and other Bioperl modules. Send
your comments and suggestions preferably to one
of the Bioperl mailing lists. Your participation
is much appreciated.

  bioperl-l@lists.open-bio.org               - General discussion
  http://www.bioperl.org/wiki/Mailing_lists  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to
help us keep track the bugs and their resolution.
Bug reports can be submitted via the web.

  http://bugzilla.open-bio.org/

=head1 AUTHOR 

Email cjfields at uiuc dot edu

=head1 APPENDIX

The rest of the documentation details each of the
object methods. Internal methods are usually
preceded with a _

=cut

# Let the code begin...

package Bio::Tools::EUtilities::Query;
use strict;
use warnings;
use Bio::Tools::EUtilities::Query::GlobalQuery;
use Bio::Tools::EUtilities::History;
use Data::Dumper;

use base qw(Bio::Tools::EUtilities);

=head1 Bio::Tools::EUtilities::Query methods

=cut

# private EUtilDataI method

{
my %TYPE = (
    'espell'    => 'spelling',
    'esearch'   => 'singledbquery',
    'egquery'   => 'multidbquery',
    'epost'     => 'history'
    );

sub _add_data {
    my ($self, $qdata) = @_;
    my $eutil = $self->eutil;
    if (!$qdata || ref($qdata) !~ /HASH/i) {
        $self->throw("Bad $eutil data");
    }
    if (exists $qdata->{WebEnv}) {
        my $cookie = Bio::Tools::EUtilities::History->new(-eutil => $eutil,
                            -verbose => $self->verbose);
        $cookie->_add_data($qdata);
        push @{$self->{'_histories'}}, $cookie;
    }
    my $type = exists $TYPE{$eutil} ? $TYPE{$eutil} :
        $self->throw("Unrecognized eutil $eutil");
    $self->datatype($type); # reset type based on what's present
    for my $key (sort keys %$qdata) {
        if ($key eq 'eGQueryResult' && exists $qdata->{$key}->{ResultItem}) {
            for my $gquery (@{ $qdata->{eGQueryResult}->{ResultItem} }) {
                $self->{'_term'} = $gquery->{Term} = $qdata->{Term};
                my $qd = Bio::Tools::EUtilities::Query::GlobalQuery->new(-eutil => 'egquery',
                                                            -datatype => 'globalquery',
                                                            -verbose => $self->verbose);
                $qd->_add_data($gquery);
                push @{ $self->{'_globalqueries'} }, $qd;
            }
        }
        if ($key eq 'IdList' &&
            exists $qdata->{IdList}->{Id}) {
            $self->{'_id'} = $qdata->{IdList}->{Id};
            delete $qdata->{IdList};
        }
        if ($key eq 'TranslationSet' &&
            exists $qdata->{TranslationSet}->{Translation}) {
            $self->{'_translation'} = $qdata->{TranslationSet}->{Translation};
            delete $qdata->{TranslationSet};
        }
        next if (ref $qdata->{$key} eq 'HASH' && !keys %{$qdata->{$key}});
        $self->{'_'.lc $key} = $qdata->{$key};
    }
}

}

1;

