#
# BioPerl module for Bio::Tools::EUtilities::Query
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
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

Bio::Tools::EUtilities::Query - parse and collect esearch, epost, espell,
egquery information

=head1 SYNOPSIS

  ### should not create instance directly; Bio::Tools::EUtilities does this ###

  # can also use '-response' (for HTTP::Response objects) or '-fh' (for
  # filehandles)
  
  my $info = Bio::Tools::EUtilities->new(-eutil => 'esearch',
                                         -file => 'esearch.xml');

  # esearch

  # esearch with history

  # egquery

  # espell (just for completeness, really)

=head1 DESCRIPTION

Pluggable module for handling query-related data returned from eutils.  

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the
evolution of this and other Bioperl modules. Send
your comments and suggestions preferably to one
of the Bioperl mailing lists. Your participation
is much appreciated.

  bioperl-l@lists.open-bio.org               - General discussion
  http://www.bioperl.org/wiki/Mailing_lists  - About the mailing lists

=head2 Support 

Please direct usage questions or support issues to the mailing list:

I<bioperl-l@bioperl.org>

rather than to the module maintainer directly. Many experienced and 
reponsive experts will be able look at the problem and quickly 
address it. Please include a thorough description of the problem 
with code and data examples if at all possible.

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to
help us keep track the bugs and their resolution.
Bug reports can be submitted via the web.

  http://bugzilla.open-bio.org/

=head1 AUTHOR 

Email cjfields at bioperl dot org

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

=head2 to_string

 Title    : to_string
 Usage    : $foo->to_string()
 Function : converts current object to string
 Returns  : none
 Args     : (optional) simple data for text formatting
 Note     : Used generally for debugging and for the print_* methods

=cut

sub to_string {
    my $self = shift;
    my %data = (
        'DB'    => [1, join(', ',$self->get_databases) || ''],
        'Query' => [2, $self->get_term || ''],
        'IDs'   => [4, join(', ',$self->get_ids) || ''],
    );
    my $string = $self->SUPER::to_string;
    if ($self->eutil eq 'esearch') {
        $data{'Count'} = [3, $self->get_count ];
        $data{'Translation From'} = [5, $self->get_translation_from || ''];
        $data{'Translation To'} = [6, $self->get_translation_to || ''];
        $data{'RetStart'} = [7, $self->get_retstart];
        $data{'RetMax'} = [8, $self->get_retmax];
        $data{'Translation'} = [9, $self->get_query_translation || ''];
    }
    if ($self->eutil eq 'espell') {
        $data{'Corrected'} = [3, $self->get_corrected_query || ''];
        $data{'Replaced'} = [4, join(',',$self->get_replaced_terms) || ''];
    }
    for my $k (sort {$data{$a}->[0] <=> $data{$b}->[0]} keys %data) {
        $string .= sprintf("%-20s:%s\n",$k, $self->_text_wrap('',' 'x 20 .':', $data{$k}->[1]));
    }
    while (my $h = $self->next_History) {
        $string .= $h->to_string;
    }
    while (my $gq = $self->next_GlobalQuery) {
        $string .= $gq->to_string;
    }
    return $string;
}

1;

