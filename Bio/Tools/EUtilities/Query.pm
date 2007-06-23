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

use base qw(Bio::Tools::EUtilities Bio::Tools::EUtilities::HistoryI);

use Bio::Tools::EUtilities::Query::GlobalQuery;

=head1 Bio::Tools::EUtilities::EUtilDataI methods

See Bio::Tools::EUtilities::EUtilDataI for more details

=cut

{
my %TYPE = (
    'espell'    => 'spelling',
    'esearch'   => 'simple',
    'egquery'   => 'base',
    );

sub _add_data {
    my ($self, $qdata) = @_;
    my $eutil = $self->eutil;
    if (!$qdata || ref($qdata) !~ /HASH/i) {
        $self->throw("Bad $eutil data");
    }
    my $type = exists $qdata->{WebEnv} ? 'history' : $TYPE{$eutil};
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
    $self->{'_parsed'} = 1;
}

}

=head1 Bio::Tools::EUtilities::HistoryI methods

=head2 history

 Title    : history
 Usage    : my ($webenv, $qk) = $hist->history
 Function : returns two-element list of webenv() and query_key()
 Returns  : array
 Args     : none

=cut

=head2 get_webenv

 Title    : get_webenv
 Usage    : my $webenv = $hist->get_webenv
 Function : returns web environment key needed to retrieve results from
            NCBI server
 Returns  : string (encoded key)
 Args     : none

=cut

=head2 get_query_key

 Title    : get_query_key
 Usage    : my $qk = $hist->get_query_key
 Function : returns query key (integer) for the history number for this session
 Returns  : integer
 Args     : none

=cut

=head2 has_history

 Title    : has_history
 Usage    : if ($hist->has_history) {...}
 Function : returns TRUE if full history (webenv, query_key) is present 
 Returns  : BOOLEAN, value eval'ing to TRUE or FALUE
 Args     : none

=cut

=head1 Bio::Tools::EUtilities::Query methods

=head1 Query-related methods

=head2 get_count

 Title    : get_count
 Usage    : my $ct = $parser->get_count
 Function : returns the count (hits for a search)
 Returns  : integer
 Args     : [CONDITIONAL] not needed in most cases; required when using egquery
            to retrieve the count for a particular database

=cut

# egquery and espell

sub get_count {
    my ($self, @args) = @_;
    $self->parse_data unless $self->data_parsed;
    return $self->{'_count'};
}

=head2 get_ids

 Title    : get_ids
 Usage    : my @ids = $parser->get_ids
 Function : returns array or array ref of IDs
 Returns  : array or array ref (based on wantarray)
 Args     : none

=cut

# egquery and espell

sub get_ids {
    my ($self, @args) = @_;
    $self->parse_data unless $self->data_parsed;
    return wantarray && $self->{'_id'} ? @{ $self->{'_id'} } : $self->{'_id'} ;
}

=head2 get_term

 Title   : get_term
 Usage   : $st = $qd->get_term;
 Function: retrieve the term for the global search
 Returns : string
 Args    : none

=cut

sub get_term {
    my ($self, @args) = @_;
    $self->parse_data unless $self->data_parsed;
    return $self->{'_term'};
}

=head2 get_database

 Title   : get_database
 Usage   : $ct = $qd->get_database;
 Function: retrieve the database
 Returns : string
 Args    : none

=cut

sub get_database {
    my ($self) = @_;
    $self->parse_data unless $self->data_parsed;
    $self->{'_database'};
}

=head2 get_translation_from

 Title   : get_translation_from
 Usage   : $string = $qd->get_translation_from();
 Function: portion of the original query replaced with translated_to()
 Returns : string
 Args    : none

=cut

# esearch

sub get_translation_from {
    my $self = shift;
    $self->parse_data unless $self->data_parsed;
    return $self->{'_translation'}->{From};
}

=head2 get_translation_to

 Title   : get_translation_to
 Usage   : $string = $qd->get_translation_to();
 Function: replaced string used in place of the original query term in translation_from()
 Returns : string
 Args    : none

=cut

sub get_translation_to {
    my $self = shift;
    $self->parse_data unless $self->data_parsed;
    return $self->{'_translation'}->{To};
}

=head2 get_retstart

 Title   : get_retstart
 Usage   : $start = $qd->get_retstart();
 Function: retstart setting for the query (either set or NCBI default)
 Returns : Integer
 Args    : none

=cut

sub get_retstart {
    my $self = shift;
    $self->parse_data unless $self->data_parsed;    
    return $self->{'_retstart'};
}

=head2 get_retmax

 Title   : get_retmax
 Usage   : $max = $qd->get_retmax();
 Function: retmax setting for the query (either set or NCBI default)
 Returns : Integer
 Args    : none

=cut

sub get_retmax {
    my $self = shift;
    $self->parse_data unless $self->data_parsed;    
    return $self->{'_retmax'};
}

=head2 get_query_translation

 Title   : get_query_translation
 Usage   : $string = $qd->get_query_translation();
 Function: returns the translated query used for the search (if any)
 Returns : string
 Args    : none
 Note    : this differs from the original term.

=cut

sub get_query_translation {
    my $self = shift;
    $self->parse_data unless $self->data_parsed;
    return $self->{'_querytranslation'};
}

=head2 get_corrected_query

 Title    : get_corrected_query
 Usage    : my $cor = $eutil->get_corrected_query;
 Function : retrieves the corrected query when using espell
 Returns  : string 
 Args     : none

=cut

sub get_corrected_query {
    my $self = shift;
    $self->parse_data unless $self->data_parsed;
    return $self->{'_correctedquery'};
}

=head2 get_replaced_terms

 Title    : get_replaced_terms
 Usage    : my $term = $eutil->get_replaced_term
 Function : returns array of strings replaced in the query
 Returns  : string 
 Args     : none

=cut

sub get_replaced_terms {
    my $self = shift;
    $self->parse_data unless $self->data_parsed;
    if ($self->{'_spelledquery'} && $self->{'_spelledquery'}->{Replaced}) {
        ref $self->{'_spelledquery'}->{Replaced} ?
        return @{ $self->{'_spelledquery'}->{Replaced} } : return;
    }
}

=head2 next_GlobalQuery

 Title    : next_GlobalQuery
 Usage    : while (my $query = $eutil->next_GlobalQuery) {...}
 Function : iterates through the queries returned from an egquery search
 Returns  : GlobalQuery object
 Args     : none

=cut

sub next_GlobalQuery {
    my $self = shift;
    unless ($self->{'_globalqueries_it'}) {
        my $cb = $self->{"_globalqueries_cb"};
        $self->throw("Callback must be a code reference")
            if $cb && ref $cb ne 'CODE';
        my $qcount = $self->get_GlobalQueries;
        my $current = 0;
        $self->{"_globalqueries_it"} = sub {
            while ($current < $qcount) {
                if ($cb) {
                    $cb->($self->{'_globalqueries'}->[$current++]) ?
                    return $self->{'_globalqueries'}->[$current] :
                    next;
                } else {
                    return $self->{'_globalqueries'}->[$current++]
                }
            }
        }
    }    
    $self->{'_globalqueries_it'}->();
}

=head2 get_GlobalQueries

 Title    : get_GlobalQueries
 Usage    : @queries = $eutil->get_GlobalQueries
 Function : returns list of GlobalQuery objects
 Returns  : array of GlobalQuery objects
 Args     : none

=cut

sub get_GlobalQueries {
    my $self = shift;
    $self->parse_data unless $self->data_parsed;
    ref $self->{'_globalqueries'} ? return @{ $self->{'_globalqueries'} } : return ();
}

=head2 rewind

 Title    : rewind
 Usage    : $eutil->rewind;
 Function : rewinds the globalquery iterator if present
 Returns  : none
 Args     : none

=cut

sub rewind {
    my $self = shift;
    delete $self->{'_globalqueries_it'};
}



1;

__END__