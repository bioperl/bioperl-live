# $Id$
#
# BioPerl module for Bio::DB::EUtilities
#
# Cared for by Chris Fields <cjfields at uiuc dot edu>
#
# Copyright Chris Fields
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code
# 
# EUtil-based extension of GenericWebDBI interface 

=head1 NAME

...

=head1 SYNOPSIS

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

package Bio::DB::EUtilities;
use strict;
use Bio::DB::EUtilParameters;
use Bio::Tools::EUtilities;

use base qw(Bio::DB::GenericWebDBI);

sub new {
    my($class,@args) = @_;
    my $self = $class->SUPER::new(@args);
    my $params = Bio::DB::EUtilParameters->new(-verbose => $self->verbose,
                                               @args);
    # cache parameters
    $self->parameter_base($params);
    return $self;
}

=head1 Bio::DB::GenericWebDBI method

=head2 delay_policy

  Title   : delay_policy
  Usage   : $secs = $self->delay_policy
  Function: return number of seconds to delay between calls to remote db
  Returns : number of seconds to delay
  Args    : none

  NOTE: NCBI requests a delay of 3 seconds between requests.  This method
        implements that policy.  This may change to check time of day for lengthening delays if needed

=cut

sub delay_policy {
  my $self = shift;
  return 3;
}

=head2 get_Parser

 Title   : get_Parser
 Usage   : $agent->get_Parser;
 Function: Parse HTTP::Response content/file/fh using defined parser
 Returns : The Bio::Tools::EUtilities parser used to parse the HTTP::Response
           content
 Args    : None
 Note    : Abstract method; defined by implementation

=cut

sub get_Parser {
    my ($self) = @_;
    my $pobj = $self->parameter_base;
    if ($pobj->parameters_changed || !$self->{'_parser'}) {
        my $eutil = $pobj->eutil ;
        if ($eutil eq 'efetch') {
            $self->throw("No parser defined for efetch; use get_response() directly");
        };
        # if we are to add pipe/tempfile support this would probably be the
        # place to add it....
        $self->{'_parser'} = Bio::Tools::EUtilities->new(
                            -eutil => $eutil,
                            -response => $self->get_Response,
                            -verbose => $self->verbose);
    }
    return $self->{'_parser'};
}

=head2 set_history

 Title   : set_history
 Usage   : $agent->set_history;
           $agent->set_history($myhistory);
 Function: sets the history() parameter using the internal parser or (if
           supplied) a passed Bio::Tools::EUtilities::HistoryI object
 Returns : TRUE upon success
 Args    : [optional] Bio::Tools::EUtilities::HistoryI object; uses internal
           parser by default.

=cut

sub set_history {
    my ($self, $hist) = @_;
    $hist ||= $self->get_Parser;
    if ($hist->has_history) {
        $self->parameter_base->history($hist);
        return 1;
    } else {
        $self->warn('No history returned');
        return 0;
    }
}

=head1 Bio::DB::EUtilParameters-delegating methods

This is only a subset of parameters available from Bio::DB::EUtilParameters.

=cut

sub set_parameters {
    my ($self, @args) = @_;
    delete $self->{'_parser'};
    $self->parameter_base->set_parameters(@args);
}

sub reset_parameters {
    my ($self, @args) = @_;
    delete $self->{'_parser'};
    $self->parameter_base->reset_parameters(@args);
}

sub available_parameters {
    my ($self, @args) = @_;
    return $self->parameter_base->available_parameters(@args);
}

=head2 get_parameters

 Title   : get_parameters
 Usage   : @params = $pobj->get_parameters;
           %params = $pobj->get_parameters;
 Function: Returns list of key/value pairs, parameter => value
 Returns : Flattened list of key-value pairs. All key-value pairs returned,
           though subsets can be returned based on the '-type' parameter.  
           Data passed as an array ref are returned based on whether the
           '-join_id' flag is set (default is the same array ref). 
 Args    : -type : the eutil name or 'history', for returning a subset of
                parameters (Default: returns all)
           -join_ids : Boolean; join IDs based on correspondence (Default: no join)

=cut

sub get_parameters {
    my ($self, @args) = @_;
    return $self->parameter_base->available_parameters(@args);
}

=head1 Bio::Tools::EUtilities-delegating methods

=cut

=head1 Bio::Tools::EUtilities::EUtilDataI methods

=head2 eutil

 Title    : eutil
 Usage    : $eutil->$foo->eutil
 Function : Get/Set eutil
 Returns  : string
 Args     : string (eutil)
 Throws   : on invalid eutil
 
=cut

sub eutil {
    my ($self, @args) = @_;
    return $self->get_Parser->eutil(@args);
}

=head2 datatype

 Title    : datatype
 Usage    : $type = $foo->datatype;
 Function : Get/Set data object type
 Returns  : string
 Args     : string

=cut

sub datatype {
    my ($self, @args) = @_;
    return $self->get_Parser->datatype(@args);
}

=head1 Bio::Tools::EUtilities::HistoryI methods

These are defined in the HistoryI interface

=head2 history

 Title    : history
 Usage    : my ($webenv, $qk) = $hist->history
 Function : returns two-element list of webenv() and query_key()
 Returns  : array
 Args     : none

=cut

sub history {
    my ($self, @args) = @_;
    return $self->get_Parser->history(@args);
}

=head2 get_webenv

 Title    : get_webenv
 Usage    : my $webenv = $hist->get_webenv
 Function : returns web environment key needed to retrieve results from
            NCBI server
 Returns  : string (encoded key)
 Args     : none

=cut

sub get_webenv {
    my ($self, @args) = @_;
    return $self->get_Parser->get_webenv(@args);
}

=head2 get_query_key

 Title    : get_query_key
 Usage    : my $qk = $hist->get_query_key
 Function : returns query key (integer) for the history number for this session
 Returns  : integer
 Args     : none

=cut

sub get_query_key {
    my ($self, @args) = @_;
    return $self->get_Parser->get_query_key(@args);
}

=head2 has_history

 Title    : has_history
 Usage    : if ($hist->has_history) {...}
 Function : returns TRUE if full history (webenv, query_key) is present 
 Returns  : BOOLEAN, value eval'ing to TRUE or FALUE
 Args     : none

=cut

sub has_history {
    my ($self, @args) = @_;
    return $self->get_Parser->has_history(@args);
}

=head1 Methods useful for multiple eutils

=head2 get_ids

 Title    : get_ids
 Usage    : my @ids = $parser->get_ids
 Function : returns array or array ref of requestes IDs
 Returns  : array or array ref (based on wantarray)
 Args     : [conditional] not required except when running elink queries against
            multiple databases. In case of the latter, the database name is
            optional (but recommended) when retrieving IDs as the ID list will
            be globbed together. If a db name isn't provided a warning is issued
            as a reminder.

=cut

sub get_ids {
    my ($self, @args) = @_;
    return $self->get_Parser->get_ids(@args);
}

=head2 get_database

 Title    : get_database
 Usage    : my $db = $info->get_database;
 Function : returns database name (eutil-compatible)
 Returns  : string
 Args     : none
 Note     : implemented for einfo and espell
 
=cut

sub get_database {
    my ($self, @args) = @_;
    return $self->get_Parser->get_database(@args);
}

=head2 get_db (alias for get_database)

=cut

sub get_db {
    my ($self, @args) = @_;
    return $self->get_Parser->get_db(@args);
}

=head1 Query-related methods

=head2 get_count

 Title    : get_count
 Usage    : my $ct = $parser->get_count
 Function : returns the count (hits for a search)
 Returns  : integer
 Args     : [CONDITIONAL] string with database name - used to retrieve
            count from specific database when using egquery

=cut

sub get_count {
    my ($self, @args) = @_;
    return $self->get_Parser->get_count(@args);
}
=head2 get_queried_databases

 Title    : get_queried_databases
 Usage    : my @dbs = $parser->get_queried_databases
 Function : returns list of databases searched with global query
 Returns  : array of strings
 Args     : none
 Note     : predominately used for egquery; if used with other eutils will
            return a list with the single database

=cut

sub get_queried_databases {
    my ($self, @args) = @_;
    return $self->get_Parser->get_queried_databases(@args);
}
=head2 get_term

 Title   : get_term
 Usage   : $st = $qd->get_term;
 Function: retrieve the term for the global search
 Returns : string
 Args    : none

=cut

# egquery and espell

sub get_term {
    my ($self, @args) = @_;
    return $self->get_Parser->get_term(@args);
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
    my ($self, @args) = @_;
    return $self->get_Parser->get_translation_from(@args);
}

=head2 get_translation_to

 Title   : get_translation_to
 Usage   : $string = $qd->get_translation_to();
 Function: replaced string used in place of the original query term in translation_from()
 Returns : string
 Args    : none

=cut

sub get_translation_to {
    my ($self, @args) = @_;
    return $self->get_Parser->get_translation_to(@args);
}

=head2 get_retstart

 Title   : get_retstart
 Usage   : $start = $qd->get_retstart();
 Function: retstart setting for the query (either set or NCBI default)
 Returns : Integer
 Args    : none

=cut

sub get_retstart {
    my ($self, @args) = @_;
    return $self->get_Parser->get_retstart(@args);
}

=head2 get_retmax

 Title   : get_retmax
 Usage   : $max = $qd->get_retmax();
 Function: retmax setting for the query (either set or NCBI default)
 Returns : Integer
 Args    : none

=cut

sub get_retmax {
    my ($self, @args) = @_;
    return $self->get_Parser->get_retmax(@args);
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
    my ($self, @args) = @_;
    return $self->get_Parser->get_query_translation(@args);
}

=head2 get_corrected_query

 Title    : get_corrected_query
 Usage    : my $cor = $eutil->get_corrected_query;
 Function : retrieves the corrected query when using espell
 Returns  : string 
 Args     : none

=cut

sub get_corrected_query {
    my ($self, @args) = @_;
    return $self->get_Parser->get_corrected_query(@args);
}

=head2 get_replaced_terms

 Title    : get_replaced_terms
 Usage    : my $term = $eutil->get_replaced_term
 Function : returns array of strings replaced in the query
 Returns  : string 
 Args     : none

=cut

sub get_replaced_terms {
    my ($self, @args) = @_;
    return $self->get_Parser->get_replaced_terms(@args);
}

=head2 next_GlobalQuery

 Title    : next_GlobalQuery
 Usage    : while (my $query = $eutil->next_GlobalQuery) {...}
 Function : iterates through the queries returned from an egquery search
 Returns  : GlobalQuery object
 Args     : none

=cut

sub next_GlobalQuery {
    my ($self, @args) = @_;
    return $self->get_Parser->next_GlobalQuery(@args);
}

=head2 get_GlobalQueries

 Title    : get_GlobalQueries
 Usage    : @queries = $eutil->get_GlobalQueries
 Function : returns list of GlobalQuery objects
 Returns  : array of GlobalQuery objects
 Args     : none

=cut

sub get_GlobalQueries {
    my ($self, @args) = @_;
    return $self->get_Parser->get_GlobalQueries(@args);
}

=head1 Summary-related methods

=head2 next_DocSum

 Title    : next_DocSum
 Usage    : while (my $ds = $esum->next_DocSum) {...}
 Function : iterate through DocSum instances
 Returns  : single Bio::Tools::EUtilities::Summary::DocSum
 Args     : none yet

=cut

sub next_DocSum {
    my ($self, @args) = @_;
    return $self->get_Parser->next_DocSum(@args);
}

=head2 get_DocSums

 Title    : get_DocSums
 Usage    : my @docsums = $esum->get_DocSums
 Function : retrieve a list of DocSum instances
 Returns  : array of Bio::Tools::EUtilities::Summary::DocSum
 Args     : none

=cut

sub get_DocSums {
    my ($self, @args) = @_;
    return $self->get_Parser->get_DocSums(@args);
}

=head1 Info-related methods

=head2 get_available_databases

 Title    : get_available_databases
 Usage    : my @dbs = $info->get_available_databases
 Function : returns list of available eutil-compatible database names
 Returns  : Array of strings 
 Args     : none

=cut

sub get_available_databases {
    my ($self, @args) = @_;
    return $self->get_Parser->get_available_databases(@args);
}

=head2 get_record_count

 Title    : get_record_count
 Usage    : my $ct = $eutil->get_record_count;
 Function : returns database record count
 Returns  : integer
 Args     : none

=cut

sub get_record_count {
    my ($self, @args) = @_;
    return $self->get_Parser->get_record_count(@args);
}

=head2 get_last_update

 Title    : get_last_update
 Usage    : my $time = $info->get_last_update;
 Function : returns string containing time/date stamp for last database update
 Returns  : integer
 Args     : none

=cut

sub get_last_update {
    my ($self, @args) = @_;
    return $self->get_Parser->get_last_update(@args);
} 
=head2 get_menu_name

 Title    : get_menu_name
 Usage    : my $nm = $info->get_menu_name;
 Function : returns string of database menu name
 Returns  : string
 Args     : none
 
=cut

sub get_menu_name {
    my ($self, @args) = @_;
    return $self->get_Parser->get_menu_name(@args);
}

=head2 get_description

 Title    : get_description
 Usage    : my $desc = $info->get_description;
 Function : returns database description
 Returns  : string
 Args     : none

=cut

sub get_description {
    my ($self, @args) = @_;
    return $self->get_Parser->get_description(@args);
}

=head2 next_FieldInfo

 Title    : next_FieldInfo
 Usage    : while (my $field = $info->next_FieldInfo) {...}
 Function : iterate through FieldInfo objects
 Returns  : Field object
 Args     : none
 Note     : uses callback() for filtering if defined for 'fields'
 
=cut

sub next_FieldInfo {
    my ($self, @args) = @_;
    return $self->get_Parser->next_FieldInfo(@args);
}

=head2 get_FieldInfo

 Title    : get_FieldInfo
 Usage    : my @fields = $info->get_FieldInfo;
 Function : returns list of FieldInfo objects
 Returns  : array (FieldInfo objects)
 Args     : none

=cut

sub get_FieldInfo {
    my ($self, @args) = @_;
    return $self->get_Parser->get_FieldInfo(@args);
}

*get_FieldInfos = \&get_FieldInfo;

=head2 next_LinkInfo

 Title    : next_LinkInfo
 Usage    : while (my $link = $info->next_LinkInfo) {...}
 Function : iterate through LinkInfo objects
 Returns  : LinkInfo object
 Args     : none
 
=cut

sub next_LinkInfo {
    my ($self, @args) = @_;
    return $self->get_Parser->next_LinkInfo(@args);
}

=head2 get_LinkInfo

 Title    : get_LinkInfo
 Usage    : my @links = $info->get_LinkInfo;
 Function : returns list of LinkInfo objects
 Returns  : array (LinkInfo objects)
 Args     : none

=cut

sub get_LinkInfo {
    my ($self, @args) = @_;
    return $self->get_Parser->get_LinkInfo(@args);
}

*get_LinkInfos = \&get_LinkInfo;

=head1 Bio::Tools::EUtilities::Link-related methods

=head2 next_LinkSet

 Title    : next_LinkSet
 Usage    : while (my $ls = $eutil->next_LinkSet {...}
 Function : 
 Returns  : 
 Args     : 

=cut

sub next_LinkSet {
    my ($self, @args) = @_;
    return $self->get_Parser->next_LinkSet(@args);
}

=head2 get_LinkSets

 Title    : get_LinkSets
 Usage    : 
 Function : 
 Returns  : 
 Args     : 

=cut

# add support for retrieval of data if lazy parsing is enacted

sub get_LinkSets {
    my ($self, @args) = @_;
    return $self->get_Parser->get_LinkSets(@args);
}

=head2 get_linked_databases

 Title    : get_linked_databases
 Usage    : my @dbs = $eutil->get_linked_databases
 Function : returns list of databases linked to in linksets
 Returns  : array of databases
 Args     : none

=cut

sub get_linked_databases {
    my ($self, @args) = @_;
    return $self->get_Parser->get_linked_databases(@args);
}

=head1 Iterator- and callback-related methods

=cut

=head2 rewind

 Title    : rewind
 Usage    : $esum->rewind()
            $esum->rewind('recursive')
 Function : retrieve a list of DocSum instances
 Returns  : array of Bio::Tools::EUtilities::Summary::DocSum
 Args     : [optional] Scalar; string ('all') to reset all iterators, or string 
            describing the specific main object iterator to reset. The following
            are recognized (case-insensitive):
            
            'all' - rewind all objects and also recursively resets nested object interators
                    (such as LinkSets and DocSums).
            'globalqueries'
            'fieldinfo' or 'fieldinfos'
            'linkinfo' or 'linkinfos'
            'linksets'
            'docsums'
            

=cut

sub rewind {
    my ($self, @args) = @_;
    return $self->get_Parser->rewind(@args);
}

=head2 generate_iterator

 Title    : generate_iterator
 Usage    : my $coderef = $esum->generate_iterator('linkinfo')
 Function : generates an iterator (code reference) which iterates through
            the relevant object indicated by the args
 Returns  : code reference
 Args     : [REQUIRED] Scalar; string describing the specific object to iterate.
            The following are currently recognized (case-insensitive):
            
            'globalqueries'
            'fieldinfo' or 'fieldinfos'
            'linkinfo' or 'linkinfos'
            'linksets'
            'docsums'
            
            A second argument can also be passed to generate a 'lazy' iterator,
            which loops through and returns objects as they are created (instead
            of creating all data instances up front, then iterating through,
            which is the default). Use of these iterators precludes use of
            rewind() for the time being as we can't guarantee you can rewind(),
            as this depends on whether the data source is seek()able and thus
            'rewindable'. We will add rewind() support at a later time which
            will work for 'seekable' data.
            
            A callback specified using callback() will be used to filter objects
            for any generated iterator. This behaviour is implemented for both
            normal and lazy iterator types and is the default. If you don't want
            this, make sure to reset any previously set callbacks via
            reset_callback() (which just deletes the code ref).
            
=cut

sub generate_iterator {
    my ($self, @args) = @_;
    return $self->get_Parser->generate_iterator(@args);
}

=head2 callback

 Title    : callback
 Usage    : $parser->callback(sub {$_[0]->get_database eq 'protein'});
 Function : Get/set callback code ref used to filter returned data objects
 Returns  : code ref if previously set
 Args     : single argument:
            code ref - evaluates a passed object and returns true or false value
                       (used in iterators)
            'reset' - string, resets the iterator.
            returns upon any other args
=cut

sub callback {
    my ($self, @args) = @_;
    return $self->get_Parser->callback(@args);
}

1;
__END__
