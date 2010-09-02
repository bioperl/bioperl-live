#
# BioPerl module for Bio::Tools::EUtilities
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

=head1 NAME

Bio::Tools::EUtilities - NCBI eutil XML parsers

=head1 SYNOPSIS

  # from file or fh
  my $parser = Bio::Tools::EUtilities->new(
                                       -eutil    => 'einfo',
                                       -file     => 'output.xml'
                                        );
  
  # or HTTP::Response object...
  my $parser = Bio::Tools::EUtilities->new(
                                       -eutil => 'esearch',
                                       -response => $response
                                        );
  # esearch, esummary, elink
  @ids = $parser->get_ids(); # returns array or array ref of IDs

  # egquery, espell
  
  $term = $parser->get_term(); # returns array or array ref of IDs
  
  # elink, einfo
  
  $db = $parser->get_database(); # returns database
  
  # Query-related methods (esearch, egquery, espell data)
  # eutil data centered on use of search terms
  
  my $ct = $parser->get_count; # uses optional database for egquery count
  my $translation = $parser->get_count;
  
  my $corrected = $parser->get_corrected_query; # espell

  while (my $gquery = $parser->next_GlobalQuery) {
     # iterates through egquery data
  }
  
  # Info-related methods (einfo data)
  # database-related information
  
  my $desc = $parser->get_description;
  my $update = $parser->get_last_update;
  my $nm = $parser->get_menu_name;
  my $ct = $parser->get_record_count;
  
  while (my $field = $parser->next_FieldInfo) {
      # ...
  }
  while (my $field = $parser->next_LinkInfo) {
      # ...
  }
  
  # History methods (epost data, some data returned from elink)
  # data which enables one to retrieve and query against user-stored
  # information on the NCBI server
  
  while (my $cookie = $parser->next_History) {
      # ...
  }
  
  my @hists = $parser->get_Histories;
  
  # Bio::Tools::EUtilities::Summary (esummary data)
  # information on a specific database record
  
  # retrieve nested docsum data
  while (my $docsum = $parser->next_DocSum) {
      print "ID:",$docsum->get_ids,"\n";
      while (my $item = $docsum->next_Item) {
          # do stuff here...
          while (my $listitem = $docsum->next_ListItem) {
              # do stuff here...
              while (my $listitem = $docsum->next_Structure) {
                  # do stuff here...
              }
          }
      }
  }
  
  # retrieve flattened item list per DocSum
  while (my $docsum = $parser->next_DocSum) {
     my @items = $docsum->get_all_DocSum_Items;
  }

=head1 DESCRIPTION

Parses NCBI eutils XML output for retrieving IDs and other information. Part of
the BioPerl EUtilities system.

This is a general parser for eutils XML; data from efetch is NOT parsed (this
requires separate format-dependent parsers). All other XML for eutils is parsed.
These modules can be used independently of Bio::DB::EUtilities and
Bio::Tools::EUtilities::EUtilParameters; if used in this way, only data present
in the XML will be parsed out (other bits are retrieved from a passed-in
Bio::Tools::EUtilities::EUtilParameters instance used while querying the
database)

=head1 TODO

This module is largely complete. However there are a few holes which will
eventually be filled in. TranslationSets from esearch are not currently parsed,
for instance.

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

package Bio::Tools::EUtilities;
use strict;
use warnings;

use base qw(Bio::Root::IO Bio::Tools::EUtilities::EUtilDataI);
use XML::Simple;

=head2 Constructor methods

=cut

=head2 new

 Title    : new
 Usage    : my $parser = Bio::Tools::EUtilities->new(-file => 'my.xml',
                                                    -eutil => 'esearch');
 Function : create Bio::Tools::EUtilities instance
 Returns  : new Bio::Tools::EUtilities instance
 Args     : -file/-fh - File or filehandle
            -eutil    - eutil parser to use (supports all but efetch)
            -response - HTTP::Response object (optional)

=cut

{

my %DATA_MODULE = (
    'esearch'   => 'Query',
    'egquery'   => 'Query',
    'espell'    => 'Query',
    'epost'     => 'Query',
    'elink'     => 'Link',
    'einfo'     => 'Info',
    'esummary'  => 'Summary',
    );

sub new {
    my($caller,@args) = @_;
    my $class = ref $caller || $caller;
    if ($class =~ m{Bio::Tools::EUtilities::(\S+)}) {
        my ($self) = $class->SUPER::new(@args);
        $self->_initialize(@args);
        return $self;
    } else {
        my %param = @args;
		@param{ map { lc $_ } keys %param } = values %param; # lowercase keys
        my $eutil = $param{'-eutil'} || $class->throw("Need eutil to make instance");
        return unless( $class->_load_eutil_module( $DATA_MODULE{$eutil}) );
        return "Bio::Tools::EUtilities::$DATA_MODULE{$eutil}"->new(-datatype => lc $DATA_MODULE{$eutil},
                                                                   -eutil => $eutil,
                                                                   @args);
    }
}

sub _initialize {
    my ($self, @args) = @_;
    my ($response, $pobj, $type, $eutil, $cache, $lazy) =
    $self->_rearrange([qw(RESPONSE
                       PARAMETERS
                       DATATYPE
                       EUTIL
                       CACHE_RESPONSE
                       LAZY)], @args);
    $lazy ||= 0;
    $cache ||= 0;
    $self->datatype($type);
    $self->eutil($eutil);
    # lazy parsing only implemented for elink and esummary (where returned data
    # can be quite long).  Also, no point to parsing lazily when the data is
    # already in memory in an HTTP::Response object, so turn it off and chunk
    # the Response object after parsing.
    $response  && $self->response($response);
    $pobj && $self->parameter_base($pobj);
    $self->cache_response($cache);
    $lazy = 0 if ($response) || ($eutil ne 'elink' && $eutil ne 'esummary');
    # setting parser to 'lazy' mode is permanent (can't reset later)
    $self->{'_lazy'} = $lazy;
    $self->{'_parsed'} = 0;
}

}

=head1 Bio::Tools::EUtilities methods

=head2 cache_response

 Title    : cache_response
 Usage    : $parser->cache_response(1)
 Function : sets flag to cache response object (off by default)
 Returns  : value eval'ing to TRUE or FALSE
 Args     : value eval'ing to TRUE or FALSE
 Note     : must be set prior to any parsing run

=cut

sub cache_response {
    my ($self, $cache) = @_;
    if (defined $cache) {
        $self->{'_cache_response'} = ($cache) ? 1 : 0;
    }
    return $self->{'_cache_response'};
}

=head2 response

 Title    : response
 Usage    : my $response = $parser->response;
 Function : Get/Set HTTP::Response object
 Returns  : HTTP::Response
 Args     : HTTP::Response
 Note     : to prevent object from destruction set cache_response() to TRUE

=cut

sub response {
    my ($self, $response) = @_;
    if ($response) {
        $self->throw('Not an HTTP::Response object') unless (ref $response && $response->isa('HTTP::Response'));
        $self->{'_response'} = $response; 
    }
    return $self->{'_response'};
}

=head2 parameter_base

 Title    : parameter_base
 Usage    : my $response = $parser->parameter_base;
 Function : Get/Set Bio::ParameterBaseI object (should be Bio::Tools::EUtilities::EUtilParameters)
 Returns  : Bio::Tools::EUtilities::EUtilParameters || undef
 Args     : (optional) Bio::Tools::EUtilities::EUtilParameters
 Note     : If this object is present, it may be used as a last resort for
            some data values if parsed XML does not contain said values (for
            instance, database, term, IDs, etc).

=cut

sub parameter_base {
    my ($self, $pb) = @_;
    if ($pb) {
        $self->throw('Not an Bio::ParameterBaseI object') unless (ref $pb && $pb->isa('Bio::ParameterBaseI'));
        $self->warn('Not an Bio::Tools::EUtilities::EUtilParameters object; may experience some turbulence...') unless (ref $pb && $pb->isa('Bio::Tools::EUtilities::EUtilParameters'));
        $self->{'_parameter_base'} = $pb; 
    }
    return $self->{'_parameter_base'};
}

=head2 data_parsed

 Title    : data_parsed
 Usage    : if ($parser->data_parsed) {...}
 Function : returns TRUE if data has been parsed
 Returns  : value eval'ing to TRUE or FALSE
 Args     : none (set within parser)
 Note     : mainly internal method (set in case user wants to check
            whether parser is exhausted).

=cut

sub data_parsed {
    return shift->{'_parsed'};
}

=head2 is_lazy

 Title    : is_lazy
 Usage    : if ($parser->is_lazy) {...}
 Function : returns TRUE if parser is set to lazy parsing mode
            (only affects elink/esummary)
 Returns  : Boolean
 Args     : none
 Note     : Permanently set in constructor.  Still highly experimental.
            Don't stare directly at happy fun ball...

=cut

sub is_lazy {
    return shift->{'_lazy'};
}

=head2 parse_data

 Title    : parse_data
 Usage    : $parser->parse_data
 Function : direct call to parse data; normally implicitly called
 Returns  : none
 Args     : none

=cut

{
my %EUTIL_DATA = (
    'esummary'  => [qw(DocSum Item)],
    'epost'     => [],
    'egquery'   => [],
    'einfo'     => [qw(Field Link)],
    'elink'     => [qw(LinkSet LinkSetDb LinkSetDbHistory IdUrlSet 
                        Id IdLinkSet ObjUrl Link LinkInfo)],
    'espell'    => [qw(Original Replaced)],
    'esearch'   => [qw(Id ErrorList WarningList)],
    );

sub parse_data {
    my $self = shift;
    my $eutil = $self->eutil;
    my $xs = XML::Simple->new();
    my $response = $self->response ? $self->response :
                   $self->_fh      ? $self->_fh      :
        $self->throw('No response or stream specified');
    my $simple = ($eutil eq 'espell') ?
            $xs->XMLin($self->_fix_espell($response), forcearray => $EUTIL_DATA{$eutil}) :
        ($response && $response->isa("HTTP::Response")) ?
            $xs->XMLin($response->content, forcearray => $EUTIL_DATA{$eutil}) :
            $xs->XMLin($response, forcearray => $EUTIL_DATA{$eutil});
    # check for errors
    if ($simple->{ERROR}) {
        my $error = $simple->{ERROR};
        $self->throw("NCBI $eutil fatal error: ".$error) unless ref $error;
    }
    if ($simple->{InvalidIdList}) {
        $self->warn("NCBI $eutil error: Invalid ID List".$simple->{InvalidIdList});
        return;
    }    
    if ($simple->{ErrorList} || $simple->{WarningList}) {
        my @errorlist = @{ $simple->{ErrorList} } if $simple->{ErrorList};
        my @warninglist = @{ $simple->{WarningList} } if $simple->{WarningList};
        my ($err_warn);
        for my $error (@errorlist) {
            my $messages = join("\n\t",map {"$_  [".$error->{$_}.']'}
                                grep {!ref $error->{$_}} keys %$error);
            $err_warn .= "Error : $messages";
        }    
        for my $warn (@warninglist) {
            my $messages = join("\n\t",map {"$_  [".$warn->{$_}.']'}
                                grep {!ref $warn->{$_}} keys %$warn);
            $err_warn .= "Warnings : $messages";
        }
        chomp($err_warn);
        $self->warn("NCBI $eutil Errors/Warnings:\n".$err_warn)
        # don't return as some data may still be useful
    }
    delete $self->{'_response'} unless $self->cache_response;
    $self->{'_parsed'} = 1;    
    $self->_add_data($simple);
}

# implemented only for elink/esummary, still experimental

sub parse_chunk {
    my $self = shift;
    my $eutil = $self->eutil;
    my $tag = $eutil eq 'elink'    ? 'LinkSet' :
              $eutil eq 'esummary' ? 'DocSum'  :
              $self->throw("Only eutil elink/esummary use parse_chunk()");
    my $xs = XML::Simple->new();
    if ($self->response) {
        $self->throw("Lazy parsing not implemented for HTTP::Response data yet");
        delete $self->{'_response'} if !$self->cache_response && $self->data_parsed;
    } else { # has to be a file/filehandle
        my $fh = $self->_fh;
        my ($chunk, $seendoc, $line);
        CHUNK:
        while ($line = <$fh>) {
            next unless $seendoc || $line =~ m{^<$tag>};
            $seendoc = 1;
            $chunk .= $line;
            last if $line =~ m{^</$tag>};
        }
        if (!defined $line) {
            $self->{'_parsed'} = 1;
            return;
        }
        $self->_add_data(
            $xs->XMLin($chunk, forcearray => $EUTIL_DATA{$eutil}, KeepRoot => 1)
            );
    }
}

}

=head2 to_string

 Title    : to_string
 Usage    : $foo->to_string()
 Function : converts current object to string
 Returns  : none
 Args     : (optional) simple data for text formatting
 Note     : Implemented in plugins

=cut

sub to_string {
    my $self = shift;
    $self->parse_data if ($self->can('parse_data') && !$self->data_parsed);
    return sprintf("%-20s:%s\n\n", 'EUtil', $self->eutil);
}

=head2 print_all

 Title    : print_all
 Usage    : $info->print_all();
            $info->print_all(-fh => $fh, -cb => $coderef);
 Function : prints (dumps) all data in parser.  Unless a coderef is supplied,
            this just dumps the parser-specific to_string method to either a
            file/fh or STDOUT
 Returns  : none
 Args     : [optional]
           -file : file to print to
           -fh   : filehandle to print to (cannot be used concurrently with file)
           -cb   : coderef to use in place of default print method.  This is
                   passed in the parser object 
           -wrap : number of columns to wrap default text output to (def = 80)
 Notes    : only applicable for einfo.  If -file or -fh are not defined,
            prints to STDOUT

=cut

sub print_all {
    my ($self, @args) = @_;
    $self->_print_handler(@args);
}

=head1 Bio::Tools::EUtilities::EUtilDataI methods

=head2 eutil

 Title    : eutil
 Usage    : $eutil->$foo->eutil
 Function : Get/Set eutil
 Returns  : string
 Args     : string (eutil)
 Throws   : on invalid eutil

=cut

=head2 datatype

 Title    : datatype
 Usage    : $type = $foo->datatype;
 Function : Get/Set data object type
 Returns  : string
 Args     : string

=cut

=head1 Methods useful for multiple eutils

=head2 get_ids

 Title    : get_ids
 Usage    : my @ids = $parser->get_ids
 Function : returns array of requested IDs (see Notes for more specifics)
 Returns  : array
 Args     : [conditional] not required except when running elink queries against
            multiple databases. In case of the latter, the database name is
            optional but recommended when retrieving IDs as the ID list will
            be globbed together. In such cases, if a db name isn't provided a
            warning is issued as a reminder.
 Notes    : esearch    : returned ID list
            elink      : returned ID list (see Args above for caveats)
            all others : from parameter_base->id or undef

=cut

sub get_ids {
    my ($self, $request) = @_;
    my $eutil = $self->eutil;
    if ($self->is_lazy) {
        $self->warn('get_ids() not implemented when using lazy mode');
        return;
    }
    $self->parse_data unless $self->data_parsed;
    if ($eutil eq 'esearch') {
        return $self->{'_id'} ? @{ $self->{'_id'} } : ();
    } elsif ($eutil eq 'elink')  {
        my @ids;
        if ($request) {
            if (ref $request eq 'CODE') {
                push @ids, map {$_->get_ids }
                    grep { $request->($_) } $self->get_LinkSets;
            } else {
                push @ids,
                    map { @{$_->[0]} }
                    grep {grep { $_ eq $request } @{$_->[1]}}
                    map {[[$_->get_ids], [$_->get_databases]]} $self->get_LinkSets;
            }
        } else {
            $self->warn('Multiple database present, IDs will be globbed together')
                if $self->get_linked_databases > 1;
            push @ids, map {$_->get_ids } $self->get_LinkSets;
        }
        return @ids;
    } elsif ($eutil eq 'esummary') {
        unless (exists $self->{'_id'}) {
            push @{$self->{'_id'}}, map {$_->get_id } $self->get_DocSums;
        }
        return @{$self->{'_id'}};
    } elsif (my $pb = $self->parameter_base) {
        my $ids = $pb->id;
        return $ids ? @{$ids} : ();
    } else {
        return ()
    }
}

=head2 get_database

 Title    : get_database
 Usage    : my $db = $info->get_database;
 Function : returns single database name (eutil-compatible).  This is the
            queried database. For most eutils this is straightforward. For
            elinks (which have 'db' and 'dbfrom') this is db/dbto, for egquery,
            it is the first db in the list (you probably want get_databases
            instead)
 Returns  : string
 Args     : none
 Notes    : egquery    : first db in the query (you probably want get_databases)
            einfo      : the queried database
            espell     : the queried database
            all others : from parameter_base->db or undef

=cut

sub get_database {
    return ($_[0]->get_databases)[0];
}

=head2 get_db (alias for get_database)

=cut

sub get_db {
    return shift->get_database;
}

=head2 get_databases

 Title    : get_databases
 Usage    : my @dbs = $parser->get_databases
 Function : returns list of databases 
 Returns  : array of strings
 Args     : none
 Notes    : This is guaranteed to return a list of databases. For a single
            database use the convenience method get_db/get_database
            
            egquery    : list of all databases in the query
            einfo      : the queried database, or the available databases
            espell     : the queried database
            elink      : collected from each LinkSet
            all others : from parameter_base->db or undef

=cut

sub get_databases {
    my ($self, $db) = @_;
    $self->parse_data unless $self->data_parsed;
    my $eutil = $self->eutil;
    my @dbs;
    if ($eutil eq 'einfo' || $eutil eq 'espell') {
        @dbs = $self->{'_dbname'} ||
        $self->{'_database'} ||
        $self->get_available_databases;
    } elsif ($eutil eq 'egquery') {
        @dbs = map {$_->get_database} ($self->get_GlobalQueries);
    } elsif ($eutil eq 'elink') {
        # only unique dbs
        my %tmp;
        @dbs = sort grep {!$tmp{$_}++} 
            map {($_->get_databases)} $self->get_LinkSets;
    } elsif ($self->parameter_base) {
        if ($self->parameter_base->eutil eq 'elink') {
            @dbs = $self->parameter_base->dbfrom;
        } else {
            @dbs = $self->parameter_base->db;
        }
    }
    return @dbs;
}

=head2 get_dbs (alias for get_databases)

=cut

sub get_dbs {
    return shift->get_databases;
}

=head2 next_History

 Title    : next_History
 Usage    : while (my $hist=$parser->next_History) {...}
 Function : returns next HistoryI (if present).
 Returns  : Bio::Tools::EUtilities::HistoryI (Cookie or LinkSet)
 Args     : none
 Note     : esearch, epost, and elink are all capable of returning data which
            indicates search results (in the form of UIDs) is stored on the
            remote server. Access to this data is wrapped up in simple interface
            (HistoryI), which is implemented in two classes:
            Bio::DB::EUtilities::History (the simplest) and
            Bio::DB::EUtilities::LinkSet. In general, calls to epost and esearch
            will only return a single HistoryI object (formerly known as a
            Cookie), but calls to elink can generate many depending on the
            number of IDs, the correspondence, etc. Hence this iterator, which
            allows one to retrieve said data one piece at a time.

=cut

sub next_History {
    my $self = shift;
    $self->parse_data unless $self->data_parsed;    
    $self->{'_histories_it'} = $self->generate_iterator('histories')
        if (!exists $self->{'_histories_it'});
    my $hist =  $self->{'_histories_it'}->();
}

=head2 next_cookie (alias for next_History)

=cut 

sub next_cookie {
    return shift->next_History;
}

=head2 get_Histories

 Title    : get_Histories
 Usage    : my @hists = $parser->get_Histories
 Function : returns list of HistoryI objects.
 Returns  : list of Bio::Tools::EUtilities::HistoryI (History or LinkSet)
 Args     : none

=cut

sub get_Histories {
    my $self = shift;
    $self->parse_data unless $self->data_parsed;
    ref $self->{'_histories'} ? return @{ $self->{'_histories'} } : return ();
}

=head1 Query-related methods

=head2 get_count

 Title    : get_count
 Usage    : my $ct = $parser->get_count
 Function : returns the count (hits for a search)
 Returns  : integer
 Args     : [CONDITIONAL] string with database name - used to retrieve
            count from specific database when using egquery
 Notes    : egquery    : count for specified database (specified above)
            esearch    : count for last search
            all others : undef

=cut

sub get_count {
    my ($self, $db) = @_;
    $self->parse_data unless $self->data_parsed;
    # egquery
    if ($self->datatype eq 'multidbquery') {
        if (!$db) {
            $self->warn('Must specify database to get count from');
            return;
        }
        my ($gq) = grep {$_->get_database eq $db} $self->get_GlobalQueries;
        $gq && return $gq->get_count;
        $self->warn("Unknown database $db");
        return;
    } else {
        return $self->{'_count'};
    }
}

=head2 get_term

 Title    : get_term
 Usage    : $st = $qd->get_term;
 Function : retrieve the term for the global search
 Returns  : string
 Args     : none
 Notes    : egquery    : search term
            espell     : search term
            esearch    : from parameter_base->term or undef
            all others : undef

=cut

sub get_term {
    my ($self, @args) = @_;
    $self->parse_data unless $self->data_parsed;
    $self->{'_term'}  ? $self->{'_term'}  :
    $self->{'_query'} ? $self->{'_query'} :
    $self->parameter_base ? $self->parameter_base->term :
    return;
}

=head2 get_translation_from

 Title   : get_translation_from
 Usage   : $string = $qd->get_translation_from();
 Function: portion of the original query replaced with translated_to()
 Returns : string
 Args    : none
 Note    : only applicable for esearch

=cut

sub get_translation_from {
    my $self = shift;
    $self->parse_data unless $self->data_parsed;
    return $self->{'_translation'}->{'From'};
}

=head2 get_translation_to

 Title   : get_translation_to
 Usage   : $string = $qd->get_translation_to();
 Function: replaced string used in place of the original query term in translation_from()
 Returns : string
 Args    : none
 Note    : only applicable for esearch 

=cut

sub get_translation_to {
    my $self = shift;
    $self->parse_data unless $self->data_parsed;
    return $self->{'_translation'}->{'To'};
}

=head2 get_retstart

 Title    : get_retstart
 Usage    : $start = $qd->get_retstart();
 Function : retstart setting for the query (either set or NCBI default)
 Returns  : Integer
 Args     : none
 Notes    : esearch    : retstart
            esummary   : retstart
            all others : from parameter_base->retstart or undef

=cut

sub get_retstart {
    my $self = shift;
    $self->parse_data unless $self->data_parsed;    
    return $self->{'_retstart'};
}

=head2 get_retmax

 Title    : get_retmax
 Usage    : $max = $qd->get_retmax();
 Function : retmax setting for the query (either set or NCBI default)
 Returns  : Integer
 Args     : none
 Notes    : esearch    : retmax
            esummary   : retmax
            all others : from parameter_base->retmax or undef

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
 Notes   : only applicable for esearch.  This is the actual term used for
           esearch.

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
 Notes    : only applicable for espell.

=cut

sub get_corrected_query {
    my $self = shift;
    $self->parse_data unless $self->data_parsed;
    return $self->{'_correctedquery'};
}

=head2 get_replaced_terms

 Title    : get_replaced_terms
 Usage    : my $term = $eutil->get_replaced_terms
 Function : returns array of strings replaced in the query
 Returns  : string 
 Args     : none
 Notes    : only applicable for espell

=cut

sub get_replaced_terms {
    my $self = shift;
    $self->parse_data unless $self->data_parsed;
    if ($self->{'_spelledquery'} && $self->{'_spelledquery'}->{Replaced}) {
        ref $self->{'_spelledquery'}->{Replaced} ?
        return @{ $self->{'_spelledquery'}->{Replaced} } : return ();
    }
}

=head2 next_GlobalQuery

 Title    : next_GlobalQuery
 Usage    : while (my $query = $eutil->next_GlobalQuery) {...}
 Function : iterates through the queries returned from an egquery search
 Returns  : GlobalQuery object
 Args     : none
 Notes    : only applicable for egquery

=cut

sub next_GlobalQuery {
    my $self = shift;
    $self->parse_data unless $self->data_parsed;    
    $self->{'_globalqueries_it'} = $self->generate_iterator('globalqueries')
        if (!exists $self->{'_globalqueries_it'});
    $self->{'_globalqueries_it'}->();
}

=head2 get_GlobalQueries

 Title    : get_GlobalQueries
 Usage    : @queries = $eutil->get_GlobalQueries
 Function : returns list of GlobalQuery objects
 Returns  : array of GlobalQuery objects
 Args     : none
 Notes    : only applicable for egquery

=cut

sub get_GlobalQueries {
    my $self = shift;
    $self->parse_data unless $self->data_parsed;
    ref $self->{'_globalqueries'} ? return @{ $self->{'_globalqueries'} } : return ();
}

=head2 print_GlobalQueries

 Title    : print_GlobalQueries
 Usage    : $docsum->print_GlobalQueries();
            $docsum->print_GlobalQueries(-fh => $fh, -callback => $coderef);
 Function : prints item data for all global queries.  The default printing
            method is each item per DocSum is printed with relevant values if
            present in a simple table using Text::Wrap. 
 Returns  : none
 Args     : [optional]
           -file : file to print to
           -fh   : filehandle to print to (cannot be used concurrently with file)
           -cb   : coderef to use in place of default print method.  This is passed
                   in a GlobalQuery object;
           -wrap : number of columns to wrap default text output to (def = 80)
 Notes    : only applicable for esummary.  If -file or -fh are not defined,
            prints to STDOUT

=cut

sub print_GlobalQueries {
    my ($self, @args) = @_;
    $self->_print_handler(@args, -type => 'GlobalQuery');
}

=head1 Summary-related methods

=head2 next_DocSum

 Title    : next_DocSum
 Usage    : while (my $ds = $esum->next_DocSum) {...}
 Function : iterate through DocSum instances
 Returns  : single Bio::Tools::EUtilities::Summary::DocSum
 Args     : none yet
 Notes    : only applicable for esummary

=cut

sub next_DocSum {
    my $self = shift;
    if(!$self->data_parsed && !$self->is_lazy) {
        $self->parse_data;
    }
    $self->{'_docsums_it'} = $self->generate_iterator('docsums')
        if (!exists $self->{'_docsums_it'});
    $self->{'_docsums_it'}->();
}

=head2 get_DocSums

 Title    : get_DocSums
 Usage    : my @docsums = $esum->get_DocSums
 Function : retrieve a list of DocSum instances
 Returns  : array of Bio::Tools::EUtilities::Summary::DocSum
 Args     : none
 Notes    : only applicable for esummary

=cut

sub get_DocSums {
    my $self = shift;
    if ($self->is_lazy) {
        $self->warn('get_DocSums() not implemented when using lazy mode');
        return ();
    }
    $self->parse_data unless $self->data_parsed;
    return ref $self->{'_docsums'} ? @{ $self->{'_docsums'} } : return ();
}

=head2 print_DocSums

 Title    : print_DocSums
 Usage    : $docsum->print_DocSums();
            $docsum->print_DocSums(-fh => $fh, -cb => $coderef);
 Function : prints item data for all docsums.  The default data is generated
            via DocSum::to_string
 Returns  : none
 Args     : [optional]
           -file : file to print to
           -fh   : filehandle to print to (cannot be used concurrently with file)
           -cb   : coderef to use in place of default print method.  This is passed
                   in a DocSum object
           -wrap : number of columns to wrap default text output to (def = 80)
 Notes    : only applicable for esummary.  If -file or -fh are not defined,
            prints to STDOUT

=cut

sub print_DocSums {
    my ($self, @args) = @_;
    $self->_print_handler(@args, -type => 'DocSum');
}

=head1 Info-related methods

=head2 get_available_databases

 Title    : get_available_databases
 Usage    : my @dbs = $info->get_available_databases
 Function : returns list of available eutil-compatible database names
 Returns  : Array of strings 
 Args     : none
 Notes    : only applicable for einfo. 

=cut

sub get_available_databases {
    my $self = shift;
    $self->parse_data unless $self->data_parsed;
    ($self->{'_available_databases'}) ?
        return @{($self->{'_available_databases'})} :
        return ();
}

=head2 get_record_count

 Title    : get_record_count
 Usage    : my $ct = $eutil->get_record_count;
 Function : returns database record count
 Returns  : integer
 Args     : none
 Notes    : only applicable for einfo.  

=cut

sub get_record_count {
    my $self = shift;
    $self->parse_data unless $self->data_parsed;
    return $self->{'_count'}
}

=head2 get_last_update

 Title    : get_last_update
 Usage    : my $time = $info->get_last_update;
 Function : returns string containing time/date stamp for last database update
 Returns  : integer
 Args     : none
 Notes    : only applicable for einfo. 

=cut

sub get_last_update {
    my $self = shift;
    $self->parse_data unless $self->data_parsed;
    return $self->{'_lastupdate'}
}

=head2 get_menu_name

 Title    : get_menu_name
 Usage    : my $nm = $info->get_menu_name;
 Function : returns string of database menu name
 Returns  : string
 Args     : none
 Notes    : only applicable for einfo. 

=cut

sub get_menu_name {
    my $self = shift;
    $self->parse_data unless $self->data_parsed;    
    exists $self->{'_menuname'} ? return $self->{'_menuname'} :
    exists $self->{'_menu'} ? return $self->{'_menu'} :
    return;
}

=head2 get_description

 Title    : get_description
 Usage    : my $desc = $info->get_description;
 Function : returns database description
 Returns  : string
 Args     : none
 Notes    : only applicable for einfo. 

=cut

sub get_description {
    my $self = shift;
    $self->parse_data unless $self->data_parsed;
    return $self->{'_description'};
}

=head2 next_FieldInfo

 Title    : next_FieldInfo
 Usage    : while (my $field = $info->next_FieldInfo) {...}
 Function : iterate through FieldInfo objects
 Returns  : Field object
 Args     : none
 Notes    : only applicable for einfo. Uses callback() for filtering if defined
            for 'fields'

=cut

sub next_FieldInfo {
    my $self = shift;
    $self->parse_data unless $self->data_parsed;    
    $self->{'_fieldinfo_it'} = $self->generate_iterator('fieldinfo')
        if (!exists $self->{'_fieldinfo_it'});
    $self->{'_fieldinfo_it'}->();
}

=head2 get_FieldInfo

 Title    : get_FieldInfo
 Usage    : my @fields = $info->get_FieldInfo;
 Function : returns list of FieldInfo objects
 Returns  : array (FieldInfo objects)
 Args     : none
 Notes    : only applicable for einfo. 

=cut

sub get_FieldInfo {
    my $self = shift;
    $self->parse_data unless $self->data_parsed;        
    return ref $self->{'_fieldinfo'} ? @{ $self->{'_fieldinfo'} } : return ();
}

*get_FieldInfos = \&get_FieldInfo;

=head2 next_LinkInfo

 Title    : next_LinkInfo
 Usage    : while (my $link = $info->next_LinkInfo) {...}
 Function : iterate through LinkInfo objects
 Returns  : LinkInfo object
 Args     : none
 Notes    : only applicable for einfo.  Uses callback() for filtering if defined
            for 'linkinfo'

=cut

sub next_LinkInfo {
    my $self = shift;
    $self->parse_data unless $self->data_parsed;    
    $self->{'_linkinfo_it'} = $self->generate_iterator('linkinfo')
        if (!exists $self->{'_linkinfo_it'});
    $self->{'_linkinfo_it'}->();
}

=head2 get_LinkInfo

 Title    : get_LinkInfo
 Usage    : my @links = $info->get_LinkInfo;
 Function : returns list of LinkInfo objects
 Returns  : array (LinkInfo objects)
 Args     : none
 Notes    : only applicable for einfo.  

=cut

sub get_LinkInfo {
    my $self = shift;
    $self->parse_data unless $self->data_parsed;        
    return ref $self->{'_linkinfo'} ? @{ $self->{'_linkinfo'} } : return ();
}

*get_LinkInfos = \&get_LinkInfo;

=head2 print_FieldInfo

 Title    : print_FieldInfo
 Usage    : $info->print_FieldInfo();
            $info->print_FieldInfo(-fh => $fh, -cb => $coderef);
 Function : prints link data for each FieldInfo object. The default is generated
            via FieldInfo::to_string
 Returns  : none
 Args     : [optional]
           -file : file to print to
           -fh   : filehandle to print to (cannot be used concurrently with file)
           -cb   : coderef to use in place of default print method.  This is
                   passed in a FieldInfo object
           -wrap : number of columns to wrap default text output to (def = 80)
 Notes    : only applicable for einfo.  If -file or -fh are not defined,
            prints to STDOUT

=cut

sub print_FieldInfo {
    my ($self, @args) = @_;
    $self->_print_handler(@args, -type => 'FieldInfo');
}

=head2 print_LinkInfo

 Title    : print_LinkInfo
 Usage    : $info->print_LinkInfo();
            $info->print_LinkInfo(-fh => $fh, -cb => $coderef);
 Function : prints link data for each LinkInfo object. The default is generated
            via LinkInfo::to_string
 Returns  : none
 Args     : [optional]
           -file : file to print to
           -fh   : filehandle to print to (cannot be used concurrently with file)
           -cb   : coderef to use in place of default print method.  This is passed
                   in a LinkInfo object
           -wrap : number of columns to wrap default text output to (def = 80)
 Notes    : only applicable for einfo.  If -file or -fh are not defined,
            prints to STDOUT

=cut

sub print_LinkInfo {
    my ($self, @args) = @_;
    $self->_print_handler(@args, -type => 'LinkInfo');
}

=head1 Bio::Tools::EUtilities::Link-related methods

=head2 next_LinkSet

 Title    : next_LinkSet
 Usage    : while (my $ls = $eutil->next_LinkSet {...}
 Function : iterate through LinkSet objects
 Returns  : LinkSet object
 Args     : none
 Notes    : only applicable for elink.  Uses callback() for filtering if defined
            for 'linksets'

=cut

sub next_LinkSet {
    my $self = shift;
    #$self->parse_data unless $self->data_parsed;
    if(!$self->data_parsed && !$self->is_lazy) {
        $self->parse_data;
    }
    $self->{'_linksets_it'} = $self->generate_iterator('linksets')
        if (!exists $self->{'_linksets_it'});
    $self->{'_linksets_it'}->();
}

=head2 get_LinkSets

 Title    : get_LinkSets
 Usage    : my @links = $info->get_LinkSets;
 Function : returns list of LinkSets objects
 Returns  : array (LinkSet objects)
 Args     : none
 Notes    : only applicable for elink.  

=cut

# add support for retrieval of data if lazy parsing is enacted

sub get_LinkSets {
    my $self = shift;
    if ($self->is_lazy) {
        $self->warn('get_LinkSets() not implemented when using lazy mode');
        return ();
    }
    $self->parse_data unless $self->data_parsed;
    return ref $self->{'_linksets'} ? @{ $self->{'_linksets'} } : return ();
}

=head2 print_LinkSets

 Title    : print_LinkSets
 Usage    : $info->print_LinkSets();
            $info->print_LinkSets(-fh => $fh, -cb => $coderef);
 Function : prints link data for each LinkSet object. The default is generated
            via LinkSet::to_string
 Returns  : none
 Args     : [optional]
           -file : file to print to
           -fh   : filehandle to print to (cannot be used concurrently with file)
           -cb   : coderef to use in place of default print method.  This is passed
                   in a LinkSet object
           -wrap : number of columns to wrap default text output to (def = 80)
 Notes    : only applicable for einfo.  If -file or -fh are not defined,
            prints to STDOUT

=cut

sub print_LinkSets {
    my ($self, @args) = @_;
    $self->_print_handler(@args, -type => 'LinkSet');
}

=head2 get_linked_databases

 Title    : get_linked_databases
 Usage    : my @dbs = $eutil->get_linked_databases
 Function : returns list of databases linked to in linksets
 Returns  : array of databases
 Args     : none
 Notes    : only applicable for elink.  Now defers to get_databases.

=cut

sub get_linked_databases {
    my $self = shift;
    return $self->get_databases if $self->eutil eq 'elink';
    return ();
}

=head1 Iterator- and callback-related methods

=cut

{
    my %VALID_ITERATORS = (
        'globalqueries' => 'globalqueries',
        'fieldinfo' =>  'fieldinfo',
        'fieldinfos' => 'fieldinfo',
        'linkinfo' =>  'linkinfo',
        'linkinfos' => 'linkinfo',
        'linksets' => 'linksets',
        'docsums' => 'docsums',
        'histories' => 'histories'
        );

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
            'globalqueries' - GlobalQuery objects
            'fieldinfo' or 'fieldinfos' - FieldInfo objects
            'linkinfo' or 'linkinfos' - LinkInfo objects in this layer
            'linksets' - LinkSet objects
            'docsums' - DocSum objects
            'histories' - HistoryI objects (Cookies, LinkSets)

=cut

sub rewind {
    my ($self, $arg) = ($_[0], lc $_[1]);
    my $eutil = $self->eutil;
    if ($self->is_lazy) {
        $self->warn('rewind() not implemented yet when running in lazy mode');
        return;
    }
    $arg ||= 'all';
    if (exists $VALID_ITERATORS{$arg}) {
        delete $self->{'_'.$arg.'_it'};
    } elsif ($arg eq 'all') {
        for my $it (values %VALID_ITERATORS){
            delete $self->{'_'.$it.'_it'} if
                exists $self->{'_'.$it.'_it'};
            map {$_->rewind('all')} $self->get_LinkSets;
            map {$_->rewind('all')} $self->get_DocSums;
        }
    }
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
            'fieldinfo' or 'fieldinfos' (the latter sounds clumsy, but I alias it JIC)
            'linkinfo' or 'linkinfos' (the latter sounds clumsy, but I alias it JIC)
            'linksets'
            'docsums'
            'histories'

 Note     : This function generates a simple coderef that one can use
            independently of the various next_* functions (in fact, the next_*
            functions use lazily created iterators generated via this method,
            while rewind() merely deletes them so they can be regenerated on the
            next call).

            A callback specified using callback() will be used to filter objects
            for any generated iterator. This behaviour is implemented for both
            normal and lazy iterator types and is the default. If you don't want
            this, make sure to reset any previously set callbacks via
            reset_callback() (which just deletes the code ref).  Note that setting
            callback() also changes the behavior of the next_* functions as the
            iterators are generated here (as described above); this is a feature
            and not a bug.

            'Lazy' iterators are considered an experimental feature and may be
            modified in the future. A 'lazy' iterator, which loops through and
            returns objects as they are created (instead of creating all data
            instances up front, then iterating through) is returned if the
            parser is set to 'lazy' mode. This mode is only present for elink
            and esummary output as they are the two formats parsed which can
            generate potentially thousands of individual objects (note efetch
            isn't parsed, so isn't counted). Use of rewind() with these
            iterators is not supported for the time being as we can't guarantee
            you can rewind(), as this depends on whether the data source is
            seek()able and thus 'rewindable'. We will add rewind() support at a
            later time which will work for 'seekable' data or possibly cached
            objects via Storable or BDB.

=cut

sub generate_iterator {
    my ($self, $obj) = @_;
    if (!$obj) {
        $self->throw('Must provide object type to iterate');
    } elsif (!exists $VALID_ITERATORS{$obj}) {
        $self->throw("Unknown object type [$obj]");
    }
    my $cb = $self->callback;
    if ($self->is_lazy) {
        my $type = $self->eutil eq 'esummary' ? '_docsums' : '_linksets';
        $self->{$type} = [];
        return sub {
            if (!@{$self->{$type}}) {
                $self->parse_chunk; # fill the queue
            }
            while (my $obj = shift @{$self->{$type}}) {
                if ($cb) {
                    ($cb->($obj)) ? return $obj : next;
                } else {
                    return $obj;
                }
            }
            undef;
        }
    } else {
        my $loc = '_'.$VALID_ITERATORS{$obj};
        my $index = $#{$self->{$loc}};
        my $current = 0;
        return sub {
            while ($current <= $index) {
                if ($cb) {
                    if (my $d = $cb->($self->{$loc}->[$current])) {
                        return $self->{$loc}->[$current++] }
                    else {
                        $current++;
                        next;
                    }
                } else {
                    return $self->{$loc}->[$current++]
                }
            }
            undef;
        }
    }
}

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
    my ($self, $cb) = @_;
    if ($cb) {
        delete $self->{'_cb'} if ($cb eq 'reset');
        return if ref $cb ne 'CODE';
        $self->{'_cb'} = $cb;
    }
    return $self->{'_cb'};
}

# Object printing methods

{
    my $DEF_HANDLER = sub {
        my $obj = shift;
        return $obj->to_string."\n";
    };
    
    my %HANDLER = (
        'DocSum'        => 1,
        'FieldInfo'     => 1,
        'LinkInfo'      => 1,
        'GlobalQuery'   => 1,
        'LinkSet'       => 1,
        'all'           => 1,
                   );
    
    sub _print_handler {
        my $self = shift;
        my ($file, $fh, $cb, $wrap, $type, $all) = $self->_rearrange([qw(FILE FH CB WRAP TYPE ALL)], @_);
        $type ||= 'all';
        
        # default formatting delegates to_string
        if (!$cb) {
            $self->throw("Type $type not registered with print handler, exiting...")
                if !exists($HANDLER{$type});
            $cb = $DEF_HANDLER;
        } else {
            $self->throw("Callback must be a code reference") if ref $cb ne 'CODE';
        }
        
        $file ||= $fh;
        $self->throw("Have defined both file and filehandle; only use one!") if $file && $fh;
        my $io = ($file) ? Bio::Root::IO->new(-input => $file, -flush => 1) :
                 Bio::Root::IO->new(-flush => 1); # defaults to STDOUT
                 
        if ($type eq 'all') {
            my $string = $cb->($self);
            $io->_print($string) if $string;
        } else {
            # set up iterator
            my $it = "next_$type";
            $self->throw("Unknown iterator method $it") unless $self->can($it);
            while (my $obj = $self->$it) {
                my $string = $cb->($obj);
                $io->_print($string) if $string;
            }
        }
        $io->close;
    }
}

# Private methods

sub _seekable {
    return shift->{'_seekable'}
}

# fixes odd bad XML issue espell data (still present 6-24-07)

sub _fix_espell {
    my ($self, $response) = @_;
    my $temp;
    my $type = ref($response);
    if ($type eq 'GLOB') {
        $temp .= $_ for <$response>;
    } elsif ($type eq 'HTTP::Response') {
        $temp = $response->content;
    } else {
        $self->throw("Unrecognized ref type $type");
    }
    if ($temp =~ m{^<html>}) {
        $self->throw("NCBI espell nonrecoverable error: HTML content returned")
    }
    $temp =~ s{<ERROR>(.*?)<ERROR>}{<ERROR>$1</ERROR>};
    return $temp;
}

sub _load_eutil_module {
    my ($self, $class) = @_;
    my $ok;
    my $module = "Bio::Tools::EUtilities::" . $class;

    eval {
        $ok = $self->_load_module($module);
    };
    if ( $@ ) {
        print STDERR <<END;
$self: data module $module cannot be found
Exception $@
For more information about the EUtilities system please see the EUtilities docs. 
END
       ;
    }
    return $ok;
}

1;
