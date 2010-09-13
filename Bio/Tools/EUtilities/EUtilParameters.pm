#
# BioPerl module for Bio::Tools::EUtilities::EUtilParameters
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Chris Fields <cjfields at bioperl dot org>
#
# Copyright Chris Fields
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=head1 NAME

Bio::Tools::EUtilities::EUtilParameters - Manipulation of NCBI eutil-based
parameters for remote database requests.

=head1 SYNOPSIS

 # Bio::Tools::EUtilities::EUtilParameters implements Bio::ParameterBaseI

 my @params = (-eutil => 'efetch',
              db => 'nucleotide',
              id => \@ids,
              email => 'me@foo.bar',
              retmode => 'xml');

 my $p = Bio::Tools::EUtilities::EUtilParameters->new(@params);

 if ($p->parameters_changed) {
                              # ...
                             } # state information

 $p->set_parameters(@extra_params); # set new NCBI parameters, leaves others preset

 $p->reset_parameters(@new_params); # reset NCBI parameters to original state

 $p->to_string(); # get a URI-encoded string representation of the URL address

 $p->to_request(); # get an HTTP::Request object (to pass on to LWP::UserAgent)

=head1 DESCRIPTION

Bio::Tools::EUtilities::EUtilParameters is-a Bio::ParameterBaseI implementation
that allows simple manipulation of NCBI eutil parameters for CGI-based queries.
SOAP-based methods may be added in the future.

For simplicity parameters do not require dashes when passed and do not need URI
encoding (spaces are converted to '+', symbols encoded, etc). Also, the
following extra parameters can be passed to the new() constructor or via
set_parameters() or reset_parameters():

  eutil - the eutil to be used. The default is 'efetch' if not set.
  correspondence - Flag for how IDs are treated. Default is undef (none).
  history - a Bio::Tools::EUtilities::HistoryI object. Default is undef (none).

At this point minimal checking is done for potential errors in parameter
passing, though these should be easily added in the future when necessary.

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

package Bio::Tools::EUtilities::EUtilParameters;
use strict;
use warnings;

use base qw(Bio::Root::Root Bio::ParameterBaseI);
use URI;
use HTTP::Request;
use Bio::Root::IO;

# eutils only has one hostbase URL

# mode : GET or POST (HTTP::Request)
# location : CGI location
# params : allowed parameters for that eutil
my %MODE = (
    'einfo'     => {
        'mode'     => ['GET'],
        'location' => 'einfo.fcgi',
        'params'   => [qw(db tool email)],
                   },
    'epost'     => {
        'mode'     => ['POST','GET'],
        'location' => 'epost.fcgi',
        'params'   => [qw(db retmode id tool email WebEnv query_key)],
                   },
    'efetch'    => {
        'mode'     => ['GET','POST'],
        'location' => 'efetch.fcgi',
        'params'   => [qw(db retmode id retmax retstart rettype strand seq_start
                       seq_stop complexity report tool email WebEnv query_key)],
                   },
    'esearch'   => {
        'mode'     => ['GET','POST'],
        'location' => 'esearch.fcgi',
        'params'   => [qw(db retmode usehistory term field reldate mindate
                       maxdate datetype retmax retstart rettype sort tool email
                       WebEnv query_key)],
                   },
    'esummary'  => {
        'mode'     => ['GET','POST'],
        'location' => 'esummary.fcgi',
        'params'   => [qw(db retmode id retmax retstart rettype tool email
                       WebEnv query_key)],
                   },
    'elink'     => {
        'mode'     => ['GET','POST'],
        'location' => 'elink.fcgi',
        'params'   => [qw(db retmode id reldate mindate maxdate datetype term 
                    dbfrom holding cmd version tool email linkname WebEnv
                    query_key)],
                   },
    'egquery'   => {
        'mode'     => ['GET','POST'],
        'location' => 'egquery.fcgi',
        'params'   => [qw(term retmode tool email)],
                   },
    'espell'    => {
        'mode'     => ['GET','POST'],
        'location' => 'espell.fcgi',
        'params'   => [qw(db retmode term tool email )],
                   }
);

my @PARAMS;

# generate getter/setters (will move this into individual ones at some point)

BEGIN {
    @PARAMS = qw(db id email retmode rettype usehistory term field tool
    reldate mindate maxdate datetype retstart retmax sort seq_start seq_stop
    strand complexity report dbfrom cmd holding version linkname WebEnv
    query_key);
    for my $method (@PARAMS) {
        eval <<END;
sub $method {
    my (\$self, \$val) = \@_;
    if (defined \$val) {
        if ((!defined \$self->{'_$method'}) ||
        (defined \$self->{'_$method'} && \$self->{'_$method'} ne \$val)) {
            \$self->{'_statechange'} = 1;
            \$self->{'_$method'} = \$val;
        }
    }
    return \$self->{'_$method'};
}
END
    }
}

sub new {
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(@args);
    my ($retmode) = $self->_rearrange(["RETMODE"],@args);
    # order is important here, eutil must be set first so that proper error
    # checking occurs for the later attributes
    $self->_set_from_args(\@args,
        -methods => [@PARAMS, qw(eutil history correspondence id_file request_mode)]);
    $self->eutil() || $self->eutil('efetch');
    $self->tool() || $self->tool('BioPerl');
    # set default retmode if not explicitly set    
    $self->set_default_retmode if (!$retmode);
    $self->{'_statechange'} = 1;
    return $self;
}

=head1 Bio::ParameterBaseI implemented methods

=head2 set_parameters

 Title   : set_parameters
 Usage   : $pobj->set_parameters(@params);
 Function: sets the NCBI parameters listed in the hash or array
 Returns : None
 Args    : [optional] hash or array of parameter/values.  
 Note    : This sets any parameter passed but leaves previously set data alone.
           In addition to regular eutil-specific parameters, you can set the
           following:

           -eutil    - the eUtil to be used (default 'efetch')
           -history  - pass a HistoryI-implementing object, which
                       sets the WebEnv, query_key, and possibly db and linkname
                       (the latter two only for LinkSets)
           -correspondence - Boolean flag, set to TRUE or FALSE; indicates how
                       IDs are to be added together for elink request where
                       ID correspondence might be needed
                       (default 0)

=cut

sub set_parameters {
    my ($self, @args) = @_;
    # allow automated resetting; must check to ensure that retmode isn't explicitly passed
    my ($newmode,$file) = $self->_rearrange([qw(RETMODE ID_FILE)],@args);
    $self->_set_from_args(\@args, -methods => [@PARAMS, qw(eutil correspondence history)]);
    # set default retmode if not explicitly passed
    $self->set_default_retmode unless $newmode;
    $file && $self->id_file($file);
    return;
}

=head2 reset_parameters

 Title   : reset_parameters
 Usage   : resets values
 Function: resets parameters to either undef or value in passed hash
 Returns : none
 Args    : [optional] hash of parameter-value pairs
 Note    : This sets any parameter passed, but resets all others (deletes them).
           In addition to regular eutil-specific parameters, you can set the
           following:

           -eutil    - the eUtil to be used (default 'efetch')
           -history  - pass a HistoryI-implementing object, which
                       sets the WebEnv, query_key, and possibly db and linkname
                       (the latter two only for LinkSets)
           -correspondence - Boolean flag, set to TRUE or FALSE; indicates how
                       IDs are to be added together for elink request where
                       ID correspondence might be needed
                       (default 0)

=cut

sub reset_parameters {
    my ($self, @args) = @_;
    # is there a better way of doing this?  probably, but this works...
    my ($retmode,$file) = $self->_rearrange([qw(RETMODE ID_FILE)],@args);
    map { defined $self->{"_$_"} && undef $self->{"_$_"} } (@PARAMS, qw(eutil correspondence history_cache request_cache));
    $self->_set_from_args(\@args, -methods => [@PARAMS, qw(eutil correspondence history)]);
    $self->eutil() || $self->eutil('efetch');
    $self->set_default_retmode unless $retmode;
    $file && $self->id_file($file);
    $self->{'_statechange'} = 1;
}

=head2 carryover

 Title    : carryover
 Usage    : $obj->carryover(qw(email tool db))
 Function : Carries over the designated parameters when using reset_parameters()
 Returns  : a list of carried-over parameters
 Args     : An array reference of parameters to carry over, followed optionally
            by the mode ('add' or 'delete', indicating whether to append to or
            remove the specified values passed in). To clear all values, pass in
            an empty array reference (the mode in this case doesn't matter).
            
            In addition to the normal eUtil-specific parameters, the following
            additional parameters are allowed:
            
            -eutil    - the eUtil to be used (default 'efetch')
            -history  - pass a HistoryI-implementing object, which
                       sets the WebEnv, query_key, and possibly db and linkname
                       (the latter two only for LinkSets)
            -correspondence - Boolean flag, set to TRUE or FALSE; indicates how
                       IDs are to be added together for elink request where
                       ID correspondence might be needed
                       (default 0)
 Default  : None (no carried over parameters)
 Status   : NYI (dev in progress, carry on, nothing to see here)
 
=cut

sub carryover {
    my ($self, $params, $mode) = @_;
    my %allowed = map {$_ => 1} (@PARAMS, qw(eutil history correspondence));
    if ($params) {
        $self->throw("Must pass in an array ref of parameters") unless
            ref($params) eq 'ARRAY';
        my $mode ||= 'add';
        $self->throw("Mode must be 'add' or 'delete'") unless $mode eq 'add' || $mode eq 'delete';
        if (!scalar(@$params)) { # empty array ref
            $self->{_carryover} = {};
        } else {
            for my $p (@$params) {
                if (!exists $allowed{$p}) {
                    $self->warn("$p is not a recognized eutil parameter");
                    next;
                }
                if ($mode eq 'add') {
                    $self->{_carryover}->{$p} = 1;
                } else {
                    delete $self->{_carryover}->{$p} if exists
                        $self->{_carryover}->{$p};
                }
            }
        }
    }
    sort keys %{$self->{_carryover}} || ();
}

sub _reset_except_carryover {
    my $self = shift;
    #for my $p (@PARAMS, qw(eutil correspondence history_cache request_cache)) {
    #    undef $self->{"_$p"} if defined $self->{"_$p"};
    #}
}

=head2 request_mode

 Title    : request_mode
 Usage    : $obj->request_mode
 Function : get/set the mode for the user agent to use for generating a request
 Returns  : either a preset mode (checked against the eutil) or a best-possible
            option based upon the currently-set parameters
 Args     : 
 Status   :
 
=cut

sub request_mode {
    my ($self, $mode) = @_;
    $mode = uc $mode if defined $mode;
    my $eutil = $self->eutil;
    if ($mode) {
        my %valid = map {$_ => 1} @{$MODE{$eutil}{mode}};
        $self->throw("Mode $mode not supported for $eutil") unless
            exists $valid{$mode};
        $self->{_request_mode} = $mode;
    }
    return $self->{_request_mode} if $self->{_request_mode};
    # let's try to make this a bit smarter...
    
    # If not explicitly set, in cases where
    # the number of IDs is greater than 200, or the search term is longer than
    # 200, use POST when available
    
    if (scalar(@{$MODE{$eutil}{mode}}) > 1) { # allows both GET and POST
        my ($id, $term) = ($self->id || [], $self->term || '');
        if (scalar(@$id) > 200 || CORE::length($term) > 300) {
            return 'POST'
        }
    }
    # otherwise, fallback to default
    $MODE{$eutil}{mode}[0]; # first is default    
}

=head2 parameters_changed

 Title   : parameters_changed
 Usage   : if ($pobj->parameters_changed) {...}
 Function: Returns TRUE if parameters have changed
 Returns : Boolean (0 or 1)
 Args    : [optional] Boolean

=cut

sub parameters_changed {
    my ($self) = @_;
    $self->{'_statechange'};
}

=head2 available_parameters

 Title   : available_parameters
 Usage   : @params = $pobj->available_parameters()
 Function: Returns a list of the available parameters
 Returns : Array of available parameters (no values)
 Args    : [optional] A string with the eutil name (for returning eutil-specific
           parameters)

=cut

sub available_parameters {
    my ($self, $type) = @_;
    $type ||= 'all';
    if ($type eq 'all') {
        return @PARAMS;
    } else {
        $self->throw("$type parameters not supported") if !exists $MODE{$type};
        return @{$MODE{$type}->{params}};
    }
}

=head2 get_parameters

 Title   : get_parameters
 Usage   : @params = $pobj->get_parameters;
           %params = $pobj->get_parameters;
 Function: Returns list of key/value pairs, parameter => value
 Returns : Flattened list of key-value pairs. All key-value pairs returned,
           though subsets can be returned based on the '-type' parameter. Data
           originally set as an array ref are returned based on whether the
           '-join_id' flag is set (default is the same array ref).
 Args    : -type : the eutil name (Default: returns all).  Use of '-list'
                    supercedes this
           -list : array ref of specific parameters
           -join_ids : Boolean; join IDs based on correspondence (Default: no join)

=cut

sub get_parameters {
    my ($self, @args) = @_;
    my ($type, $list, $join) = $self->_rearrange([qw(TYPE LIST JOIN_IDS)], @args);
    $self->throw("Parameter list not an array ref") if $list && ref $list ne 'ARRAY';
    $type ||= '';
    my @final = $list ? grep {$self->can($_)} @{$list} : $self->available_parameters($type);
    my @p;
    for my $param (@final) {
        if ($param eq 'id' && $self->id && $join) {
            my $id = $self->id;
            if ($self->correspondence && $self->eutil eq 'elink') {
                for my $id_group (@{ $id }) {
                    if (ref($id_group) eq 'ARRAY') {
                        push @p, ('id' => join(q(,), @{ $id_group }));
                    }
                    elsif (!ref($id_group)) {
                        push @p, ('id'  => $id_group);
                    }
                    else {
                        $self->throw("Unknown ID type: $id_group");
                    }
                }
            } else {
                # add a check for undef
                push @p, ref $id eq 'ARRAY' ?
                ($param => join(',', grep {defined($_)} @{ $id })):
                ($param => $id);
            }
        }
        elsif ($param eq 'db' && $self->db && $join) {
            my $db = $self->db;
            push @p, (ref $db eq 'ARRAY') ? 
                ($param => join(',', @{ $db })) :
                ($param => $db) ;
        }
        else {
            push @p, ($param => $self->{"_$param"}) if defined $self->{"_$param"};
        }
    }
    return @p;
}

=head1 Implementation-specific to_* methods

=head2 to_string

 Title   : to_string
 Usage   : $string = $pobj->to_string;
 Function: Returns string (URL only in this case)
 Returns : String (URL only for now)
 Args    : [optional] 'all'; build URI::http using all parameters
           Default : Builds based on allowed parameters (presence of history data
           or eutil type in %MODE).
 Note    : Changes state of object.  Absolute string

=cut

sub to_string {
    my ($self, @args) = @_;
    # calling to_uri changes the state
    if ($self->parameters_changed || !defined $self->{'_string_cache'}) {
        my $string = $self->to_request(@args)->uri->as_string;
        $self->{'_statechange'} = 0;
        $self->{'_string_cache'} = $string;
    }
    return $self->{'_string_cache'};
}

=head2 to_request

 Title   : to_request
 Usage   : $uri = $pobj->to_request;
 Function: Returns HTTP::Request object
 Returns : HTTP::Request
 Args    : [optional] 'all'; builds request using all parameters
           Default : Builds based on allowed parameters (presence of history data
           or eutil type in %MODE).
 Note    : Changes state of object (to boolean FALSE).  Used for CGI-based GET/POST
 TODO    : esearch, esummary, elink now accept POST for batch submission
           (something NCBI apparently allowed but didn't advertise). Should we
           switch most of these to utilize POST instead, or make it dep on the
           number of submitted IDs?

=cut

sub to_request {
    my ($self, $type) = @_;
    if ($self->parameters_changed || !defined $self->{'_request_cache'}) {
        my $eutil = $self->eutil;
        $self->throw("No eutil set") if !$eutil;
        #set default retmode
        $type ||= $eutil;
        my ($location, $mode) = ($MODE{$eutil}->{location}, $self->request_mode);
        my $request;
        my $uri = URI->new($self->url_base_address . $location);
        if ($mode eq 'GET') {
            $uri->query_form($self->get_parameters(-type => $type, -join_ids => 1) );
            $request = HTTP::Request->new($mode => $uri);
            $self->{'_request_cache'} = $request;
        } elsif ($mode eq 'POST') {
            $request = HTTP::Request->new($mode => $uri->as_string);
            $uri->query_form($self->get_parameters(-type => $type, -join_ids => 1) );
            $request->content_type('application/x-www-form-urlencoded');
            $request->content($uri->query);
            $self->{'_request_cache'} = $request;
        } else {
            $self->throw("Unrecognized request mode: $mode");
        }
        $self->{'_statechange'} = 0;
        $self->{'_request_cache'} = $request;
    }
    return $self->{'_request_cache'};
}

=head1 Implementation specific-methods

=head2 eutil

 Title   : eutil
 Usage   : $p->eutil('efetch')
 Function: gets/sets the eutil for this set of parameters
 Returns : string (eutil)
 Args    : [optional] string (eutil)
 Throws  : '$eutil not supported' if eutil not present
 Note    : This does not reset retmode to the default if called directly.

=cut

sub eutil {
    my ($self, $eutil) = @_;
    if ($eutil) {
        $self->throw("$eutil not supported") if !exists $MODE{$eutil};
        if (!defined $self->{'_eutil'} || ($self->{'_eutil'} && $self->{'_eutil'} ne $eutil)) {
            $self->{'_eutil'} = $eutil;
            $self->{'_statechange'} = 1;
        }
    }
    return $self->{'_eutil'};
}

=head2 history

 Title   : history
 Usage   : $p->history($history);
 Function: gets/sets the history object to be used for these parameters
 Returns : Bio::Tools::EUtilities::HistoryI (if set)
 Args    : [optional] Bio::Tools::EUtilities::HistoryI 
 Throws  : Passed something other than a Bio::Tools::EUtilities::HistoryI 
 Note    : This overrides WebEnv() and query_key() settings when set.  This
           caches the last history object passed and returns like a Get/Set

=cut

sub history {
    my ($self, $history) = @_;
    if ($history) {
        $self->throw('Not a Bio::Tools::EUtilities::HistoryI object!') if
            !$history->isa('Bio::Tools::EUtilities::HistoryI');
        my ($webenv, $qkey) = $history->history;
        $self->WebEnv($webenv);
        $self->query_key($qkey);
        $self->{'_statechange'} = 1;
        $self->{'_history_cache'} = $history;
    }
    return $self->{'_history_cache'};
}

=head2 correspondence

 Title   : correspondence
 Usage   : $p->correspondence(1);
 Function: Sets flag for posting IDs for one-to-one correspondence
 Returns : Boolean
 Args    : [optional] boolean value

=cut

sub correspondence {
    my ($self, $corr) = @_;
    if (defined $corr) {
        $self->{'_correspondence'} = $corr;
        $self->{'_statechange'} = 1;
    }
    return $self->{'_correspondence'};
}

=head2 id_file

 Title   : id_file
 Usage   : $p->id_file('<foo');
 Function: convenience method; passes in file containing a list of IDs for
           searches (one per line), sets id() to list
 Returns : none
 Args    : either string indicating file to use, a file handle, or an IO::Handle
           object
 Note    : use of this overrides concurrent use of the '-id' parameter when both
           are passed.  The filename is not retained, merely parsed for IDs.

=cut

sub id_file {
    my ($self, $file) = @_;
    if ($file) {
        # do this in a way that allows file, fh, IO::Handle
        my $io = $self->_io;
        $io->_initialize_io(-input => $file);
        my @ids;
        while (my $line = $io->_readline) {
            chomp $line;
            push @ids, $line;
        }
        $self->_io->close;
        $self->id(\@ids);
    }
}

=head2 url_base_address

 Title   : url_base_address
 Usage   : $address = $p->url_base_address();
 Function: Get URL base address
 Returns : String
 Args    : None in this implementation; the URL is fixed

=cut

{
    my $HOSTBASE = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/';
    
    sub url_base_address {
        my ($self, $address) = @_;
        return $HOSTBASE;
    }
}

=head2 set_default_retmode

 Title   : set_default_retmode
 Usage   : $p->set_default_retmode();
 Function: sets retmode to default value specified by the eutil() and the value
           in %NCBI_DATABASE (for efetch only) if called
 Returns : none
 Args    : none

=cut

{
    # default retmode if one is not supplied
    my %NCBI_DATABASE = (
        'protein'          => 'text',
        'nucleotide'       => 'text',
        'nuccore'          => 'text',
        'nucgss'           => 'text',
        'nucest'           => 'text',
        'structure'        => 'text',
        'genome'           => 'text',
        'gene'             => 'asn1',
        'journals'         => 'text',
    );

    sub set_default_retmode {
        my $self = shift;
        if ($self->eutil eq 'efetch') {
            my $db = $self->db || return; # assume retmode will be set along with db
            my $mode = exists $NCBI_DATABASE{$db} ? $NCBI_DATABASE{$db} : 'xml';
            $self->retmode($mode);
        } else {
            $self->retmode('xml');
        }
    }
}

sub _io {
    my $self = shift;
    if (!defined $self->{'_io'}) {
        $self->{'_io'} = Bio::Root::IO->new();
    }
    return $self->{'_io'};
}

1;
