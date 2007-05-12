# $Id$
#
# BioPerl module for Bio::DB::EUtilParameters
#
# Cared for by Chris Fields <cjfields at uiuc dot edu>
#
# Copyright Chris Fields
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=head1 NAME

Bio::DB::EUtilParameters - Manipulation of NCBI eutil-based parameters for
remote database requests.

=head1 SYNOPSIS

 # Bio::DB::EUtilParameters implements Bio::ParameterBaseI
 
 my @params = (-eutil => 'efetch',
              db => 'nucleotide',
              id => \@ids,
              email => 'me@foo.bar',
              retmode => 'xml');
 
 my $p = Bio::DB::EUtilParameters->new(@params);
 
 if ($p->parameters_changed) {...} # state information
 
 $p->set_parameters(@extra_params); # set new NCBI parameters, leaves others preset
 
 $p->reset_parameters(@new_params); # reset NCBI parameters to original state
 
 $p->to_string(); # get a URI-encoded string representation of the URL address
 
 $p->to_request(); # get an HTTP::Request object (to pass on to LWP::UserAgent)
 
=head1 DESCRIPTION

Bio::DB::EUtilParameters is-a Bio::ParameterBaseI implementation that allows
simple manipulation of NCBI eutil parameters for CGI-based queries. SOAP-based
methods may be added in the future.

For simplicity parameters do not require dashes when passed and do not need URI
encoding (spaces are converted to '+', symbols encoded, etc). Also, the
following extra parameters can be passed to the new() constructor or via
set_parameters() or reset_parameters():

  eutil - the eutil to be used. The default is 'efetch' if not set.
  correspondence - Flag for how IDs are treated. Default is undef (none).
  cookie - a Bio::DB::EUtilities::Cookie object. Default is undef (none).

At this point minimal checking is done for potential errors in parameter
passing, though these should be easily added in the future when necessary.

=head1 TODO

Possibly integrate SOAP-compliant methods. SOAP::Lite may be undergoing an
complete rewrite so I'm hesitant about adding this in immediately.

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

package Bio::DB::EUtilParameters;
use strict;
use warnings;

use base qw(Bio::Root::Root Bio::ParameterBaseI);
use URI;
use HTTP::Request;

# eutils only has one hostbase URL
my $HOSTBASE = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/';

# mode : GET or POST (HTTP::Request)
# location : CGI location
# params : allowed parameters for that eutil
my %MODE = (
    'einfo'     => {
        'mode'     => 'get',
        'location' => 'einfo.fcgi',
        'params'   => [qw(db retmode tool email)],
                   },
    'epost'     => {
        'mode'     => 'post',
        'location' => 'epost.fcgi',
        'params'   => [qw(db retmode id tool email)],
                   },
    'efetch'    => {
        'mode'     => 'get',
        'location' => 'efetch.fcgi',
        'params'   => [qw(db retmode id retmax retstart rettype strand seq_start
                       seq_stop complexity report tool email )],
                   },
    'esearch'   => {
        'mode'     => 'get',
        'location' => 'esearch.fcgi',
        'params'   => [qw(db retmode usehistory term field reldate mindate
                       maxdate datetype retmax retstart rettype sort tool email)],
                   },
    'esummary'  => {
        'mode'     => 'get',
        'location' => 'esummary.fcgi',
        'params'   => [qw(db retmode id retmax retstart rettype tool email )],
                   },
    'elink'     => {
        'mode'     => 'get',
        'location' => 'elink.fcgi',
        'params'   => [qw(db retmode id reldate mindate maxdate datetype term 
                    dbfrom holding cmd version tool email)],
                   },
    'egquery'   => {
        'mode'     => 'get',
        'location' => 'egquery.fcgi',
        'params'   => [qw(term retmode tool email)],
                   },
    'espell'    => {
        'mode'     => 'get',
        'location' => 'espell.fcgi',
        'params'   => [qw(db retmode term tool email )],
                   }
);

# used only if cookie is present
my @COOKIE_PARAMS = qw(db sort seq_start seq_stop strand complexity rettype
    retstart retmax cmd linkname retmode WebEnv query_key);            

# default retmode if one is not supplied
my %NCBI_DATABASE = (
    'pubmed'           => 'xml',
    'protein'          => 'text',
    'nucleotide'       => 'text',
    'nuccore'          => 'text',
    'nucgss'           => 'text',
    'nucest'           => 'text',
    'structure'        => 'text',
    'genome'           => 'text',
    'books'            => 'xml',
    'cancerchromosomes'=> 'xml',
    'cdd'              => 'xml',
    'domains'          => 'xml',
    'gene'             => 'asn1',
    'genomeprj'        => 'xml',
    'gensat'           => 'xml',
    'geo'              => 'xml',
    'gds'              => 'xml',
    'homologene'       => 'xml',
    'journals'         => 'text',
    'mesh'             => 'xml',
    'ncbisearch'       => 'xml',
    'nlmcatalog'       => 'xml',
    'omia'             => 'xml',
    'omim'             => 'xml',
    'pmc'              => 'xml',
    'popset'           => 'xml',
    'probe'            => 'xml',
    'pcassay'          => 'xml',
    'pccompound'       => 'xml',
    'pcsubstance'      => 'xml',
    'snp'              => 'xml',
    'taxonomy'         => 'xml',
    'unigene'          => 'xml',
    'unists'           => 'xml',
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
        \$self->{'_statechange'} = 1 if (!defined \$self->{'_$method'}) ||
            (defined \$self->{'_$method'} && \$self->{'_$method'} ne \$val);
        \$self->{'_$method'} = \$val;
    }
    return \$self->{'_$method'};
}
END
    }
}

sub new {
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(@args);
    $self->_set_from_args(\@args,
        -methods => [@PARAMS, qw(eutil cookie correspondence)]);
    # set default retmode if not explicitly set
    $self->eutil() || $self->eutil('efetch');
    $self->_set_default_retmode if (!$self->retmode);
    $self->{'_statechange'} = 1;  
    return $self;
}

=head1 Bio::ParameterBaseI implemented methods

=head2 

 Title   : set_parameters
 Usage   : $pobj->set_parameters(%params);
 Function: sets the NCBI parameters listed in the hash or array
 Returns : None
 Args    : [optional] hash or array of parameter/values.  
 Note    : This sets any parameter (i.e. doesn't screen them using $MODE or via
           set cookies).  

=cut

sub set_parameters {
    my ($self, @args) = @_;
    $self->_set_from_args(\@args, -methods => [@PARAMS]);
}

=head2 

 Title   : reset_parameters
 Usage   : resets values
 Function: resets parameters to either undef or value in passed hash
 Returns : none
 Args    : [optional] hash of parameter-value pairs
 Note    : this also resets eutil(), correspondence(), and the cookie and request
           cache

=cut

sub reset_parameters {
    my ($self, @args) = @_;
    # is there a better way of doing this?  probably, but this works
    for my $param (@PARAMS, qw(eutil correspondence cookie_cache request_cache)) {
        defined $self->{"_$param"} && undef $self->{"_$param"};
    }
    $self->_set_from_args(\@args, -methods => [@PARAMS, qw(eutil correspondence cookie)]);
    $self->eutil() || $self->eutil('efetch');
    $self->_set_default_retmode if (!$self->retmode);
    $self->{'_statechange'} = 1;
}

=head2 

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

=head2 

 Title   : available_parameters
 Usage   : @params = $pobj->available_parameters()
 Function: Returns a list of the available parameters
 Returns : Array of available parameters (no values)
 Args    : [optional] A string; either eutil name (for returning eutil-specific
           parameters) or 'cookie' (for those parameters allowed when retrieving
           data stored on the remote server using a 'Cookie').  

=cut

sub available_parameters {
    my ($self, $type) = @_;
    $type ||= 'all';
    if ($type eq 'all') {
        return @PARAMS;
    } elsif ($type eq 'cookie') {
        return @COOKIE_PARAMS;
    } else {
        $self->throw("$type parameters not supported") if !exists $MODE{$type};
        return @{$MODE{$type}->{params}};
    }
}

=head2 

 Title   : get_parameters
 Usage   : @params = $pobj->get_parameters;
           %params = $pobj->get_parameters;
 Function: Returns list of key/value pairs, parameter => value
 Returns : Flattened list of key-value pairs. IDs are returned based on the
           correspondence value (a string joined by commas or as an array ref).
 Args    : -type : the eutil name or 'cookie', for returning a subset of
                parameters (Default: returns all)
                
           -join_ids : Boolean; join IDs based on correspondence (Default: no join)

=cut

sub get_parameters {
    my ($self, @args) = @_;
    my ($type, $join) = $self->_rearrange([qw(TYPE JOIN_IDS)], @args);
    $type ||= '';
    my @final = $self->available_parameters($type);
    my @p;
    for my $param (@final) {
        if ($param eq 'id' && $join) {
            if ($self->correspondence && $self->eutil eq 'elink') {
                for my $id_group (@{ $self->id }) {
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
                push @p, ($param => join(',', @{ $self->id }));
            }
        } elsif ($param eq 'retmode' && !$self->retmode) {
            
        } else {
            push @p, ($param => $self->{"_$param"}) if defined $self->{"_$param"};
        }
    }
    return @p;
}

=head2 

 Title   : to_string
 Usage   : $string = $pobj->to_string;
 Function: Returns string (URL only in this case)
 Returns : String (URL only for now)
 Args    : [optional] 'all'; build URI::http using all parameters
           Default : Builds based on allowed parameters (presence of cookie data
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

=head2 

 Title   : to_request
 Usage   : $uri = $pobj->to_request;
 Function: Returns HTTP::Request object
 Returns : HTTP::Request
 Args    : [optional] 'all'; builds request using all parameters
           Default : Builds based on allowed parameters (presence of cookie data
           or eutil type in %MODE).
 Note    : Changes state of object.  Used for CGI-based GET/POST
 
=cut

sub to_request {
    my ($self, $type) = @_;
    if ($self->parameters_changed || !defined $self->{'_uri_cache'}) {
        my $eutil = $self->eutil;
        $self->throw("No eutil set") if !$eutil;
        #set default retmode
        my $cookie = ($self->cookie) ? 1 : 0;
        $type ||= ($cookie) ? 'cookie' : $eutil;
        my $uri = URI->new($HOSTBASE . $MODE{$eutil}->{location});
        $uri->query_form($self->get_parameters(-type => $type, -join_ids => 1) );
        my $method = ($eutil eq 'epost') ? 'POST' : 'GET';
        my $request = HTTP::Request->new($method => $uri);
        $self->{'_statechange'} = 0;
        $self->{'_request_cache'} = $request;
    }
    return $self->{'_request_cache'};
}

=head1 Implementation specific-methods

=head2 

 Title   : eutil
 Usage   : $p->eutil('efetch')
 Function: gets/sets the eutil for this set of parameters
 Returns : string (eutil)
 Args    : [optional] string (eutil)
 Throws  : '$eutil not supported' if eutil not present

=cut

sub eutil {
    my ($self, $eutil) = @_;
    if ($eutil) {
        $self->throw("$eutil not supported") if !exists $MODE{$eutil};
        $self->{'_eutil'} = $eutil;
        $self->{'_statechange'} = 1;
    }
    return $self->{'_eutil'};
}

=head2 

 Title   : cookie
 Usage   : $p->cookie($cookie);
 Function: gets/sets the cookie (history) to be used for these parameters
 Returns : Bio::DB::EUtilities::Cookie (if set)
 Args    : [optional] Bio::DB::EUtilities::Cookie 
 Throws  : Passed something other than a Bio::DB::EUtilities::Cookie
 Note    : This overrides WebEnv() and query_key() if set

=cut

# cookie not changed over to ParameterBaseI yet...

sub cookie {
    my ($self, $cookie) = @_;
    if ($cookie) {
        $self->throw('Not a Bio::DB::EUtilities::Cookie object!') if
            !$cookie->isa('Bio::DB::EUtilities::Cookie');
        my ($webenv, $qkey) = @{$cookie->cookie};
        $webenv     && $self->WebEnv($webenv);
        $qkey       && $self->query_key($qkey);

        #TODO: set db(), dbfrom() based on eutil

        $self->{'_statechange'} = 1;
        $self->{'_cookie_cache'} = $cookie;
    }
    return $self->{'_cookie_cache'};
}

=head2 

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

# Title   : _set_default_retmode
# Usage   : $p->_set_default_retmode();
# Function: sets retmode to default value if called
# Returns : none
# Args    : none

sub _set_default_retmode {
    my $self = shift;
    if ($self->eutil eq 'efetch') {
        my $db = $self->db || $self->throw('No database defined for efetch!');
        $self->throw('Database $db not recognized') if !exists $NCBI_DATABASE{$db};
        # set efetch-based retmode
        $self->retmode($NCBI_DATABASE{$db});
    } else {
        $self->retmode('xml');
    }
}

1;

