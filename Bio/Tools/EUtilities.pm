# $Id$
#
# BioPerl module for Bio::Tools::EUtilities
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

  # esearch, egquery
  
    $term = $parser->get_term(); # returns array or array ref of IDs
  
  # elink, egquery
  
    $db = $parser->get_database(); # returns database
  
  # Other methods are focused into four subclasses:
  
  # Bio::Tools::EUtilities::Query (esearch, egquery, espell data)
  # eutil data centered on use of search terms
    
    my $ct = $parser->get_count; # uses optional database for egquery count
    my $translation = $parser->get_count;
    
    my $corrected = $parser->get_corrected_query; # espell
    
    while (my $gquery = $parser->next_GlobalQuery) {
       # iterates through egquery data
    }
  
  # Bio::Tools::EUtilities::Info (einfo data)
  # database-related information
  
    my $desc = $parser->get_description;  
    my $update = $parser->get_last_update;
    my $nm = $parser->get_menu_name;
    my $ct = $parser->get_record_count;
  
    while (my $field = $parser->next_InfoField) {...}
    while (my $field = $parser->next_InfoLink) {...}
  
  # Bio::Tools::EUtilities::History (epost data)
  # data which enables one to retrieve and query against user-stored information on the NCBI server
  
    if ($parser->has_history) {
        # use parser as direct input in future queries
    }
    
    my ($webenv, $querykey) = $parser->history;
  
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

  # Bio::Tools::EUtilities::Link (elink data)
  # data retrieved using links between related information in databases
  
    # still working on new API
    
=head1 DESCRIPTION

Parses NCBI eutils XML output for retrieving IDs and other information. Part of
the BioPerl EUtilities system.

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

package Bio::Tools::EUtilities;
use strict;
use warnings;

use base qw(Bio::Root::Root);

use XML::Simple;
use Bio::Root::IO;
use Data::Dumper;

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
    'epost'     => 'History',
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
    my ($response, $type, $eutil, $cache) = $self->_rearrange([qw(RESPONSE DATATYPE EUTIL CACHE_RESPONSE)], @args);
    $cache ||= 0;
    $self->datatype($type);
    $self->eutil($eutil);
    $response  && $self->response($response);
    $self->cache_response($cache);
    my $io = Bio::Root::IO->new(@args);
    #$io->_initialize_io(@args);
    $self->io($io);
    $self->{'_parsed'} = 0;
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

=head2 lazy

 *** NOT IMPLEMENTED YET ***
 Title    : lazy
 Usage    : $parser->lazy(1)
 Function : sets flag to set lazy parsing iterator for elink/esummary data
 Returns  : value eval'ing to TRUE or FALSE
 Args     : value eval'ing to TRUE or FALSE
 Note     : must be set prior to any parsing run
 
=cut

sub lazy {
    my ($self, $flag) = @_;
    if (defined $flag) {
        $self->{'_lazy'} = ($flag) ? 1 : 0;
    }
    return $self->{'_lazy'};
}

=head2 io

 Title    : io
 Usage    : my $io = $parser->io;
 Function : Get/Set Bio::Root::IO object
 Returns  : Bio::Root::IO
 Args     : Bio::Root::IO

=cut

sub io {
    my ($self, $io) = @_;
    if ($io) {
        $self->throw('Not a Bio::Root::IO') if !$io->isa('Bio::Root::IO');
        $self->{'_io'} = $io;
    }
    return $self->{'_io'};
}

=head2 data_parsed

 Title    : data_parsed
 Usage    : if ($parser->data_parsed) {...}
 Function : returns TRUE if data has been parsed
 Returns  : value eval'ing to TRUE or FALSE
 Args     : none (set within parser)

=cut

sub data_parsed {
    return shift->{'_parsed'};
}

{
my %EUTIL_DATA = (
    'esummary'  => [qw(DocSum Item)],
    'epost'     => [],
    'egquery'   => [],
    'einfo'     => [qw(Field Link)],
    'elink'     => [qw(LinkSet LinkSetDb LinkSetDbHistory IdUrlSet 
                        Id IdLinkSet ObjUrl Link LinkInfo)],
    'espell'    => [qw(Original Replaced)],
    'esearch'   => [qw(Id)],
    );

sub parse_data {
    my $self = shift;
    my $eutil = $self->eutil;
    my $xs = XML::Simple->new();
    my $response = $self->response ? $self->response :
                   $self->io       ? $self->io->_fh  :
        $self->throw('No response or stream specified');
    my $simple = ($eutil eq 'espell') ?
            $xs->XMLin($self->_fix_espell($response), forcearray => $EUTIL_DATA{$eutil}) :
        ($response && $response->isa("HTTP::Response")) ?
            $xs->XMLin($response->content, forcearray => $EUTIL_DATA{$eutil}) :
            $xs->XMLin($response, forcearray => $EUTIL_DATA{$eutil});
    # check for errors
    if ($simple->{ERROR}) {
        my $error = $simple->{ERROR};
        $self->throw("NCBI $eutil nonrecoverable error: ".$error) unless ref $error;
    }
    if ($simple->{InvalidIdList}) {
        $self->warn("NCBI $eutil error: Invalid ID List".$simple->{InvalidIdList});
        return;
    }    
    if ($simple->{ErrorList} || $simple->{WarningList}) {
        my %errorlist = %{ $simple->{ErrorList} } if $simple->{ErrorList};
        my %warninglist = %{ $simple->{WarningList} } if $simple->{WarningList};
        my ($err_warn);
        for my $key (sort keys %errorlist) {
            my $messages = join("\n",grep {!ref $_} @{$errorlist{$key}});
            $err_warn .= "Error : $key = $messages";
        }    
        for my $key (sort keys %warninglist) {
            my $messages = join("\n",grep {!ref $_} @{$warninglist{$key}});
            $err_warn .= "Warning : $key = $messages";
        }
        chomp($err_warn);
        $self->warn("NCBI $eutil Errors/Warnings:\n".$err_warn)
    }
    delete $self->{'_response'} unless $self->cache_response;
    $self->_add_data($simple);    
}

}

=head2 Private methods

=cut

# Private data handlers for various eutils

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

1;


