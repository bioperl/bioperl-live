#
# BioPerl module for Bio::DB::NCBIHelper
#
# Please direct questions and support issues to <bioperl-l@bioperl.org>
#
# Cared for by Jason Stajich
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code
#
# Interfaces with new WebDBSeqI interface

=head1 NAME

Bio::DB::NCBIHelper - A collection of routines useful for queries to
NCBI databases.

=head1 SYNOPSIS

 # Do not use this module directly.

 # get a Bio::DB::NCBIHelper object somehow
 my $seqio = $db->get_Stream_by_acc(['J00522']);
 foreach my $seq ( $seqio->next_seq ) {
     # process seq
 }

=head1 DESCRIPTION

Provides a single place to setup some common methods for querying NCBI
web databases.  This module just centralizes the methods for
constructing a URL for querying NCBI GenBank and NCBI GenPept and the
common HTML stripping done in L<postprocess_data>().

The base NCBI query URL used is:
http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the
evolution of this and other Bioperl modules. Send
your comments and suggestions preferably to one
of the Bioperl mailing lists. Your participation
is much appreciated.

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

Report bugs to the Bioperl bug tracking system to
help us keep track the bugs and their resolution.
Bug reports can be submitted via the web.

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Jason Stajich

Email jason@bioperl.org

=head1 APPENDIX

The rest of the documentation details each of the
object methods. Internal methods are usually
preceded with a _

=cut

# Let the code begin...

package Bio::DB::NCBIHelper;
use strict;

use Bio::DB::Query::GenBank;
use HTTP::Request::Common;
use URI;
use Bio::Root::IO;
use Bio::DB::RefSeq;
use URI::Escape qw(uri_unescape);

use base qw(Bio::DB::WebDBSeqI Bio::Root::Root);

our $HOSTBASE = 'http://eutils.ncbi.nlm.nih.gov';
our $MAX_ENTRIES = 19000;
our $REQUEST_DELAY = 3;
our %CGILOCATION = (
        'batch'   => [ 'post' => '/entrez/eutils/epost.fcgi' ],
        'query'   => [ 'get'  => '/entrez/eutils/efetch.fcgi' ],
        'single'  => [ 'get'  => '/entrez/eutils/efetch.fcgi' ],
        'version' => [ 'get'  => '/entrez/eutils/efetch.fcgi' ],
        'gi'      => [ 'get'  => '/entrez/eutils/efetch.fcgi' ],
        'webenv'  => [ 'get'  => '/entrez/eutils/efetch.fcgi' ]
    );
our %FORMATMAP = (
        'gb'          => 'genbank',
        'gp'          => 'genbank',
        'fasta'       => 'fasta',
        'asn.1'       => 'entrezgene',
        'gbwithparts' => 'genbank',
    );
our $DEFAULTFORMAT = 'gb';

=head2 new

 Title   : new
 Usage   :
 Function: the new way to make modules a little more lightweight
 Returns : 
 Args    : 

=cut

sub new {
    my ( $class, @args ) = @_;
    my $self = $class->SUPER::new(@args);
    my ($seq_start, $seq_stop,   $no_redirect,
        $redirect,  $complexity, $strand
        )
        = $self->_rearrange(
        [ qw(SEQ_START SEQ_STOP NO_REDIRECT REDIRECT_REFSEQ COMPLEXITY STRAND) ],
        @args
        );
    $seq_start   && $self->seq_start($seq_start);
    $seq_stop    && $self->seq_stop($seq_stop);
    $no_redirect && $self->no_redirect($no_redirect);
    $redirect    && $self->redirect_refseq($redirect);
    $strand      && $self->strand($strand);

    # adjust statement to accept zero value
    defined $complexity
        && ( $complexity >= 0 && $complexity <= 4 )
        && $self->complexity($complexity);
    return $self;
}


=head2 get_params

 Title   : get_params
 Usage   : my %params = $self->get_params($mode)
 Function: returns key,value pairs to be passed to NCBI database
           for either 'batch' or 'single' sequence retrieval method
 Returns : a key,value pair hash
 Args    : 'single' or 'batch' mode for retrieval

=cut

sub get_params {
    my ($self, $mode) = @_;
    $self->throw("subclass did not implement get_params");
}

=head2 default_format

 Title   : default_format
 Usage   : my $format = $self->default_format
 Function: returns default sequence format for this module
 Returns : string
 Args    : none

=cut

sub default_format {
    return $DEFAULTFORMAT;
}

=head2 get_request

 Title   : get_request
 Usage   : my $url = $self->get_request
 Function: HTTP::Request
 Returns :
 Args    : %qualifiers = a hash of qualifiers (ids, format, etc)

=cut

sub get_request {
    my ( $self, @qualifiers ) = @_;
    my ( $mode, $uids, $format, $query, $seq_start, $seq_stop, $strand,
        $complexity )
        = $self->_rearrange(
        [qw(MODE UIDS FORMAT QUERY SEQ_START SEQ_STOP STRAND COMPLEXITY)],
        @qualifiers );
    $mode = lc $mode;
    ($format) = $self->request_format() unless ( defined $format );
    if ( !defined $mode || $mode eq '' ) { $mode = 'single'; }
    my %params = $self->get_params($mode);
    if ( !%params ) {
        $self->throw(
            "must specify a valid retrieval mode 'single' or 'batch' not '$mode'"
        );
    }
    my $url = URI->new( $HOSTBASE . $CGILOCATION{$mode}[1] );
    unless ( $mode eq 'webenv' || defined $uids || defined $query ) {
        $self->throw("Must specify a query or list of uids to fetch");
    }
    if ( $query && $query->can('cookie') ) {
        @params{ 'WebEnv', 'query_key' } = $query->cookie;
        $params{'db'} = $query->db;
    }
    elsif ($query) {
        $params{'id'} = join ',', $query->ids;
    }

    # for batch retrieval, non-query style
    elsif ( $mode eq 'webenv' && $self->can('cookie') ) {
        @params{ 'WebEnv', 'query_key' } = $self->cookie;
    }
    elsif ($uids) {
        if ( ref($uids) =~ /array/i ) {
            $uids = join( ",", @$uids );
        }
        $params{'id'} = $uids;
    }
    $seq_start && ( $params{'seq_start'} = $seq_start );
    $seq_stop  && ( $params{'seq_stop'}  = $seq_stop );
    $strand    && ( $params{'strand'}    = $strand );
    if ( defined $complexity && ( $seq_start || $seq_stop || $strand ) ) {
        $self->warn(
            "Complexity set to $complexity; seq_start and seq_stop may not work!"
        ) if ( $complexity != 1 && ( $seq_start || $seq_stop ) );
        $self->warn(
            "Complexity set to 0; expect strange results with strand set to 2"
        ) if ( $complexity == 0 && $strand == 2 && $format eq 'fasta' );
    }
    defined $complexity && ( $params{'complexity'} = $complexity );
    $params{'rettype'} = $format unless $mode eq 'batch';

    # for now, 'post' is batch retrieval
    if ( $CGILOCATION{$mode}[0] eq 'post' ) {
        my $response = $self->ua->request( POST $url, [%params] );
        $response->proxy_authorization_basic( $self->authentication )
            if ( $self->authentication );
        $self->_parse_response( $response->content );
        my ( $cookie, $querykey ) = $self->cookie;
        my %qualifiers = (
            '-mode'       => 'webenv',
            '-seq_start'  => $seq_start,
            '-seq_stop'   => $seq_stop,
            '-strand'     => $strand,
            '-complexity' => $complexity,
            '-format'     => $format
        );
        return $self->get_request(%qualifiers);
    }
    else {
        $url->query_form(%params);
        return GET $url;
    }
}


=head2 get_Stream_by_batch

  Title   : get_Stream_by_batch
  Usage   : $seq = $db->get_Stream_by_batch($ref);
  Function: Retrieves Seq objects from Entrez 'en masse', rather than one
            at a time.  For large numbers of sequences, this is far superior
            than get_Stream_by_id or get_Stream_by_acc.
  Example :
  Returns : a Bio::SeqIO stream object
  Args    : $ref : either an array reference, a filename, or a filehandle
            from which to get the list of unique ids/accession numbers.

            NOTE: deprecated API.  Use get_Stream_by_id() instead.

=cut

*get_Stream_by_batch = sub {
   my $self = shift;
   $self->deprecated('get_Stream_by_batch() is deprecated; use get_Stream_by_id() instead');
   $self->get_Stream_by_id(@_)
};

=head2 get_Stream_by_query

  Title   : get_Stream_by_query
  Usage   : $seq = $db->get_Stream_by_query($query);
  Function: Retrieves Seq objects from Entrez 'en masse', rather than one
            at a time.  For large numbers of sequences, this is far superior
            to get_Stream_by_id and get_Stream_by_acc.
  Example :
  Returns : a Bio::SeqIO stream object
  Args    : An Entrez query string or a Bio::DB::Query::GenBank object.
            It is suggested that you create a Bio::DB::Query::GenBank object and get 
            the entry count before you fetch a potentially large stream.

=cut

sub get_Stream_by_query {
    my ($self, $query) = @_;
    unless (ref $query && $query->can('query')) {
       $query = Bio::DB::Query::GenBank->new($query);
    }
    return $self->get_seq_stream('-query' => $query, '-mode'=>'query');
}

=head2 postprocess_data

 Title   : postprocess_data
 Usage   : $self->postprocess_data ( 'type' => 'string',
				                             'location' => \$datastr );
 Function: Process downloaded data before loading into a Bio::SeqIO. This
           works for Genbank and Genpept, other classes should override
           it with their own method.
 Returns : void
 Args    : hash with two keys:

           'type' can be 'string' or 'file'
           'location' either file location or string reference containing data

=cut

sub postprocess_data {
	# retain this in case postprocessing is needed at a future date
}


=head2 request_format

 Title   : request_format
 Usage   : my ($req_format, $ioformat) = $self->request_format;
           $self->request_format("genbank");
           $self->request_format("fasta");
 Function: Get/Set sequence format retrieval. The get-form will normally not
           be used outside of this and derived modules.
 Returns : Array of two strings, the first representing the format for
           retrieval, and the second specifying the corresponding SeqIO format.
 Args    : $format = sequence format

=cut

sub request_format {
    my ( $self, $value ) = @_;
    if ( defined $value ) {
        $value = lc $value;
        if ( defined $FORMATMAP{$value} ) {
            $self->{'_format'} = [ $value, $FORMATMAP{$value} ];
        }
        else {
            # Try to fall back to a default. Alternatively, we could throw
            # an exception
            $self->{'_format'} = [ $value, $value ];
        }
    }
    return @{ $self->{'_format'} };
}


=head2 redirect_refseq

 Title   : redirect_refseq
 Usage   : $db->redirect_refseq(1)
 Function: simple getter/setter which redirects RefSeqs to use Bio::DB::RefSeq
 Returns : Boolean value
 Args    : Boolean value (optional)
 Throws  : 'unparseable output exception'
 Note    : This replaces 'no_redirect' as a more straightforward flag to
           redirect possible RefSeqs to use Bio::DB::RefSeq (EBI interface)
           instead of retrieving the NCBI records

=cut

sub redirect_refseq {
    my $self = shift;
    return $self->{'_redirect_refseq'} = shift if @_;
    return $self->{'_redirect_refseq'};
}

=head2 complexity

 Title   : complexity
 Usage   : $db->complexity(3)
 Function: get/set complexity value
 Returns : value from 0-4 indicating level of complexity
 Args    : value from 0-4 (optional); if unset server assumes 1
 Throws  : if arg is not an integer or falls outside of noted range above
 Note    : From efetch docs, the complexity regulates the display:

           0 - get the whole blob
           1 - get the bioseq for gi of interest (default in Entrez)
           2 - get the minimal bioseq-set containing the gi of interest
           3 - get the minimal nuc-prot containing the gi of interest
           4 - get the minimal pub-set containing the gi of interest

=cut

sub complexity {
    my ( $self, $comp ) = @_;
    if ( defined $comp ) {
        $self->throw("Complexity value must be integer between 0 and 4")
            if $comp !~ /^\d+$/ || $comp < 0 || $comp > 4;
        $self->{'_complexity'} = $comp;
    }
    return $self->{'_complexity'};
}

=head2 strand

 Title   : strand
 Usage   : $db->strand(1)
 Function: get/set strand value
 Returns : strand value if set
 Args    : value of 1 (plus) or 2 (minus); if unset server assumes 1
 Throws  : if arg is not an integer or is not 1 or 2
 Note    : This differs from BioPerl's use of strand: 1 = plus, -1 = minus 0 = not relevant.
           We should probably add in some functionality to convert over in the future.

=cut

sub strand {
    my ($self, $str) = @_;
    if ($str) {
        $self->throw("strand() must be integer value of 1 (plus strand) or 2 (minus strand) if set") if
            $str !~ /^\d+$/ || $str < 1 || $str > 2;
        $self->{'_strand'} = $str;
    }
    return $self->{'_strand'};
}

=head2 seq_start

 Title   : seq_start
 Usage   : $db->seq_start(123)
 Function: get/set sequence start location
 Returns : sequence start value if set
 Args    : integer; if unset server assumes 1
 Throws  : if arg is not an integer

=cut

sub seq_start {
    my ($self, $start) = @_;
    if ($start) {
        $self->throw("seq_start() must be integer value if set") if
            $start !~ /^\d+$/;
        $self->{'_seq_start'} = $start;
    }
    return $self->{'_seq_start'};
}

=head2 seq_stop

 Title   : seq_stop
 Usage   : $db->seq_stop(456)
 Function: get/set sequence stop (end) location
 Returns : sequence stop (end) value if set
 Args    : integer; if unset server assumes 1
 Throws  : if arg is not an integer

=cut

sub seq_stop {
    my ($self, $stop) = @_;
    if ($stop) {
        $self->throw("seq_stop() must be integer if set") if
            $stop !~ /^\d+$/;
        $self->{'_seq_stop'} = $stop;
    }
    return $self->{'_seq_stop'};
}

=head2 Bio::DB::WebDBSeqI methods

Overriding WebDBSeqI method to help newbies to retrieve sequences

=head2 get_Stream_by_acc

  Title   : get_Stream_by_acc
  Usage   : $seq = $db->get_Stream_by_acc([$acc1, $acc2]);
  Function: gets a series of Seq objects by accession numbers
  Returns : a Bio::SeqIO stream object
  Args    : $ref : a reference to an array of accession numbers for
                   the desired sequence entries
  Note    : For GenBank, this just calls the same code for get_Stream_by_id()

=cut

sub get_Stream_by_acc {
    my ( $self, $ids ) = @_;
    my $newdb = $self->_check_id($ids);
    if ( defined $newdb && ref($newdb) && $newdb->isa('Bio::DB::RefSeq') ) {
        return $newdb->get_seq_stream( '-uids' => $ids, '-mode' => 'single' );
    }
    else {
        return $self->get_seq_stream( '-uids' => $ids, '-mode' => 'single' );
    }
}

=head2 _check_id

  Title   : _check_id
  Usage   :
  Function:
  Returns : a Bio::DB::RefSeq reference or throws
  Args    : $id(s), $string

=cut

sub _check_id {
    my ( $self, $ids ) = @_;

    # NT contigs can not be retrieved
    $self->throw("NT_ contigs are whole chromosome files which are not part of regular"
            . "database distributions. Go to ftp://ftp.ncbi.nih.gov/genomes/.") 
      if $ids =~ /NT_/;

    # Asking for a RefSeq from EMBL/GenBank
    if ( $self->redirect_refseq ) {
        if ( $ids =~ /N._/ ) {
            $self->warn(
                "[$ids] is not a normal sequence database but a RefSeq entry."
                    . " Redirecting the request.\n" )
                if $self->verbose >= 0;
            return Bio::DB::RefSeq->new();
        }
    }
}


=head2 delay_policy

  Title   : delay_policy
  Usage   : $secs = $self->delay_policy
  Function: NCBI requests a delay of 3 seconds between requests. This method
            implements that policy.
  Returns : number of seconds to delay
  Args    : none

=cut

sub delay_policy {
    my $self = shift;
    return $REQUEST_DELAY;
}

=head2 cookie

 Title   : cookie
 Usage   : ($cookie,$querynum) = $db->cookie
 Function: return the NCBI query cookie, this information is used by 
           Bio::DB::GenBank in conjunction with efetch, ripped from 
           Bio::DB::Query::GenBank
 Returns : list of (cookie,querynum)
 Args    : none

=cut

sub cookie {
    my $self = shift;
    if (@_) {
        $self->{'_cookie'}   = shift;
        $self->{'_querynum'} = shift;
    }
    else {
        return @{$self}{qw(_cookie _querynum)};
    }
}

=head2 _parse_response

 Title   : _parse_response
 Usage   : $db->_parse_response($content)
 Function: parse out response for cookie, this is a trimmed-down version 
           of _parse_response from Bio::DB::Query::GenBank
 Returns : empty
 Args    : none
 Throws  : 'unparseable output exception'

=cut

sub _parse_response {
    my $self    = shift;
    my $content = shift;
    if ( my ($warning) = $content =~ m!<ErrorList>(.+)</ErrorList>!s ) {
        $self->warn("Warning(s) from GenBank: $warning\n");
    }
    if ( my ($error) = $content =~ /<OutputMessage>([^<]+)/ ) {
        $self->throw("Error from Genbank: $error");
    }
    my ($cookie)   = $content =~ m!<WebEnv>(\S+)</WebEnv>!;
    my ($querykey) = $content =~ m!<QueryKey>(\d+)!;
    $self->cookie( uri_unescape($cookie), $querykey );
}

=head2 no_redirect

 Title   : no_redirect
 Usage   : $db->no_redirect($content)
 Function: DEPRECATED - Used to indicate that Bio::DB::GenBank instance retrieves
           possible RefSeqs from EBI instead; default behavior is now to
           retrieve directly from NCBI
 Returns : None
 Args    : None
 Throws  : Method is deprecated in favor of positive flag method 'redirect_refseq'

=cut

sub no_redirect {
    shift->throw(
    "Use of no_redirect() is deprecated.  Bio::DB::GenBank default is to always\n".
    "retrieve from NCBI.  In order to redirect possible RefSeqs to EBI, set\n".
    "redirect_refseq flag to 1");
}

1;

__END__
