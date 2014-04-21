# $Id: HIV.pm 232 2008-12-11 14:51:51Z maj $
#
# BioPerl module for Bio::DB::HIV
#
# Please direct questions and support issues to <bioperl-l@bioperl.org>
#
# Cared for by Mark A. Jensen <maj@fortinbras.us>
#
# Copyright Mark A. Jensen
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::DB::HIV - Database object interface to the Los Alamos HIV Sequence Database

=head1 SYNOPSIS

    $db = new Bio::DB::HIV;

    $seq = $db->get_Seq_by_id('94284');                                 # LANL sequence id
    $seq = $db->get_Seq_by_acc('EF432710');                             # GenBank accession

    $q = new Bio::DB::Query::HIVQuery( " (C D)[subtype] SI[phenotype] (symptomatic AIDS)[patient_health] " );

    $seqio = $db->get_Stream_by_query($q);
    $seq = $seqio->next_seq();
    ($seq->annotation->get_Annotations('Virus'))[0]->{subtype}          # returns 'D'
    ($seq->annotation->get_Annotations('Patient'))[0]->{patient_health} # returns 'AIDS'
    ($seq->annotation->get_Annotations('accession'))[0]->{value}        # returns 'K03454'

=head1 DESCRIPTION

Bio::DB::HIV, along with L<Bio::DB::Query::HIVQuery>, provides an
interface for obtaining annotated HIV and SIV sequences from the Los
Alamos National Laboratory (LANL) HIV Sequence Database (
L<http://www.hiv.lanl.gov/content/sequence/HIV/mainpage.html>
). Unannotated sequences can be retrieved directly from the database
object, using either LANL ids or GenBank accessions. Annotations are
obtained via a query object, and are attached to the correct C<Bio::Seq>
objects when the query is handled by C<Bio::DB::HIV::get_Seq_by_query>
or C<Bio::DB::HIV::get_Stream_by_query>.

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
of the bugs and their resolution. Bug reports can be submitted via
the web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Mark A. Jensen

Email maj@fortinbras.us

=head1 CONTRIBUTORS

Mark A. Jensen

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

# Let the code begin...


package Bio::DB::HIV;
use strict;
use warnings;
use vars qw( $LANL_BASE $LANL_MAP_DB $LANL_MAKE_SEARCH_IF $LANL_SEARCH );

# Object preamble - inherits from Bio::DB::WebDBSeqI

use Bio::Root::Root;
use HTTP::Request::Common;
use Bio::DB::HIV::HIVAnnotProcessor;

use base qw(Bio::DB::WebDBSeqI);


BEGIN {
    # base change of 01/14/09
    $LANL_BASE = "http://www.hiv.lanl.gov/components/sequence/HIV/asearch";
    $LANL_MAP_DB = "map_db.comp";
    $LANL_MAKE_SEARCH_IF = "make_search_if.comp";
    $LANL_SEARCH = "search.comp";
    @Bio::ResponseProblem::Exception::ISA = qw( Bio::Root::Exception );
    @Bio::HIVSorry::Exception::ISA = qw ( Bio::Root::Exception );
    @Bio::WebError::Exception::ISA = qw( Bio::Root::Exception );
}

=head1 Constructor

=head2 new

 Title   : new
 Usage   : my $obj = new Bio::DB::HIV();
 Function: Builds a new Bio::DB::HIV object
 Returns : an instance of Bio::DB::HIV
 Args    :

=cut

sub new {
  my($class,@args) = @_;

  my $self = $class->SUPER::new(@args);
  my ($lanl_base, $lanl_map_db, $lanl_make_search_if, $lanl_search) =
      $self->_rearrange([qw(
                           LANL_BASE
                           LANL_MAP_DB
                           LANL_MAKE_SEARCH_IF
                           LANL_SEARCH
                           )], @args);

  $lanl_base                  && $self->lanl_base($lanl_base);
  $lanl_map_db                && $self->map_db($lanl_map_db);
  $lanl_make_search_if        && $self->make_search_if($lanl_make_search_if);
  $lanl_search                && $self->search_($lanl_search);
  # defaults
  $self->lanl_base            || $self->lanl_base($LANL_BASE);
  $self->map_db               || $self->map_db($LANL_MAP_DB);
  $self->make_search_if       || $self->make_search_if($LANL_MAKE_SEARCH_IF);
  $self->search_              || $self->search_($LANL_SEARCH);
  $self->url_base_address     || $self->url_base_address($self->lanl_base);

  $self->request_format("fasta");

  return $self;
}

=head1 WebDBSeqI compliance

=head2 get_request

 Title   : get_request
 Usage   : my $url = $self->get_request
 Function: returns a HTTP::Request object
 Returns :
 Args    : %qualifiers = a hash of qualifiers with keys in
            (-ids, -format, -mode, -query)
 Note    : Several layers of requests are performed to get to the sequence;
           see Bio::DB::Query::HIVQuery.

=cut

sub get_request {
    my $self = shift;
    my %quals = @_;
    my ($resp);
    my (@ids, $mode, @interface, @query_parms, $query);

    # html parsing regexps
    my $tags_re = qr{(?:\s*<[^>]+>\s*)};
    my $session_id_re = qr{<input.*name="id".*value="([0-9a-f]+)"}m;
    my $search_form_re = qr{<form[^>]*action=".*/search.comp"};
    my $seqs_found_re = qr{Displaying$tags_re*(?:\s*[0-9-]*\s*)*$tags_re*of$tags_re*\s*([0-9]+)$tags_re*sequences found};
    my $no_seqs_found_re = qr{Sorry.*no sequences found};
    my $too_many_re = qr{too many records: $tags_re*([0-9]+)};
    # find something like:
    #  <strong>tables without join:</strong><br>SequenceAccessions<br>
    my $tbl_no_join_re = qr{tables without join}i;
#    my $sorry_bud_re = qr{};

    # handle "qualifiers"
    foreach (keys %quals) {
	m/mode/ && do {
	    $mode = $quals{$_};
	    next;
	};
	m/uids/ && do {
	    $self->throw(-class=>"Bio::Root::BadParameter",
			 -text=>"Arrayref required for qualifier \"$_\"",
			 -value=>$quals{$_}) unless ref($quals{$_}) eq 'ARRAY';
	    @ids = @{$quals{$_}};
	    next;
	};
	m/query/ && do {
	    $self->throw(-class=>"Bio::Root::BadParameter",
			 -text=>"Bio::DB::Query::HIVQuery required for qualifier \"$_\"",
			 -value=>$quals{$_}) unless $quals{$_}->isa("Bio::DB::Query::HIVQuery");
	    $query = $quals{$_};
	    next;
	};
	do {
	    1; #else stub
	};
    }
	# what kind of request?
    for my $m ($mode) {
	($m =~ m/single/) && do {
	    @interface = (
		'sequenceentry' => 'se_sequence',
		'sequenceentry' => 'se_id',
		'action' => 'Search Interface'
		);
	    @query_parms = map { ('sequenceentry.se_id' => $_ ) } @ids;
	    push @query_parms, (
		'sequenceentry.se_sequence'=>'Any',
		'order' => 'sequenceentry.se_id',
		'sort_dir' => 'ASC',
		'action' => 'Search'
	    );
	};
	($mode =~ m/acc/) && do {
	    @interface = (
		'sequenceentry' => 'se_sequence',
		'sequenceentry' => 'se_id',
		'sequenceaccessions' => 'sa_genbankaccession',
		'sequenceaccessions' => 'sa_se_id',
		'action' => 'Search Interface'
		);
	    @query_parms = map {('sequenceaccessions.sa_genbankaccession' => $_)} @ids;
	    push @query_parms, (
		'sequenceentry.se_sequence' => 'Any',
		'order' => 'sequenceaccessions.sa_genbankaccession',
		'sort_dir' => 'ASC',
		'action' => 'Search'
	    );
	};
	($mode =~ m/gi/) && do {
	    $self->_sorry("-mode=>gi");
	};
	($mode =~ m/version/) && do {
	    $self->_sorry("-mode=>version");
	};
	($mode =~ m/query/) && do {
	    $self->throw(-class=>"Bio::Root::BadParameter",
			 -text=>"Query ".($query->{'_RUN_LEVEL'} ? "has been run only at run level ".$query->{'_RUN_LEVEL'} : "has not been run").", run at level 2 with _do_query(2)",
			 -value=>$query->{'_RUN_LEVEL'}) unless $query->{'_RUN_LEVEL'} == 2;
	    @interface = (
		'sequenceentry' => 'se_sequence',
		'sequenceentry' => 'se_id',
		'action' => 'Search Interface'
		);
	    @query_parms = ("sequenceentry.se_id" =>sprintf("'%s'",join("\t", $query->ids)));
#	    @query_parms = map { ( "sequenceentry.se_id" => $_ ) } $query->ids;
	    push @query_parms, (
		'sequenceentry.se_sequence' => 'Any',
		'order' => 'sequenceentry.se_id',
		'sort_dir' => 'ASC',
		'action' => 'Search'
	    );
	};
	do {
	    1; # else stub
	};
    }
    # web work
    eval { # capture web errors; throw below...
	# negotiate a session with lanl db
	if (!$self->_session_id) {
	    $resp = $self->ua->get($self->_map_db_uri);
	    $resp->is_success || die "Connect failed";
	    # get the session id
	    if (!$self->_session_id) {
		($self->{'_session_id'}) = ($resp->content =~ /$session_id_re/);
		$self->_session_id || die "Session not established";
	    }
	}

	# establish correct "interface" for this session id
	$resp = $self->ua->post($self->_make_search_if_uri, [@interface, id=>$self->_session_id]);
	$resp->is_success || die "Interface request failed (1)";
	$self->_response($resp);
	$resp->content =~ /$search_form_re/ || die "Interface request failed (2)";

	# interface successful, do the "pre-search"
	$resp = $self->ua()->post($self->_search_uri, [(@query_parms, 'id' => $self->_session_id)] );
	unless ($resp->is_success) {
	     die "Search post failed";
	}
	$self->_response($resp);
	# check for error conditions
	for ($resp->content) {
	    /$no_seqs_found_re/ && do {
		die "No sequences found";
		last;
	    };
	    /$too_many_re/ && do {
		die "Too many records ($1): must be <10000";
		last;
	    };
	    /$tbl_no_join_re/ && do {
		die "Some required tables went unjoined to query";
		last;
	    };
	    /$seqs_found_re/ && do {
		last;
	    };
	    do {
		die "Unparsed failure";
		last;
	    };
	}

    };
    $self->throw(-class=>'Bio::WebError::Exception',
		 -text=>$@,
		 -value=>$resp->content) if $@;

    # "pre-search" successful, return request
###  check this post update
    return POST $self->_search_uri,
    ['action Download.x' => 1,
     'action Download.y'=>1,
     'id'=>$self->_session_id
    ];


}

=head2 postprocess_data

 Title   : postprocess_data
 Usage   : $self->postprocess_data ( 'type' => 'string',
				     'location' => \$datastr);
 Function: process downloaded data before loading into a Bio::SeqIO
 Returns : void
 Args    : hash with two keys - 'type' can be 'string' or 'file'
                              - 'location' either file location or string
                                           reference containing data

=cut

sub postprocess_data {
    # parse tab-separated value content from LANL db
	my ( $self, %args) = @_;
	my ($type, $loc) = ($args{type}, $args{location});
	my (@data, @cols, %rec, $idkey, @flines);
	$self->throw(-class=>'Bio::Root::BadParameter',
		     -text=>"Argument hash requires values for keys \"type\" and \"location\"",
		     -value=>\%args) unless ($type && $loc);
	for ($type) {
	    m/string/ && do {
		@data = split(/\n|\r/, ${$loc});
		last;
	    };
	    m/file/ && do {
		local $/ = undef;
		open my $F, '<', $loc or
		    $self->throw(
		    -class => 'Bio::Root::FileOpenException',
		    -text  => "Error opening tempfile '$loc' for reading",
		    -value => $!
		    );
		@data = split( /\n|\r/, <$F>);
		close $F;
		last;
	    };
	    do {
		1; # else stub
	    };
	}
	$self->throw(-class=>'Bio::Root::BadParameter',
		     -text=>'No data found in repsonse',
		     -value=>%args) unless (@data);
	my $l;
	do {
	    $l = shift @data;
	} while  ( defined $l && $l !~ /Number/ ); # number-returned line
	@cols = split( /\t/, shift @data);

	# if Accession column is present, get_Stream_by_acc was called
	# otherwise, return lanl ids
	($idkey) = grep /SE.id/i, @cols unless ($idkey) = grep /Accession/i, @cols;
	$self->throw(-class=>"Bio::ResponseProblem::Exception",
		     -text=>"Trouble with column headers in LANL response",
		     -value=>join(' ',@cols)) unless $idkey;

 	foreach (@data) {
	    chop;
	    @rec{@cols} = split /\t/;
	    push @flines, ">$rec{$idkey}\n".$rec{'Sequence'}."\n";
	}
	for ($type) {
	    m/string/ && do {
		${$loc} = join("", @flines);
		last;
	    };
	    m/file/ && do {
		open my $F, '>', $loc or $self->throw(-class=>'Bio::Root::FileOpenException',
					     -text=>"Error opening tempfile '$loc' for writing",
					     -value=>$!);
		print $F join("", @flines);
		close $F;
		last;
	    };
	    do {
		1; #else stub
	    };
	}
	return;
}

=head1 WebDBSeqI overrides

=head2 get_seq_stream

 Title   : get_seq_stream
 Usage   : my $seqio = $self->get_seq_stream(%qualifiers)
 Function: builds a url and queries a web db
 Returns : a Bio::SeqIO stream capable of producing sequence
 Args    : %qualifiers = a hash qualifiers that the implementing class
           will process to make a url suitable for web querying
 Note    : Some tightening up of the baseclass version

=cut

sub get_seq_stream {
    my ($self, %qualifiers) = @_;
    my ($rformat, $ioformat) = $self->request_format();

    my ($key) = grep /format$/, keys %qualifiers;
    $qualifiers{'-format'} = ($key ? $qualifiers{$key} : $rformat);
    ($rformat, $ioformat) = $self->request_format($qualifiers{'format'});

# web work is here/maj
    my $request = $self->get_request(%qualifiers);

# authorization is here/maj
    $request->proxy_authorization_basic($self->authentication)
	if ( $self->authentication);
    $self->debug("request is ". $request->as_string(). "\n");

# workaround for MSWin systems (no forking available/maj)
    $self->retrieval_type('io_string') if $self->retrieval_type =~ /pipeline/ && $^O =~ /^MSWin/;

    if ($self->retrieval_type =~ /pipeline/) {
	# Try to create a stream using POSIX fork-and-pipe facility.
	# this is a *big* win when fetching thousands of sequences from
	# a web database because we can return the first entry while
	# transmission is still in progress.
	# Also, no need to keep sequence in memory or in a temporary file.
	# If this fails (Windows, MacOS 9), we fall back to non-pipelined access.

	# fork and pipe: _stream_request()=><STREAM>
	my ($result,$stream) = $self->_open_pipe();

	if (defined $result) {
	    $DB::fork_TTY = File::Spec->devnull; # prevents complaints from debugge
	    if (!$result) { # in child process
		$self->_stream_request($request,$stream);
		POSIX::_exit(0); #prevent END blocks from executing in this forked child
	    }
	    else {
		return Bio::SeqIO->new('-verbose' => $self->verbose,
				       '-format'  => $ioformat,
				       '-fh'      => $stream);
	    }
	}
	else {
	    $self->retrieval_type('io_string');
	}
    }

    if ($self->retrieval_type =~ /temp/i) {
	my $dir = $self->io->tempdir( CLEANUP => 1);
	my ( $fh, $tmpfile) = $self->io()->tempfile( DIR => $dir );
	close $fh;
	my $resp = $self->_request($request, $tmpfile);
	if( ! -e $tmpfile || -z $tmpfile || ! $resp->is_success() ) {
	    $self->throw("WebDBSeqI Error - check query sequences!\n");
	}
	$self->postprocess_data('type' => 'file','location' => $tmpfile);
	# this may get reset when requesting batch mode
	($rformat,$ioformat) = $self->request_format();
	if( $self->verbose > 0 ) {
	    open my $ERR, '<', $tmpfile or $self->throw("Could not read file '$tmpfile': $!");
	    while(<$ERR>) { $self->debug($_);}
	    close $ERR;
	}

	return Bio::SeqIO->new('-verbose' => $self->verbose,
			       '-format' => $ioformat,
			       '-file'   => $tmpfile);
    }

    if ($self->retrieval_type =~ /io_string/i ) {
	my $resp = $self->_request($request);
	my $content = $resp->content_ref;
	$self->debug( "content is $$content\n");
	if (!$resp->is_success() || length($$content) == 0) {
	    $self->throw("WebDBSeqI Error - check query sequences!\n");
	}
	($rformat,$ioformat) = $self->request_format();
	$self->postprocess_data('type'=> 'string',
				'location' => $content);
	$self->debug( "str is $$content\n");
	return Bio::SeqIO->new('-verbose' => $self->verbose,
			       '-format' => $ioformat,
			       '-fh'   => new IO::String($$content));
    }

    # if we got here, we don't know how to handle the retrieval type
    $self->throw("retrieval type " .
		 $self->retrieval_type .
		 " unsupported\n");
}

=head2 get_Stream_by_acc

  Title   : get_Stream_by_acc
  Usage   : $seq = $db->get_Stream_by_acc([$acc1, $acc2]);
  Function: Gets a series of Seq objects by GenBank accession numbers
  Returns : a Bio::SeqIO stream object
  Args    : an arrayref of accession numbers for
            the desired sequence entries
  Note    : For LANL DB, alternative to LANL seqids

=cut

sub get_Stream_by_acc {
    my ($self, $ids ) = @_;
    return $self->get_seq_stream('-uids' => [$ids], '-mode' => 'acc');
}

=head2 get_Stream_by_query

  Title   : get_Stream_by_query
  Usage   : $stream = $db->get_Stream_by_query($query);
  Function: Gets a series of Seq objects by way of a query string or oject
  Returns : a Bio::SeqIO stream object
  Args    : $query : Currently, only a Bio::DB::Query::HIVQuery object.
            It's a good idea to create the query object first and interrogate
            it for the entry count before you fetch a potentially large stream.

=cut

sub get_Stream_by_query {
    my ($self, $query ) = @_;
    my $stream = $self->get_seq_stream('-query' => $query, '-mode'=>'query');
    return new Bio::DB::HIV::HIVAnnotProcessor( -hiv_query=>$query, -source_stream=>$stream );
}

sub _request {
	my ($self, $request,$tmpfile) = @_;
	my ($resp);

	if( defined $tmpfile && $tmpfile ne '' ) {
		$resp =  $self->ua->request($request, $tmpfile);
	} else {
		$resp =  $self->ua->request($request);
	}

	if( $resp->is_error  ) {
		$self->throw("WebDBSeqI Request Error:\n".$resp->as_string);
	}
	return $resp;
}

=head1 Internals

=head2 lanl_base

 Title   : lanl_base
 Usage   : $obj->lanl_base($newval)
 Function: get/set the base url of the LANL HIV database
 Example :
 Returns : value of lanl_base (a scalar)
 Args    : on set, new value (a scalar or undef, optional)

=cut

sub lanl_base{
    my $self = shift;

    return $self->{'lanl_base'} = shift if @_;
    return $self->{'lanl_base'};
}

=head2 map_db

 Title   : map_db
 Usage   : $obj->map_db($newval)
 Function: get/set the cgi filename for map_db ("Database Map")
 Example :
 Returns : value of map_db (a scalar)
 Args    : on set, new value (a scalar or undef, optional)

=cut

sub map_db{
    my $self = shift;

    return $self->{'map_db'} = shift if @_;
    return $self->{'map_db'};
}

=head2 make_search_if

 Title   : make_search_if
 Usage   : $obj->make_search_if($newval)
 Function: get/set the cgi filename for make_search_if ("Make Search Interface")
 Example :
 Returns : value of make_search_if (a scalar)
 Args    : on set, new value (a scalar or undef, optional)

=cut

sub make_search_if{
    my $self = shift;

    return $self->{'make_search_if'} = shift if @_;
    return $self->{'make_search_if'};
}

=head2 search_

 Title   : search_
 Usage   : $obj->search_($newval)
 Function: get/set the cgi filename for the search query page
           ("Search Database")
 Example :
 Returns : value of search_ (a scalar)
 Args    : on set, new value (a scalar or undef, optional)

=cut

sub search_{
    my $self = shift;

    return $self->{'search_'} = shift if @_;
    return $self->{'search_'};
}

=head2 _map_db_uri

 Title   : _map_db_uri
 Usage   :
 Function: return the full map_db uri ("Database Map")
 Example :
 Returns : scalar string
 Args    : none

=cut

sub _map_db_uri{
    my $self = shift;
    return $self->url_base_address."/".$self->map_db;
}


=head2 _make_search_if_uri

 Title   : _make_search_if_uri
 Usage   :
 Function: return the full make_search_if uri ("Make Search Interface")
 Example :
 Returns : scalar string
 Args    : none

=cut

sub _make_search_if_uri{
    my $self = shift;
    return $self->url_base_address."/".$self->make_search_if;
}

=head2 _search_uri

 Title   : _search_uri
 Usage   :
 Function: return the full search cgi uri ("Search Database")
 Example :
 Returns : scalar string
 Args    : none

=cut

sub _search_uri{
    my $self = shift;
    return $self->url_base_address."/".$self->search_;
}

=head2 _session_id

 Title   : _session_id
 Usage   : $obj->_session_id($newval)
 Function: Contains HIV db session id (initialized in _do_lanl_request)
 Example :
 Returns : value of _session_id (a scalar)
 Args    : on set, new value (a scalar or undef, optional)

=cut

sub _session_id{
    my $self = shift;

    return $self->{'_session_id'} = shift if @_;
    return $self->{'_session_id'};
}

=head2 _response

 Title   : _response
 Usage   : $obj->_response($newval)
 Function: hold the response to search post
 Example :
 Returns : value of _response (a scalar)
 Args    : on set, new value (a scalar or undef, optional)

=cut

sub _response{
    my $self = shift;

    return $self->{'_response'} = shift if @_;
    return $self->{'_response'};
}

=head2 Dude, sorry

 Title   : _sorry
 Usage   : $hiv->_sorry
 Function: Throws an exception for unsupported option or parameter
 Example :
 Returns :
 Args    : scalar string

=cut

sub _sorry{
    my $self = shift;
    my $parm = shift;
    $self->throw(-class=>"Bio::HIVSorry::Exception",
		 -text=>"Sorry, option/parameter \"$parm\" not (yet) supported. See manpage to complain.",
		 -value=>$parm);
    return;
}


1;
