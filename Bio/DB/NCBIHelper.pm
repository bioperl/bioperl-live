# $Id$
#
# BioPerl module for Bio::DB::NCBIHelper
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

 Do not use this module directly.
 # get a Bio::DB::NCBIHelper object somehow
 my $seqio = $db->get_Stream_by_acc(['MUSIGHBA1']);
 foreach my $seq ( $seqio->next_seq ) {
  # process seq
 }

=head1 DESCRIPTION

Provides a single place to setup some common methods for querying NCBI
web databases.  This module just centralizes the methods for
constructing a URL for querying NCBI GenBank and NCBI GenPept and the
common HTML stripping done in L<postprocess_data>().

The NCBI query URLs used are http://www.ncbi.nlm.nih.gov as the base URL,
/cgi-bin/Entrez/qserver.cgi as the query interface for batch mode, and 
/entrez/utils/qmap.cgi for single-query mode.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the
evolution of this and other Bioperl modules. Send
your comments and suggestions preferably to one
of the Bioperl mailing lists. Your participation
is much appreciated.

  bioperl-l@bioperl.org              - General discussion
  http://bioperl.org/MailList.shtml  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to
help us keep track the bugs and their resolution.
Bug reports can be submitted via email or the
web:

  bioperl-bugs@bio.perl.org
  http://bugzilla.bioperl.org/

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
use vars qw(@ISA $HOSTBASE %CGILOCATION %FORMATMAP 
	    $DEFAULTFORMAT $MAX_ENTRIES $VERSION);

use Bio::DB::WebDBSeqI;
use Bio::DB::Query::GenBank;
use HTTP::Request::Common;
use URI;
use Bio::Root::IO;
use Bio::DB::RefSeq;
use Bio::Root::Root;

@ISA = qw(Bio::DB::WebDBSeqI Bio::Root::Root);
$VERSION = '0.8';

BEGIN {
    $MAX_ENTRIES = 19000;
    $HOSTBASE = 'http://www.ncbi.nih.gov';
    %CGILOCATION = (
		    'batch'  => ['post' => '/entrez/eutils/efetch.fcgi'],
		    'query'  => ['get'  => '/entrez/eutils/efetch.fcgi'],
		    'single' => ['get'  => '/entrez/eutils/efetch.fcgi'],
		    'version'=> ['get'  => '/entrez/eutils/efetch.fcgi'],
		    'gi'   =>   ['get'  => '/entrez/eutils/efetch.fcgi'],
		     );

    %FORMATMAP = ( 'gb' => 'genbank',
		   'gp' => 'genbank',
		   'fasta'   => 'fasta',
		   );

    $DEFAULTFORMAT = 'gb';
}

# the new way to make modules a little more lightweight

sub new {
    my ($class, @args ) = @_;
    my $self = $class->SUPER::new(@args);

    return $self;
}


=head2 get_params

 Title   : get_params
 Usage   : my %params = $self->get_params($mode)
 Function: Returns key,value pairs to be passed to NCBI database
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
 Function: Returns default sequence format for this module
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
    my ($self, @qualifiers) = @_;
    my ($mode, $uids, $format, $query) = $self->_rearrange([qw(MODE UIDS FORMAT QUERY)],
							   @qualifiers);

    $mode = lc $mode;
    ($format) = $self->request_format() if( !defined $format);
    if( !defined $mode || $mode eq '' ) { $mode = 'single'; }
    my %params = $self->get_params($mode);
    if( ! %params ) {
	$self->throw("must specify a valid retrieval mode 'single' or 'batch' not '$mode'") 
    }
    my $url = URI->new($HOSTBASE . $CGILOCATION{$mode}[1]);

    unless( defined $uids or defined $query) {
	$self->throw("Must specify a query or list of uids to fetch");
    }

    if ($uids) {
      if( ref($uids) =~ /array/i ) {
	$uids = join(",", @$uids);
      }
      $params{'id'}      = $uids;
    }

    elsif ($query && $query->can('cookie')) {
      @params{'WebEnv','query_key'} = $query->cookie;
      $params{'db'}                 = $query->db;
    }

    elsif ($query) {
      $params{'id'} = join ',',$query->ids;
    }

    $params{'rettype'} = $format;
    if ($CGILOCATION{$mode}[0] eq 'post') {
      return POST $url,[%params];
    } else {
      $url->query_form(%params);
      $self->debug("url is $url \n");
      return GET $url;
    }
}

=head2 get_Stream_by_batch

  Title   : get_Stream_by_batch
  Usage   : $seq = $db->get_Stream_by_batch($ref);
  Function: Retrieves Seq objects from Entrez 'en masse', rather than one
            at a time.  For large numbers of sequences, this is far superior
            than get_Stream_by_[id/acc]().
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
            than get_Stream_by_[id/acc]().
  Example :
  Returns : a Bio::SeqIO stream object
  Args    : $query :   An Entrez query string or a
            Bio::DB::Query::GenBank object.  It is suggested that you
            create a Bio::DB::Query::GenBank object and get the entry
            count before you fetch a potentially large stream.

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
				     'location' => \$datastr);
 Function: process downloaded data before loading into a Bio::SeqIO
 Returns : void
 Args    : hash with two keys - 'type' can be 'string' or 'file'
                              - 'location' either file location or string 
                                           reference containing data

=cut

# the default method, works for genbank/genpept, other classes should
# override it with their own method.

sub postprocess_data {
    my ($self, %args) = @_;
    my $data;
    my $type = uc $args{'type'};
    my $location = $args{'location'};
    if( !defined $type || $type eq '' || !defined $location) {
	return;
    } elsif( $type eq 'STRING' ) {
	$data = $$location; 
    } elsif ( $type eq 'FILE' ) {
	open(TMP, $location) or $self->throw("could not open file $location");
	my @in = <TMP>;
	close TMP;
	$data = join("", @in);
    }
    # transform links to appropriate descriptions
    if ($data =~ /\nCONTIG\s+/) {
	$self->warn("CONTIG found. GenBank get_Stream_by_batch about to run."); 
    	my(@batch,@accession,%accessions,@location,$id,
	   $contig,$stream,$aCount,$cCount,$gCount,$tCount);
    	my $gb = new Bio::DB::GenBank();
	
    	# process GenBank CONTIG join(...) into two arrays
    	$data =~ /(?:CONTIG\s+join\()((?:.+\n)+)(?:\/\/)/;
    	$contig = $1;
    	$contig =~ s/\n|\)//g;
	foreach (split /,/,$contig){
	    if (/>(.+)<.+>:(.+)/) {
		($id) = split /\./, $1;
		if (!$accessions{$id}) { push @batch, $id; }
		push @accession, $id;
		push @location, $2;
		$accessions{$id}->{'count'}++;
	    }
	}

	# grab multiple sequences by batch and join based location variable
	#$stream = $gb->get_Stream_by_batch(\@accession);

	$stream = $gb->get_Stream_by_batch(\@batch);
	$contig = "";
	
	for (my $i = 0; $i < @accession; $i++) {
	    my $seq;
	    if ($accessions{$accession[$i]}->{'seq'} ne '') {
				# retrieve stored sequence
				#my $seq =  $accessions{$accession[$i]}->{'seq'}   ;
		$seq = Bio::Seq::RichSeq->new(-seq => $accessions{$accession[$i]}->{'seq'});
	    } else {
				# seq not cached, get next sequence
		$seq = $stream->next_seq();
		if( defined $seq ) {
		    if ($accessions{$accession[$i]}->{'count'} > 1) {
			# cache sequence for later use		    
			$accessions{$accession[$i]}->{'seq'} = $seq->seq();
		    }
		} else { 
		    $self->warn("No Sequence available on stream");
		    return undef;
		}
	    }
	    my($start,$end) = split(/\.\./, $location[$i]);
	    $contig .= $seq->subseq($start,$end);
	}

	# count number of each letter in sequence
	$aCount = () = $contig =~ /a/ig;
	$cCount = () = $contig =~ /c/ig;
	$gCount = () = $contig =~ /g/ig;
	$tCount = () = $contig =~ /t/ig;

	# remove everything after and including CONTIG
	$data =~ s/(CONTIG[\s\S]+$)//i;

		    # build ORIGIN part of data file using sequence and counts
		    $data .= "BASE COUNT     $aCount a   $cCount c   $gCount g   $tCount t\n";
		    $data .= "ORIGIN      \n";
		    $data .= "$contig\n//";
		}
	else {
	    $data =~ s/<a\s+href\s*=.+>\s*(\S+)\s*<\s*\/a\s*\>/$1/ig;
	}

	# fix gt and lt
	$data =~ s/&gt;/>/ig;
	$data =~ s/&lt;/</ig;
	if( $type eq 'FILE'  ) {
	    open(TMP, ">$location") or $self->throw("couldn't overwrite file $location");
	    print TMP $data;
	    close TMP;
	} elsif ( $type eq 'STRING' ) {
	    ${$args{'location'}} = $data;
    }
    $self->debug("format is ". join(',',$self->request_format()). " data is $data\n");
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
    my ($self, $value) = @_;    
    if( defined $value ) {
	$value = lc $value;	
	if( defined $FORMATMAP{$value} ) {
	    $self->{'_format'} = [ $value, $FORMATMAP{$value}];
	} else {
	    # Try to fall back to a default. Alternatively, we could throw
	    # an exception
	    $self->{'_format'} = [ $value, $value ];
	}
    }
    return @{$self->{'_format'}};
}

=head2 Bio::DB::WebDBSeqI methods

Overriding WebDBSeqI method to help newbies to retrieve sequences

=head2 get_Stream_by_acc

  Title   : get_Stream_by_acc
  Usage   : $seq = $db->get_Seq_by_acc([$acc1, $acc2]);
  Function: Gets a series of Seq objects by accession numbers
  Returns : a Bio::SeqIO stream object
  Args    : $ref : a reference to an array of accession numbers for
                   the desired sequence entries
  Note    : For GenBank, this just calls the same code for get_Stream_by_id()

=cut

sub get_Stream_by_acc {
    my ($self, $ids ) = @_;
    my $newdb = $self->_check_id($ids);
    if (defined $newdb && ref($newdb) && $newdb->isa('Bio::DB::RefSeq')) {
	return $newdb->get_seq_stream('-uids' => $ids, '-mode' => 'single');
    } else {
	return $self->get_seq_stream('-uids' => $ids, '-mode' => 'single');
    }
}


=head2 _check_id

  Title   : _check_id
  Usage   : 
  Function: 
  Returns : A Bio::DB::RefSeq reference or throws
  Args    : $id(s), $string
=cut

sub _check_id {
    my ($self, $ids) = @_;

    # NT contigs can not be retrieved
    $self->throw("NT_ contigs are whole chromosome files which are not part of regular".
		 "database distributions. Go to ftp://ftp.ncbi.nih.gov/genomes/.") 
	if $ids =~ /NT_/;

    # Asking for a RefSeq from EMBL/GenBank

    if ($ids =~ /N._/) {
	$self->warn("[$ids] is not a normal sequence database but a RefSeq entry.".
		   " Redirecting the request.\n")
	    if $self->verbose >= 0;
	return  new Bio::DB::RefSeq;
    }
}

=head2 delay_policy

 Title   : delay_policy
 Usage   : $secs = $self->delay_policy
 Function: return number of seconds to delay between calls to remote db
 Returns : number of seconds to delay
 Args    : none

NOTE: NCBI requests a delay of 3s between requests.  This method
implements that policy.

=cut

sub delay_policy {
  my $self = shift;
  return 3;
}

1;
__END__
