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
  http://bio.perl.org/bioperl-bugs/

=head1 AUTHOR - Jason Stajich

Email jason@chg.mc.duke.edu

=head1 APPENDIX

The rest of the documentation details each of the
object methods. Internal methods are usually
preceded with a _

=cut

# Let the code begin...

package Bio::DB::NCBIHelper;
use strict;
use vars qw(@ISA $HOSTBASE %CGILOCATION %FORMATMAP $DEFAULTFORMAT);

use Bio::DB::WebDBSeqI;
use HTTP::Request::Common;

@ISA = qw(Bio::DB::WebDBSeqI);

BEGIN { 	    
    $HOSTBASE = 'http://www.ncbi.nlm.nih.gov';
    %CGILOCATION = ( 'batch' => '/cgi-bin/Entrez/qserver.cgi',
		     'single'=> '/entrez/utils/qmap.cgi',
		     'version'=> '/htbin-post/Entrez/girevhist',
		     'gi'=> '/entrez/query.fcg');

    %FORMATMAP = ( 'genbank' => 'genbank',
		   'genpept' => 'genbank',
		   'fasta'   => 'fasta' );
    
    $DEFAULTFORMAT = 'genbank';
}

# the new way to make modules a little more lightweight

sub new {
    my ($class, @args ) = @_;
    return $class->SUPER::new(@args);
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
    my ($mode, $uids, $format) = $self->_rearrange([qw(MODE UIDS FORMAT)],
							 @qualifiers);
    
    $mode = lc $mode;
    ($format) = $self->request_format() if( !defined $format);
    if( !defined $mode || $mode eq '' ) { $mode = 'single'; }
    my %params = $self->get_params($mode);    
    if( ! %params ) {
	$self->throw("must specify a valid retrival mode 'single' or 'batch' not '$mode'") 
    }
    my $url = $HOSTBASE . $CGILOCATION{$mode};
#    print $url, "\t($mode)\n"; exit;
    if( !defined $uids ) {
	$self->throw("Must specify a value for uids to query");
    }

    if ($mode eq 'version') {
	$params{'val'} = $uids;
    } else {
	if( ref($uids) =~ /array/i ) {
	    $uids = join(",", @$uids);
	}
	$params{'uid'} = $uids;
    }

    if( $mode eq 'batch' ) {
	# has to be genbank at this point in time
	my $sformat = $format;
	if( $self->default_format !~ /$format/i ) {
	    $self->warn("must reset format to ". $self->default_format. 
			" for batch retrieval mode\n".
			"the only format supported by NCBI batch mode");
	    ($format) = $self->request_format($self->default_format);
	}
	$params{'format'} = $format;
	my $querystr = '?' . join("&", map { "$_=$params{$_}" } keys %params);
	print STDERR "url is ", $HOSTBASE . 
	    $CGILOCATION{$mode} . $querystr, "\n"
	    unless ( $self->verbose == 0 ); 	
	return POST ( $url, \%params );
    } elsif( $mode eq 'single' || $mode eq 'gi') {
	$params{'dopt'} = $format;
	my $querystr = '?' . join("&", map { "$_=$params{$_}" } keys %params);
	print STDERR "url is ", $url . $querystr, "\n"
	    unless ( $self->verbose == 0 );	
	return GET $url . $querystr;
    } elsif( $mode eq 'version') {
	my $querystr = '?' . join("&", map { "$_=$params{$_}" } keys %params);
	print STDERR "url is ", $url . $querystr, "\n"
	    unless ( $self->verbose == 0 );	
	return GET $url . $querystr;
    }  else { 
	return undef;
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

=cut

sub get_Stream_by_batch {
    my ($self, $ids) = @_;
    return $self->get_seq_stream('-uids' => $ids, '-mode'=>'batch');
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
		$data = ${$location}; 
    } elsif ( $type eq 'FILE' ) {
		open(TMP, $location) or $self->throw("could not open file $location");
		my @in = <TMP>;
		close TMP;
		$data = join("", @in);
    }
    my ($format) = $self->request_format();    
    if( $format =~ /fasta/i ) {
		$data =~ s/^[\s\S]+<dd>([\s\S]+)<pre>/$1/i;
    } else { 
       # remove everything before <PRE>
       $data =~ s/^[\s\S]+<pre>//i;
    }
    # remove everything after </PRE>
    $data =~ s/<\/pre>[\s\S]+$//i;

    # transform links to appropriate descriptions
    if ($data =~ /\nCONTIG\s+/) {
		$self->warn("CONTIG found. GenBank get_Stream_by_batch about to run."); 
    	my(@batch,@accession,%accessions,@location,$id,$contig,$stream,$aCount,$cCount,$gCount,$tCount);
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
			
		# grab multiple sequnces by batch and join based location variable
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
				if ($accessions{$accession[$i]}->{'count'} > 1) {
					# cache sequence for later use
					$accessions{$accession[$i]}->{'seq'} = $seq->seq();
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
    	$data =~ s/<a href=.+>(\S+)<\/a\>/$1/ig;
    }
    
    # fix gt and lt
    $data =~ s/&gt;/>/ig;
    $data =~ s/&lt;/</ig;
    if( $type eq 'FILE'  ) {
		open(TMP, ">$location") or $self->throw("could overwrite file $location");
		print TMP $data;
		close TMP;
    } elsif ( $type eq 'STRING' ) {
		${$args{'location'}} = $data;
    }
    $self->debug("format is ". $self->request_format(). "data is $data\n");
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


=head2 get_Seq_by_version

 Title   : get_Seq_by_version
 Usage   : $seq = $db->get_Seq_by_version('X77802.1');
 Function: Gets a Bio::Seq object by sequence version
 Returns : A Bio::Seq object
 Args    : accession.version (as a string)
 Throws  : "acc.version does not exist" exception

=cut 

sub get_Seq_by_version {
    my ($self,$seqid) = @_;  
    my ($acc, $version) =  $seqid =~ /(\w+).(\d+)/; 
    $self->throw("Use accesion.version notation, not[$seqid]") if( !defined $version );
    my $request = $self->get_Stream_by_version($acc);
    $self->throw("accession [$acc] does not exist") if( !defined $request );
    my $res = $self->ua->request($request);
    
    my $data  = $res->content;
    $data =~ s/<.*?>/ /gs;
    my($gi) = $data =~ /\s+(\d+)\s+$version\s+[A-Z][a-z]/;
    $self->throw("Version number [$version] does not exist for sequence [$acc]") unless $gi;
    return $self->get_Seq_by_gi($gi);
}

=head2 get_Stream_by_version

  Title   : get_Stream_by_version
  Usage   : 
  Function: DO NOT USE. HACK.
            Reuses the method defined by the interface file to retrieve
            a HTML table with all GIs (versions) for a accession number.
  Returns : a HTTP::Request object
  Args    : $ref : a reference to an array of accession.version strings for
                   the desired sequence entries

=cut

sub get_Stream_by_version {
    my ($self, $ids ) = @_;
    return $self->_get_version_request('-uids' => $ids, '-mode' => 'version');
}


sub _get_version_request{	# internal method to format a request 
                                # for a sequence version table
    my ($self, %qualifiers) = @_;
    my ($rformat, $ioformat) = $self->request_format();
    my $seen = 0;
    foreach my $key ( keys %qualifiers ) {
	if( $key =~ /format/i ) {
	    $rformat = $qualifiers{$key};
	    $seen = 1;
	}
    }
    $qualifiers{'-format'} = $rformat if( !$seen);
    ($rformat, $ioformat) = $self->request_format($rformat);
    
    my $request = $self->get_request(%qualifiers);
}

1;
__END__
