#
# BioPerl module for Bio::DB::DBFetch
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Heikki Lehvaslaiho <heikki-at-bioperl-dot-org>
#
# Copyright Heikki Lehvaslaiho
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::DB::DBFetch - Database object for retrieving using the dbfetch script

=head1 SYNOPSIS

  #do not use this module directly

=head1 DESCRIPTION

Allows the dynamic retrieval of entries from databases using the
dbfetch script at EBI:
L<http:E<sol>E<sol>www.ebi.ac.ukE<sol>cgi-binE<sol>dbfetch>.

In order to make changes transparent we have host type (currently only
ebi) and location (defaults to ebi) separated out.  This allows later
additions of more servers in different geographical locations.

This is a superclass which is called by instantiable subclasses with
correct parameters.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists.  Your participation is much appreciated.

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
the bugs and their resolution.  Bug reports can be submitted via the
web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Heikki Lehvaslaiho

Email Heikki Lehvaslaiho E<lt>heikki-at-bioperl-dot-orgE<gt>

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::DB::DBFetch;
use strict;
use vars qw($MODVERSION $DEFAULTFORMAT $DEFAULTLOCATION
	         $DEFAULTSERVERTYPE);

$MODVERSION = '0.1';
use HTTP::Request::Common;

use base qw(Bio::DB::WebDBSeqI);

# the new way to make modules a little more lightweight

BEGIN { 	
    # global vars
    $DEFAULTSERVERTYPE = 'dbfetch';
    $DEFAULTLOCATION = 'ebi';
}


=head1 Routines from Bio::DB::WebDBSeqI

=head2 get_request

 Title   : get_request
 Usage   : my $url = $self->get_request
 Function: returns a HTTP::Request object
 Returns :
 Args    : %qualifiers = a hash of qualifiers (ids, format, etc)

=cut

sub get_request {
	my ($self, @qualifiers) = @_;
	my ($uids, $format) = $self->_rearrange([qw(UIDS FORMAT)],
														 @qualifiers);

	$self->throw("Must specify a value for UIDs to fetch")
	  unless defined $uids;
	my $tmp;
	my $format_string = '';
	$format ||= $self->default_format;
	($format, $tmp) = $self->request_format($format);
	$format_string = "&format=$format"; 
	my $url = $self->location_url();
	my $uid;
	if( ref($uids) =~ /ARRAY/i ) {
		$uid = join (',', @$uids);
		$self->warn ('The server will accept maximum of 50 entries in a request. The rest are ignored.')
		  if scalar @$uids >50;
	} else {
		$uid = $uids;
	}

	return GET $url. $format_string. '&id='. $uid;
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

# remove occasional blank lines at top of web output
sub postprocess_data {
  my ($self, %args) = @_;
  if ($args{type} eq 'string') {
    ${$args{location}} =~ s/^\s+//;  # get rid of leading whitespace
  }
  elsif ($args{type} eq 'file') {
    my $F;
    open $F,"<", $args{location} or $self->throw("Cannot open $args{location}: $!");
    my @data = <$F>;
    for (@data) {
      last unless /^\s+$/;
      shift @data;
    }
    open $F,">", $args{location} or $self->throw("Cannot write to $args{location}: $!");
    print $F @data;
    close $F;
  }
}

=head2 default_format

 Title   : default_format
 Usage   : my $format = $self->default_format
 Function: Returns default sequence format for this module
 Returns : string
 Args    : none

=cut

sub default_format {
    my ($self) = @_;
    return $self->{'_default_format'};
}

=head1 Bio::DB::DBFetch specific routines

=head2 get_Stream_by_id

  Title   : get_Stream_by_id
  Usage   : $seq = $db->get_Stream_by_id($ref);
  Function: Retrieves Seq objects from the server 'en masse', rather than one
            at a time.  For large numbers of sequences, this is far superior
            than get_Stream_by_[id/acc]().
  Example :
  Returns : a Bio::SeqIO stream object
  Args    : $ref : either an array reference, a filename, or a filehandle
            from which to get the list of unique ids/accession numbers.

NOTE: for backward compatibility, this method is also called
get_Stream_by_batch.

=cut

sub get_Stream_by_id {
    my ($self, $ids) = @_;
    return $self->get_seq_stream('-uids' => $ids, '-mode' => 'batch');
}

=head2 get_Seq_by_version

 Title   : get_Seq_by_version
 Usage   : $seq = $db->get_Seq_by_version('X77802.1');
 Function: Gets a Bio::Seq object by accession number
 Returns : A Bio::Seq object
 Args    : version number (as a string)
 Throws  : "version does not exist" exception

=cut

sub get_Seq_by_version {
   my ($self,$seqid) = @_;
   my $seqio = $self->get_Stream_by_acc([$seqid]);
   $self->throw("version does not exist") if( !defined $seqio );
   return $seqio->next_seq();
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
	$self->{'_format'} = $value;
	return ($value, $value);
    }
    $value = $self->{'_format'};
    if( $value and defined $self->formatmap->{$value} ) {
	return ($value, $self->formatmap->{$value});
    } else {
	# Try to fall back to a default.
	return ($self->default_format, $self->default_format );
    }
}


=head2 servertype

 Title   : servertype
 Usage   : my $servertype = $self->servertype
	    $self->servertype($servertype);
 Function: Get/Set server type
 Returns : string
 Args    : server type string [optional]

=cut

sub servertype {
    my ($self, $servertype) = @_;
    if( defined $servertype && $servertype ne '') {		
	 $self->throw("You gave an invalid server type ($servertype)".
			  " - available types are ".
			  keys %{$self->hosts}) unless( $self->hosts->{$servertype} );
	 $self->{'_servertype'} = $servertype;
    }
    $self->{'_servertype'} = $DEFAULTSERVERTYPE unless $self->{'_servertype'};
    return $self->{'_servertype'};
}

=head2 hostlocation

 Title   : hostlocation
 Usage   : my $location = $self->hostlocation()
          $self->hostlocation($location)
 Function: Set/Get Hostlocation
 Returns : string representing hostlocation
 Args    : string specifying hostlocation [optional]

=cut

sub hostlocation {
    my ($self, $location ) = @_;
    my $servertype = $self->servertype;
    $self->throw("Must have a valid servertype defined not $servertype")
	unless defined $servertype; 
    my %hosts = %{$self->hosts->{$servertype}->{'hosts'}};
    if( defined $location && $location ne '' ) {
    $location = lc $location;
	if( ! $hosts{$location} ) {
	    $self->throw("Must specify a known host, not $location,".
			 " possible values (".
			 join(",", sort keys %hosts ). ")");
	}
	$self->{'_hostlocation'} = $location;
    }
    $self->{'_hostlocation'} = $DEFAULTLOCATION unless $self->{'_hostlocation'};
    return $self->{'_hostlocation'};
}

=head2 location_url

 Title   : location
 Usage   : my $url = $self->location_url()
 Function: Get host url
 Returns : string representing url
 Args    : none

=cut

sub location_url {
    my ($self) = @_;
    my $servertype = $self->servertype();
    my $location = $self->hostlocation();
    if( ! defined $location || !defined $servertype )  {	
	$self->throw("must have a valid hostlocation and servertype set before calling location_url");
    }
    return sprintf($self->hosts->{$servertype}->{'baseurl'},
		   $self->hosts->{$servertype}->{'hosts'}->{$location});
}		

=head1 Bio::DB::DBFetch routines

These methods allow subclasses to pass parameters.

=head2 hosts

 Title   : hosts
 Usage   : 
 Function: get/set for host hash 
 Returns : 
 Args    : optional hash

=cut

sub hosts {
    my ($self, $value) = @_;
    if (defined $value) {
	$self->{'_hosts'} = $value;
    }
    unless (exists $self->{'_hosts'}) {
	return ('');
    } else {
	return $self->{'_hosts'};
    }
}		

=head2 formatmap

 Title   : formatmap
 Usage   : 
 Function: get/set for format hash
 Returns : 
 Args    : optional hash

=cut

sub formatmap {
    my ($self, $value) = @_;
    if (defined $value) {
	$self->{'_formatmap'} = $value;
    }
    unless (exists $self->{'_formatmap'}) {
	return ('');
    } else {
	return $self->{'_formatmap'};
    }
}		


1;
__END__
