#
# $Id$
#
# BioPerl module for Bio::DB::EMBL
#
# Cared for by Heikki Lehvaslaiho <Heikki@ebi.ac.uk>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::DB::EMBL - Database object interface to EMBL retrieval

=head1 SYNOPSIS

  use Bio::DB::EMBL;

  $embl = new Bio::DB::EMBL;

  # remember that EMBL_ID does not equal GenBank_ID!
  $seq = $embl->get_Seq_by_id('BUM'); # EMBL ID
 	print "cloneid is ", $seq->id, "\n";

  # or changeing to accession number and Fasta format ...
  $embl->request_format('fasta');
  $seq = $embl->get_Seq_by_acc('J02231'); # EMBL ACC
 	print "cloneid is ", $seq->id, "\n";

  # especially when using versions, you better be prepared
  # in not getting what what want
  eval {
      $seq = $embl->get_Seq_by_version('J02231.1'); # EMBL VERSION
  }
  print "cloneid is ", $seq->id, "\n" unless $@;

  # or ... best when downloading very large files, prevents
  # keeping all of the file in memory

  # also don't want features, just sequence so let's save bandwith
  # and request Fasta sequence
  $embl = new Bio::DB::EMBL(-retrievaltype => 'tempfile' ,
 			       -format => 'fasta');
   my $seqio = $embl->get_Stream_by_batch(['AC013798', 'AC021953'] );
  while( my $clone =  $seqio->next_seq ) {
 	print "cloneid is ", $clone->id, "\n";
  }

=head1 DESCRIPTION

Allows the dynamic retrieval of sequence objects L<Bio::Seq> from the
EMBL database using the emblfetch script at EBI:
L<http://www.ebi.ac.uk/cgi-bin/emblfetch>.

In order to make changes transparent we have host type (currently only
ebi) and location (defaults to ebi) separated out.  This allows later
additions of more servers in different geographical locations.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this
and other Bioperl modules. Send your comments and suggestions preferably
 to one of the Bioperl mailing lists.
Your participation is much appreciated.


  bioperl-l@bioperl.org                         - General discussion
  http://bio.perl.org/MailList.html             - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
 the bugs and their resolution.
 Bug reports can be submitted via email or the web:

  bioperl-bugs@bio.perl.org
  http://bio.perl.org/bioperl-bugs/

=head1 AUTHOR - Heikki Lehvaslaiho

Email Heikki Lehvaslaiho E<lt>Heikki@ebi.ac.ukE<gt>
=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::DB::EMBL;
use strict;
use vars qw(@ISA $MODVERSION %HOSTS $DEFAULTFORMAT $DEFAULTLOCATION
	    $DEFAULTSERVERTYPE %FORMATMAP);

$MODVERSION = '0.1';
use HTTP::Request::Common;
use Bio::DB::WebDBSeqI;

@ISA = qw(Bio::DB::WebDBSeqI);


BEGIN { 	
    # global vars
    $DEFAULTSERVERTYPE = 'ebi';
    $DEFAULTFORMAT = 'embl';
    $DEFAULTLOCATION = 'ebi';
    # you can add your own here theoretically.
    %HOSTS = (
	       'ebi' => {
		   baseurl => 'http://%s/cgi-bin/dbfetch?db=embl&style=raw&format=',
		   hosts   => {
		       'ebi'  => 'www.ebi.ac.uk'
		       }
	       });
    %FORMATMAP = ( 'embl' => 'embl',
		   'fasta'   => 'fasta'
		   );
}

# the new way to make modules a little more lightweight

sub new {
    my ($class, @args ) = @_;
    return $class->SUPER::new(@args);
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

    $self->throw("Must specify a value for uids to query")
	unless defined $uids;
    my $tmp;
    ($format, $tmp) = $self->request_format($format);

    my $url = $self->location_url();
    my $uid;
    if( ref($uids) =~ /ARRAY/i ) {
	$uid = join (',', @$uids);
	$self->warn ('The server will accept maximum of 50 entries in a request. The rest are ignored.')
	    if scalar @$uids >50;
    } else {
	$uid = $uids;
    }

    return GET $url. $format. '&id='. $uid;
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

# don't need to do anything

sub postprocess_data {
    my ($self, %args) = @_;
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

=head1 Bio::DB::EMBL specific routines

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
    return $self->get_seq_stream('-uids' => $ids, '-mode' => 'single');
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
			  keys %HOSTS) unless( $HOSTS{$servertype} );
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
    $location = lc $location;
    my $servertype = $self->servertype;
    $self->throw("Must have a valid servertype defined not $servertype")
	unless defined $servertype;
    my %hosts = %{$HOSTS{$servertype}->{'hosts'}};
    if( defined $location && $location ne '' ) {
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
    return sprintf($HOSTS{$servertype}->{'baseurl'},
		   $HOSTS{$servertype}->{'hosts'}->{$location});
}		

1;
__END__
