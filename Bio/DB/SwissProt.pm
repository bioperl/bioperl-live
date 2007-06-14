#
# $Id$
#
# BioPerl module for Bio::DB::SwissProt
#
# Cared for by Jason Stajich <jason@bioperl.org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code
# Reworked to use Bio::DB::WebDBSeqI 2000-12-11

=head1 NAME

Bio::DB::SwissProt - Database object interface to SwissProt retrieval

=head1 SYNOPSIS

    use Bio::DB::SwissProt;

    $sp = Bio::DB::SwissProt->new();

    $seq = $sp->get_Seq_by_id('KPY1_ECOLI'); # SwissProt ID
    # <4-letter-identifier>_<species 5-letter code>
    # or ...
    $seq = $sp->get_Seq_by_acc('P43780'); # SwissProt AC      
    # [OPQ]xxxxx


    # In fact in this implementation 
    # these methods call the same webscript so you can use 
    # then interchangeably

    # choose a different server to query
    $sp = Bio::DB::SwissProt->new('-servertype' => 'expasy',
				 '-hostlocation' => 'us');

    $seq = $sp->get_Seq_by_id('BOLA_HAEIN'); # SwissProtID

=head1 DESCRIPTION

SwissProt is a curated database of proteins managed by the Swiss
Bioinformatics Institute. Additional tools for
parsing and manipulating swissprot files can be found at
ftp://ftp.ebi.ac.uk/pub/software/swissprot/Swissknife/.

Allows the dynamic retrieval of Sequence objects (Bio::Seq) from the
SwissProt database via an Expasy retrieval.

In order to make changes transparent we have host type (currently only
expasy) and location (default to Switzerland) separated out.  This
allows the user to pick the closest Expasy mirror for running their
queries.


=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists.  Your participation is much appreciated.


  bioperl-l@bioperl.org                  - General discussion
  http://bioperl.org/wiki/Mailing_lists  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via the
web:

  http://bugzilla.open-bio.org/

=head1 AUTHOR - Jason Stajich

Email Jason Stajich  E<lt>jason@bioperl.org E<lt>

Thanks go to Alexandre Gattiker E<lt>gattiker@isb-sib.chE<gt> of Swiss
Institute of Bioinformatics for helping point us in the direction of
the correct expasy scripts and for swissknife references.

Also thanks to Heikki Lehvaslaiho E<lt>heikki-at-bioperl-dot-orgE<gt> 
for help with adding EBI swall server.

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::DB::SwissProt;
use strict;
use vars qw($MODVERSION %HOSTS $DEFAULTFORMAT $DEFAULTSERVERTYPE);

$MODVERSION = '0.8.1';
use HTTP::Request::Common;

use base qw(Bio::DB::WebDBSeqI);

# global vars
$DEFAULTSERVERTYPE = 'ebi';
$DEFAULTFORMAT = 'swissprot';

# you can add your own here theoretically.
%HOSTS = ( 
	   'expasy' => { 
	       'default' => 'us',
	       'baseurl' => 'http://%s/cgi-bin/sprot-retrieve-list.pl',
	       'hosts'   => 	       
	       { 
		   'switzerland'  => 'ch.expasy.org',
		   'canada' => 'ca.expasy.org',
		   'china'  => 'cn.expasy.org',
		   'taiwan' => 'tw.expasy.org',
		   'australia' => 'au.expasy.org',
		   'korea'  => 'kr.expasy.org',
		   'us'     => 'us.expasy.org',
	       },
	       # ick, CGI variables
	       'jointype' => ' ',
	       'idvar'    => 'list',
	       'basevars' => [ ],	       
	   },
	   'ebi'    => {
	       'default' => 'uk',
	       'baseurl' => 'http://%s/cgi-bin/dbfetch',
	       'hosts' => { 
		   'uk'   => 'www.ebi.ac.uk',
	       },
	       'jointype' => ',',
	       'idvar'    => 'id',
	       'basevars' => [ 'db'    => 'UniProtKB',
			       'style' => 'raw' ],
	   }
	   );

# new modules should be a little more lightweight and
# should use Bio::Root::Root
sub new {
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(@args);

    my ($format, $hostlocation,$servertype) = 
	$self->_rearrange([qw(FORMAT HOSTLOCATION SERVERTYPE)],
			  @args);    

    if( $format && $format !~ /(swiss)|(fasta)/i ) {
	$self->warn("Requested Format $format is ignored because only SwissProt and Fasta formats are currently supported");
	$format = $self->default_format;
    } 
    $servertype = $DEFAULTSERVERTYPE unless $servertype;
    $servertype = lc $servertype;
    $self->servertype($servertype);
    if (  $hostlocation ) {
	$self->hostlocation(lc $hostlocation);
    }

    $self->request_format($format); # let's always override the format, as it must be swiss or fasta
    return $self;
}

=head2 Routines from Bio::DB::RandomAccessI

=cut

=head2 get_Seq_by_id

 Title   : get_Seq_by_id
 Usage   : $seq = $db->get_Seq_by_id('ROA1_HUMAN')
 Function: Gets a Bio::Seq object by its name
 Returns : a Bio::Seq object
 Args    : the id (as a string) of a sequence
 Throws  : "id does not exist" exception

=cut

=head2 get_Seq_by_acc

 Title   : get_Seq_by_acc
 Usage   : $seq = $db->get_Seq_by_acc('X77802');
 Function: Gets a Bio::Seq object by accession number
 Returns : A Bio::Seq object
 Args    : accession number (as a string)
 Throws  : "acc does not exist" exception

=cut

=head2 get_Stream_by_id

  Title   : get_Stream_by_id
  Usage   : $stream = $db->get_Stream_by_id( [$uid1, $uid2] );
  Function: Gets a series of Seq objects by unique identifiers
  Returns : a Bio::SeqIO stream object
  Args    : $ref : a reference to an array of unique identifiers for
                   the desired sequence entries

=cut

=head2 get_Stream_by_acc

  Title   : get_Stream_by_acc
  Usage   : $seq = $db->get_Seq_by_acc([$acc1, $acc2]);
  Function: Gets a series of Seq objects by accession numbers
  Returns : a Bio::SeqIO stream object
  Args    : $ref : a reference to an array of accession numbers for
                   the desired sequence entries
  Note    : For GenBank, this just calls the same code for get_Stream_by_id()

=cut

=head2 get_Stream_by_batch

  Title   : get_Stream_by_batch
  Usage   : $seq = $db->get_Stream_by_batch($ref);
  Function: Retrieves Seq objects from SwissProt 'en masse', rather than one
            at a time.  This is implemented the same way as get_Stream_by_id, 
            but is provided here in keeping with access methods of NCBI 
            modules.
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

=head2 Implemented Routines from Bio::DB::WebDBSeqI interface

=cut

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

    if( !defined $uids ) {
	$self->throw("Must specify a value for uids to query");
    }
    my ($f,undef) = $self->request_format($format);
    
    my %vars = ( 
		 @{$HOSTS{$self->servertype}->{'basevars'}}, 
		 ( 'format' => $f )
		 );
    
    my $url = $self->location_url;
    
    my $uid;
    my $jointype = $HOSTS{$self->servertype}->{'jointype'} || ' ';
    my $idvar = $HOSTS{$self->servertype}->{'idvar'} || 'id';
    
    if( ref($uids) =~ /ARRAY/i ) {	
	# HTTP::Request automagically converts the ' ' to %20
	$uid = join($jointype, @$uids);	
    } else {
	$uid = $uids;
    }
    $vars{$idvar} = $uid;

    return POST $url, \%vars;
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
    return;
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

=head2 Bio::DB::SwissProt specific routines

=cut

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
	$self->{'_hostlocation'} = $HOSTS{$servertype}->{'default'};
	
	# make sure format is reset properly in that different
	# servers have different syntaxes
	my ($existingformat,$seqioformat) = $self->request_format;
	$self->request_format($existingformat);		
    }
    return $self->{'_servertype'} || $DEFAULTSERVERTYPE;
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

=head2 request_format

 Title   : request_format
 Usage   : my ($req_format, $ioformat) = $self->request_format;
           $self->request_format("genbank");
           $self->request_format("fasta");
 Function: Get/Set sequence format retrieval. The get-form will normally
           not be used outside of this and derived modules.
 Returns : Array of two strings, the first representing the format for
           retrieval, and the second specifying the corresponding SeqIO 
           format.
 Args    : $format = sequence format

=cut

sub request_format {
    my ($self, $value) = @_;
    if( defined $value ) {
	if( $self->servertype =~ /expasy/ ) {
	    if( $value =~ /sprot/ || $value =~ /swiss/ ) {
		$self->{'_format'} = [ 'sprot', 'swiss'];	    
	    } elsif( $value =~ /^fa/ ) {
		$self->{'_format'} = [ 'fasta', 'fasta'];
	    } else {
		$self->warn("Unrecognized format $value requested");
		$self->{'_format'} = [ 'fasta', 'fasta'];
	    }
	} elsif( $self->servertype =~ /ebi/ ) {
	    if( $value =~ /sprot/ || $value =~ /swiss/ ) {		
		$self->{'_format'} = [ 'swissprot', 'swiss' ];
	    } elsif( $value =~ /^fa/ ) {
		$self->{'_format'} = [ 'fasta', 'fasta'];
	    } else { 
		$self->warn("Unrecognized format $value requested");
		$self->{'_format'} = [ 'swissprot', 'swiss'];
	    }
	}
    }
    return @{$self->{'_format'}};
}

1;

__END__
