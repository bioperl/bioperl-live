#
# $Id$
#
# BioPerl module for Bio::DB::SwissProt
#
# Cared for by Jason Stajich <jason@chg.mc.duke.edu>
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

    $sp = new Bio::DB::SwissProt;

    $seq = $sp->get_Seq_by_id('KPY1_ECOLI'); # SwissProt ID
    # <4-letter-identifier>_<species 5-letter code>
    # or ...
    $seq = $sp->get_Seq_by_acc('P43780'); # SwissProt AC      
    # [OPQ]xxxxx


    # In fact in this implementation 
    # these methods call the same webscript so you can use 
    # then interchangeably

    # choose a different server to query
    $sp = new Bio::DB::SwissProt('-hostlocation' => 'canada');

    $seq = $sp->get_Seq_by_id('BOLA_HAEIN'); # SwissProtID

=head1 DESCRIPTION

SwissProt is a curated database of proteins managed by the Swiss
Bioinformatics Institute.  This is in contrast to EMBL/GenBank/DDBJ which are archives of protein information.  Additional tools for parsing and manipulating swissprot files can be found at ftp://ftp.ebi.ac.uk/pub/software/swissprot/Swissknife/.

Allows the dynamic retrieval of Sequence objects (Bio::Seq) from the SwissProt
database via an expasy retrieval.  Perhaps through SRS later.

In order to make changes transparent we have host type (currently only
expasy) and location (default to switzerland) separated out.  This
allows the user to pick the closest expasy mirror for running their
queries.


=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists.  Your participation is much appreciated.


  bioperl-l@bioperl.org                         - General discussion
  http://bio.perl.org/MailList.html             - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via email
or the web:

  bioperl-bugs@bio.perl.org
  http://bio.perl.org/bioperl-bugs/

=head1 AUTHOR - Jason Stajich

Email Jason Stajich  E<lt>jason@chg.mc.duke.edu E<lt>

Thanks go to Alexandre Gattiker E<lt>gattiker@isb-sib.chE<gt> of Swiss
Institute of Bioinformatics for helping point us in the direction of
the correct expasy scripts and for swissknife references.

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::DB::SwissProt;
use strict;
use vars qw(@ISA $MODVERSION %HOSTS $DEFAULTFORMAT $DEFAULTLOCATION 
	    $DEFAULTSERVERTYPE);

$MODVERSION = '0.7.1';
use HTTP::Request::Common;
use Bio::DB::WebDBSeqI;

@ISA = qw(Bio::DB::WebDBSeqI);

# global vars
$DEFAULTSERVERTYPE = 'expasy';
$DEFAULTFORMAT = 'sprot';
$DEFAULTLOCATION = 'switzerland';
# you can add your own here theoretically.
%HOSTS = ( 
	   'expasy' => { 
	       baseurl => 'http://%s/cgi-bin/sprot-retrieve-list.pl',
	       hosts   => 
	       { 'switzerland'  => 'ch.expasy.org',
		 'canada' => 'ca.expasy.org',
		 'china'  => 'cn.expasy.org',
		 'taiwan' => 'tw.expasy.org',
		 'australia' => 'au.expasy.org',
		 'korea'  => 'kr.expasy.org'
	     }
	   });


# new modules should be a little more lightweight and
# should use Bio::Root::RootI
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
    $hostlocation = $DEFAULTLOCATION unless( $hostlocation );    

    $self->request_format($format); # let's always override the format, as it must be swiss from this location

    $hostlocation = lc $hostlocation;
    $servertype = lc $servertype;
    $self->servertype($servertype);
    $self->hostlocation($hostlocation);
    return $self;
}

=head2 Routines fro Bio::DB::WebDBSeqI

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
    $self->request_format($format) if( defined $format );
    
    my %vars = ( 'format' => $format );
    my $url = $self->location_url;
    my $uid;
    if( ref($uids) =~ /ARRAY/i ) {	
	# HTTP::Request automagically converts the ' ' to %20
	$uid = join(' ', @$uids);
    } else {
	$uid = $uids;
    }
    $vars{'list'} = $uid;
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
 Function: Get/Set sequence format retrieval. The get-form will normally not
           be used outside of this and derived modules.
 Returns : Array of two strings, the first representing the format for
           retrieval, and the second specifying the corresponding SeqIO format.
 Args    : $format = sequence format

=cut

sub request_format {
    my ($self, $value) = @_;
    if( defined $value ) {
	if( $value =~ /sprot/ || $value =~ /swiss/ ) {
	    $self->{'_format'} = [ 'sprot', 'swiss'];	    
	} else {
	    $self->{'_format'} = [ $value, $value];
	}
    }
    return @{$self->{'_format'}};
}

1;
__END__
















