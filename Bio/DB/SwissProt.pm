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

    $seq = $sp->get_Seq_by_id('P43780'); # SwissProtID
    # or ...
    $seq = $sp->get_Seq_by_acc('P43780'); # SwissProtID     
    # can only query on SwissProtID at expasy right now

    # choose a different server to query
    $sp = new Bio::DB::SwissProt('-hostlocation' => 'canada');

    $seq = $sp->get_Seq_by_id('P43780'); # SwissProtID
    

=head1 DESCRIPTION

Allows the dynamic retrieval of Sequence objects (Bio::Seq) from the SwissProt
database via an expasy retrieval.  Perhaps through SRS later.

In order to make changes transparent we have host type (currently only
expasy) and location (default to switzerland) separated out.  This
allows the user to pick the closest expasy mirror for running their
queries.



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

=head1 AUTHOR - Jason Stajich

Email Jason Stajich <jason@chg.mc.duke.edu>

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::DB::SwissProt;
use strict;
use vars qw(@ISA $MODVERSION %HOSTS $DEFAULTFORMAT $DEFAULTLOCATION 
	    $DEFAULTSERVERTYPE);

$MODVERSION = '0.8';
use HTTP::Request::Common;
use Bio::DB::WebDBSeqI;

@ISA = qw(Bio::DB::WebDBSeqI);

# global vars
$DEFAULTSERVERTYPE = 'expasy';
$DEFAULTFORMAT = 'swiss';
$DEFAULTLOCATION = 'switzerland';
# you can add your own here theoretically.
%HOSTS = ( 
	   'expasy' => { 
	       baseurl => 'http://%s/cgi-bin/get-sprot-raw.pl?',
	       hosts   => 
	       { 'switzerland'  => 'www.expasy.ch',
		 'canada' => 'expasy.cbr.nrc.ca',
		 'china'  => 'expasy.pku.edu.cn',
		 'taiwan' => 'expasy.nhri.org.tw',
		 'australia' => 'expasy.proteome.org.au', 
	     }
	   });


# new modules should be a little more lightweight and
# should use Bio::Root::RootI
sub new {
    my ($class, @args) = @_;
    my $self = bless {}, $class;
    my $make = $self->_initialize(@args);
    return $self;
}

sub _initialize {
    my ($self, @args) = @_;
    $self->SUPER::_initialize(@args);

    my ($format, $hostlocation,$servertype) = 
	$self->_rearrange([qw(fFORMAT HOSTLOCATION SERVERTYPE)],
			  @args);    

    if( $format && $format !~ /swiss/i ) {
	$self->warn("Requested Format $format is ignored because only SwissPort format is currently supported");
    } 
    $servertype = $DEFAULTSERVERTYPE unless $servertype;
    $hostlocation = $DEFAULTLOCATION unless( $hostlocation );    

    $self->request_format($self->default_format); # let's always override the format, as it must be swiss from this location

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

    my $url = $self->location_url();
    my $uid;
    if( ref($uids) =~ /ARRAY/i ) {
	$uid = $uids->[0];
	$self->warn("Currently can only process 1 swiss prot request at a time -- only processing $uid");
    } else {
	$uid = $uids;
    }
    return GET $url . $uid;
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
    my %hosts = %{$HOSTS{$servertype}->{hosts}};
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
    return sprintf($HOSTS{$servertype}->{baseurl}, 
		   $HOSTS{$servertype}->{hosts}->{$location});
}		   

1;
__END__
















