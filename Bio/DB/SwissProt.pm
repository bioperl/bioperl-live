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

=head1 NAME

Bio::DB::SwissProt - Database object interface to SwissProt retrieval

=head1 SYNOPSIS

    $sp = new Bio::DB::SwissProt;

    $seq = $sp->get_Seq_by_id('P43780'); # Unique ID

    # or ...

    $seq = $sp->get_Seq_by_acc('P43780'); # Accession Number

=head1 DESCRIPTION

Allows the dynamic retrieval of Sequence objects (Bio::Seq) from the SwissProt
database via an expasy retrieval.  Perhaps through SRS later.

In order to make changes transparent we have host type (currently only expasy) 
and location (default to switzerland) separated out.

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
use HTTP::Request;
use LWP::UserAgent;
use vars qw(@ISA $MODVERSION %HOSTS $DEFAULTFORMAT $DEFAULTLOCATION 
	    $DEFAULTSTYPE);

$MODVERSION = '0.8';
# Object preamble - inherits from Bio::DB::BioSeqI

use Bio::DB::RandomAccessI;
use Bio::Root::RootI;
use Bio::SeqIO;
use IO::File;
use IO::String;
use LWP::UserAgent;
use HTTP::Request::Common;
use HTTP::Response;

@ISA = qw(Bio::Root::RootI Bio::DB::RandomAccessI);

# global vars
$DEFAULTSTYPE = 'expasy';
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
    my $self = bless {
	_stype         => $DEFAULTSTYPE,
	_location      => $DEFAULTLOCATION,	    
	_url           => '',
	_ua            => undef,
    }, $class;
    my $host;
    my ($format,$hostlocation,$stype) = 
	$self->_rearrange([qw(FORMAT HOSTLOCATION  SERVERTYPE)],
			  @args);
    $format = $DEFAULTFORMAT unless $format;
    $stype = $DEFAULTSTYPE unless $stype;
    $self->throw("You gave an invalid server type ($stype) - available types are ".  
		 keys %HOSTS) unless( $HOSTS{$stype} );
    $self->throw('Must specify a location ('.
		 join(',', keys %{$HOSTS{$stype}->{hosts}}) .')')
	if( $hostlocation);    

    $hostlocation = $DEFAULTLOCATION unless( $hostlocation );    
    $self->stype($stype);
    $self->streamfmt($DEFAULTFORMAT);
    $self->_request_host($stype,$hostlocation);

    my $ua = new LWP::UserAgent;
    $ua->agent("$class/$MODVERSION");
    $self->ua($ua);
    return $self;
}

sub ua {
    my ($self, $ua) = @_;
    if( defined $ua && $ua->isa("LWP::UserAgent") ) {
	$self->{_ua} = $ua;
    }
    return $self->{_ua};
}

sub streamfmt {
    my ($self,$fmt) = @_;
    if( defined $fmt ) {
	$self->{_format} = $fmt;	
    }
    return $self->{_format};
}

sub stype {
    my ($self, $stype) = @_;
    if( defined $stype ) {
	$self->{_stype} = $stype;
    }
    return $self->{_stype};
}

=head2 proxy

 Title   : proxy
 Usage   : $httpproxy = $db->proxy('http')  or 
           $db->proxy(['http','ftp'], 'http://myproxy' )
 Function: Get/Set a proxy for use
 Returns : a string indicating the proxy
 Args    : $protocol : an array ref of the protocol(s) to set/get
           $proxyurl : url of the proxy to use for the specified protocol

=cut


sub proxy {
    my ($self,$protocol,$proxy) = @_;
    return undef if( !defined $self->ua );
    return $self->ua->proxy($protocol,$proxy);
}


=head2 get_url

 Title   : get_url
 Usage   : my $url = $self->get_url
 Function: returns the url, urls are set through the request_host method
 Returns : 
 Args    : 
=cut

sub get_url {
    my ($self) = @_;
    return $self->{_url} || $self->throw("Did not call _request_host, or url is set invalidly");
}

=head2 get_Seq_by_id

 Title   : get_Seq_by_id
 Usage   : $seq = $db->get_Seq_by_id($uid);
 Function: Gets a Bio::Seq object by its unique identifier/name
 Returns : a Bio::Seq object
 Args    : $uid : the id (as a string) of the desired sequence entry

=cut

sub get_Seq_by_id {

  my $self = shift;
  my $uid = shift or $self->throw("Must supply an identifier!\n");
  $self->throw("Must have requested a host") 
      unless ( defined $self->get_url()); 
  my $stream = $self->_get_stream([$uid]);
  my $seq = $stream->next_seq();
  $self->throw("Unable to get seq for id $uid, is it really a swissprot id?\n")
      if( !defined $seq );
  return $seq;
}

=head2 get_Seq_by_acc

  Title   : get_Seq_by_acc
  Usage   : $seq = $db->get_Seq_by_acc($acc);
  Function: Gets a Bio::Seq object by its accession number
  Returns : a Bio::Seq object
  Args    : $acc : the accession number of the desired sequence entry
  Note    : For GenPept, this just calls the same code for get_Seq_by_id()

=cut

sub get_Seq_by_acc {

  my $self = shift;
  my $acc = shift or $self->throw("Must supply an accesion number!\n");
  
  return $self->get_Seq_by_id($acc);
}

=head2 get_Stream_by_id

  Title   : get_Stream_by_id
  Usage   : $stream = $db->get_Stream_by_id( [$uid1, $uid2] );
  Function: Gets a series of Seq objects by unique identifiers
  Returns : a Bio::SeqIO stream object
  Args    : $ref : a reference to an array of unique identifiers for
                   the desired sequence entries


=cut

sub get_Stream_by_id {

  my $self = shift;
  my $id = shift or $self->throw("Must supply a unique identifier!\n");
  ref($id) eq "ARRAY" or $self->throw("Must supply an array ref!\n");

  return $self->_get_stream([$id]);

}

=head2 get_Stream_by_acc

  Title   : get_Stream_by_acc
  Usage   : $seq = $db->get_Seq_by_acc($acc);
  Function: Gets a series of Seq objects by accession numbers
  Returns : a Bio::SeqIO stream object
  Args    : $ref : a reference to an array of accession numbers for
                   the desired sequence entries
  Note    : For GenPept, this just calls the same code for get_Stream_by_id()

=cut

sub get_Stream_by_acc {

  my $self = shift;
  my $acc = shift or $self->throw("Must supply an accesion number!\n");

  return $self->get_Seq_by_id($acc);
}

sub _get_stream {

  my($self, $uid, $streamfmt) = @_;
  my @all;
  my $url = $self->get_url;  
  foreach my $id ( @$uid ) {
      my $resp;      
      eval { 
	  $resp = $self->ua->request(GET $url.$id);
	  push @all, $resp->content if( $resp->is_success);
      };
      if( $@ || $resp->is_error ) {
	  $self->throw($@);
      }
  }
  my $stream = new IO::String(join("\n", @all));
  return Bio::SeqIO->new('-fh' => $stream, '-format' => $self->streamfmt);
}

sub _request_host {
    my ($self,$stype,$location) = @_;    
    if( defined $stype && defined $location )  {
	$self->{_url} = sprintf($HOSTS{$stype}->{baseurl}, 
				$HOSTS{$stype}->{hosts}->{$location});
    } elsif( defined $stype ) {
	return keys %{$HOSTS{$stype}->{hosts}};
    } else {
	return $self->{_url};
    }
}

1;
__END__
















