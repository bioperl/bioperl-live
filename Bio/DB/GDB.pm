# $Id$
#
# BioPerl module for Bio::DB::GDB
#
# Cared for by Jason Stajich <jason@bioperl.org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code
# 

=head1 NAME

Bio::DB::GDB - Database object interface to GDB HTTP query

=head1 SYNOPSIS

    use Bio::DB::GDB;

    $gdb = Bio::DB::GDB->new();

    $info = $gdb->get_info(-type => 'marker',
			                  -id => 'D1S243'); # Marker name

   print "genbank id is ", $info->{'gdbid'},
    "\nprimers are (fwd, rev) ", join(",", @{$info->{'primers'}}), 
    "\nproduct length is ", $info->{'length'}, "\n";

=head1 DESCRIPTION

This class allows connections to the Genome Database (GDB) and queries
to retrieve any database objects. See http://www.gdb.org/ or any
mirror for details.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the
evolution of this and other Bioperl modules. Send
your comments and suggestions preferably to one
of the Bioperl mailing lists. Your participation
is much appreciated.

  bioperl-l@bioperl.org                  - General discussion
  http://bioperl.org/wiki/Mailing_lists  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via the
web:

  http://bugzilla.open-bio.org/

=head1 AUTHOR - Jason Stajich

Email jason@bioperl.org

=head1 APPENDIX

The rest of the documentation details each of the
object methods. Internal methods are usually
preceded with a _

=cut

# Let the code begin...

package Bio::DB::GDB;
use strict;
use LWP::UserAgent;
use HTTP::Request::Common;
use HTML::Parser;

use vars qw($BASEADDRESS %PARAMSTRING $MODVERSION);

use base qw(Bio::Root::Root);

$MODVERSION = '0.01';
$BASEADDRESS = 'http://www.gdb.org/gdb-bin/genera/genera/hgd/GenomicSegment';
%PARAMSTRING = ( 
		 gene   => { '!action' => 'query' }, 
		 marker => { '!action' => 'query' },
		 );

# the new way to make modules a little more lightweight
sub new {
    my($class,@args) = @_;
    my $self = $class->SUPER::new(@args);
    
    my $ua = new LWP::UserAgent(env_proxy => 1);
    $ua->agent(ref($self) ."/$MODVERSION");
    $self->ua($ua);    

    return $self;
}

=head2 ua

 Title   : ua
 Usage   : my $ua = $self->ua or 
           $self->ua($ua)
 Function: Get/Set a LWP::UserAgent for use
 Returns : reference to LWP::UserAgent Object
 Args    : $ua - must be a LWP::UserAgent

=cut

sub ua {
    my ($self, $ua) = @_;
    if( defined $ua && $ua->isa("LWP::UserAgent") ) {
	$self->{_ua} = $ua;
    }
    return $self->{_ua};
}

# helper method to get specific options

=head2 get_params

 Title   : get_params
 Usage   : my %params = $self->get_params($mode)
 Function: Returns key,value pairs to be passed to query
            for mode ('marker', 'gene')
 Returns : a key,value pair hash
 Args    : 'marker' or 'gene' mode for retrieval

=cut

sub get_params {
    my ($self, $mode) = @_;
    return %{$PARAMSTRING{$mode}};
}

=head2 get_info

 Title   : get_info
 Usage   : my $info = $self->get_info(-type => 'marker',
				      -id   => 'D1S234'); 
 Function: Returns key,value pairs specific
 Returns : a key,value pair hash
 Args    : -type => 'marker' or 'gene' mode for retrieval
           -id   => unique id to query for

=cut

sub get_info {
    my ($self, @args) = @_;
    my ( $type, $id) = $self->_rearrange([qw(TYPE ID)], @args);
    if( !defined $type ) {
	$self->throw("Must specify a type you are querying for");
    } elsif( !defined $id ) {
	$self->throw("Must specify a id to query for");
    }
    my %params = $self->get_params($type);

    $params{'displayName'} = $id;

    if( $type eq 'marker' ) {
	# do more specific stuff?
    } elsif( $type eq 'gene' ) {
	# do more specific stuff?
    }
    my $url = $self->get_request(%params);    
    
    my ($resp) = $self->_request($url);
    if( ! defined $resp || ! ref($resp) ) {
	$self->warn("Did not get any data for url ". $url->uri);
	return;
    }
    my $content = $resp->content;	
    if( $content =~ /ERROR/ || length($resp->content) == 0 ) {
	$self->warn("Error getting for url " . $url->uri . "!\n");
	return;
    }
    my (@primers, $length, $markerurl, $realname);
    my $state = 0;
    my $title = 0;
    my $p;
    $p = new HTML::Parser( api_version => 3,
			   start_h => [ sub { 
			       return if( $title == 2 || $state == 3);
			       my($tag,$attr,$text) = @_;
			       return if( !defined $tag);
			       if( $tag eq 'table' ) {
				   $state = 1;
			       } elsif( $tag eq 'title' ) {
				   $title = 1;
			       } elsif( $state == 2 && 
					$tag eq 'a' &&
					$attr->{'href'} ) {
				   $state = 3; 
				   if( $text =~ m(href="?(http://.+)"?\s*>) ) { 
				       $markerurl = $1;
				   }
			       } 
			   }, "tagname, attr, text" ],
			   end_h   => [ sub { 
			       return if ($title == 2 || $state == 3);
			       my ( $tag ) = @_;
			       $title = 0 if( $tag eq 'title' );
			   }, "tagname" ],
			   text_h  => [ sub { 
			       return if( $title == 2 || $state == 3);
			       my($text) = @_;
			       if( $title && $text =~ /Amplimer/ ) {
				   $markerurl = 'this';
				   $title = 2;
			       }
			       $state = 2 if( $state == 1 && $text =~ /Amplimer/);
			   }, "text" ],
			   marked_sections =>1);
    $p->parse($content) or $self->throw("Can't open: $!");        
    if( ! defined $markerurl ) {
	@primers = ('notfound','notfound', '?');
    } elsif( $markerurl eq 'this' ) {

    }
    else { 
	my $resp = $self->_request(GET $markerurl);
        return if ( !defined $resp );
	$content = $resp->content();
    }
    $state = 0;
    $realname = 'unknown';
    my $lasttag = '';
    $p = HTML::Parser->new(api_version => 3,			      
			   start_h => [ sub { my ($tag) = @_;
					      $tag = lc $tag;
					      $lasttag = $tag;
					      if( $state == 3 && $tag eq 'dd' ) {
						  $state = 4;
					      }
					  } , 'tagname'],			   
			   text_h  => [ sub { 
			       my($text) = @_; 
			       if( $text =~ /Primer Sequence/ ) {
				   $state =1;
			       } elsif( $state == 1 ) {
				   foreach my $l ( split(/\n+/,$text) ) {
				       $l =~ s/\s+(\S+)/$1/;
				       my ($name,$primer) = split(/\s+/,$l);
				       next if( !defined $name);
				       push @primers, $primer;
				       $state = 2;
				   }
			       } elsif( $state == 2 && 
					($text =~ /Seq Min Len/i ||
					 $text =~ /Seq Max Len/i) ) {
				   $state = 3;
			       }  elsif ( $state == 4 ) {
				   my ($len) = ( $text =~ /(\d+\.\d+)/
);
				   $length = $len;
				   $length *= 1000 if( $len < 1 );
				   $state = 0;
			       } elsif( $lasttag eq 'dd' && 
					$text =~ /(GDB:\d+)/i ) {
				   $realname = $1;
			       }
			   }  , "text" ],
			   marked_sections =>1,
			   );
    $p->parse($content) || $self->throw("Can't open: $!");

    return { 'gdbid' => $realname, 'length' => $length, 'primers' => \@primers };
}

=head2 get_request

 Title   : get_request
 Usage   : my $url = $self->get_request
 Function: HTTP::Request
 Returns : 
 Args    : %qualifiers = a hash of qualifiers (ids, format, etc)

=cut

sub get_request {
    my ($self, %params) = @_;
    if( ! %params ) {
	$self->throw("must provide parameters with which to query");
    }
    my $url = $BASEADDRESS;    
    my $querystr = '?' . join("&", map { "$_=$params{$_}" } keys %params);
    return GET $url . $querystr;
}

# private methods
sub _request {

    my ($self, $url,$tmpfile) = @_;
    my ($resp);
    if( defined $tmpfile && $tmpfile ne '' ) { 
	$resp =  $self->ua->request($url, $tmpfile);
    } else { $resp =  $self->ua->request($url); } 

    if( $resp->is_error  ) {
	$self->throw($resp->as_string() . "\nError getting for url " .
		     $url->uri . "!\n");
	return;
    }
    return $resp;
}

sub _gdb_search_tag_start {

}

1;
__END__
