# $Id$
#
# BioPerl module for Bio::DB::Taxonomy::entrez
#
# Cared for by Jason Stajich <jason-at-bioperl.org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::DB::Taxonomy::entrez - Taxonomy Entrez driver

=head1 SYNOPSIS

# Do not use this object directly, rather through the Bio::DB::Taxonomy
# interface

  use Bio::DB::Taxonomy;

  my $db = new Bio::DB::Taxonomy(-source => 'entrez');

  my $taxaid = $db->get_taxaid('Homo sapiens');

=head1 DESCRIPTION

A driver for querying NCBI Entrez Taxonomy database.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

  bioperl-l@bioperl.org              - General discussion
  http://bioperl.org/MailList.shtml  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution. Bug reports can be submitted via
the web:

  http://bugzilla.bioperl.org/

=head1 AUTHOR - Jason Stajich

Email jason-at-bioperl.org

=head1 CONTRIBUTORS

Additional contributors names and emails here

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::DB::Taxonomy::entrez;
use vars qw(@ISA $EntrezLocation $UrlParamSeparatorValue %EntrezParams
	    $EntrezGet $EntrezSummary
	    $XMLTWIG);
use strict;

use Bio::DB::Taxonomy;
use Bio::Root::HTTPget;
use Bio::Taxonomy::Node;

eval {
    require XML::Twig;
    $XMLTWIG = 1;
};
if( $@ ) {
    $XMLTWIG = 0;
}
@ISA = qw(Bio::DB::Taxonomy Bio::Root::HTTPget);

$EntrezLocation = 'http://www.ncbi.nih.gov/entrez/eutils/';
$EntrezGet      = 'esearch.fcgi';
$EntrezSummary  = 'esummary.fcgi';

%EntrezParams = ( 'db' => 'taxonomy');
$UrlParamSeparatorValue = '&';

=head2 new

 Title   : new
 Usage   : my $obj = new Bio::DB::Taxonomy::entrez();
 Function: Builds a new Bio::DB::Taxonomy::entrez object
 Returns : an instance of Bio::DB::Taxonomy::entrez
 Args    : -location => URL to Entrez (if you want to override the default)
           -params   => Hashref of URL params if you want to override the
                        default

=cut

sub _initialize {
  my($self) = shift;
  if( ! $XMLTWIG ) {
      $self->throw("Need to have installed XML::Twig");
  }

  $self->SUPER::_initialize(@_);

  my ($location,$params) = $self->_rearrange([qw(LOCATION PARAMS)],@_);

  $self->entrez_url($location || $EntrezLocation );
  if( $params ) {
      if( ref($params) !~ /HASH/i ) {
	  $self->warn("Must have provided a valid HASHref for -params");
	  $params = \%EntrezParams;
      }
  } else {
      $params = \%EntrezParams;
  }
  $self->entrez_params($params);
  $self->entrez_url($location || $EntrezLocation );
}


=head2 get_Taxonomy_Node

 Title   : get_Taxonomy_Node
 Usage   : my $species = $db->get_Taxonomy_Node(-taxonid => $taxonid)
 Function: Get a Bio::Taxonomy::Taxon object
 Returns : Bio::Taxonomy::Taxon object(s) [more than one
 Args    : -taxonid => taxonomy id (to query by taxonid)
            OR
           -name   => string (to query by a taxonomy name: common name,
                              species, genus, etc)
           or just a single value which is the taxid.

=cut

sub get_Taxonomy_Node{
   my ($self) = shift;
   my %p = $self->entrez_params;
   my $taxonid;
   if( @_ > 1 ) {
       my %params = @_;
       if( $params{'-taxonid'} ) {
	   $taxonid = $params{'-taxonid'};
       } elsif( $params{'-name'} ) {
	   my @taxaids = $self->get_taxonid($params{'-name'});
	   if( @taxaids > 1 ) { 
	       $self->warn("Got > 1 taxid for ".$params{'-name'}. " only using the first one");
	   }
	   $taxonid = shift @taxaids;
       } else { 
	   $self->warn("Need to have provided either a -taxonid or -name value to get_Taxonomy_Node");
       } 
   } else { 
       $taxonid= shift;
   }
   $p{'id'}      = $taxonid;

   my $params = join($UrlParamSeparatorValue, map { "$_=".$p{$_} } keys %p);
   my $url = sprintf("%s%s?%s",$self->entrez_url,$EntrezSummary,$params);
   $self->debug("url is $url\n")  if( $self->verbose > 0);
   my $response;
   eval {
       $response = $self->get($url);
   };
   if( $@ ) {
       $self->warn("Can't query website: $@");
       return;
   }
   my $twig = new XML::Twig;
   $self->debug( "resp is $response\n") if( $self->verbose > 0);
   $twig->parse($response);
   my $root = $twig->root;
   my $list = $root->first_child('DocSum');
   if( ! $list ) { 
       $self->warn("Could not find any value for $taxonid");
       return undef;
   }
   my ($id) = map { $_->text } $list->children('Id');
   my (%item) = map {  uc($_->{'att'}->{'Name'}) => $_->text } $list->children('Item');
   if( $item{'RANK'} eq 'species') {
       my $node = Bio::Taxonomy::Node->new(-name      => $item{'COMMONNAME'},
					   -object_id => $item{'TAXID'},
					   -parent_id => undef,
					   -rank      => $item{'RANK'},
					   -division  => $item{'DIVISION'},
					   -dbh       => $self);
       $node->classification(reverse split(/\s+/,$item{'SCIENTIFICNAME'}));
       return $node;
   } else {
       $self->warn(sprintf("can't create a species object for %s (%s) because it isn't a species but is a '%s' instead",$item{'SCIENTIFICNAME'},
			   $item{'COMMONNAME'}, $item{'RANK'}));
   }
   \%item;
}


=head2 get_taxonid

 Title   : get_taxonid
 Usage   : my $taxonid = $db->get_taxonid('Homo sapiens');
 Function: Searches for a taxonid (typically ncbi_taxon_id)
           based on a query string
 Returns : Integer ID
 Args    : Array of Strings representing species/node name


=cut

sub get_taxonid {
   my ($self,$query) = @_;
   my %p = $self->entrez_params;
   $query        =~ s/\s/\+/g;
   $p{'term'}      = $query;
   my $params = join($UrlParamSeparatorValue, map { "$_=".$p{$_} } keys %p);
   my $url = sprintf("%s%s?%s",$self->entrez_url,$EntrezGet,$params);
   my $response;
   eval {
       $response = $self->get($url);
   };
   if( $@ ) {
       $self->warn("Can't query website: $@");
       return;
   }
   $self->debug( "response is $response\n") if( $self->verbose > 0);

   my $twig = new XML::Twig;
   $twig->parse($response);
   my $root = $twig->root;
   my $list = $root->first_child('IdList');
   my @data = map { $_->text } $list->children('Id');
   ( wantarray ) ? @data : shift @data;
}



=head2 Some Get/Setter methods

=cut

=head2 entrez_url

 Title   : entrez_url
 Usage   : $obj->entrez_url($newval)
 Function: Get/set entrez URL
 Returns : value of entrez url (a scalar)
 Args    : on set, new value (a scalar or undef, optional)


=cut

sub entrez_url{
    my $self = shift;

    return $self->{'_entrez_url'} = shift if @_;
    return $self->{'_entrez_url'};
}

=head2 entrez_params

 Title   : entrez_params
 Usage   : $obj->entrez_params($newval)
 Function: Get/set entrez params
 Returns : value of entrez_params (a hashref)
 Args    : on set, new value Hashref


=cut

sub entrez_params{
    my $self = shift;
    my $f;
    if( @_ ) {
	$f = $self->{'_entrez_params'} = shift;
    } else {
	$f = $self->{'_entrez_params'};
    }
    return %$f;
}

=head2 Bio::DB::WebBase methods

=cut

=head2 proxy_string

 Title   : proxy_string
 Usage   : my $proxy_string = $self->proxy_string($protocol)
 Function: Get the proxy string (plus user/pass )
 Returns : string
 Args    : protocol ('http' or 'ftp'), default 'http'


=head2 proxy

 Title   : proxy
 Usage   : $httpproxy = $db->proxy('http')  or
           $db->proxy(['http','ftp'], 'http://myproxy' )
 Function: Get/Set a proxy for use of proxy
 Returns : a string indicating the proxy
 Args    : $protocol : an array ref of the protocol(s) to set/get
           $proxyurl : url of the proxy to use for the specified protocol
           $username : username (if proxy requires authentication)
           $password : password (if proxy requires authentication)


=head2 authentication

 Title   : authentication
 Usage   : $db->authentication($user,$pass)
 Function: Get/Set authentication credentials
 Returns : Array of user/pass
 Args    : Array or user/pass


=cut

1;
