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

  my $taxonid = $db->get_taxonid('Homo sapiens');
  my $node   = $db->get_Taxonomy_Node(-taxonid => $taxonid);

  my $gi = 71836523;
  my $node = $db->get_Taxonomy_Node(-gi => $gi, -db => 'protein');
  print $node->binomial, "\n";
  my ($species,$genus,$family) =  $node->classification;
  print "family is $family\n";

  # Can also go up 4 levels
  my $p = $node;  
  for ( 1..4 ) { 
    $p = $db->get_Taxonomy_Node(-taxonid => $p->parent_id);
  }
  print $p->rank, " ", ($p->classification)[0], "\n";

  # could then classify a set of BLAST hits based on their GI numbers
  # into taxonomic categories.


It is not currently possibly to query a node for its children so we
cannot completely replace the advantage of the flatfile
Bio::DB::Taxonomy::flatfile module.


=head1 DESCRIPTION

A driver for querying NCBI Entrez Taxonomy database.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

  bioperl-l@bioperl.org                  - General discussion
  http://bioperl.org/wiki/Mailing_lists  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution. Bug reports can be submitted via
the web:

  http://bugzilla.open-bio.org/

=head1 AUTHOR - Jason Stajich

Email jason-at-bioperl.org

=head1 CONTRIBUTORS

Sendu Bala: bix@sendu.me.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::DB::Taxonomy::entrez;
use vars qw(@ISA $EntrezLocation $UrlParamSeparatorValue %EntrezParams
	    $EntrezGet $EntrezSummary $EntrezFetch %SequenceParams
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
$EntrezFetch    = 'efetch.fcgi';
$EntrezSummary  = 'esummary.fcgi';

%EntrezParams = ( 'db'     => 'taxonomy', 
                  'report' => 'xml',
                  'retmode'=> 'xml',
                  'tool'   => 'Bioperl');

%SequenceParams = ( 'db'      => 'nucleotide', # or protein
		            'retmode' => 'xml',
		            'tool'    => 'Bioperl');

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

  $self->SUPER::_initialize(@_);

  my ($location,$params) = $self->_rearrange([qw(LOCATION PARAMS)],@_);

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
 Function: Get a Bio::Taxonomy::Node object
 Returns : Bio::Taxonomy::Node object(s) - in list context will return an array,
           in scalar context will return one of the Nodes.
 Args    : -taxonid => taxonomy id (to query by taxonid)
            OR
           -name   => string (to query by a taxonomy name: common name,
                              species, genus, etc)
           or just a single value which is the taxid.

           To retrieve a taxonomy node for a GI number
           provide the -gi option with the gi number
            and -db with either 'nucleotide' or 'protein'
            to define the db 

=cut

sub get_Taxonomy_Node {
   my ($self) = shift;
   if( ! $XMLTWIG ) {
       $self->throw("Need to have installed XML::Twig");
   }

   my %p = $self->entrez_params;
   my $taxonid;
   if( @_ > 1 ) {
       my %params = @_;
       if( $params{'-taxonid'} ) {
	   $taxonid = $params{'-taxonid'};
       } elsif( $params{'-gi'} ) {
	   my $db = $params{'-db'};
	   # we're going to do all the work here and then redirect
	   # the call based on the TaxId
	   my %p = %SequenceParams;
	   my %items;
	   if( ref($params{'-gi'}) =~ /ARRAY/i ) {	       
	       $p{'id'} = join(',', @{$params{'-gi'}});
	   } else { 
	       $p{'id'} = $params{'-gi'}; 
	   }
	   $p{'db'} = $db if defined $db;
	   my $params = join($UrlParamSeparatorValue, map { "$_=".$p{$_} } keys %p);
	   my $url = sprintf("%s%s?%s",$self->entrez_url,$EntrezSummary,$params);
	   $self->debug("url is $url\n");
	   my $response;
	   eval {
	       $response = $self->get($url);
	   };
	   if( $@ ) {
	       $self->throw("Can't query website: $@");
	   }
	   $self->debug("resp is $response\n");
	   my $twig = XML::Twig->new;
	   $twig->parse($response);
	   my $root = $twig->root;
	   my @ids;
	   for my $topnode ( $root->children('DocSum') ) {
	       for my $child ( $topnode->children('Item') ) {
		   if( uc($child->{att}->{'Name'}) eq 'TAXID' ) {
		       push @ids, $child->text;
		   }
	       }
	   }
	   return $self->get_Taxonomy_Node(-taxonid => \@ids);
       } elsif( $params{'-name'} ) {
	   my @taxaids = $self->get_taxonid($params{'-name'});
	   $taxonid = join(",", @taxaids);
       } else { 
	   $self->warn("Need to have provided either a -taxonid or -name value to get_Taxonomy_Node");
       } 
   } else { 
       $taxonid= shift;
   }
   if( ref($taxonid) =~ /ARRAY/i ) {
       $taxonid = join(',', @$taxonid);
   }
   $p{'id'}      = $taxonid;
   $self->debug("id is $taxonid\n");
   my $params = join($UrlParamSeparatorValue, map { "$_=".$p{$_} } keys %p);
   my $url = sprintf("%s%s?%s",$self->entrez_url,$EntrezFetch,$params);
   $self->debug("url is $url\n");
   my $response;
   eval {
       $response = $self->get($url);
   };
   if( $@ ) {
       $self->throw("Can't query website: $@");
   }
   $self->debug("resp is $response\n");
   my $twig = XML::Twig->new;
   $twig->parse($response);
   my $root = $twig->root;        
   my @results;
   for my $taxon ( $root->children('Taxon') ) {
       my %init = ('-dbh'            => $self,
		   '-rank'           => '',
		   '-classification' => []);
       
       for my $child ( $taxon->children ) {
	   my $name = uc ($child->gi);
	   if( $name eq 'OTHERNAMES' ) {
	       $init{'-common_names'} = [$child->children_text()];
	   } elsif( $name eq 'TAXID' ) {
	       $init{'-object_id'} = $child->text;
	   } elsif( $name eq 'PARENTTAXID' ) {
	       $init{'-parent_id'} = $child->text;
	   } elsif( $name eq 'DIVISION' ) {
	       $init{'-division'} = $child->text;
	   } elsif( $name eq 'RANK' ) {
	       $init{'-rank'} = $child->text;
	   } elsif( $name eq 'GENETICCODE' ) {
	       ($init{'-genetic_code'}) = ( map { $_->text } 
					    $child->children('GCId'));
	   } elsif( $name eq 'MITOGENETICCODE' ) {
	       ($init{'-mito_genetic_code'}) = ( map { $_->text } 
						 $child->children('MGCId'));
	   } elsif( $name eq 'CREATEDATE' ) {
	       $init{'-create_date'} = $child->text;
	   } elsif( $name eq 'UPDATEDATE' ) {
	       $init{'-update_date'} = $child->text;
	   } elsif( $name eq 'PUBDATE' ) {
	       $init{'-pub_date'} = $child->text;
	   } elsif( $name eq 'LINEAGEEX' ) {
	       for my $taxon ( $child->children('Taxon') ) {
		   my ($id) = map { $_->text } $taxon->children('TaxId');
		   my ($sname) = map { $_->text } $taxon->children('ScientificName');
		   my ($rank) = map { $_->text } $taxon->children('Rank');
		   next if( $rank eq 'no rank' );
		   unshift @{$init{'-classification'}}, $sname;
	       }
	   } elsif( $name eq 'SCIENTIFICNAME' ) {
           my $sci_name = $child->text;
           my $orig_sci_name = $sci_name;
           $sci_name =~ s/ \(class\)$//;
           push(@{$init{'-common_names'}}, $orig_sci_name) if $orig_sci_name ne $sci_name;
           $init{'-name'} = $sci_name;
	   } 
       }
       
       if( ! $init{'-object_id'} ) { 
	   $self->warn("Could not find any value for $taxonid");
	   return;
       }
       
       unless ($init{'-rank'} eq 'no rank') {
            unshift @{$init{'-classification'}}, $init{'-name'};	   
       }
       
       push @results, Bio::Taxonomy::Node->new(%init);
   }
   ( wantarray ) ? @results : shift @results;
}

=head2 get_taxonids

 Title   : get_taxonids
 Usage   : my $taxonid = $db->get_taxonids('Homo sapiens');
 Function: Searches for a taxonid (typically ncbi_taxon_id) based on a query
           string. Note that multiple taxonids can match to the same supplied
           name.
 Returns : array of integer ids in list context, one of these in scalar context
 Args    : string representing taxanomic (node) name

=cut

sub get_taxonids {
    my ($self,$query) = @_;
    my %p = $self->entrez_params;
    
    # queries don't work correctly with special characters, so get rid of them.
    if ($query =~ /<.+>/) {
        # queries with <something> will fail, so workaround by removing, doing
        # the query, getting multiple taxonids, then picking the one id that
        # has a parent node with a scientific_name() or common_names()
        # case-insensitive matching to the word(s) within <>
        $query =~ s/ <(.+?)>//;
        my $desired_parent_name = lc($1);
        
        ID: foreach my $start_id ($self->get_taxonids($query)) {
            my $node = $self->get_Taxonomy_Node($start_id) || next ID;
            
            # walk up the parents until we hit a node with a named rank
            while (1) {
                my $parent_id = $node->parent_id || next ID;
                my $parent_node = $self->get_Taxonomy_Node($parent_id) || next ID;
                my $parent_sci_name = $parent_node->scientific_name || next ID;
                my @parent_common_names = $parent_node->common_names;
                
                foreach my $name ($parent_sci_name, @parent_common_names) {
                    if (lc($name) eq $desired_parent_name) {
                        return wantarray() ? ($start_id) : $start_id;
                    }
                }
                
                my $parent_rank = $parent_node->rank || '';
                if ($parent_rank && $parent_rank ne 'no rank') {
                    last;
                }
                else {
                    $node = $parent_node;
                }
            }
        }
        return;
    }
    $query =~ s/[\"\(\)]//g; # not an exhaustive list; these are just the ones I know cause problems
    $query =~ s/\s/+/g;
    
    $p{'term'} = $query;
    my $params = join($UrlParamSeparatorValue, map { "$_=".$p{$_} } keys %p);
    my $url = sprintf("%s%s?%s",$self->entrez_url,$EntrezGet,$params);
    my $response;
    eval {
        $response = $self->get($url);
    };
    if( $@ ) {
        $self->throw("Can't query website: $@");
    }
    $self->debug("response is $response\n");
    my $twig = XML::Twig->new;
    $twig->parse($response);
    my $root = $twig->root;
    my $list = $root->first_child('IdList');
    my @data = map { $_->text } $list->children('Id');
    ( wantarray ) ? @data : shift @data;
}

*get_taxonid = \&get_taxonids;

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
