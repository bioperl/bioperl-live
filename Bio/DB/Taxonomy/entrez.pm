#
# BioPerl module for Bio::DB::Taxonomy::entrez
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
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

  my $db = Bio::DB::Taxonomy->new(-source => 'entrez');

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

=head2 Support 

Please direct usage questions or support issues to the mailing list:

I<bioperl-l@bioperl.org>

rather than to the module maintainer directly. Many experienced and 
reponsive experts will be able look at the problem and quickly 
address it. Please include a thorough description of the problem 
with code and data examples if at all possible.

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution. Bug reports can be submitted via
the web:

  https://github.com/bioperl/bioperl-live/issues

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
use vars qw($EntrezLocation $UrlParamSeparatorValue %EntrezParams
            $EntrezGet $EntrezSummary $EntrezFetch %SequenceParams
            $XMLTWIG $DATA_CACHE $RELATIONS);
use strict;
use Bio::Taxon;
use Bio::DB::NCBIHelper;

eval {
    require XML::Twig;
    $XMLTWIG = 1;
};
if( $@ ) {
    $XMLTWIG = 0;
}

use base qw(Bio::WebAgent Bio::DB::Taxonomy);

$EntrezLocation = $Bio::DB::NCBIHelper::HOSTBASE . '/entrez/eutils/';
$EntrezGet      = 'esearch.fcgi';
$EntrezFetch    = 'efetch.fcgi';
$EntrezSummary  = 'esummary.fcgi';

$DATA_CACHE = {};
$RELATIONS  = {};

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
 Usage   : my $obj = Bio::DB::Taxonomy::entrez->new();
 Function: Builds a new Bio::DB::Taxonomy::entrez object
 Returns : an instance of Bio::DB::Taxonomy::entrez
 Args    : -location => URL to Entrez (if you want to override the default)
           -params   => Hashref of URL params if you want to override the
                        default

=cut

sub new {
    my ($class, @args) = @_;
    
    # need to initialise Bio::WebAgent...
    my ($self) = $class->SUPER::new(@args);
    
    # ... as well as our normal Bio::DB::Taxonomy selves:
    $self->_initialize(@args);
    return $self;
}


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


=head2 get_num_taxa

 Title   : get_num_taxa
 Usage   : my $num = $db->get_num_taxa();
 Function: Get the number of taxa stored in the database.
 Returns : A number
 Args    : None

=cut

sub get_num_taxa {
    my ($self) = @_;
    # Use this URL query to get the ID of all the taxa in the NCBI Taxonomy database:
    # http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=taxonomy&term=all[Filter]
    # Only the first 20 taxa IDs are returned (good because the list is long),
    # and the total number is reported as well (which is what we are interested
    # in).
    my %p = $self->entrez_params;
    $p{'term'} = 'all[Filter]'; 
    my $twig = $self->_run_query($self->_build_url($EntrezGet, \%p));
    my $count = $twig->root->first_child('Count')->first_child->text;
    return $count;
}


=head2 get_taxon

 Title   : get_taxon
 Usage   : my $taxon = $db->get_taxon(-taxonid => $taxonid)
 Function: Get a Bio::Taxon object from the database.
 Returns : Bio::Taxon object
 Args    : just a single value which is the database id, OR named args:
           -taxonid => taxonomy id (to query by taxonid)
            OR
           -name    => string (to query by a taxonomy name: common name, 
                               scientific name, etc)
            OR
           To retrieve a taxonomy node for a GI number provide the -gi option
           with the gi number and -db with either 'nucleotide' or 'protein' to
           define the db.
            AND optionally,
           -full    => 1 (to force retrieval of full information - sometimes
                          minimal information about your taxon may have been
                          cached, which is normally used to save database
                          accesses)

=cut

sub get_taxon {
    my $self = shift;
    if (! $XMLTWIG) {
        eval { require XML::Twig };
        $self->throw("Could not load XML::Twig for get_taxon(): $@") if $@;
    }

    my %p = $self->entrez_params;

    # convert input request to one or more ids
    my (@taxonids, $taxonid, $want_full);
    if (@_ > 1) {
        my %params = @_;
        if ($params{'-taxonid'}) {
            $taxonid = $params{'-taxonid'};
        }
        elsif ($params{'-gi'}) {
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
            my $url = $self->_build_url($EntrezSummary, \%p);
            my @ids;
            if (exists $DATA_CACHE->{gi_to_ids}->{$url}) {
                @ids = @{$DATA_CACHE->{gi_to_ids}->{$url}};
            }
            else {
                my $twig = $self->_run_query($url);
                my $root = $twig->root;                
                for my $topnode ( $root->children('DocSum') ) {
                    for my $child ( $topnode->children('Item') ) {
                        if( uc($child->{att}->{'Name'}) eq 'TAXID' ) {
                            push @ids, $child->text;
                        }
                    }
                }
                
                $DATA_CACHE->{gi_to_ids}->{$url} = \@ids;
            }
            
            return $self->get_taxon(-taxonid => \@ids);
        }
        elsif ($params{'-name'}) {
            @taxonids = $self->get_taxonid($params{'-name'});
        }
        else { 
            $self->warn("Need to have provided either a -taxonid or -name value to get_taxon");
        }
        
        if ($params{'-full'}) {
            $want_full = 1;
        }
    }
    else {
        $taxonid = shift;
    }
    
    if (ref($taxonid) =~ /ARRAY/i ) {
        @taxonids = @{$taxonid};
    }
    else {
        push(@taxonids, $taxonid) if $taxonid;
    }
    
    # return answer(s) from the cache if possible
    my @results;
    my @uncached;
    foreach my $taxonid (@taxonids) {
        $taxonid || $self->throw("In taxonids list one was undef! '@taxonids'\n");
        if (defined $DATA_CACHE->{full_info}->{$taxonid}) {
            push(@results, $self->_make_taxon($DATA_CACHE->{full_info}->{$taxonid}));
        }
        elsif (! $want_full && defined $DATA_CACHE->{minimal_info}->{$taxonid}) {
            push(@results, $self->_make_taxon($DATA_CACHE->{minimal_info}->{$taxonid}));
        }
        else {
            push(@uncached, $taxonid);
        }
    }
    
    if (@uncached > 0) {
        $taxonid = join(',', @uncached);
        
        $p{'id'} = $taxonid;
        $self->debug("id is $taxonid\n");

        my $twig = $self->_run_query($self->_build_url($EntrezFetch, \%p));
        my $root = $twig->root;
        for my $taxon ( $root->children('Taxon') ) {
            my $taxid = $taxon->first_child_text('TaxId');
            $self->throw("Got a result with no TaxId!") unless $taxid;
            
            my $data = {};
            if (exists $DATA_CACHE->{minimal_info}->{$taxid}) {
                $data = $DATA_CACHE->{minimal_info}->{$taxid};
            }
            
            $data->{id} = $taxid;
            $data->{rank} = $taxon->first_child_text('Rank');
            
            my $other_names = $taxon->first_child('OtherNames');
            my @other_names = $other_names->children_text() if $other_names;
            my $sci_name = $taxon->first_child_text('ScientificName');
            my $orig_sci_name = $sci_name;
            $sci_name =~ s/ \(class\)$//;
            push(@other_names, $orig_sci_name) if $orig_sci_name ne $sci_name;
            $data->{scientific_name} = $sci_name;
            $data->{common_names} = \@other_names;
            
            $data->{division} = $taxon->first_child_text('Division');
            $data->{genetic_code} = $taxon->first_child('GeneticCode')->first_child_text('GCId');
            $data->{mitochondrial_genetic_code} = $taxon->first_child('MitoGeneticCode')->first_child_text('MGCId');
            $data->{create_date} = $taxon->first_child_text('CreateDate');
            $data->{update_date} = $taxon->first_child_text('UpdateDate');
            $data->{pub_date} = $taxon->first_child_text('PubDate');
            
            # since we have some information about all the ancestors of our
            # requested node, we may as well cache data for the ancestors to
            # reduce the number of accesses to website in future
            my $lineage_ex = $taxon->first_child('LineageEx');
            my ($ancestor, $lineage_data, @taxa);
            foreach my $lineage_taxon ($lineage_ex->children) {
                my $lineage_taxid = $lineage_taxon->first_child_text('TaxId');
                
                if (exists $DATA_CACHE->{minimal_info}->{$lineage_taxid} || exists $DATA_CACHE->{full_info}->{$lineage_taxid}) {
                    $lineage_data = $DATA_CACHE->{minimal_info}->{$lineage_taxid} || $DATA_CACHE->{full_info}->{$lineage_taxid};
                    next;
                }
                else {
                    $lineage_data = {};
                }
                
                $lineage_data->{id} = $lineage_taxid;
                $lineage_data->{scientific_name} = $lineage_taxon->first_child_text('ScientificName');
                $lineage_data->{rank} = $lineage_taxon->first_child_text('Rank');
                
                $RELATIONS->{ancestors}->{$lineage_taxid} = $ancestor->{id} if $ancestor;
                
                $DATA_CACHE->{minimal_info}->{$lineage_taxid} = $lineage_data;
            } continue { $ancestor = $lineage_data; unshift(@taxa, $lineage_data); }
            
            $RELATIONS->{ancestors}->{$taxid} = $ancestor->{id} if $ancestor;
            
            # go through the lineage in reverse so we can remember the children
            my $child = $data;
            foreach my $lineage_data (@taxa) {
                $RELATIONS->{children}->{$lineage_data->{id}}->{$child->{id}} = 1;
            } continue { $child = $lineage_data; }
            
            delete $DATA_CACHE->{minimal_info}->{$taxid};
            $DATA_CACHE->{full_info}->{$taxid} = $data;
            push(@results, $self->_make_taxon($data));
        }
    }
    
    wantarray() ? @results : shift @results;
}

*get_Taxonomy_Node = \&get_taxon;


=head2 get_taxonids

 Title   : get_taxonids
 Usage   : my $taxonid = $db->get_taxonids('Homo sapiens');
 Function: Searches for a taxonid (typically ncbi_taxon_id) based on a query
           string. Note that multiple taxonids can match to the same supplied
           name.
 Returns : array of integer ids in list context, one of these in scalar context
 Args    : string representing taxon's name

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
            my $node = $self->get_taxon($start_id) || next ID;
            
            # walk up the parents until we hit a node with a named rank
            while (1) {
                my $parent_node = $self->ancestor($node) || next ID;
                my $parent_sci_name = $parent_node->scientific_name || next ID;
                my @parent_common_names = $parent_node->common_names;
                unless (@parent_common_names) {
                    # ensure we're not using a minimal-info cached version
                    $parent_node = $self->get_taxon(-taxonid => $parent_node->id, -full => 1);
                    @parent_common_names = $parent_node->common_names;
                }
                
                foreach my $name ($parent_sci_name, @parent_common_names) {
                    if (lc($name) eq $desired_parent_name) {
                        return wantarray() ? ($start_id) : $start_id;
                    }
                }
                
                my $parent_rank = $parent_node->rank || 'no rank';
                if ($parent_rank ne 'no rank') {
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
    
    my @data;
    if (defined $DATA_CACHE->{name_to_id}->{$query}) {
        @data = @{$DATA_CACHE->{name_to_id}->{$query}};
    }
    else {
        $p{'term'} = $query;

        my $twig = $self->_run_query($self->_build_url($EntrezGet, \%p));
        my $root = $twig->root;
        my $list = $root->first_child('IdList');
        @data = map { $_->text } $list->children('Id');
        
        $DATA_CACHE->{name_to_id}->{$query} = [@data];
    }
    
    return wantarray() ? @data : shift @data;
}

*get_taxonid = \&get_taxonids;


=head2 ancestor

 Title   : ancestor
 Usage   : my $ancestor_taxon = $db->ancestor($taxon)
 Function: Retrieve the ancestor taxon of a supplied Taxon from the database.

           Note that unless the ancestor has previously been directly
           requested with get_taxon(), the returned Taxon object will only have
           a minimal amount of information.

 Returns : Bio::Taxon
 Args    : Bio::Taxon (that was retrieved from this database)

=cut

sub ancestor {
    my ($self, $taxon) = @_;
    $self->throw("Must supply a Bio::Taxon") unless ref($taxon) && $taxon->isa('Bio::Taxon');
    $self->throw("The supplied Taxon must belong to this database") unless $taxon->db_handle && $taxon->db_handle eq $self;
    my $id = $taxon->id || $self->throw("The supplied Taxon is missing its id!");
    
    my $ancestor_id = $RELATIONS->{ancestors}->{$id} || return;
    return $self->_make_taxon($DATA_CACHE->{full_info}->{$ancestor_id} || $DATA_CACHE->{minimal_info}->{$ancestor_id});
}


=head2 each_Descendent

 Title   : each_Descendent
 Usage   : my @taxa = $db->each_Descendent($taxon);
 Function: Get all the descendents of the supplied Taxon (but not their
           descendents, ie. not a recursive fetchall).

           Note that this implementation is unable to return a taxon that
           hasn't previously been directly fetched with get_taxon(), or wasn't
           an ancestor of such a fetch.

 Returns : Array of Bio::Taxon objects
 Args    : Bio::Taxon (that was retrieved from this database)

=cut

sub each_Descendent {
    my ($self, $taxon) = @_;
    $self->throw("Must supply a Bio::Taxon") unless ref($taxon) && $taxon->isa('Bio::Taxon');
    $self->throw("The supplied Taxon must belong to this database") unless $taxon->db_handle && $taxon->db_handle eq $self;
    my $id = $taxon->id || $self->throw("The supplied Taxon is missing its id!");
    
    my @children_ids = keys %{$RELATIONS->{children}->{$id} || {}};
    my @children;
    foreach my $child_id (@children_ids) {
        push(@children, $self->_make_taxon($DATA_CACHE->{full_info}->{$child_id} || $DATA_CACHE->{minimal_info}->{$child_id}));
    }
    
    return @children;
}


=head2 Some Get/Setter methods

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


# make a Taxon object from data hash ref
sub _make_taxon {
    my ($self, $data) = @_;
    
    my $taxon = Bio::Taxon->new();
    
    my $taxid;
    while (my ($method, $value) = each %{$data}) {
        if ($method eq 'id') {
            $method = 'ncbi_taxid'; # since this is a real ncbi taxid, explicitly set it as one
            $taxid = $value;
        }
        $taxon->$method(ref($value) eq 'ARRAY' ? @{$value} : $value);
    }
    
    # we can't use -dbh or the db_handle() method ourselves or we'll go
    # infinite on the merge attempt
    $taxon->{'db_handle'} = $self;
    
    $self->_handle_internal_id($taxon);
    
    return $taxon;
}


sub _build_url {
    # Given a eutility (esearch.fcgi, efetch.fcgi or esummary.fcgi) and a
    # hashref or parameters, build a url suitable for eutil query
    my ($self, $eutility, $p) = @_;
    my $params = join($UrlParamSeparatorValue, map { $_.'='.$p->{$_} } keys %$p);
    my $url = $self->entrez_url.$eutility.'?'.$params;
    $self->debug("url is $url\n");
    return $url;
}


sub _run_query {
    # Given an eutil url, run the eutil query and parse the response into an
    # XML Twig object
    my ($self, $url) = @_;
    my $response = $self->get($url);
    if ($response->is_success) {
        $response = $response->content;
    }else {
        $self->throw("Can't query website: ".$response->status_line);
    }
    $self->debug("response is $response\n");
    my $twig = XML::Twig->new;
    $twig->parse($response);
    return $twig;
}


1;
