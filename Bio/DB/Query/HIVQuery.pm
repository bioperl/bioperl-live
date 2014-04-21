# to do: support for comment, reference annotations

# $Id: HIVQuery.pm 232 2008-12-11 14:51:51Z maj $
#
# BioPerl module for Bio::DB::Query::LANLQuery
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Mark A. Jensen <maj@fortinbras.us>
#
# Copyright Mark A. Jensen
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::DB::Query::HIVQuery - Query interface to the Los Alamos HIV Sequence Database

=head1 SYNOPSIS

    $q = new Bio::DB::Query::HIVQuery(" C[subtype] ZA[country] CXCR4[coreceptor] ");
    $q = new Bio::DB::Query::HIVQuery(
         -query=>{'subtype'=>'C', 
                  'country'=>'ZA', 
                  'coreceptor'=>'CXCR4'});

    $ac = $q->get_annotations_by_id(($q->ids)[0]);
    $ac->get_value('Geo', 'country')                    # returns 'SOUTH AFRICA'

    $db = new Bio::DB::HIV();
    $seqio = $db->get_Stream_by_query($q);              # returns annotated Bio::Seqs 

    # get subtype C sequences from South Africa and Brazil, 
    # with associated info on patient health, coreceptor use, and 
    # infection period:

    $q = new Bio::DB::Query::HIVQuery(
         -query => {
                    'query' => {'subtype'=>'C',
		    'country'=>['ZA', 'BR']},
                    'annot' => ['patient_health', 
                                'coreceptor', 
                                'days_post_infection']
                    });
	

=head1 DESCRIPTION

Bio::DB::Query::HIVQuery provides a query-like interface to the
cgi-based Los Alamos National Laboratory (LANL) HIV Sequence
Database. It uses Bioperl facilities to capture both sequences and
annotations in batch in an automated and computable way. Use with
L<Bio::DB::HIV> to create C<Bio::Seq> objects and annotated C<Bio::SeqIO>
streams.

=head2 Query format

The interface implements a simple query language emulation that understands AND,
OR, and parenthetical nesting. The basic query unit is

 (match1 match2 ...)[fieldname]

Sequences are returned for which C<fieldname> equals C<match1 OR match2 OR ...>.
These units can be combined with AND, OR and parentheses. For example:

 (B, C)[subtype] AND (2000, 2001, 2002, 2003)[year] AND ((CN)[country] OR (ZA)[country])

which can be shortened to

 (B C)[subtype] (2000 2001 2002 2003)[year] (CN ZA)[country]

The user can specify annotation fields, that do not restrict the query, but
arrange for the return of the associated field data for each sequence returned.
Specify annotation fields between curly braces, as in:

 (B C)[subtype] 2000[year] {country cd4_count cd8_count}

Annotations can be accessed off the query using methods described in APPENDIX.

=head2 Hash specifications for query construction

Single query specifications can be made as hash references provided to the
C<-query> argument of the constructor. There are two forms:

 -query => { 'country'=>'BR', 'phenotype'=>'NSI', 'cd4_count'=>'Any' }

equivalent to

 -query => [ 'country'=>'BR', 'phenotype'=>'NSI', 'cd4_count'=>'Any' ]

or

 -query => { 'query' => {'country'=>'BR', 'phenotype'=>'NSI'},
             'annot' => ['cd4_count'] }

In both cases, the CD4 count is included in the annotations returned, but does
not restrict the rest of the query.

To 'OR' multiple values of a field, use an anonymous array ref:

 -query => { 'country'=>['ZA','BR','NL'], 'subtype'=>['A', 'C', 'D'] }

=head2 Valid query field names

An attempt was made to make the query field names natural and easy to
remember. Aliases are specified in an XML file (C<lanl-schema.xml>) that is part
of the distribution. Custom field aliases can be set up by modifying this file.

An HTML cheatsheet with valid field names, aliases, and match data can be
generated from the XML by using C<hiv_object-E<gt>help('help.html')>. A query
can also be validated locally before it is unleashed on the server; see below.

=head2 Annotations

LANL DB annotations have been organized into a number of natural
groupings, tagged C<Geo>, C<Patient>, C<Virus>, and C<StdMap>.  After a
successful query, each id is associated with a tree of
L<Bio::Annotation::SimpleValue> objects. These can be accessed with
methods C<get_value> and C<put_value> described in APPENDIX.

=head2 Delayed/partial query runs

Accessing the LANL DB involves multiple HTTP requests. The query can
be instructed to proceed through all (the default) or only some of
them, using the named parameter C<RUN_OPTION>.

To validate a query locally, use 

 $q = new Bio::DB::Query::HIVQuery( -query => {...}, -RUN_OPTION=>0 )

which will throw an exception if a field name or option is invalid. 

To get a query count only, you can save a server hit by using

 $q = new Bio::DB::Query::HIVQuery( -query => {...}, -RUN_OPTION=>1 )

and asking for C<$q-E<gt>count>. To finish the query, do

 $q->_do_query(2)

which picks up where you left off. 

C<-RUN_OPTION=E<gt>2>, the default, runs the full query, returning ids and
annotations.

=head2 Query re-use

You can clear the query results, retaining the same LANL session and query spec,
by doing C<$q-E<gt>_reset>. Change the query, and rerun with
C<$q-E<gt>_do_query($YOUR_RUN_OPTION)>. 

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

=head1 AUTHOR - Mark A. Jensen

Email maj@fortinbras.us

=head1 CONTRIBUTORS

Mark A. Jensen

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::DB::Query::HIVQuery;
use strict;
use vars qw( $LANL_BASE $LANL_MAP_DB $LANL_MAKE_SEARCH_IF $LANL_SEARCH $SCHEMA_FILE $RUN_OPTION );

# Object preamble - inherits from Bio::DB::QueryI
use Bio::Root::Root;
use Bio::Annotation::Collection;
use Bio::Annotation::Comment;
use Bio::Annotation::Reference;
use Bio::WebAgent;
use XML::Simple;
use CGI;

use Bio::DB::HIV::HIVQueryHelper;

use base qw(Bio::Root::Root Bio::DB::QueryI);

# globals
BEGIN {
    # change base to new search page 01/14/09 /maj
    $LANL_BASE = "http://www.hiv.lanl.gov/components/sequence/HIV/asearch";
    $LANL_MAP_DB = "map_db.comp";
    $LANL_MAKE_SEARCH_IF = "make_search_if.comp";
    $LANL_SEARCH = "search.comp";
    $SCHEMA_FILE = Bio::Root::IO->catfile(qw(Bio DB HIV lanl-schema.xml));
    $RUN_OPTION = 2; # execute query
# exceptions
    @Bio::SchemaNotInit::Exception::ISA = qw( Bio::Root::Exception );
    @Bio::WebError::Exception::ISA = qw( Bio::Root::Exception );
    @Bio::QueryNotMade::Exception::ISA = qw( Bio::Root::Exception );
    @Bio::QueryStringException::Exception::ISA = qw( Bio::Root::Exception );
    @Bio::HIVSorry::Exception::ISA = qw ( Bio::Root::Exception );

}

=head1 Constructor

=head2 new

 Title   : new
 Usage   : my $hiv_query = new Bio::DB::Query::HIVQuery();
 Function: Builds a new Bio::DB::Query::HIVQuery object,
           running a sequence query against the Los Alamos
           HIV sequence database
 Returns : an instance of Bio::DB::Query::HIVQuery
 Args    :

=cut

sub new {
  my($class,@args) = @_;
  my $self = $class->SUPER::new(@args);
  # constructor option for web agent parameter spec: added 01/14/09 /maj
  my ($query, $ids, $lanl_base, $lanl_map_db, $lanl_make_search_if, $lanl_search, $schema_file,$run_option, $uahash) =
      $self->_rearrange([ qw(QUERY
                             IDS
                             LANL_BASE
                             LANL_MAP_DB
                             LANL_MAKE_SEARCH_IF
                             LANL_SEARCH
                             SCHEMA_FILE
                             RUN_OPTION
                             USER_AGENT_HASH
                            )], @args);

  # default globals
  $lanl_base||= $LANL_BASE;
  $lanl_map_db||=$LANL_MAP_DB;
  $lanl_make_search_if||=$LANL_MAKE_SEARCH_IF;
  $lanl_search||=$LANL_SEARCH;
  $schema_file||=$SCHEMA_FILE;
  $uahash ||= {timeout => 90};
  defined $run_option                || ($run_option = $RUN_OPTION);

  $self->lanl_base($lanl_base);
  $self->map_db($lanl_map_db);
  $self->make_search_if($lanl_make_search_if);
  $self->search_($lanl_search);
  $self->_run_option($run_option);
  $self->_ua_hash($uahash);
  
  # catch this at the top
  if (-e $schema_file) {
      $self->_schema_file($schema_file);
  }
  else { # look around
      my ($p) = $self->_schema_file( [grep {$_} map {
	  my $p = Bio::Root::IO->catfile($_, $schema_file);
	  $p if -e $p
				      } (@INC,"")]->[0]);
      $self->throw(-class=>"Bio::Root::NoSuchThing",
		   -text=>"Schema file \"".$self->_schema_file."\" cannot be found",
		   -value=>$self->_schema_file) unless -e $self->_schema_file;
  }

  $self->count(0);
  $self->{_schema} =  HIVSchema->new($self->_schema_file);

  # internal storage and flags
  $self->{'_lanl_query'} = [];
  $self->{'_lanl_response'} = []; 
  $self->{'_annotations'} = {}; # container for annotation collections assoc. with ids
  $self->{'_RUN_LEVEL'} = undef; # set in _do_query()
  
  # work
  defined $query               && $self->query($query);
  defined $ids                 && $self->ids($ids);
  
  # exec query

  $self->_do_query($self->_run_option) if $self->query;

  return $self;
}

=head1 QueryI compliance

=head2 count

 Title   : count
 Usage   : $hiv_query->count($newval)
 Function: return number of sequences found
 Example : 
 Returns : value of count (a scalar)
 Args    : on set, new value (a scalar or undef, optional)
 Note    : count warns if it is accessed for reading before query
           has been executed to at least level 1

=cut

sub count{
    my $self = shift;
    return $self->{'count'} = shift if @_;
    if (!$self->{'_RUN_LEVEL'} || ($self->{'_RUN_LEVEL'} < 1)) {
	$self->warn('Query not yet run at > level 1');
    }
    return $self->{'count'};
}

=head2 ids

 Title   : ids
 Usage   : $hiv_query->ids($newval)
 Function: LANL ids of returned sequences 
 Example : 
 Returns : value of ids (an arrayref of sequence accessions/ids)
 Args    : on set, new value (an arrayref or undef, optional)

=cut

sub ids{
    my $self = shift;
    if (@_) {
	my $a = shift;
	$self->throw(-class=>'Bio::Root::BadParameter',
		     -text=>'Arrayref required',
		     -value=> ref $a) unless ref($a) eq 'ARRAY';
	if (@$a) {
	    @{$self->{'ids'}}{@$a} = (1) x @$a;
	    return $a;
	}
	else { #with empty arrayref, clear the hash
	    $self->{'ids'} = {};
	}
    }
    return keys %{$self->{'ids'}} if $self->{'ids'};
}

=head2 query

 Title   : query
 Usage   : $hiv_query->query
 Function: Get/set the submitted query hash or string
 Example :
 Returns : hashref or string
 Args    : query in hash or string form (see DESCRIPTION)

=cut

sub query {
    my $self = shift;
    return $self->{'query'} = shift if @_;
    return $self->{'query'};
}

=head1 Bio::DB::Query::HIVQuery specific methods

=head2 help

 Title   : help
 Usage   : $hiv_query->help("help.html")
 Function: get html-formatted listing of valid fields/aliases/options
           based on current schema xml
 Example : perl -MBio::DB::Query::HIVQuery -e "new Bio::DB::Query::HIVQuery()->help" | lynx -stdin
 Returns : HTML
 Args    : optional filename; otherwise prints to stdout

=cut

sub help{
   my ($self, $fname) = @_;
   my (@ret, @tok);
   my $schema = $self->_schema;
   my $h = CGI->new();

   my (@tbls, @flds, @als, @opts, $fh);
   if ($fname) {
       open $fh, '>', $fname or $self->throw(-class => 'Bio::Root::IOException',
                                             -text  => "Error opening help html file $fname for writing",
                                             -value => $!);
   }
   else {
       open $fh, ">&1";
   }
   @tbls = $schema->tables;
   @tbls = ('COMMAND', grep !/COMMAND/,@tbls);
   print $fh (
       $h->start_html(-title=>"HIVQuery Help")
       );
   print $fh $h->a({-id=>'TOP'}, $h->h2("Valid <span style='font-variant:small-caps'>HIVQuery</span> query fields and match data"));
   print $fh "Fields are organized below according to their Los Alamos HIV database tables. Use aliases in place of full field names in queries; for example:<br/>";
   print $fh "<blockquote><code> (CCR5 CXCR4)[coreceptor]</code></blockquote>";
   print $fh "rather than";
   print $fh "<blockquote><code>(CCR5 CXCR4)[seq_sample.ssam_second_receptor] </code></blockquote>";
   print $fh "(which does work, however). Click hyperlinks to see valid search options within the field. The token <code><b>Any</b></code> is the wildcard for all fields.<br/><br/>";
   print $fh $h->start_table({-style=>"font-family:sans-serif;"}) ;
   foreach my $tbl (@tbls) {
       @flds = grep /^$tbl/, $schema->fields;
       @flds = grep !/_id/, @flds;
       print $fh (
           $h->start_Tr({-style=>"background-color: lightblue;"}), 
           $h->td([$h->a({-id=>$tbl},$tbl), $h->span({-style=>"font-style:italic"},"fields"), $h->span({-style=>"font-style:italic"}, "aliases")]),
           $h->end_Tr
       );
       foreach my $fld (@flds) {
           @als = reverse $schema->aliases($fld);
           print $fh (
               # note that aliases can sometimes be empty
               $h->Tr( $h->td( ["", $h->a({-href=>"#opt$fld"}, shift @als || '???'), $h->code(join(',',@als))] ))
           );
           my @tmp = grep {$_} $schema->options($fld);
           @tmp = sort {(($a =~ /^[0-9]+$/) && $b =~ /^[0-9]+$/) ? $a<=>$b : $a cmp $b} @tmp;
           if (grep /Any/,@tmp) {
               @tmp = grep !/Any/, @tmp;
               unshift @tmp, 'Any';
           }
           #print STDERR join(', ',@tmp)."\n";
           push @opts, $h->div(
               {-style=>"font-family:sans-serif;font-size:small"},
               $h->hr,
               $h->a(
                   {-id=>"opt$fld"},
                   "<i>Valid options for</i> <b>$fld</b>: "
               ),
               $h->blockquote(
                   @tmp ? $h->code(join(", ", @tmp)) : $h->i("free text")
               ),
               $h->span(
                   "<i>Other aliases</i>: "
               ),
               $h->blockquote(
                   @als ? $h->code(join(",",@als)) : "<i>none</i>"
               ),
               " ", 
               $h->table(
                   $h->Tr(
                       $h->td([
                           $h->a({-href=>"#$tbl"}, $h->small('BACK')), 
                           $h->a({-href=>"#TOP"}, $h->small('TOP'))
                       ])
                   )
               )
           );
   
       }
   }
   print $fh $h->end_table;
   print $fh @opts;
   print $fh $h->end_html;
   close($fh);
   return 1;
}

=head1 Annotation manipulation methods    

=head2 get_annotations_by_ids

 Title   : get_annotations_by_ids (or ..._by_id)
 Usage   : $ac = $hiv_query->get_annotations_by_ids(@ids)
 Function: Get the Bio::Annotation::Collection for these sequence ids
 Example :
 Returns : A Bio::Annotation::Collection object
 Args    : an array of sequence ids

=cut

sub get_annotations_by_ids{
    my $self = shift;
    my @ids = @_;
    my @ret;
    if (!$self->{'_RUN_LEVEL'} || ($self->{'_RUN_LEVEL'} < 2)) {
	$self->warn('Requires query run at level 2');
	return ();
    }
    @ret = map {$self->{'_annotations'}->{$_}} @ids if exists($self->{'_annotations'});

    return (wantarray ? @ret : $ret[0]) if @ret;
    return {};
}

# singular alias
sub get_annotations_by_id {
    shift->get_annotations_by_ids(@_);
}

=head2 add_annotations_for_id

 Title   : add_annotations_for_id
 Usage   : $hiv_query->add_annotations_for_id( $id ) to create a new 
            empty collection for $id
           $hiv_query->add_annotations_for_id( $id, $ac ) to associate 
           $ac with $id
 Function: Associate a Bio::Annotation::Collection with this sequence id
 Example :
 Returns : a Bio::Annotation::Collection object
 Args    : sequence id [, Bio::Annotation::Collection object]

=cut

sub add_annotations_for_id{
    my $self = shift;
    my ($id, $ac) = @_;
    $id = "" unless defined $id; # avoid warnings
    $ac = Bio::Annotation::Collection->new() unless defined $ac;
    $self->throw(-class=>'Bio::Root::BadParameter'
		 -text=>'Bio::Annotation::Collection required at arg 2',
		 -value=>"") unless ref($ac) eq 'Bio::Annotation::Collection';
    
    $self->{'_annotations'}->{$id} = $ac unless exists($self->{'_annotations'}->{$id});
    return $ac;
}

=head2 remove_annotations_for_ids

 Title   : remove_annotations_for_ids (or ..._for_id)
 Usage   : $hiv_query->remove_annotations_for_ids( @ids)
 Function: Remove annotation collection for this sequence id
 Example :
 Returns : An array of the previous annotation collections for these ids
 Args    : an array of sequence ids

=cut

sub remove_annotations_for_ids {
    my $self = shift;
    my @ids = @_;
    my @ac;
    foreach (@ids) {
	push @ac, delete $self->{'_annotations'}->{$_};
    }
    return @ac;
}

# singular alias
sub remove_annotations_for_id {
    shift->remove_annotations_for_ids(@_);
}

=head2 remove_annotations

 Title   : remove_annotations
 Usage   : $hiv_query->remove_annotations()
 Function: Remove all annotation collections for this object
 Example :
 Returns : The previous annotation collection hash for this object
 Args    : none

=cut

sub remove_annotations {
    my $self = shift;

    my $ach = $self->{'_annotations'};
    $self->{'_annotations'} = {};
    return $ach;
}

=head2 get_value

 Title   : get_value
 Usage   : $ac->get_value($tagname) -or-
           $ac->get_value( $tag_level1, $tag_level2,... )
 Function: access the annotation value assocated with the given tags
 Example :
 Returns : a scalar
 Args    : an array of tagnames that descend into the annotation tree
 Note    : this is a L<Bio::AnnotationCollectionI> method added in 
           L<Bio::DB::HIV::HIVQueryHelper>

=cut

=head2 put_value

 Title   : put_value
 Usage   : $ac->put_value($tagname, $value) -or-
           $ac->put_value([$tag_level1, $tag_level2, ...], $value) -or-
           $ac->put_value( [$tag_level1, $tag_level2, ...] )
 Function: create a node in an annotation tree, and assign a scalar value to it
           if a value is specified
 Example :
 Returns : scalar or a Bio::AnnotationCollection object
 Args    : $tagname, $value scalars (can be specified as -KEYS=>$tagname,
           -VALUE=>$value) -or- 
           \@tagnames, $value (or as -KEYS=>\@tagnames, -VALUE=>$value )
 Notes   : This is a L<Bio::AnnotationCollectionI> method added in 
           L<Bio::DB::HIV::HIVQueryHelper>.
           If intervening nodes do not exist, put_value creates them, replacing 
           existing nodes. So if $ac->put_value('x', 10) was done, then later,
           $ac->put_value(['x', 'y'], 20), the original value of 'x' is trashed,
           and $ac->get_value('x') will now return the annotation collection 
           with tagname 'y'. 

=cut

=head2 get_keys

 Title   : get_keys
 Usage   : $ac->get_keys($tagname_level_1, $tagname_level_2,...)
 Function: Get an array of tagnames underneath the named tag nodes
 Example : # prints the values of the members of Category 1...
           print map { $ac->get_value($_) } $ac->get_keys('Category 1') ;
 Returns : array of tagnames or empty list if the arguments represent a leaf
 Args    : [array of] tagname[s]

=cut

=head1 GenBank accession manipulation methods

=head2 get_accessions

 Title   : get_accessions
 Usage   : $hiv_query->get_accessions()
 Function: Return an array of GenBank accessions associated with these 
           sequences (available only after a query is subjected to a 
           full run (i.e., when $RUN_OPTION == 2)
 Example :
 Returns : array of gb accession numbers, or () if none found for this query
 Args    : none

=cut

sub get_accessions{
    my $self = shift; 
    my @ret;
    if (!$self->{'_RUN_LEVEL'} || ($self->{'_RUN_LEVEL'} < 2)) {
	$self->warn('Requires query run at level 2');
	return ();
    }
    my @ac = $self->get_annotations_by_ids($self->ids);
    foreach (@ac) {
	push @ret, $_->get_value('Special','accession');
    };
    return @ret;
}

=head2 get_accessions_by_ids

 Title   : get_accessions_by_ids (or ..._by_id)
 Usage   : $hiv_query->get_accessions_by_ids(@ids)
 Function: Return an array of GenBank accessions associated with these 
           LANL ids (available only after a query is subjected to a 
           full run (i.e., when $RUN_OPTION == 2)
 Example :
 Returns : array of gb accession numbers, or () if none found for this query
 Args    : none

=cut

sub get_accessions_by_ids {
    my $self = shift; 
    my @ids = @_;
    my @ret;
    if (!$self->{'_RUN_LEVEL'} || ($self->{'_RUN_LEVEL'} < 2)) {
	$self->warn('Requires query run at level 2');
	return ();
    }
    my @ac = $self->get_annotations_by_ids(@ids);
    foreach (@ac) {
	push @ret, $_->get_value('Special', 'accession');
    };
    return wantarray ? @ret : $ret[0];
}

# singular alias
sub get_accessions_by_id {
    shift->get_accessions_by_ids(@_);
}

##########    

=head1 Query control methods

=head2 _do_query

 Title   : _do_query
 Usage   : $hiv_query->_do_query or $hiv_query->_do_query($run_level)
 Function: Execute the query according to argument or $RUN_OPTION
           and set _RUN_LEVEL
           extent of query reflects the value of argument
            0 : validate only (no HTTP action)
            1 : return sequence count only
            2 : return sequence ids (full query, returns with annotations)
           noop if current _RUN_LEVEL of query is >= argument or $RUN_OPTION,
 Example :
 Returns : actual _RUN_LEVEL (0, 1, or 2) achieved
 Args    : desired run level (optional, global $RUN_OPTION is default)

=cut

sub _do_query{
   my ($self,$rl) = @_;
   $rl = $RUN_OPTION unless defined $rl;
   $self->throw(-class=>"Bio::Root::BadParameter",
		-text=>"Invalid run option \"$RUN_OPTION\"",
		-value=>$RUN_OPTION) unless grep /^$RUN_OPTION$/, (0, 1, 2);
   (!defined($self->{'_RUN_LEVEL'})) && do {
       $self->_create_lanl_query();
       $self->{'_RUN_LEVEL'} = 0;
   };
   ($rl > 0) && (!defined($self->{'_RUN_LEVEL'}) || ($self->{'_RUN_LEVEL'} <= 0)) && do {
       $self->_do_lanl_request();
       $self->{'_RUN_LEVEL'} = 1;
   };
   ($rl > 1) && (!defined($self->{'_RUN_LEVEL'}) || ($self->{'_RUN_LEVEL'} <= 1)) && do {
       $self->_parse_lanl_response();
       $self->{'_RUN_LEVEL'} = 2;
   };
   return $self->{'_RUN_LEVEL'};
}

=head2 _reset

 Title   : _reset
 Usage   : $hiv_query->_reset
 Function: Resets query storage, count, and ids, while retaining session id, 
           original query string, and db schema
 Example : 
 Returns : void
 Args    : none

=cut

sub _reset{
    my $self = shift;
    $self->ids([]);
    $self->count(0);
    $self->{'_annotations'} = {};
    $self->{'_lanl_response'} = [];
    $self->{'_lanl_query'} = [];
    $self->{'_RUN_LEVEL'} = undef;
    return;
}

=head2 _session_id

 Title   : _session_id
 Usage   : $hiv_query->_session_id($newval)
 Function: Get/set HIV db session id (initialized in _do_lanl_request)
 Example : 
 Returns : value of _session_id (a scalar)
 Args    : on set, new value (a scalar or undef, optional)

=cut

sub _session_id{
    my $self = shift;

    return $self->{'_session_id'} = shift if @_;
    return $self->{'_session_id'};
}
=head2 _run_level

 Title   : _run_level
 Usage   : $obj->_run_level($newval)
 Function: returns the level at which the query has so far been run
 Example : 
 Returns : value of _run_level (a scalar)
 Args    : on set, new value (a scalar or undef, optional)

=cut

sub _run_level{
    my $self = shift;

    return $self->{'_RUN_LEVEL'} = shift if @_;
    return $self->{'_RUN_LEVEL'};
}

=head2 _run_option

 Title   : _run_option
 Usage   : $hiv_query->_run_option($newval)
 Function: Get/set HIV db query run option (see _do_query for values)
 Example : 
 Returns : value of _run_option (a scalar)
 Args    : on set, new value (a scalar or undef, optional)

=cut

sub _run_option{
    my $self = shift;

    return $self->{'_run_option'} = shift if @_;
    return $self->{'_run_option'};
}

=head2 _ua_hash

 Title   : _ua_hash
 Usage   : $obj->_ua_hash($newval)
 Function: 
 Example : 
 Returns : value of _ua_hash (a scalar)
 Args    : on set, new value (a scalar or undef, optional)

=cut

sub _ua_hash{
    my $self = shift;
    if (@_) {
	for (ref $_[0]) {
	    $_ eq 'HASH' && do {
		$self->{'_ua_hash'} = $_[0];
		last;
	    };
	    !$_ && do {
		$self->{'_ua_hash'} = {@_};
		last;
	    };
	    do {
		$self->throw("Type ".ref($_)." unsupported as arg in _ua_hash");
	    };
	    
	}
    }
    return %{$self->{'_ua_hash'}};
}


#######

=head1 Internals

=head2 add_id

 Title   : add_id
 Usage   : $hiv_query->add_id($id)
 Function: Add new id to ids
 Example : 
 Returns : the new id
 Args    : a sequence id

=cut

sub add_id {
    my $self = shift;
    my $id = shift;
    $id = "" unless defined $id; # avoid warnings
    ${$self->{'ids'}}{$id}++;
    return $id;
}


sub lanl_base{
    my $self = shift;
    return $self->{'lanl_base'} = shift if @_;
    return $self->{'lanl_base'};
}

=head2 map_db

 Title   : map_db
 Usage   : $obj->map_db($newval)
 Function: 
 Example : 
 Returns : value of map_db (a scalar)
 Args    : on set, new value (a scalar or undef, optional)

=cut

sub map_db{
    my $self = shift;
    return $self->{'map_db'} = shift if @_;
    return $self->{'map_db'};
}

=head2 make_search_if

 Title   : make_search_if
 Usage   : $obj->make_search_if($newval)
 Function: 
 Example : 
 Returns : value of make_search_if (a scalar)
 Args    : on set, new value (a scalar or undef, optional)

=cut

sub make_search_if{
    my $self = shift;
    return $self->{'make_search_if'} = shift if @_;
    return $self->{'make_search_if'};
}

=head2 search_

 Title   : search_
 Usage   : $obj->search_($newval)
 Function: 
 Example : 
 Returns : value of search_ (a scalar)
 Args    : on set, new value (a scalar or undef, optional)

=cut

sub search_{
    my $self = shift;
    return $self->{'search_'} = shift if @_;
    return $self->{'search_'};
}

=head2 _map_db_uri

 Title   : _map_db_uri
 Usage   :
 Function: return the full map_db uri ("Database Map")
 Example :
 Returns : scalar string
 Args    : none

=cut

sub _map_db_uri{
    my $self = shift;
    return $self->lanl_base."/".$self->map_db;
}
  

=head2 _make_search_if_uri

 Title   : _make_search_if_uri
 Usage   :
 Function: return the full make_search_if uri ("Make Search Interface")
 Example :
 Returns : scalar string
 Args    : none

=cut

sub _make_search_if_uri{
    my $self = shift;
    return $self->lanl_base."/".$self->make_search_if;
}

=head2 _search_uri

 Title   : _search_uri
 Usage   :
 Function: return the full search cgi uri ("Search Database")
 Example :
 Returns : scalar string
 Args    : none

=cut

sub _search_uri{
    my $self = shift;
    return $self->lanl_base."/".$self->search_;
}

=head2 _schema_file

 Title   : _schema_file
 Usage   : $hiv_query->_schema_file($newval)
 Function: 
 Example : 
 Returns : value of _schema_file (an XML string or filename)
 Args    : on set, new value (an XML string or filename, or undef, optional)

=cut

sub _schema_file {
    my $self = shift;

    return $self->{'_schema_file'} = shift if @_;
    return $self->{'_schema_file'};
}

=head2 _schema

 Title   : _schema
 Usage   : $hiv_query->_schema($newVal)
 Function: 
 Example : 
 Returns : value of _schema (an HIVSchema object in package 
           L<Bio::DB::HIV::HIVQueryHelper>)
 Args    : none (field set directly in new())

=cut

sub _schema{
    my $self = shift;
    
    $self->{'_schema'} ?
	return $self->{'_schema'} :
	$self->throw(-class=>'Bio::SchemaNotInit::Exception', 
		     -text=>"DB schema not initialized",
		     -value=>"");
	
}

=head2 _lanl_query

 Title   : _lanl_query
 Usage   : $hiv_query->_lanl_query(\@query_parms)
 Function: pushes \@query_parms onto @{$self->{'_lanl_query'}
 Example : 
 Returns : value of _lanl_query (an arrayref)
 Args    : on set, new value (an arrayref or undef, optional)

=cut

sub _lanl_query{
    my $self = shift;
    my $a = shift;
    return $self->{'_lanl_query'} unless $a;
    if (ref $a eq 'ARRAY') {
	push @{$self->{'_lanl_query'}}, $a;
	return $a;
    }
    else {
	$self->throw(-class=>'Bio::Root::BadParameter',
		     -text=>'Array ref required for argument.',
		     -value=>$a);
    }

}

=head2 _lanl_response

 Title   : _lanl_response
 Usage   : $hiv_query->_lanl_response($response)
 Function: pushes $response onto @{$hiv_query->{'_lanl_response'}}
 Example : 
 Returns : value of _lanl_response (an arrayref of HTTP::Response objects)
 Args    : on set, new value (an HTTP::Response object or undef, optional)

=cut

sub _lanl_response{
    my $self = shift;
    if (@_) {
	my $r = shift;
	$self->throw(-class=>'Bio::Root::BadParameter',
		     -text=>'Requires an HTTP::Response object',
		     -value=> ref $r) unless ref($r) eq 'HTTP::Response';
	push @{$self->{'_lanl_response'}}, $r;
	return $r;
    }
    return $self->{'_lanl_response'};
}

=head2 _create_lanl_query

 Title   : _create_lanl_query
 Usage   : $hiv_query->_create_lanl_query()
 Function: validate query hash or string, prepare for _do_lanl_request
 Example : 
 Returns : 1 if successful; throws exception on invalid query
 Args    :

=cut

sub _create_lanl_query {
    my $self = shift;
    my (%inhash, @query, @qhashes);
    my ($schema, @validFields, @validAliases);
    
    for ($self->query) {
	!defined && do {
	    $self->throw(-class=>'Bio::Root::NoSuchThing',
			 -text=>'Query not specified',
			 -value=>'');
	    last;
	};
	ref eq 'HASH' && do {
	    %inhash = %$_;
	    if ( grep /HASH/, map {ref} values %inhash ) {
		# check for {query=>{},annot=>[]} style
		$self->throw(-class=>'Bio::Root::BadParameter',
			     -text=>'Query style unrecognized',
			     -value=>"") unless defined $inhash{query};
		push @qhashes, $_;
	    }
	    last;
	};
	ref eq 'ARRAY' && do {
	    $inhash{'query'} = {@$_};
	    push @qhashes, \%inhash;
	    last;
	};
	#else
	do {
	    @qhashes = $self->_parse_query_string($_);
	};
    }
    $schema = $self->_schema;
    @validFields = $schema->fields;
    @validAliases = $schema->aliases;

    # validate args based on the xml specification file
    # only checks blanks and fields with explicitly specified options
    # text fields can put anything, and the query will be run before 
    # an error is caught in these 
    foreach my $qh (@qhashes) {
	@query=();
	foreach my $k (keys %{$$qh{'query'}}) {
	    my $fld;
	    # validate field
	    if (grep /^$k$/, @validFields) {
		$fld = $k;
	    }
	    elsif (grep /^$k$/, @validAliases) {
		foreach (@validFields) {
		    if (grep (/^$k$/, $schema->aliases($_))) {
			$fld = $_;
			last;
		    }
		    # $fld contains the field corresp. to the alias
		}
	    }
	    else {
		$self->throw(-class=>'Bio::Root::BadParameter',
			     -text=>"Invalid field or alias \"$k\"",
			     -value=>$qh);
	    }
	    # validate matchdata
	    my $vf = $schema->_sfieldh($fld);
	    my @md = (ref($qh->{'query'}{$k}) eq 'ARRAY') ? @{$qh->{'query'}{$k}} : $qh->{'query'}{$k};
	    if ($$vf{type} eq 'text') {
		foreach (@md) {
		    $self->throw(-class=>'Bio::Root::BadParameter',
				 -text=>'Value for field \"$k\" cannot be empty',
				 -value=>$qh)
			if ($_ eq "") && ($$vf{blank_ok} eq 'false');
		}
	    }
	    elsif ($$vf{type} eq 'option') {
		foreach my $md (@md) {
		    $self->throw(-class=>'Bio::Root::BadParameter',
				 -text=>"Invalid value \"".$md."\" for field \"$fld\"",
				 -value=>$md)
			unless $$vf{option} && grep {defined $_ && /^$md$/} @{$$vf{option}};
		}
	    }
	    # validated; add to query
	    foreach (@md) {
		push @query, ($fld => $_);
	    }
	}
	if ($qh->{'annot'}) {
	    # validate the column names to be included in the query
	    # to obtain annotations
	    my @annot_cols = @{$qh->{'annot'}};
	    foreach my $k (@annot_cols) {
		my $fld;
		# validate field
		if (grep /^$k$/, @validFields) {
		    $fld = $k;
		}
		elsif (grep /^$k$/, @validAliases) {
		    foreach (@validFields) {
			if (grep (/^$k$/, $schema->aliases($_))) {
			    $fld = $_;
			    last;
			}
			# $fld should contain the field corresp. to the alias
		    }
		}
		else {
		    $self->throw(-class=>'Bio::Root::NoSuchThing',
				 -text=>"Invalid field or alias \"$k\"",
				 -value=>$k);
		}
		# lazy: 'Any' may not be the right default (but appears to 
		# be, based on the lanl html)
		push @query, ($fld => 'Any');
	    }
	}

	# insure that LANL and GenBank ids are retrieved
	push @query, ('sequenceentry.se_id' => 'Any') unless grep /SequenceEntry\.SE_id/, @query;
	push @query, ('sequenceaccessions.sa_genbankaccession' => 'Any')
	    unless grep /SequenceAccessions\.SA_GenBankAccession/, @query;

	# an "order" field is required by the LANL CGI
	# if not specified, default to SE_id

	push @query, ('order'=>'sequenceentry.se_id') unless grep /order/, @query;

	# @query now contains sfield=>matchdata pairs, as specified by user
	# include appropriate indexes to create correct automatic joins
	# established by the LANL CGI
	my (@qtbl, @qpk, @qfk);

	# the tables represented in query:
	my %q = @query; # squish the tables in the current query into hash keys
	@qtbl = $schema->tbl('-s', keys %q);

	if (@qtbl > 1) {
	    # more than one table, see if they can be connected
	    # get primary keys of query tables
	    @qpk = $schema->pk(@qtbl);

	    # we need to get each query table to join to 
	    # SequenceEntry.
	    #
	    # The schema is a graph with tables as nodes and 
	    # foreign keys<->primary keys as branches. To get a 
	    # join that works, need to include in the query
	    # all branches along a path from SequenceEntry 
	    # to each query table.
	    #
	    # find_join does it...
	    my @joink = map {
		my @k = $schema->find_join($_,'sequenceentry');
		map {$_ || ()} @k
	    } @qtbl;
	    # squish the keys in @joink 
	    my %j;
	    @j{@joink} = (1) x @joink;
	    @joink = keys %j;
	    # add the fields not currently in the query
	    foreach (@qpk, @joink) {
		my $fld = $_;
		if (!grep(/^$fld$/,keys %q)) {
		    # lazy: 'Any' may not be the right default (but appears to 
		    # be, based on the lanl html)
		    push @query, ($_ => 'Any');
		}
	    }

	}
	
	# set object property
	$self->_lanl_query([@query]);
    }
    return 1;
}

# _do_lanl_request : post the queries created by _create_lanl_query
# 
# @args (or {@args}) should be unaliased Table.Column=>Matchdata
# pairs (these will be used directly in the POSTs)

=head2 _do_lanl_request

 Title   : _do_lanl_request
 Usage   : $hiv_query->_do_lanl_request()
 Function: Perform search request on _create_lanl_query-validated query
 Example : 
 Returns : 1 if successful
 Args    : 

=cut

sub _do_lanl_request {
    my $self = shift;
    my (@queries, @query, @interface,$interfGet,$searchGet,$response);
    my ($numseqs, $count);

    # handle args
    if (!$self->_lanl_query) {
	$self->throw(-class=>"Bio::Root::BadParameter",
		     -text=>"_lanl_query empty, run _create_lanl_request first",
		     -value=>"");
    }
    else {
	@queries = @{$self->_lanl_query};
    }
    
    ## utility vars
    ## search site specific CGI parms
    my @search_pms = ('action'=>'Search');
    my @searchif_pms = ('action'=>'Search Interface');
    # don't get the actual sequence data here (i.e., the cgi parm 
    # 'incl_seq' remains undefined...
    my @download_pms = ('action Download.x'=>1, 'action Download.y'=>1);

    ## HTML-testing regexps
    my $tags_re = qr{(?:\s*<[^>]+>\s*)};
    my $session_id_re = qr{<input.*name="id".*value="([0-9a-f]+)"}m;
    my $search_form_re = qr{<form[^>]*action=".*/search.comp"};
    my $seqs_found_re = qr{Displaying$tags_re*(?:\s*[0-9-]*\s*)*$tags_re*of$tags_re*\s*([0-9]+)$tags_re*sequences found};
    my $no_seqs_found_re = qr{Sorry.*no sequences found};
    my $too_many_re = qr{too many records: $tags_re*([0-9]+)};
    my $sys_error_re = qr{[Ss]ystem error};
    my $sys_error_extract_re = qr{${tags_re}error:.*?<td[^>]+>${tags_re}(.*?)<br>};
    # find something like:
    #  <strong>tables without join:</strong><br>SequenceAccessions<br>
    my $tbl_no_join_re = qr{tables without join}i;
#    my $sorry_bud_re = qr{};


    foreach my $q (@queries) {
	@query = @$q;
	# default query control parameters
	my %qctrl = (
	    max_rec=>100,
	    sort_dir=>'ASC',
	    translate=>'FALSE' # nucleotides 
	    );
	
	# do work...

	# pull out commands, designated by the COMMAND pseudo-table...
	my @commands = map { $query[$_] =~ s/^COMMAND\.// ? @query[$_..$_+1] : () } (0..$#query-1);
	@query = map { $query[$_] =~ /^COMMAND/ ? () : @query[2*$_..2*$_+1] } (0..($#query-1)/2);

	
	# set control parameters explicitly made in query
	foreach my $cp (keys %qctrl) {
	    if (!grep( /^$cp$/, @query)) {
		push @query, ($cp, $qctrl{$cp});
	    }
	}

	# note that @interface must be an array, since a single 'key' (the table)
	# can be associated with multiple 'values' (the columns) in the POST 

	# squish fieldnames into hash keys
	my %q = @query;
	@interface = grep {defined} map {my ($tbl,$col) = /^(.*)\.(.*)$/} keys %q;
	my $err_val = ""; # to contain informative (ha!) value if error is parsed

	eval { # encapsulate communication errors here, defer biothrows...
        
        #mark the useragent should be setable from outside (so we can modify timeouts, etc)
	    my $ua = Bio::WebAgent->new($self->_ua_hash);
	    my $idPing = $ua->get($self->_map_db_uri);
	    $idPing->is_success || do {
		$response=$idPing; 
		die "Connect failed";
	    };
	    # get the session id
	    if (!$self->_session_id) {
		($self->{'_session_id'}) = ($idPing->content =~ /$session_id_re/);
		$self->_session_id || do {
		    $response=$idPing; 
		    die "Session not established";
		};
	    }
	    # 10/07/08:
	    # strange bug: if action=>'Search+Interface' below (note "+"), 
	    # the response to the search (in $searchGet) shows the correct 
	    # >number< of sequences found, but also an error "No sequences 
	    # match" and an SQL barf. Changing the "+" to a " " sets up the 
	    # interface to lead to the actual sequences being delivered as 
	    # expected. maj
	    $interfGet = $ua->post($self->_make_search_if_uri, [@interface, @searchif_pms, id=>$self->_session_id]);
	    $interfGet->is_success || do {
		$response=$interfGet;
		die "Interface request failed";
	    };
	    # see if a search form was returned...
	    
	    $interfGet->content =~ /$search_form_re/ || do {
		$response=$interfGet; 
		die "Interface request failed";
	    };
	    
	    $searchGet = $ua->post($self->_search_uri, [@query, @commands, @search_pms, id=>$self->_session_id]);
	    $searchGet->is_success || do {
		$response = $searchGet;
		die "Search failed";
	    };
	    $response = $searchGet;
	    for ($searchGet->content) {
		/$no_seqs_found_re/ && do {
		    $err_val = 0;
		    die "No sequences found";
		    last;
		};
		/$too_many_re/ && do {
		    $err_val = $1;
		    die "Too many records ($1): must be <10000";
		    last;
		};
		/$tbl_no_join_re/ && do {
		    die "Some required tables went unjoined to query";
		    last;
		};
		/$sys_error_re/ && do {
		    /$sys_error_extract_re/;
		    $err_val = $1;
		    die "LANL system error";
		};
		/$seqs_found_re/ && do {
		    $numseqs = $1;
		    $count += $numseqs;
		    last;
		};
		# else...
		do {
		    die "Search failed (response not parsed)";
		};
	    }
	    $response = $ua->post($self->_search_uri, [@download_pms, id=>$self->_session_id]);
	    $response->is_success || die "Query failed";
	    # $response->content is a tab-separated value table of sequences 
	    # and metadata, first line starts with \# and contains fieldnames
	};
	$self->_lanl_response($response);
	# throw, if necessary
	if ($@) {
	    ($@ !~ "No sequences found") && do {
		$self->throw(-class=>'Bio::WebError::Exception',
			     -text=>$@,
			     -value=>$err_val);
	    };
	}
    }

    $self->warn("No sequences found for this query") unless $count;
    $self->count($count);
    return 1; # made it.

}

=head2 _parse_lanl_response

 Title   : _parse_lanl_response
 Usage   : $hiv_query->_parse_lanl_response()
 Function: Parse the tab-separated-value response obtained by _do_lanl_request
           for sequence ids, accessions, and annotations
 Example : 
 Returns : 1 if successful
 Args    : 

=cut

sub _parse_lanl_response {

### handle parsing and merging multiple responses into the query object
### (ids and annotations)
    my $self = shift;
    
    my ($seqGet) = (@_);
    my (@data, @cols, %antbl, %antype);
    my $numseq = 0;
    my ($schema, @retseqs, %rec, $ac);
    $schema = $self->_schema;
    
    $self->_lanl_response || 
	$self->throw(-class=>"Bio::QueryNotMade::Exception",
		     -text=>"Query not yet performed; call _do_lanl_request()",
		     -value=>"");
    foreach my $rsp (@{$self->_lanl_response}) {
	@data = split(/\r|\n/, $rsp->content);
	my $l;
	do {
	    $l = shift @data;
	} while ($l !~ /Number/);
	$numseq += ( $l =~ /Number.*:\s([0-9]+)/ )[0];
	@cols = split(/\t/, shift(@data));
	# mappings from column headings to annotation keys
	# squish into hash keys
	my %q = @{ shift @{$self->_lanl_query} };
	%antbl = $schema->ankh(keys %q);
	# get the category for each annotation
	map { $antype{ $_->{ankey} } = $_->{antype} } values %antbl;
	# normalize column headers
	map { tr/ /_/; $_ = lc; } @cols;
	foreach (@data) {
	    @rec{@cols} = split /\t/;
	    my $id = $rec{'se_id'};
	    $self->add_id($id);
	    $ac = Bio::Annotation::Collection->new();
	    #create annotations
	    foreach (@cols) {
                next if $_ eq '#';
		my $t = $antype{$_} || "Unclassified";
		my $d = $rec{$_}; # the data
                $ac->put_value(-KEYS=>[$t, $_], -VALUE=>$d);
	    }
	    $self->add_annotations_for_id($id, $ac);	
	}
	1;
    }
    return 1; # made it.
}
    
=head2 _parse_query_string

 Title   : _parse_query_string
 Usage   : $hiv_query->_parse_query_string($str)
 Function: Parses a query string using query language emulator QRY
         : in L<Bio::DB::Query::HIVQueryHelper>
 Example : 
 Returns : arrayref of hash structures suitable for passing to _create_lanl_query
 Args    : a string scalar

=cut

sub _parse_query_string {
    my $self = shift;
    my $qstring = shift;
    my ($ptree, @ret);
    #syntax errors thrown in QRY (in HIVQueryHelper module)
    $ptree = QRY::_parse_q( $qstring );
    @ret = QRY::_make_q($ptree);
    return @ret;
}

=head1 Dude, sorry-

=head2 _sorry

 Title   : _sorry
 Usage   : $hiv_query->_sorry("-president=>Powell")
 Function: Throws an exception for unsupported option or parameter
 Example :
 Returns : 
 Args    : scalar string

=cut

sub _sorry{
    my $self = shift;
    my $parm = shift;
    $self->throw(-class=>"Bio::HIVSorry::Exception",
		 -text=>"Sorry, option/parameter \"$parm\" not (yet) supported. See manpage to complain.",
		 -value=>$parm);
    return;
}

1;
