# $Id: HIVQueryHelper.pm 231 2008-12-11 14:32:00Z maj $
#
# BioPerl module for Bio::DB::HIV::HIVQueryHelper
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

Bio::DB::HIV::HIVQueryHelper - Routines and packages used by Bio::DB::HIV and
Bio::DB::Query::HIVQuery

=head1 SYNOPSIS

  Used in Bio::DB::Query::HIVQuery. No need to use directly.

=head1 DESCRIPTION

C<Bio::DB::HIV::HIVQueryHelper> contains a number of packages for use
by L<Bio::DB::Query::HIVQuery>. Package C<HIVSchema> parses the
C<lanl-schema.xml> file, and allows access to it in the context of the
relational database it represents (see APPENDIX for excruciating
detail). Packages C<QRY>, C<R>, and C<Q> together create the query
string parser that enables NCBI-like queries to be understood by
C<Bio::DB::Query::HIVQuery>. They provide objects and operators to
perform and simplify logical expressions involving C<AND>, C<OR>, and
C<()> and return hash structures that can be handled by
C<Bio::DB::Query::HIVQuery> routines.

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

The rest of the documentation details each of the contained packages.
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::DB::HIV::HIVQueryHelper;
use strict;
use Bio::Root::Root;

# globals
BEGIN {
#exceptions
    @Bio::QueryStringSyntax::Exception::ISA = qw( Bio::Root::Exception);
}

1;

=head2 HIVSchema -  objects/methods to manipulate a version of the LANL HIV DB schema

=head3 HIVSchema SYNOPSIS

    $schema = new HIVSchema( 'lanl-schema.xml' );
    @tables = $schema->tables;
    @validFields = $schema->fields;
    @validAliases = $schema->aliases;
    @query_aliases_for_coreceptor = $schema->aliases( 'SEQ_SAMple.SSAM_second_receptor' );
    $pk_for_SequenceEntry = $schema->primarykey('SequenceEntry');    # returns 'SequenceEntry.SE_id'
    $fk_for_SEQ_SAMple_to_SequenceEntry =
              $schema->foreignkey('SEQ_SAMple', 'SequenceEntry');    # returns 'SEQ_SAMple.SSAM_SE_id'

    $table = $schema->tablepart('SEQ_SAMple.SSAM_badseq');           # returns 'SEQ_SAMple'
    $column = $schema->columnpart('SEQ_SAMple.SSAM_badseq');         # returns 'SSAM_badseq'

=head3 HIVSchema DESCRIPTION

HIVSchema methods are used in L<Bio::DB::Query::HIVQuery> for table,
column, primary/foreign key manipulations based on the observed Los
Alamos HIV Sequence Database (LANL DB) naming conventions for their
CGI parameters. The schema is contained in an XML file
(C<lanl-schema.xml>) which is read into an HIVSchema object, in turn a
property of the HIVQuery object. HIVSchema methods are used to build
correct cgi queries in a way that attempts to preserve the context of
the relational database the query parameters represent.

=cut

package # hide from PAUSE
    HIVSchema;
# objects/methods to manipulate a version of the LANL HIV DB schema
# stored in XML
use XML::Simple;
use Bio::Root::Root;
use strict;

### constructor

=head3 HIVSchema CONSTRUCTOR

=head4 HIVSchema::new

 Title   : new
 Usage   : $schema = new HIVSchema( "lanl-schema.xml ");
 Function:
 Example :
 Returns : an HIVSchema object
 Args    : XML filename

=cut

sub new {
    my $class = shift;
    my @args = @_;
    my $self = {};
    if ($args[0]) {
	$self->{schema_ref} = loadHIVSchema($args[0]);
    }
    bless($self, $class);
    return $self;
}

### object methods

=head3 HIVSchema INSTANCE METHODS

=head4 HIVSchema tables

 Title   : tables
 Usage   : $schema->tables()
 Function: get all table names in schema
 Example :
 Returns : array of table names
 Args    : none

=cut

sub tables {
    # return array of all tables in schema
    local $_;
    my $self = shift;
    my $sref = $self->{schema_ref};
    Bio::Root::Root->throw("schema not initialized") unless $sref;
    my  @k = grep(/\./, keys %$sref);
    my %ret;
    foreach (@k) {
	s/\..*$//;
	$ret{$_}++;
    }
    @k = sort keys %ret;
    return @k;
}

=head4 HIVSchema columns

 Title   : columns
 Usage   : $schema->columns( [$tablename] );
 Function: return array of columns for specified table, or all columns in
           schema, if called w/o args
 Example :
 Returns :
 Args    : tablename or fieldname string

=cut

sub columns {
    # return array of columns for specified table
    # all columns in schema, if called w/o args
    local $_;
    my $self = shift;
    my ($tbl) = @_;
    my $sref = $self->{schema_ref};
    Bio::Root::Root->throw("schema not initialized") unless $sref;
    # trim column name
    $tbl =~ s/\..*$//;
    # check if table exists
    return () unless grep(/^$tbl$/i, $self->tables);
    my @k = sort keys %$sref;
    @k = grep (/^$tbl\./i, @k);
    foreach (@k) {
	s/^$tbl\.//;
    }
    return @k;
}

=head4 HIVSchema fields

 Title   : fields
 Usage   : $schema->fields();
 Function: return array of all fields in schema, in format "table.column"
 Example :
 Returns : array of all fields
 Args    : none

=cut

sub fields {
    # return array of all fields (Table.Column format) in schema
    my $self = shift;
    my $sref = $self->{schema_ref};
    Bio::Root::Root->throw("schema not initialized") unless $sref;
    my @k = sort keys %{$sref};
    return @k;
}

=head4 HIVSchema options

 Title   : options
 Usage   : $schema->options(@fieldnames)
 Function: get array of options (i.e., valid match data strings) available
           to specified field
 Example :
 Returns : array of match data strings
 Args    : [array of] fieldname string[s] in "table.column" format

=cut

sub options {
    # return array of options available to specified field
    my $self = shift;
    my ($sfield) = @_;
    my $sref = $self->{schema_ref};
    Bio::Root::Root->throw("schema not initialized") unless $sref;
    return $$sref{$sfield}{option} ? @{$$sref{$sfield}{option}} : ();
}

=head4 HIVSchema aliases

 Title   : aliases
 Usage   : $schema->aliases(@fieldnames)
 Function: get array of aliases to specified field[s]
 Example :
 Returns : array of valid query aliases for fields as spec'd in XML file
 Args    : [an array of] fieldname[s] in "table.column" format

=cut

sub aliases {
    # return array of aliases to specified field
    my $self = shift;
    my ($sfield) = @_;
    my $sref = $self->{schema_ref};
    my @ret;
    Bio::Root::Root->throw("schema not initialized") unless $sref;
    if ($sfield) {
	return $$sref{$sfield}{alias} ? @{$$sref{$sfield}{alias}} : ();
    }
    else { # all valid aliases
	map {push @ret, @{$$sref{$_}{alias}} if $$sref{$_}{alias}} $self->fields;
	return @ret;
    }
}

=head4 HIVSchema ankh

 Title   : ankh (annotation key hash)
 Usage   : $schema->ankh(@fieldnames)
 Function: return a hash translating fields to annotation keys for the
           spec'd fields.
           (Annotation keys are used for parsing the tab-delimited response
           to Bio::DB::Query::HIVQuery::_do_lanl_request.)
 Example :
 Returns : hash ref
 Args    : [an array of] fieldname[s] in "table.column" format

=cut

sub ankh {
    # return hash translating sfields to annotation keys for specified sfield(s)
    my $self = shift;
    my %ret = ();
    my @sfields = @_;
    my $sref = $self->{schema_ref};
    Bio::Root::Root->throw("schema not initialized") unless $sref;
    foreach (@sfields) {
	next unless $$sref{$_}{ankey};
	$ret{$_} = {'ankey'=>$$sref{$_}{ankey},'antype'=>$$sref{$_}{antype}};
    }
    return %ret;
}

=head4 HIVSchema tablepart

 Title   : tablepart (alias: tbl)
 Usage   : $schema->tbl(@fieldnames)
 Function: return the portion of the fieldname[s] that refer to the
           db table
 Example : $schema->tbl('SequenceEntry.SE_id'); # returns 'SequenceEntry'
 Returns : table name as string
 Args    : [an array of] fieldname[s] in "table.column" format

=cut

sub tablepart {
    # return the 'Table' part of the specified field(s)
    my $self = shift;
    my @sfields = @_;
    Bio::Root::Root->throw("schema not initialized") unless $self->{schema_ref};
    my ($squish,@ret, %ret);
    if ($sfields[0] eq '-s') {
	# squish : remove duplicates from the returned array
	$squish=1;
	shift @sfields;
    }
    foreach (@sfields) {
	push @ret, /^(.*)\./;
    }
    if ($squish) {
	# arg order is clobbered
	@ret{@ret} = undef;
	@ret = keys %ret;
    }
    return (wantarray ? @ret : $ret[0]);
}

sub tbl {
    # tablepart alias
    shift->tablepart(@_);
}

=head4 HIVSchema columnpart

 Title   : columnpart (alias: col)
 Usage   : $schema->col(@fieldnames)
 Function: return the portion of the fieldname[s] that refer to the
           db column
 Example : $schema->col('SequenceEntry.SE_id'); # returns 'SE_id'
 Returns : column name as string
 Args    : [an array of] fieldname[s] in "table.column" format

=cut

sub columnpart {
    # return the 'Column' part of the specified field(s)
    my $self = shift;
    my @sfields = @_;
    Bio::Root::Root->throw("schema not initialized") unless $self->{schema_ref};
    my @ret;
    foreach (@sfields) {
	push @ret, /\.(.*)$/;
    }
    return (wantarray ? @ret : $ret[0]);
}

sub col {
    # columnpart alias
    shift->columnpart(@_);
}

=head4 HIVSchema primarykey

 Title   : primarykey [alias: pk]
 Usage   : $schema->pk(@tablenames);
 Function: return the primary key of the specified table[s], as judged by
           the syntax of the table's[s'] fieldnames
 Example : $schema->pk('SequenceEntry') # returns 'SequenceEntry.SE_id'
 Returns : primary key fieldname[s] in "table.column" format, or null if
           no pk exists
 Args    : [an array of] table name[s] (fieldnames are ok, table part used)

=cut

sub primarykey {
    # return the primary key (in Table.Column format) of specified table(s)
    my $self = shift;
    my @tbl = @_;
    my @ret;
    Bio::Root::Root->throw("schema not initialized") unless $self->{schema_ref};
    foreach my $tbl (@tbl) {
	# trim column name
	$tbl =~ s/\..*$//;
	grep(/^$tbl$/i, $self->tables) ?
	    push(@ret, grep(/\.[0-9a-zA-Z]+_id/, grep(/$tbl/i,$self->fields))) :
	    push(@ret, "");
    }
    return (wantarray ? @ret : $ret[0]);
}

sub pk {
    # primarykey alias
    shift->primarykey(@_);
}

=head4 HIVSchema foreignkey

 Title   : foreignkey [alias: fk]
 Usage   : $schema->fk($intable [, $totable])
 Function: return foreign key fieldname in table $intable referring to
           table $totable, or all foreign keys in $intable if $totable
           unspec'd
 Example : $schema->fk('AUthor', 'SequenceEntry'); # returns 'AUthor_AU_SE_id'
 Returns : foreign key fieldname[s] in "table.column" format
 Args    : tablename [, optional foreign table name] (fieldnames are ok,
           table part used)

=cut

sub foreignkey {
    # return foreign key in in-table ($intbl) to to-table ($totbl)
    # or all foreign keys in in-table if to-table not specified
    # keys returned in Table.Column format
    my $self = shift;
    my ($intbl, $totbl) = @_;
    Bio::Root::Root->throw("schema not initialized") unless $self->{schema_ref};
    # trim col names
    $intbl =~ s/\..*$//;
    $totbl =~ s/\..*$// if $totbl;
    # check if in-table exists
    return () unless grep( /^$intbl/i, $self->tables);
    my @ret = grep( /$intbl\.(?:[0-9a-zA-Z]+_){2,}id/i, $self->fields);
    if ($totbl) {
	my $tpk = $self->primarykey($totbl);
	return (wantarray ? () : "") unless grep( /^$totbl/i, $self->tables) && $tpk;
	($tpk) = ($tpk =~ /\.(.*)$/);
	@ret = grep( /$tpk$/, @ret);
	return (wantarray ? @ret : $ret[0]);
    }
    else {
	# return all foreign keys in in-table
	return @ret;
    }
}

sub fk {
    # foreignkey alias
    shift->foreignkey(@_);
}

=head4 HIVSchema foreigntable

 Title   : foreigntable [alias ftbl]
 Usage   : $schema->ftbl( @foreign_key_fieldnames );
 Function: return tablename of table that foreign keys points to
 Example : $schema->ftbl( 'AUthor.AU_SE_id' ); # returns 'SequenceEntry'
 Returns : tablename
 Args    : [an array of] fieldname[s] in "table.column" format

=cut

sub foreigntable {
    # return table name that foreign key(s) point(s) to
    my $self = shift;
    my @fk = @_;
    my @ret;
    Bio::Root::Root->throw("schema not initialized") unless $self->{schema_ref};
    foreach (@fk) {
	my ($mnem, $fmnem) = /\.([0-9a-zA-Z]+)_([0-9a-zA-Z]+)_.*$/;
	next unless $mnem && $fmnem;
	# lookup based on Table.Column format of fields
	my $sf = [grep( /^[0-9a-zA-Z]+\.$fmnem\_/, $self->fields )]->[0];
	next unless $sf;
	($sf) = ($sf =~ /^([0-9a-zA-Z]+)\./);
	push @ret, $sf;
    }
    return (wantarray ? @ret : $ret[0]);
}

sub ftbl {
    # foreigntable alias
    shift->foreigntable(@_);
}

=head4 HIVSchema find_join

 Title   : find_join
 Usage   : $sch->find_join('Table1', 'Table2')
 Function: Retrieves a set of foreign and primary keys (in table.column
           format) that represents a join path from Table1 to Table2
 Example :
 Returns : an array of keys (as table.column strings) -or- an empty
           array if Table1 == Table2 -or- undef if no path exists
 Args    : two table names as strings

=cut

sub find_join {
    my $self = shift;
    my ($tgt, $tbl) = @_;
    my ($stack, $revstack, $found, $revcut) = ([],[], 0, 4);
    $self->_find_join_guts($tgt, $tbl, $stack, \$found);
    if ($found) {
	if (@$stack > $revcut) {
	    # reverse order of tables, see if a shorter path emerges
	    $found = 0;
	    $self->_find_join_guts($tgt, $tbl, $revstack, \$found, 1);
	    return (@$stack <= @$revstack ? @$stack : @$revstack);
	}
	return @$stack;
    }
    else {
	return undef;
    }
}

=head4 HIVSchema _find_join_guts

 Title   : _find_join_guts
 Usage   : $sch->_find_join_guts($table1, $table2, $stackref, \$found, $reverse)
           (call with $stackref = [], $found=0)
 Function: recursive guts of find_join
 Example :
 Returns : if a path is found, $found==1 and @$stackref contains the keys
           in table.column format representing the path; if a path is not
           found, $found == 0 and @$stackref contains garbage
 Args    : $table1, $table2 : table names as strings
           $stackref : an arrayref to an empty array
           \$found   : a scalar ref to the value 0
           $rev : if $rev==1, the arrays of table names will be reversed;
                  this can give a shorter path if cycles exist in the
                  schema graph

=cut

sub _find_join_guts {
    my $self = shift;
    my ($tbl, $tgt, $stack, $found, $rev) = @_;
    return () if $tbl eq $tgt;
    my $k = $self->pk($tbl);
    if ($k) {
	# all fks pointing to pk
	my @fk2pk = map {
	    $self->fk($_, $k) || ()
	} ($rev ? reverse $self->tables : $self->tables);
	# skip keys already on stack
	if (@$stack) {
	    (@$stack == 1) && do {
		@fk2pk = grep (!/$$stack[0]/, @fk2pk);
	    };
	    (@$stack > 1 ) && do {
		@fk2pk = map { my $f=$_; grep(/$f/, @$stack) ? () : $f } @fk2pk;
	    };
	}
	foreach my $f2p (@fk2pk) { # tables with fks pointing to pk
	    push @$stack, $f2p;
	    if ($self->tbl($f2p) eq $tgt) { # this fk's table is the target
		# found it
		$$found = 1;
		return;
	    }
	    else {
		#keep looking
		$self->_find_join_guts($self->tbl($f2p), $tgt, $stack, $found, $rev);
		return if $$found;
	    }
	}
    }
    # all fks in $tbl
    my @fks = ($rev ? reverse $self->fk($tbl) : $self->fk($tbl));
    #skip keys already on stack
    if (@$stack) {
	(@$stack == 1) && do {
	    @fks = grep(!/$$stack[0]/, @fks);
	};
	(@$stack > 1) && do {
	    @fks = map { my $f=$_; grep(/$f/, @$stack) ? () : $f } @fks;
	};
    }
    # all fks in table
    if (@fks) {
	for my $f (@fks) {
	    push @$stack, $f;
	    if ($self->ftbl($f) eq $tgt) { #found it
		$$found = 1;
		return;
	    }
	    else {
		$self->_find_join_guts($self->ftbl($f), $tgt, $stack, $found, $rev);
		$$found ? return : pop @$stack;
	    }
	}
    }
    else {
	pop @$stack;
	return;
    }
}

=head4 HIVSchema loadSchema

 Title   : loadHIVSchema [alias: loadSchema]
 Usage   : $schema->loadSchema( $XMLfilename )
 Function: read (LANL DB) schema spec from XML
 Example : $schema->loadSchema('lanl-schema.xml');
 Returns : hashref to schema data
           Keys are fieldnames in "table.column" format.
           Each value is a hashref with the following properties:
           {name}    : HIVWEB 'table.column' format fieldname,
                       can be used directly in the cgi query
           {aliases} : ref to array containing valid aliases/shortcuts for
                       {name}; can be used in routines creating the HTML query
           {options} : ref to array containing valid matchdata for this field
                       can be used directly in the HTML query
           {ankey}   : contains the annotation key for this field used with
                       Bioperl annotation objects
           {..attr..}: ..value_of_attr.. for this field (app-specific metadata)
 Args    :

=cut

sub loadHIVSchema {
    my $fn = shift;
    Bio::Root::Root->throw("loadHIVSchema: schema file not found") unless -e $fn;
    my $q = XML::Simple->new(ContentKey=>'name',NormalizeSpace=>2,ForceArray=>1);
    my %ret;
    my $ref = $q->XMLin($fn);
    my @sf = keys %{$$ref{sfield}};
    foreach (@sf) {
	my $h = $$ref{sfield}{$_};
	$ret{$_} = $h;
	foreach my $ptr ($$h{option}, $$h{alias}) {
	    if ($ptr) {
		# kludge for XMLin: appears to convert to arrays, if there
		# exists a tag without content, but to convert to hashes
		# with content as key, if all tags possess content
		if (ref($ptr) eq 'HASH') {
		    my @k = keys %{$ptr};
		    if (grep /desc/, keys %{$ptr->{$k[0]}}) {
			# slurp the desc's
			$$h{desc} = [ map { $$ptr{$_}->{desc} } @k ];
		    }
		    # now overwrite with keys (descs in same order...)
		    $ptr = [@k];
		}
		elsif (ref($ptr) eq 'ARRAY') {
		    $ptr = [map { ref eq 'HASH' ? $_->{name} : $_ } @{$ptr}]
		}
		else {
		    1; # stub : doh!
		}
	    }
	}
	for my $ptr ($$h{ankey}) {
	    # flatten
	    my $ank = [keys %{$ptr}]->[0];
	    if (!defined $ank) {
		delete $$h{ankey};
	    }
	    else {
		$h->{antype} = $ptr->{$ank}{antype};
		$ptr = $ank;
	    }
	}
    }
    return \%ret;
}

sub loadSchema {
    my $self = shift;
    $self->{schema_ref} = loadHIVSchema(shift);
}

# below, dangerous

=head4 HIVSchema _sfieldh

 Title   : _sfieldh
 Usage   : $schema->_sfieldh($fieldname)
 Function: get hashref to the specified field hash
 Example :
 Returns : hashref
 Args    : fieldname in "table.column" format

=cut

sub _sfieldh {
    # return reference to the specified field hash
    my $self = shift;
    my ($sfield) = @_;
    return ${$self->{schema_ref}}{$sfield};
}

1;

=head2 Class QRY -  a query algebra for HIVQuery

=head3 QRY SYNOPSIS

  $Q = new QRY(
               new R(
                      new Q('coreceptor', 'CXCR4'),
                      new Q('country', 'ZA')
                     )
               );
  QRY::Eq(QRY::And($Q, $Q), $Q);                     # returns 1
  QRY::Eq(QRY::Or($Q, $Q), $Q);                      # returns 1
  $Q2 = $Q1->clone;
  $Q2 = new QRY(
                new R(
                       new Q( 'coreceptor', 'CCR5' ),
                       new Q( 'country', 'ZA')
                      )
               );
  (QRY::And($Q, $Q2))->isnull;                       # returns 1
  $Q3 = QRY::Or($Q, $Q2);
  print $Q3->A;                                      # prints '(CCR5 CXCR4)[coreceptor] (ZA)[country]'

=head3 QRY DESCRIPTION

The QRY package provides a query parser for
L<Bio::DB::Query::HIVQuery>. Currently, the parser supports AND, OR,
and () operations. The structure of the LANL cgi makes it tricky to
perform NOTs, though this could be implemented if the desire were
great.

Two class methods do the work. C<QRY::_parse_q> does a first-pass
parse of the query string. C<QRY::_make_q> interprets the parse tree
as returned by C<QRY::_parse_q> and produces an array of hash
structures that can be used directly by C<Bio::DB::Query::HIVQuery>
query execution methods. Validation of query fields and options is
performed at the C<Bio::DB::Query::HIVQuery> level, not here.

C<QRY> objects are collections of C<R> (or request) objects, which are
in turn collections of C<Q> (or atomic query) objects. C<Q> objects
represent a query on a single field, with match data options C<OR>ed
together, e.g. C<(A B)[subtype]>. C<R> objects collect C<Q> objects
that could be processed in a single HTTP request; i.e., a set of
atomic queries each having different fields C<AND>ed together, such as

  (A B)[subtype] AND ('CCR5')[coreceptor] AND (US CA)[country]

The C<QRY> object collects C<R>s that cannot be reduced (through
logical operations) to a single HTTP request, e.g.

  ((C)[subtype] AND (SI)[phenotype]) OR ( (D)[subtype] AND (NSI)[phenotype] ),

which cannot be got in one go through the current LANL cgi
implementation (as far as I can tell). The parser will simplify
something like

  ((C)[subtype] AND (SI)[phenotype]) OR ((C)[subtype] AND (NSI)[phenotype])

to the single request

  (C)[subtype] AND (NSI SI)[phenotype]

however.

The operators C<&> and C<|> are overloaded to C<QRY::And> and
C<QRY::Or>, to get Perl precedence and grouping for free. C<bool> is
overloaded to get symbolic tests such as C<if ($QRY) {stuff}>. C<==>
is overloaded with C<QRY::Eq> for convenience. No overloading is done
for C<R> or C<Q>.

=cut

# a query algebra for HIVQuery
#
# Each Q object is an 'atomic' query, written as (data)[field]
# (a b ...)[X] equals (a)[X] | (b)[X] | ...
# Each R object represents a single HTTP request to the db
#  contains an array of Q (atomic) objects (q1, q2, ...)
#  the R object is interpreted as  q1 & q2 & ...
# Each QRY object represents a series of HTTP requests to the db
#  contains an array of R (request) objects (R1, R2, ...)
#  the QRY object is interpreted as R1 | R2 | ...
#
# & and | operations are specified for each type

package # hide from PAUSE
    QRY;
use strict;
$QRY::NULL = new QRY();


use overload
    "|" => \&Or,
    "&" => \&And,
    "bool" => \&Bool,
    "==" => \&Eq;


# query language emulator
# supports only AND and OR, any groupings
#
# syntax rules:
# query atom: bareword [field] OR  (bareword ...) [field]
# only single bareword allowed between []
# annotation fields in {} (only bareword lists allowed between {})
# () can group query atoms joined by operators (AND or OR)
# () containing only barewords MUST be followed by a field descriptor [field]
# empty [] not allowed
# query atoms joined with AND by default
# barewords are associated (ORed within) the next field descriptor in the line

# follow the parse tree, creating new QRY objects as needed in @q, and
# construct a logical expression using & and | symbols.
# These are overloaded for doing ands and ors on QRY objects;
# to get the final QRY object, eval the resulting expression $q_expr.
# QRY object will be translated into (possibly multiple) hashes
# conforming to HIVQuery parameter requirements.

=head4 QRY _make_q

 Title   : _make_q
 Usage   : QRY::_make_q($parsetree)
 Function: creates hash structures suitable for HIVQuery from parse tree
           returned by QRY::_parse_q
 Example :
 Returns : array of hashrefs of query specs
 Args    : a hashref

=cut

sub _make_q {
    my $ptree = shift;
    my ($q_expr, @q, @an, $query, @dbq);
    _make_q_guts($ptree, \$q_expr, \@q, \@an);
    $query = eval $q_expr;
    throw Bio::Root::Root(-class=>'Bio::Root::Exception',
			  -text=>$@,
			  -value=>$q_expr) if $@;
    return {} if $query->isnull;
    foreach my $rq ($query->requests) {
	my $h = {'query'=>{}};
	foreach ($rq->atoms) {
	    my @d = split(/\s+/, $_->dta);
	    foreach my $d (@d) {
		$d =~ s/[+]/ /g; ###! _ to [+]
		$d =~ s/'//g;
	    }
	    $h->{'query'}{$_->fld} = (@d == 1) ? $d[0] : [@d];
	}
	$h->{'annot'} = [@an] if @an;
	push @dbq, $h;
    }
    return @dbq;
}

=head4 QRY _make_q_guts

 Title   : _make_q_guts (Internal class method)
 Usage   : _make_q_guts($ptree, $q_expr, $qarry, $anarry)
 Function: traverses the parse tree returned from QRY::_parse_q, checking
           syntax and creating HIVQuery-compliant query structures
 Example :
 Returns :
 Args    : $parse_tree (hashref), $query_expression (scalar string ref),
           $query_array (array ref : stack for returning query structures),
           $annotation_array (array ref : stack for returning annotation
           fields)

=cut

sub _make_q_guts {
    my ($ptree, $q_expr, $qarry, $anarry) = @_;
    my (@words, $o);
    eval { # catch
	foreach (@{$ptree->{cont}}) {
	    m{^AND$} && do {
		$$q_expr .= "&";
		next;
	    };
	    m{^OR$} && do {
		$$q_expr .= "|";
		next;
	    };
	    m{^HASH} && do {
		for my $dl ($_->{delim}) {
		    ($dl =~ m{\(}) && do {
			if (grep /^HASH/, @{$_->{cont}}) {
			    $$q_expr .= "&" unless !$$q_expr || !length($$q_expr) || (substr($$q_expr, -1, 1) =~ /[&|(]/);
			    $$q_expr .= "(";
			    _make_q_guts($_,$q_expr,$qarry,$anarry);
			    $$q_expr .= ")";
			}
			else {
			    my @c;
			    my $c = join(' ',@{$_->{cont}});
			    $c =~ s/,/ /g;
			    Bio::Root::Root->throw("query syntax error: unmatched ['\"]") if (@c = ($c =~ /(['"])/g)) % 2;
			    @c = split(/\s*(['"])\s*/, $c);
			    do {
				$c = shift @c;
				if ($c =~ m{['"]}) {
				    $c = join('', ($c, shift @c, shift @c));
				    $c =~ s/\s+/+/g; ###! _ to +
				    push @words, $c;
				}
				else {
				    push @words, split(/\s+/,$c);
				}
			    } while @c;
			}
			last;
		    };
		    ($dl =~ m{\[}) && do {
			Bio::Root::Root->throw("syntax error: empty field descriptor") unless @{$_->{cont}};
			Bio::Root::Root->throw("syntax error: more than one field descriptor in square brackets") unless @{$_->{cont}} == 1;

			push @{$qarry}, new QRY( new R( new Q( $_->{cont}->[0], @words)));
			# add default operation if nec
			$$q_expr .= "&" unless !$$q_expr || !length($$q_expr) || (substr($$q_expr, -1, 1) =~ /[&|(]/);
			$$q_expr .= "\$q[".$#$qarry."]";
			@words = ();
			last;
		    };
		    ($dl =~ m{\{}) && do {
			foreach my $an (@{$_->{cont}}) {
			    ($an =~ /^HASH/) && do {
				if ($an->{delim} eq '[') {
				    push @$anarry, @{$an->{cont}};
				}
				else {
				    Bio::Root::Root->throw("query syntax error: only field descriptors (with or without square brackets) allowed in annotation spec");
				}
				next;
			    };
			    do { #else
				push @$anarry, $an;
				next;
			    };
			}
			last;
		    };
		    do {
			1; #else stub
		    };
		}
		next;
	    };
	    do { # else, bareword
		if ($o) {
		    $words[-1] .= "+$_"; ####! _ to +
		}
		else {
		    push @words, $_;
		}
		m/['"]/ && ($o = !$o);
	    };
	} # @{ptree->{cont}}
	Bio::Root::Root->throw("query syntax error: no search fields specified")
	    unless $$q_expr =~ /q\[[0-9]+\]/;
    };
    $@ ?
	throw Bio::Root::Root(-class=>'Bio::QueryStringSyntax::Exception',
			      -text=>$@,
			      -value=>$$q_expr)
	: return 1;
}

=head4 QRY _parse_q

 Title   : _parse_q
 Usage   : QRY::_parse_q($query_string)
 Function: perform first pass parse of a query string with some syntax
           checking, return a parse tree suitable for QRY::_make_q
 Example : QRY::_parse_q(" to[be] OR (not to)[be] ");
 Returns : hashref
 Args    : query string

=cut

# parse qry string into a branching tree structure
# each branch tagged by the opening delimiter ( key 'delim' )
# content (tokens and subbranch hashes) placed in l2r order in
# @{p->{cont}}
sub _parse_q {
    local $_;
    my $qstr = shift;
    my $illegal = qr/[^a-zA-Z0-9-_<>=,\.\(\[\{\}\]\)\s'"]/;
    my $pdlm = qr/[\{\[\(\)\]\}]/;
    my %md = ('('=>')', '['=>']','{'=>'}');
    my @tok =  grep !/^\s*$/, split /($pdlm)/, $qstr;
    return {} unless @tok;
    my @pstack = ();
    my @dstack = ();
    my ($ptree, $p);

    eval { #catch
	Bio::Root::Root->throw("query syntax error: illegal character") if $qstr =~ /$illegal/;

	$ptree = $p = {'delim'=>'*'};
	foreach (@tok) {
	    #trim whsp
	    s/^\s+//;
	    s/\s+$//;
	    m{[\(\[\{]} && do {
		my $new = {'delim'=>$_};
		$p->{cont} = [] unless $p->{cont};
		push @{$p->{cont}}, $new;
		push @pstack, $p;
		push @dstack, $_;
		$p = $new;
		next;
	    };
	    m{[\)\]\}]} && do {
		my $d = pop @dstack;
		if ($md{$d} eq $_) {
		    $p = pop @pstack;
		    Bio::Root::Root->throw("query syntax error: unmatched \"$_\"") unless $p;
		}
		else {
		    Bio::Root::Root->throw("query syntax error: saw \"$_\" before matching \"$md{$d}\"");
		}
		next;
	    };
	    do { # else
		$p->{cont} = [] unless $p->{cont};
		push @{$p->{cont}}, split(/\s+/);
	    };
	}
    };
    $@ ?
	throw Bio::Root::Root(-class=>'Bio::QueryStringSyntax::Exception',
			      -text=>$@,
			      -value=>"")
	: return $ptree;
}

## QRY constructor

=head3 QRY CONSTRUCTOR

=head4 QRY Constructor

 Title   : QRY constructor
 Usage   : $QRY = new QRY()
 Function:
 Example :
 Returns :
 Args    : array of R objects, optional

=cut

sub new {
    my $class = shift;
    my @args = @_;
    my $self = {};
    $self->{requests} = [];
    bless($self, $class);
    $self->put_requests(@args) if @args;
    return $self;
}

## QRY instance methods

=head3 QRY INSTANCE METHODS

=head4 QRY requests

 Title   : requests
 Usage   : $QRY->requests
 Function: get/set array of requests comprising this QRY object
 Example :
 Returns :
 Args    : array of class R objects

=cut

sub requests {
    my $self = shift;
    $self->put_requests(@_) if @_;
    return @{$self->{'requests'}};
}

=head4 QRY put_requests

 Title   : put_requests
 Usage   : $QRY->put_request(@R)
 Function: add object of class R to $QRY
 Example :
 Returns :
 Args    : [an array of] of class R object[s]

=cut

sub put_requests {
    my $self = shift;
    my @args = @_;
    foreach (@args) {
	Bio::Root::Root->throw('requires type R (request)') unless ref && $_->isa('R');
	push @{$self->{requests}}, $_;
    }
    return @args;
}

=head4 QRY isnull

 Title   : isnull
 Usage   : $QRY->isnull
 Function: test if QRY object is null
 Example :
 Returns : 1 if null, 0 otherwise
 Args    :

=cut

sub isnull {
    my $self = shift;
    return ($self->requests) ? 0 : 1;
}

=head4 QRY A

 Title   : A
 Usage   : print $QRY->A
 Function: get a string representation of QRY object
 Example :
 Returns : string scalar
 Args    :

=cut

sub A {
    my $self = shift;
    return join( "\n", map {$_->A} $self->requests );
}

=head4 QRY len

 Title   : len
 Usage   : $QRY->len
 Function: get number of class R objects contained by QRY object
 Example :
 Returns : scalar
 Args    :

=cut

sub len {
    my $self = shift;
    return scalar @{$self->{'requests'}};
}

=head4 QRY clone

 Title   : clone
 Usage   : $QRY2 = $QRY1->clone;
 Function: create and return a clone of the object
 Example :
 Returns : object of class QRY
 Args    :

=cut

sub clone {
    local $_;
    my $self = shift;
    my $ret = QRY->new();
    foreach ($self->requests) {
	$ret->put_requests($_->clone);
    }
    return $ret;
}

## QRY class methods

=head3 QRY CLASS METHODS

=head4 QRY Or

 Title   : Or
 Usage   : $QRY3 = QRY::Or($QRY1, $QRY2)
 Function: logical OR for QRY objects
 Example :
 Returns : a QRY object
 Args    : two class QRY objects

=cut

sub Or {
    local $_;
    my ($q, $r, $rev_f) = @_;
    Bio::Root::Root->throw('requires type QRY') unless ref($q) && $q->isa('QRY');
    Bio::Root::Root->throw('requires type QRY') unless ref($r) && $r->isa('QRY');
    if ($q->isnull) {
	return $r->clone;
    }
    elsif ($r->isnull) {
	return $q->clone;
    }
    do {my $qq = $q; $q=$r; $r=$qq} if ($q->len > $r->len);
    my @rq_r = $r->requests;
    my @rq_q = $q->requests;
    my (@cand_rq, @ret_rq);
    # search for simplifications
    my @now = @rq_q;
    my @nxt =();
    foreach (@rq_r) {
	my $found = 0;
	while (my $rq = pop @now) {
	    my @result = R::Or($rq, $_);
	    if (@result==1) {
		push @cand_rq, $result[0]->clone;
		$found = 1;
		last;
	    }
	    else {
		push @nxt, $rq;
	    }
	}
	push @cand_rq, $_->clone unless ($found);
	# @now becomes unexamined @rq_q's plus failed @rq_q's
	@now = (@now, @nxt);
    }
    push @cand_rq, map {$_->clone} @now; # add all failed @rq_q's
    # squeeze out redundant requests
    while (my $rq = pop @cand_rq) {
	push @ret_rq, $rq unless @cand_rq && grep {R::Eq($rq, $_)} @cand_rq;
    }
    return new QRY( @ret_rq );
}

=head4 QRY And

 Title   : And
 Usage   : $QRY3 = QRY::And($QRY1, $QRY2)
 Function: logical AND for QRY objects
 Example :
 Returns : a QRY object
 Args    : two class QRY objects

=cut

sub And {
    my ($q, $r, $rev_f) = @_;
    Bio::Root::Root->throw('requires type QRY') unless ref($q) && $q->isa('QRY');
    Bio::Root::Root->throw('requires type QRY') unless ref($r) && $r->isa('QRY');
    return ($QRY::NULL) if ($q->isnull || $r->isnull);
    my (@cand_rq, @ret_rq);
    foreach my $rq_r ($r->requests) {
	foreach my $rq_q ($q->requests) {
	    my ($rq) = R::And($rq_r, $rq_q);
	    push @cand_rq, $rq unless $rq->isnull;
	}
    }
    return $QRY::NULL unless @cand_rq;
    # squeeze out redundant requests
    while (my $rq = pop @cand_rq) {
	push @ret_rq, $rq unless @cand_rq && grep {R::Eq($rq, $_)} @cand_rq;
    }
    return new QRY( @ret_rq );
}

=head4 QRY Bool

 Title   : Bool
 Usage   : QRY::Bool($QRY1)
 Function: allows symbolic testing of QRY object when bool overloaded
 Example : do {stuff} if $QRY1 *same as* do {stuff} if !$QRY1->isnull
 Returns :
 Args    : a class QRY object

=cut

sub Bool {
    my $q = shift;
    Bio::Root::Root->throw('requires type QRY') unless ref($q) && $q->isa('QRY');
    return $q->isnull ? 0 : 1;
}

=head4 QRY Eq

 Title   : Eq
 Usage   : QRY::Eq($QRY1, $QRY2)
 Function: test if R objects in two QRY objects are the same
           (irrespective of order)
 Example :
 Returns : 1 if equal, 0 otherwise
 Args    : two class QRY objects

=cut

sub Eq {
    my ($q, $r, $rev_f) = @_;
    Bio::Root::Root->throw('requires type QRY') unless ref($q) && $q->isa('QRY');
    Bio::Root::Root->throw('requires type QRY') unless ref($r) && $r->isa('QRY');
    return 0 unless $q->len == $r->len;
    foreach my $rq_q ($q->requests) {
	my $found = 0;
	foreach my $rq_r ($r->requests) {
	    if (R::Eq($rq_q,$rq_r)) {
		$found = 1;
		last;
	    }
	}
	return 0 unless $found;
    }
    return 1;
}

1;

=head2 Class R - request objects for QRY algebra

=head3 R SYNOPSIS

  $R = new R( $q1, $q2 );
  $R->put_atoms($q3);
  $R->del_atoms('coreceptor', 'phenotype');
  return $R->clone;
  $R1 = new R( new Q('subtype', 'B') );
  $R2 = new R( new Q('subtype', 'B C'),
               new Q('country', 'US') );
  R::Eq( (R::And($R1, $R2))[0],
         new R( new Q('subtype', 'B' ),
                new Q('country', 'US') ));                 # returns 1
  QRY::Eq( new QRY(R::Or($R1, $R2)), new QRY($R1, $R2) );  # returns 1
  R::In( (R::And($R1, $R2))[0], $R1 );                     # returns 1

=head3 R DESCRIPTION

Class R objects contain a list of atomic queries (class Q
objects). Each class R object represents a single HTTP request to the
LANL DB. When converted to a DB query, the class Q objects contained
by an R object are effectively C<AND>ed.

=cut

package # hide from PAUSE
    R;
use strict;
$R::NULL = R->new();


## R constructor

=head3 R CONSTRUCTOR

=head4 R constructor

 Title   : R constructor
 Usage   : $R = new R()
 Function: create a new R (request) object
 Example :
 Returns : class R (request) object
 Args    : optional, array of class Q objects

=cut

sub new {
    my $class = shift;
    my @args = @_;
    my $self = {};
    $self->{atoms} = {};
    bless($self, $class);
    $self->put_atoms(@args) if @args;
    return $self;
}

## R instance methods

=head3 R INSTANCE METHODS

=head4 R len

 Title   : len
 Usage   : $R->len
 Function: get number of class Q objects contained in R object
 Example :
 Returns : scalar
 Args    :

=cut

sub len {
    my $self = shift;
    return scalar @{[keys %{$self->{'atoms'}}]};
}

=head4 R atoms

 Title   : atoms
 Usage   : $R->atoms( [optional $field])
 Function: get array of class Q (atomic query) objects in class R object
 Example : $R->atoms(); $R->atoms('coreceptor')
 Returns : array of class Q objects (all Qs or those corresponding to $field
           if present)
 Args    : optional, scalar string

=cut

sub atoms {
    local $_;
    # returns an array of atoms
    # no arg: all atoms;
    # args: atoms with specified fields
    my $self = shift;
    my @flds = (@_ ? @_ : keys %{$self->{'atoms'}});
    return wantarray ? map { $self->{'atoms'}->{$_} } @flds : $self->{'atoms'}->{$flds[0]};
}

=head4 R fields

 Title   : fields
 Usage   : $R->fields
 Function: get array of fields of all Q objects contained in $R
 Example :
 Returns : array of scalars
 Args    :

=cut

sub fields {
    my $self = shift;
    return keys %{$self->{'atoms'}};
}

=head4 R put_atoms

 Title   : put_atoms
 Usage   : $R->put_atoms( @q )
 Function: AND an atomic query (class Q object) to the class R object's list
 Example :
 Returns : void
 Args    : an [array of] class Q object[s]

=cut

sub put_atoms {
    # AND this atom to the request
    local $_;
    my $self = shift;
    my @args = @_;
    foreach (@args) {
	Bio::Root::Root->throw('requires type Q (atom)') unless ref && $_->isa('Q');
	if ($self->atoms($_->fld)) {
	    my $a = Q::qand( $self->atoms($_->fld), $_ );
	    if ($a->isnull) {
		delete $self->{'atoms'}->{$_->fld};
	    }
	    else {
		$self->{atoms}->{$_->fld} = $a->clone;
	    }
	}
	else {
	    $self->{atoms}->{$_->fld} = $_->clone;
	}
    }
    return;
}

=head4 R del_atoms

 Title   : del_atoms
 Usage   : $R->del_atoms( @qfields )
 Function: removes class Q objects from R object's list according to the
           field names given in arguments
 Example :
 Returns : the class Q objects deleted
 Args    : scalar array of field names

=cut

sub del_atoms {
    # remove atoms by field from request
    local $_;
    my $self = shift;
    my @args = @_;
    return () unless @args;
    my @ret;
    foreach (@args) {
	push @ret, delete $self->{'atoms'}->{$_};
    }
    return @ret;
}

=head4 R isnull

 Title   : isnull
 Usage   : $R->isnull
 Function: test if class R object is null
 Example :
 Returns : 1 if null, 0 otherwise
 Args    :

=cut

sub isnull {
    my $self = shift;
    return ($self->len) ? 0 : 1;
}

=head4 R A

 Title   : A
 Usage   : print $R->A
 Function: get a string representation of class R object
 Example :
 Returns : string scalar
 Args    :

=cut

sub A {
    my $self = shift;
    my @a = sort {$a->fld cmp $b->fld} $self->atoms;
    return join(" ", map {$_->A} @a);
}

=head4 R clone

 Title   : clone
 Usage   : $R2 = $R1->clone;
 Function: create and return a clone of the object
 Example :
 Returns : object of class R
 Args    :

=cut

sub clone {
    local $_;
    my $self = shift;
    my $ret = R->new();
    foreach ($self->atoms) {
	$ret->put_atoms($_->clone);
    }
    return $ret;
}

## R class methods

=head3 R CLASS METHODS

=head4 R In

 Title   : In
 Usage   : R::In($R1, $R2)
 Function: tests whether the query represented by $R1 would return a subset
           of items returned by the query represented by $R2
 Example : print "R2 gets those and more" if R::In($R1, $R2);
 Returns : 1 if R1 is subset of R2, 0 otherwise
 Args    : two class R objects

=cut

sub In {
    local $_;
    my ($s, $t) = @_;
    Bio::Root::Root->throw('requires type R (request)') unless ref($s) && $s->isa('R');
    Bio::Root::Root->throw('requires type R (request)') unless ref($t) && $t->isa('R');
    return 1 if ($s->isnull);
    # common fields
    my @cf = grep {defined} map {my $f=$_; grep /^$f$/,$s->fields} $t->fields;
    return 0 unless @cf==$t->len;
    foreach (@cf) {
	my @sd = split(/\s+/, $s->atoms($_)->dta);
	my @td = split(/\s+/, $t->atoms($_)->dta);
	my @cd = grep {defined} map {my $d=$_; grep /^$d$/, @td} @sd;
	return 0 unless @cd==@sd;
    }
    return 1;
}

=head4 R And

 Title   : And
 Usage   : @Rresult = R::And($R1, $R2)
 Function: logical AND for R objects
 Example :
 Returns : an array containing class R objects
 Args    : two class R objects

=cut

sub And {
    local $_;
    my ($s, $t) = @_;
    Bio::Root::Root->throw('requires type R (request)') unless ref($s) && $s->isa('R');
    Bio::Root::Root->throw('requires type R (request)') unless ref($t) && $t->isa('R');
    return ($R::NULL) if ($s->isnull || $t->isnull);

    do { my $ss = $s; $s = $t; $t = $ss } if ( $s->len > $t->len );
    # $t has at least as many fields defined than $s ($t is more restrictive)

    # common fields
    my @cf = grep {defined} map {my $sf = $_; grep /$sf/, $t->fields } $s->fields;
    my $ret = R->new();
    my $v = $t->clone;
    $v->del_atoms(@cf);
    my $u = $s->clone;
    $u->del_atoms(@cf);

    # And the atoms with identical fields

    foreach (@cf) {
	my ($a) = Q::qand($s->atoms($_), $t->atoms($_));
	if ($a->isnull) {
	    return $R::NULL;
	}
	else {
	    $ret->put_atoms($a);
	}
    }
    # put the private atoms
    $ret->put_atoms($u->atoms, $v->atoms);
    return ($ret);

}

=head4 R Or

 Title   : Or
 Usage   : @Rresult = R::Or($R1, $R2)
 Function: logical OR for R objects
 Example :
 Returns : an array containing class R objects
 Args    : two class R objects

=cut

sub Or {
    local $_;
    my ($s, $t) = @_;
    Bio::Root::Root->throw('requires type R (request)') unless ref($s) && $s->isa('R');
    Bio::Root::Root->throw('requires type R (request)') unless ref($t) && $t->isa('R');
    if ($s->isnull) {
	return $t->clone;
    }
    elsif ($t->isnull) {
	return $s->clone;
    }
    return $s->clone if (R::In($t, $s));
    return $t->clone if (R::In($s, $t));

    # try simplifying
    do { my $ss = $s; $s = $t; $t = $ss } if ( $s->len > $t->len );
    # common fields
    my @cf = grep {defined} map {my $sf = $_; grep /$sf/, $t->fields } $s->fields;
    #
    if ($t->len == @cf) {
    # all atoms equal within fields but one? If yes, simplify...
	my @df = grep {!Q::qeq($s->atoms($_), $t->atoms($_))} @cf;
	if (@df == 1) {
	    my ($a) = Q::qor($s->atoms($df[0]), $t->atoms($df[0]));
	    my $ret = $s->clone;
	    $ret->del_atoms($df[0]);
	    $ret->put_atoms($a);
	    return ($ret);
	}
    }

    # neither request contains the other, and the requests cannot be
    # simplified; reflect back (clones of) the input...
    return ($s->clone, $t->clone);

}

=head4 R Eq

 Title   : Eq
 Usage   : R::Eq($R1, $R2)
 Function: test if class Q objects in two R objects are the same
           (irrespective of order)
 Example :
 Returns : 1 if equal, 0 otherwise
 Args    : two class R objects

=cut

sub Eq {
    local $_;
    my ($s, $t) = @_;
    Bio::Root::Root->throw('requires type R (request)') unless ref($s) && $s->isa('R');
    Bio::Root::Root->throw('requires type R (request)') unless ref($t) && $t->isa('R');
    my @sf = $s->fields;
    my @tf = $t->fields;
    return 0 unless @sf==@tf;
    my @cf = grep {defined} map {my $f=$_; grep /^$f$/,@sf} @tf;
    return 0 unless @cf==@tf;
    foreach (@cf) {
	return 0 unless Q::qeq($s->atoms($_), $t->atoms($_));
    }
    return 1;
}
1;

=head2 Class Q -  atomic query objects for QRY algebra

=head3 Q SYNOPSIS

    $q = new Q('coreceptor', 'CXCR4 CCR5');
    $u = new Q('coreceptor', 'CXCR4');
    $q->fld;                                 # returns 'coreceptor'
    $q->dta;                                 # returns 'CXCR4 CCR5'
    print $q->A;                             # prints '(CXCR4 CCR5)[coreceptor]
    Q::qeq($q, $u);                          # returns 0
    Q::qeq( Q::qor($q, $q), $q );            # returns 1
    Q::qin($u, $q)                           # returns 1
    Q::qeq(Q::qand($u, $q), $u );            # returns 1

=head3 Q DESCRIPTION

Class Q objects represent atomic queries, that can be described by a
single LANL cgi parameter=value pair. Class R objects (requests) are
built from class Qs. The logical operations at the higher levels
(C<QRY, R>) ultimately depend on the lower level operations on Qs:
C<qeq, qin, qand, qor>.

=cut

package # hide from PAUSE
    Q;
use strict;
$Q::NULL = Q->new();

## Q constructor

=head3 Q CONSTRUCTOR

=head4 Q constructor

 Title   : Q constructor
 Usage   : $q = new Q($field, $data)
 Function: create a new Q (atomic query) object
 Example :
 Returns : class Q object
 Args    : optional $field, $data strings

=cut

sub new {
    local $_;
    my ($class,@args) = @_;
    my $self={};
    foreach (@args) { s/^\s+//; s/\s+$//; }
    my ($fld, @dta) = @args;
    $self->{fld}=$fld;
    $self->{dta}=join(" ", @dta);
    bless($self, $class);
    return $self;
}

## Q instance methods

=head3 Q INSTANCE METHODS

=head4 Q isnull

 Title   : isnull
 Usage   : $q->isnull
 Function: test if class Q object is null
 Example :
 Returns : 1 if null, 0 otherwise
 Args    :

=cut

sub isnull {
    my $self = shift;
    Bio::Root::Root->throw("requires type Q (atom)") unless ref($self) && $self->isa('Q');
    return 1 unless (($self->fld && length($self->fld)) || ($self->dta && length($self->dta)));
    return 0;
}

=head4 Q fld

 Title   : fld
 Usage   : $q->fld($field)
 Function: get/set fld (field name) property
 Example :
 Returns : scalar
 Args    : scalar

=cut

sub fld {
    my $self = shift;
    Bio::Root::Root->throw("requires type Q (atom)") unless ref($self) && $self->isa('Q');
    my $f = shift;
    if ($f) {
	$f =~ s/^\s+//;
	$f =~ s/\s+$//;
	return $self->{fld}=$f;
    }
    return $self->{fld};
}


=head4 Q dta

 Title   : dta
 Usage   : $q->dta($data)
 Function: get/set dta (whsp-separated data string) property
 Example :
 Returns : scalar
 Args    : scalar

=cut

sub dta {
    my $self = shift;
    Bio::Root::Root->throw("requires type Q (atom)") unless ref($self) && $self->isa('Q');
    my $d = join(" ", @_);
    if ($d) {
	$d =~ s/^\s+//;
	$d =~ s/\s+$//;
	return $self->{dta} = $d;
    }
    return $self->{dta};
}

=head4 Q A

 Title   : A
 Usage   : print $q->A
 Function: get a string representation of class Q object
 Example :
 Returns : string scalar
 Args    :

=cut

sub A {
    my $self = shift;
    Bio::Root::Root->throw("requires type Q (atom)") unless ref($self) && $self->isa('Q');
    my @a = split(/\s+/, $self->dta);

    return "(".join(' ', sort {$a cmp $b} @a).")[".$self->fld."]";
}

=head4 Q clone

 Title   : clone
 Usage   : $q2 = $q1->clone;
 Function: create and return a clone of the object
 Example :
 Returns : object of class Q
 Args    :

=cut

sub clone {
    my $self = shift;
    Bio::Root::Root->throw("requires type Q (atom)") unless ref($self) && $self->isa('Q');
    my $ret = Q->new($self->fld, $self->dta);
    return $ret;
}

### Q class methods

=head3 Q CLASS METHODS

=head4 Q qin

 Title   : qin
 Usage   : Q::qin($q1, $q2)
 Function: tests whether the query represented by $q1 would return a subset
           of items returned by the query represented by $q2
 Example : print "q2 gets those and more" if Q::qin($q1, $q2);
 Returns : 1 if q1 is subset of q2, 0 otherwise
 Args    : two class Q objects

=cut

sub qin {
    my ($a, $b) = @_;
    Bio::Root::Root->throw('requires type Q (atom)') unless (ref $a) && $a->isa('Q') && (ref $b) && $b->isa('Q');
    return 0 unless $a->fld eq $b->fld;
    return Q::qeq( $b, Q::qor($a, $b) );
}

=head4 Q qeq

 Title   : qeq
 Usage   : Q::qeq($q1, $q2)
 Function: test if fld and dta properties in two class Q objects are the same
           (irrespective of order)
 Example :
 Returns : 1 if equal, 0 otherwise
 Args    : two class Q objects

=cut

sub qeq {
    local $_;
    my ($a, $b) = @_;
    Bio::Root::Root->throw('requires type Q (atom)') unless (ref $a) && $a->isa('Q') && (ref $b) && $b->isa('Q');
    return 0 unless $a->fld eq $b->fld;
    my @ad = unique(split(/\s+/,$a->dta));
    my @bd = unique(split(/\s+/,$b->dta));
    return 0 unless @ad==@bd;
    my @cd = grep {defined} map {my $f = $_; grep /^$f$/, @ad} @bd;
    return @cd == @bd;
}

=head4 Q qor

 Title   : qor
 Usage   : @qresult = Q::qor($q1, $q2)
 Function: logical OR for Q objects
 Example :
 Returns : an array of class Q objects
 Args    : two class Q objects

=cut

sub qor {
    local $_;
    my @a = @_;
    foreach (@a) {
	Bio::Root::Root->throw("requires type Q (atom)") unless ref && $_->isa('Q');
    }
    my @ret;
    my (%f, @f);
    @a = grep {!$_->isnull} @a;
    return ($Q::NULL) unless @a > 0;
    # list of unique flds
    @f = unique(map {$_->fld} @a);
    foreach my $f (@f) {
	my @fobjs =  grep {$_->fld eq $f} @a;
	my @d = unique(map {split(/\s/, $_->dta)} @fobjs );
        my $r = Q->new($f, @d);
	push @ret, $r;
    }
    return @ret;
}

=head4 Q qand

 Title   : qand
 Usage   : @qresult = Q::And($q1, $q2)
 Function: logical AND for R objects
 Example :
 Returns : an array of class Q objects
 Args    : two class Q objects

=cut

sub qand {
    local $_;
    my ($a, $b) = @_;
    Bio::Root::Root->throw('requires type Q (atom)') unless (ref $a) && $a->isa('Q') && (ref $b) && $b->isa('Q');
    my @ret;
    if (ref $a eq 'ARRAY') {
	foreach my $ea (@$a) {
	    push @ret, qand( $ea, $b );
	}
	return qor(@ret); # simplify
    }
    elsif (ref $b eq 'ARRAY') {
	foreach my $eb (@$b) {
	    push @ret, qand( $a, $eb);
	    1;
	}
	return qor(@ret); # simplify
    }
    else {
	return ($Q::NULL) if ($a->isnull || $b->isnull);
	if ($a->fld eq $b->fld) {
	    # find intersection of data
	    my (%ad, @ad, @bd);
	    @ad = split(/\s+/, $a->dta);
	    @ad{@ad} = (1) x @ad;
	    @bd = split(/\s+/, $b->dta);
	    foreach (@bd) {
		$ad{$_}++;
	    }
	    my $r = Q->new($a->fld,
			  grep {$_}
			  map {$ad{$_} == 2 ? $_ : undef} keys %ad);
	    return (length($r->dta) > 0) ? ($r) : ($Q::NULL);
	}
	else {
	    return ($a, $b);
	}
    }
}

=head3 Q INTERNALS

=head4 Q unique

 Title   : unique
 Usage   : @ua = unique(@a)
 Function: return contents of @a with duplicates removed
 Example :
 Returns :
 Args    : an array

=cut

sub unique {
    my @a = @_;
    my %a;
    @a{@a} = undef;
    return keys %a;
}

1;

=head2 Additional tools for Bio::AnnotationCollectionI

=head3 Bio::AnnotationCollectionI SYNOPSIS (additional methods)

    $seq->annotation->put_value('patient_id', 1401)
    $seq->annotation->get_value('patient_ids')                   # returns 1401
    $seq->annotation->put_value('patient_group', 'MassGenH')
    $seq->annotation->put_value(['clinical', 'cd4count'], 503);
    $seq->annotation->put_value(['clinical', 'virus_load'], 150805);
    foreach ( qw( cd4count virus_load ) ) {
        $blood_readings{$_} = $seq->annonation->get_value(['clinical', $_]);
    }

=head3 Bio::AnnotationCollectionI DESCRIPTION (additional methods)

C<get_value()> and C<put_value> allow easy creation of and access to an
annotation collection tree with nodes of L<Bio::Annotation::SimpleValue>. These
methods obiviate direct accession of the SimpleValue objects.

=cut

package Bio::AnnotationCollectionI;
use strict;
use Bio::Annotation::SimpleValue;

=head2 get_value

 Title   : get_value
 Usage   : $ac->get_value($tagname) -or-
           $ac->get_value( $tag_level1, $tag_level2,... )
 Function: access the annotation value assocated with the given tags
 Example :
 Returns : a scalar
 Args    : an array of tagnames that descend into the annotation tree

=cut

sub get_value {
    local $_;
    my $self = shift;
    my @args = @_;
    my @h;
    return "" unless @_;
    while ($_ = shift @args) {
	@h = $self->get_Annotations($_);
	if (ref($h[0]->{value})) {
	    $self = $h[0]->{value}; # must be another Bio::AnnotationCollectionI
	}
	else {
	    last;
	}
    }
	return $h[0] && $h[0]->{value} ; # now the last value.
}

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
 Note    : If intervening nodes do not exist, put_value creates them, replacing
           existing nodes. So if $ac->put_value('x', 10) was done, then later,
           $ac->put_value(['x', 'y'], 20), the original value of 'x' is trashed,
           and $ac->get_value('x') will now return the annotation collection
           with tagname 'y'.

=cut

sub put_value {
    local $_;
    my $self = shift;
    my @args = @_;
    my ($keys, $value) = $self->_rearrange([qw( KEYS VALUE )], @args);
    my (@keys, $lastkey);
#    $value ||= new Bio::Annotation::Collection;
    @keys = (ref($keys) eq 'ARRAY') ? @$keys : ($keys);
    $lastkey = pop @keys;
    foreach (@keys) {
	my $a = $self->get_value($_);
	if (ref($a) && $a->isa('Bio::Annotation::Collection')) {
	    $self = $a;
	}
	else {
	    # replace an old value
	    $self->remove_Annotations($_) if $a;
	    my $ac = Bio::Annotation::Collection->new();
	    $self->add_Annotation(Bio::Annotation::SimpleValue->new(
				      -tagname => $_,
				      -value => $ac
				  )
		);
	    $self = $ac;
	}
    }
    if ($self->get_value($lastkey)) {
	# replace existing value
	($self->get_Annotations($lastkey))[0]->{value} = $value;
    }
    else {
	$self->add_Annotation(Bio::Annotation::SimpleValue->new(
				  -tagname=>$lastkey,
				  -value=>$value
			      ));
    }
    return $value;
}

=head2 get_keys

 Title   : get_keys
 Usage   : $ac->get_keys($tagname_level_1, $tagname_level_2,...)
 Function: Get an array of tagnames underneath the named tag nodes
 Example : # prints the values of the members of Category 1...
           print map { $ac->get_value($_) } $ac->get_keys('Category 1') ;
 Returns : array of tagnames or empty list if the arguments represent a leaf
 Args    : [array of] tagname[s]

=cut

sub get_keys {
    my $self = shift;
    my @keys = @_;
    foreach (@keys) {
	my $a = $self->get_value($_);
	if (ref($a) && $a->isa('Bio::Annotation::Collection')) {
	    $self = $a;
	}
	else {
	    return ();
	}
    }
    return $self->get_all_annotation_keys();
}

1;
