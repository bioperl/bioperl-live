#############################################################################
#
# BioPerl module for Bio::SeqFeature
#
# Cared for by Ian Korf <ikorf@sapiens.wustl.edu>
#
# Copyright Ian Korf
#
# You may distribute this module under the same terms as perl itself
#
##############################################################################

=head1 NAME

Bio::SeqFeature

=head1 SYNOPSIS

 use Bio::SeqFeature;
 use Bio::SeqFeatureParser qw(
	parseRepeatMasker
	parsMaskedSequence
	parseGenscan
	parseBlast2
	SortFeatures
	SelectFeatures
	ClusterFeatures
	); # class methods are importable except where noted

=head1 DESCRIPTION

B<Bio::SeqFeature> objects describe parts of a sequence. For example, in a
genomic sequence, you may have genes. A gene is a sequence feature, as are
the exons that make up the gene. In a protein sequence, features may be
active sites or regions of secondary structure. In general, a feature is
bound to some coordinates on the sequence.

=head2 Features

There are two types of features: simple and complex. Simple features must
contain a B<begin> and B<end>. Complex features contain other features, from
which B<begin>, B<end>, and other attributes are computed. A gene is a complex
feature that contains exons that are simple features. All features have
some common attribues. These include, B<begin>, B<end>, B<length>, and B<span>.
The difference between B<length> and B<span> is that B<span> is the overall
distance from B<begin> to B<end> whereas B<length> is the sum of the B<span>s
of the individual subfeatures (B<span> and B<length> are identical in
simple features).

Features derived from nucleic acid sequences have a B<strand>. This is
indicated by B<begin> and B<end>. Plus-strand features have B<begin> < B<end>
and minus-strand features have B<end> < B<begin>. Amino acid features do not
have strandedness.

=head2 Simple Features

To create a simple feature, you assign B<begin> and B<end> using the B<new>
constructor. You can also set other optional attributes (discussed in 
detail below).

 # an amino acid feature (maybe an alpha-helical region)
 my $protein_feature = new SeqFeature(-begin => 20, -end => 25);

 # a nucleic acid feature on the complement strand
 my $exon_feature = new SeqFeature(-begin => 30, -end => 1);

 # a strandless nucleic acid feature on the plus strand (eg. enhancer)
 my $enhancer = new SeqFeature(-begin => 100, -end => 120, -strand => 0);

=head2 Complex Features

To create a complex feature, you assign it some B<subfeature>s using the B<new>
constructor. The B<subfeature>s may be B<ordered>, like exons in a gene, or
just a mixed bag. The default is for B<subfeature>s to be B<ordered>.
You may want to assign other attributes to complex features (see below). Here
are some examples of creating complex features.

 # making a gene (first make exons, then the gene)
 my $exon1 = new Bio::SeqFeature(-begin => 1,   -end => 100);
 my $exon2 = new Bio::SeqFeature(-begin => 201, -end => 300);
 my $gene  = new SeqFeature(-subfeature => [$exon1, $exon2]);

 # example of setting the subfeatures as not ordered
 my $thing = new SeqFeature(-subfeature => [$exon1, $exon2], -ordered => 0);

=head2 Optional Attributes

Features by themselves may not be very interesting. To add some flavor,
there are several standard optional attributes. You may also make up your
own. Non-standard attributes are not guaranteed to work well with complex
features at this time.

Attributes can be assigned either as -attribute => value in the B<new>
constructor or at any other time with $feature->attribute($value).
You can access any attribute with $feature->attribute. For example:

 # assignment of the non-standard attribute 'foo' in new
 my $feature = new Bio::SeqFeature(-begin => 1, -end => 2, -foo => bar);

 # assignment of the non-standard attribute 'hello' later
 $feature->hello('world');

Common optional attributes include B<name>, B<score>, B<P>, and B<percent>.
Some features have B<name>s, others do not. Low-complexity regions do not call
for a B<name>, but a similarity to an Alu repeat probably should. Features are
often the result of processing the sequence with some algorithm. Algorithms
often output B<score>, B<P> value, or B<percent> similarity/identity figures.
These attributes may not mean the same thing from one algorithm to the next,
but the names are commonly used.

Often, features are regions that are similar to some other sequence. To link
up two sequences, you use the B<match> attribute. The following statement,
$feature->match->name, does what you think (gets you the B<name> of the
B<match>ing sequence). To construct such links, you must create the
B<match>ing feature and then link it up.
You may also like to store a B<sequence> corresponding to the
feature and the matching feature (especially if there are gaps present).
Here is an example of how you might create a BLAST similarity:

 # assume that some variables have been defined during parsing

 # first create the match (sbjct feature)
 my $blast_hsp = new Bio::SeqFeature(-begin => $sbjct_begin,
	-end => $sbjct_end, -name => $sbjct_def, -sequence => $sbjct_seq);

 # now create the query feature
 my $feature = new Bio::SeqFeature(-begin => $query_begin,
	-end => $query_end, -score => $score, -P => $P, -percent => $percent
	-sequence => $query_seq, -match => $blast_hsp);

There are cases when you want to store information that runs along the
feature. The B<metric> is a general attribute that measures some value such
as hydrophobicity or base quality value. A sequence may have more than one
B<metric>, in which case you should make up your own attributes.

=head2 Methods

Common methods include B<isSimple>, B<isComplex>, B<query>, and B<cmpFeature>.
B<isSimple> and B<isComplex> let you determine if a feature is simple or
complex. B<query> lets you determine if a feature has any combination of
attributes. B<cmpFeature> is used to compare one feature to another.
There are also methods for complex features: B<minimum>, B<maximum>,
B<sum>, B<product>, B<average>, and B<cat>. These methods traverse the
B<subfeature>s and operate on some attribute. You can also use these methods on
simple features, but they just return the attribute. Examples follow:

 $feature->isSimple;  # returns true if the feature is simple
 $feature->isComplex; # returns true if the feature is complex

 $feature->query('length > 50');   # returns true if length is greater than 50
 $feature->query('name =~ /Alu/'); # returns true if name matches Alu

 # here's an example that returns true if the features are nearly identical
 # see the appendix for more info on cmpFeature
 $feature->cmpFeature(-subject => $some_other_feature, -mode => 'identical',
	-boundary => 20);

 $feature->minimum('score'); # returns the minimum subfeature score
 $feature->maximum('score'); # returns the maximum subfeature score
 $feature->sum('score');     # returns the sum of the subfeature scores
 $feature->product('score'); # returns the product of the subfeature scores
 $feature->average('score'); # returns the average subfeature score
 $feature->cat('sequence');  # returns subfeature sequences concatenated

 # a rather complicated query that you might use in filtering BLASTN est hits
 # to get only the really good ones
 $feature->query(
	'P < 1e-100 or ( average("percent") > 95 and length > 30 )'
	);

=head2 Tags

Features may have B<tags> assoicated with them that identify them as being part
of some more general class. For example, a gene prediction program may have
the B<tags> "gene" and "prediction" associated with its gene features and the
B<tags> "exon" and "prediction" associated with its exon features. The
B<hasTag> method is useful for determining if a tag is present in a feature.

=head2 Tag Conventions

Although tag names can be whatever you like. There are certain conventions
which you should respect to make your ideas as clear as possible. See the
next section on DDJB/EMBL/GenBank for more information. Here are some of
the more common tags (or feature key identifiers) for nucleic acid sequence.

 Tag Name         Meaning
 --------         -------
 variation        some difference from something else
 promoter         proximal transcriptional signal
 enhancer         distal transcriptional signal
 RBS              ribosome binding site
 polyA_signal     AATAAA
 prim_transcript  the primary transcript
 mRNA             the mature transcript
 5'UTR            5' untranslated region
 3'UTR            3' untranslated region
 exon             the expressed part of the primary transcript
 intron           the spliced-out part of the primary transcript
 CDS              the coding part of the primary transcript
 tRNA             a transfer RNA sequence
 gene             a named region of biological interest
 repeat_region    some repetitive sequence
 STS              sequence tagged site (uniquely amplifiable)
 misc_feature     catch-all

Please submit protein feature conventions if you are knowledgable.


=head2 Comparison to DDJB/EMBL/GenBank Feature Table

Because the DDJB/EMBL/GenBank feature description is widely used, it is
important to note the similarities and differences between Bio::SeqFeature
and DDJB/EMBL/GenBank (I an just going to call this GenBank). One important
difference is that a GenBank feature has just one identifier, while
Bio::SeqFeatures can have any number.

A GenBank feature has three parts: an identifier, a location, and optional
qualifiers. Here is an example of some EMBL features that annotate a
prototypical eukaryotic gene.

      IDENTIFIER     LOCATION and /QUALIFIERS
      ----------     ------------------------

      5'UTR           100..200
      exon            100..300
                      /number=1
      intron          301..400
                      /number=1
      exon            401..600
                      /number=2
      intron          601..700
                      /number=2
      exon            701..800
                      /number=3
      intron          801..900
                      /number=3
      exon            901..1100
                      /number=4
      3'UTR           1002..1100
      sig_peptide     join(201..300, 401..501)
      mat_peptide     join(502..600, 701..800, 901..1001)
                      /product="prototypical protein"
      gene            100..1100
                      /gene="prototypical gene"
      CDS             join(201..300, 401..600, 701..800, 901..1001)
                      /product="prototypical protein"
      mRNA            join(100..300, 401..600, 701..800, 901..1100)
      prim_transcript 100..1100

Here is one way you can express much of the same information in
B<Bio::SeqFeature> objects.

 # these are the exons that make up the mRNA/gene
 my $exon_tags = ['exon';];
 my $exon1 = new Bio::SeqFeature(-begin => 100, -end => 300,
	-tags => $exon_tags);
 my $exon2 = new Bio::SeqFeature(-begin => 401, -end => 600
	-tags => $exon_tags);
 my $exon3 = new Bio::SeqFeature(-begin => 701, -end => 800
	-tags => $exon_tags);
 my $exon4 = new Bio::SeqFeature(-begin => 901, -end => 1100
	-tags => $exon_tags);

 # create a mRNA/gene
 my $gene  = new Bio::SeqFeature(-name => "prototypical gene",
	-subfeature => [$exon1, $exon2, $exon3, $exon4],
	-tags => ['gene', 'mRNA', 'prim_transcript');

I have not shown how to make all the features, but the process is the same.
Here is how you access the information.

 # prim_transcript/gene coordinates
 print $gene->begin, "..", $gene->end, "\n";

 # prim_transcript/gene name
 print $gene->name, "\n";

 # the coordinates of each exon in the gene
 foreach my $exon (@{$gene->subfeature}) {
	print $exon->begin, "..", $exon->end, "\n";
 }

 # look at some tags
 print "yup\n" if $gene->hasTag('gene', 'mRNA');


=head2 Class Methods

Several class methods are available to create and operate on lists of
B<Bio::SeqFeature> objects. Lists can be created using the parsers for WU-BLAST,
Genscan, RepeatMasker, and masked sequence. More parsers will be added in the
future (exactly when your favorite program is supported is partly dependent upon
need, so if you require a particular parser, speak up).

Once created, objects in the lists can be sorted, selected on various criteria,
or clustered. Importantly, these steps can occur during parsing. Because the
parsers take filehandles as arguments, you can do the following:

 my $handle = new FileHandle;
 $handle->open("blast... |");
 my $list = parseBlast2(
 	-handle         => $handle,
	-sbjctfilter    => "score > 200",
	-hspfilter      => "score > 50",
	-clustermode    => "identical",
	-clusterbound   => 10,
	-clusterattr    => "match->name");

What you get back from these three statements is the unique BLAST hits (within
10 residues of slop) whose HSP scores are greater than 50, and sum up to more
than 200. The unique hits also contain the names of other matching hits.


=head1 CONTACT

If you wish to make comments, bug reports, or just say hello, contact Ian Korf.
email: ikorf@sapiens.wustl.edu
web:   http://sapiens.wustl.edu/~ikorf


=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with an _. Attributes are accessed
by $object->attribute and may be assigned by $object->attribute($value). There
are several class methods described at the end of the appendix. These all
begin with uppercase letters and are accessed using the explicit package name
(eg. Bio::SeqFeature::Method). Most class methods are importable.

Throughout the documenation, there are various code examples. I have used
the following convention: the $exon object is a prototypical simple feature
and the $gene object is a prototypical complex feature. These are "created"
in the new constructor documentation and followed in most of the examples
throughout the appendix.

=cut


##############################################################################
# Code begins here
##############################################################################

package Bio::SeqFeature;
require 5.004;
use vars qw(@ISA @EXPORT_OK $VERSION $AUTOLOAD);
use strict;
use constant TRUE  => 1;
use constant FALSE => 0;

# Object preamble - inheriets from Bio::Root::Object
#use Bio::Root::Object;


@ISA = qw(Exporter);
@EXPORT_OK = qw(
	SortFeatures
	SelectFeatures
	ClusterFeatures
	parseRepeatMasker
	parseMaskedSequence
	parseGenscan
	parseBlast2);
$VERSION = 0.1;

#-----------------------------------------------------------------------------
# AUTOLOAD - to set/get unknown objects (from Ewan Birney)
#-----------------------------------------------------------------------------
sub AUTOLOAD {
	my $self = shift;
	my $name = $AUTOLOAD;
	$name =~ /::DESTROY/ && return;
	$name =~ s/.*://;
	if (@_) {
		return $self->{$name} = shift;
	}
	else {
		return $self->{$name};
	}
}

#-----------------------------------------------------------------------------
# Class Variables
#-----------------------------------------------------------------------------
# patterns used in query
my $FIELD  = '\S+';             # field separator
my $NOP    = '<|>|<=|>=|==';    # number operators
my $TOP    = 'eq|ne|gt|lt';     # text operators
my $POP    = '=~|!~';           # pattern operators
my $TEXT   = '\'[\w\s]+\'';     # 'text' must be surrounded by single quotes
my $MATCH  = '/\S+/';           # /match/ must use slashes and no spaces
my $NUMBER = '[\de\.\-\+]+';    # number may be scientific notation
my $ARG    = '\S+';             # an argument passed to a function


#-----------------------------------------------------------------------------
# Warnings
#-----------------------------------------------------------------------------
my $complex_warning = '
ERROR: you have attempted to assign a value to the attribute of a complex
feature. The attributes of complex features, features that contain other
features (subfeatures), are computed from the subfeatures. For example, a
gene feature contains exon subfeatures, and the length of the gene is the
sum of the exons.
';

my $assignment_impossible = '
ERROR: you cannot assign to this attribute.
';


#-----------------------------------------------------------------------------
# new constructor
#-----------------------------------------------------------------------------

=head2 new

 Function:  creates a new sequence feature.
 Arguments: either a begin and end or a subfeature. Others optional.
            You may assign any additional attributes using the syntax:
            -attribute => value.
 Returns:   Bio::SeqFeature object
 Examples:

 #----- a minus-strand feature
 $exon = new Bio::SeqFeature(
 	-begin => 100,
 	-end   => 1,
 	-name  => 'example exon',
 	-score => 10,
 	-tags  => ['exon']);

 #----- a two exon gene on the plus strand
 $exon1 = new Bio::SeqFeature(
 	-begin => 1,
 	-end   => 100,
 	-name  => 'example exon1',
 	-score => 20,
 	-tags  => ['exon']);
 $exon2 = new Bio::SeqFeature(
 	-begin => 201,
 	-end   => 300,
 	-name  => 'example exon2',
 	-score => 30,
 	-tags  => ['exon']);
 $gene  = new Bio::SeqFeature(
 	-name       => 'example gene',
 	-subfeature => [$exon1, $exon2],
	-ordered    => 1,
	-tags       => ['mRNA']);
 
 #----- a metric feature for 6 residues of a sequence
 $metric = new Bio::SeqFeature(
 	-begin  => 51,
 	-end    => 56,
	-metric => [10, 20, 14, 18, 27, 35],
	-strand => 0,
	-tags   => ['quality']);

=cut

sub new {
	my ($class, %p) = @_;
	my $f = bless {}; # f is for feature

	foreach my $arg (keys %p) {
		my ($key) = $arg =~ /^\-(.+)/;
		$f->{$key} = $p{$arg};
	}
	return $f;
}

#-----------------------------------------------------------------------------
# Standard Simple Feature Attributes
#-----------------------------------------------------------------------------

=head2 begin

 Function:  access/set the begin attribute.
 Arguments: optional to set the attribute.
 Returns:   the start of a feature. If the feaure is complex, the smallest
            begin is chosen from among the subfeatures (reverse is true on
            complement strand).
 Examples:

 print $exon->begin; # prints "100"
 print $gene->begin; # prints "1", the begin of the first subfeature

=cut

sub begin {
	my ($this, $value) = @_;
	if (defined $value) {
		if ($this->isComplex) {
			warn $complex_warning if $this->isComplex;
			return;
		}
		$this->{'begin'} = $value;
	}
	else {
		if ($this->isComplex) {
			if ($this->strand > 0)
				{return $this->minimum('begin')}
			elsif ($this->strand < 0)
				{return $this->maximum('begin')}
		}
		else {return $this->{'begin'}}
	}
}

#-----------------------------------------------------------------------------

=head2 end

 Function:  access/set the end attribute.
 Arguments: optional to set the attribute.
 Returns:   the end of a feature. If the feature is complex, the greatest
            end is chosen from among the subfeatures (reverse is true on
            complement strand).
 Examples:

 print $exon->end; # prints "1"
 print $gene->end; # prints "300", the end of the last subfeature

=cut

sub end {
	my ($this, $value) = @_;
	if (defined $value) {
		if ($this->isComplex) {
			warn $complex_warning if $this->isComplex;
			return;
		}
	}
	else {
		if ($this->isComplex) {
			if ($this->strand > 0)
				{return $this->maximum('end')}
			elsif ($this->strand() < 0)
				{return $this->minimum('end')}
		}
		else {return $this->{'end'}}
	}
}

#-----------------------------------------------------------------------------

=head2 subfeature

 Function:  links a feature to a set of subfeatures. Some subfeatures have a
            defined order (eg. exons of a gene), while others might not. See
            the ordered attribute below.
 Arguments: optional for assignment.
 Returns:   reference to a list of features.
 Examples:

 foreach my $exon (@{$gene->subfeature}) {
 	print $exon->begin;
 }

=cut

sub subfeature {
	my ($this, $value) = @_;
	if ($this->isComplex) {
		if (defined $value) {$this->{'subfeature'} = $value}
		else                {return $this->{'subfeature'}}
	}
	else {
		if (defined $value) {$this->{'subfeature'} = $value}
		else                {return $this->{'subfeature'}}
	}
}

#-----------------------------------------------------------------------------

=head2 ordered

 Function:  determines whether subfeatures are ordered.
 Arguments: optional for assignment.
 Returns:   1 or 0 (TRUE or FALSE).
 Examples:

=cut

sub ordered {
	my ($this, $value) = @_;
	if (defined $value) {$this->{ordered} = $value}
	else {
		if (defined $this->{ordered}) {
			if ($this->{ordered}) {return 1}
			else                  {return 0}
		}
		else {return 1}
	}
}

#-----------------------------------------------------------------------------

=head2 strand

 Function:  access/set the strand attribute.
 Arguments: optional for assignment to 0 (strandless).
 Returns:   strand of feature as 1, 0, or -1. If the feature is complex,
            the strand of the first subfeature is used. This could cause
            some misunderstanding if you have a mixture of plus and minus
            strand subfeatures (but who would do that).
 Examples:

 # canonical conditional statement with strands: prints "negative feature"
 if    ($exon->strand > 0) {print "positive feature"}
 elsif ($exon->strand < 0) {print "negative feature"}
 else                      {print "strandless feature"}

=cut

sub strand {
	my ($this, $value) = @_;
	if (defined $value) {
		warn $assignment_impossible if $value;
		$this->{'strand'} = 0;
		return;
	}

	if ($this->isComplex) {
		return $this->subfeature->[0]->strand;
	}
	else {
		return $this->{'strand'} if defined $this->{'strand'};
		return $this->begin > $this->end ? -1 : 1
	}
}

#-----------------------------------------------------------------------------

=head2 length

 Function:  access the length attribute.
 Arguments: none, this is a read only attribute.
 Returns:   length of feature or sum of lengths of subfeatures if the
            feature is complex.
 Examples:

 print $exon->length; # prints "100", the length of the feature
 print $gene->length; # prints "200", the sum of the subfeatures

=cut

sub length {
	my ($this, $value) = @_;
	if (defined $value) {
		warn $assignment_impossible;
		return;
	}

	if ($this->isComplex) {return $this->sum('length')}
	else                  {return abs($this->begin - $this->end) + 1}
}

#-----------------------------------------------------------------------------

=head2 span

 Function:  access the span attribute.
 Arguments: none, this is a read only attribute.
 Returns:   span of a feature, which is identical to length except for
            complex features where it is the distance from begin to end.
 Examples:

 print $exon->span; # prints "100", the span of the exon
 print $gene->span; # prints "300", the span of the gene (primary transcript)

=cut

sub span {
	my ($this, $value) = @_;
	if (defined $value) {
		warn $assignment_impossible;
		return;
	}

	return abs($this->begin - $this->end) + 1;
}

#-----------------------------------------------------------------------------
# Optional Attributes
#-----------------------------------------------------------------------------

=head2 name

 Function:  access/set the name attribute.
 Arguments: optional for assignment.
 Returns:   name of the feature or the undefined value.
 Examples:

 print $exon->name;                  # prints "example exon"
 print $gene->name;                  # prints "example gene"
 print $gene->subfeature->[0]->name; # prints "example exon1"

=cut

sub name {
	my ($this, $value) = @_;
	if (defined $value) {$this->{'name'} = $value}
	else                {return $this->{'name'}}
}

#-----------------------------------------------------------------------------

=head2 score

 Function:  access/set the score attribute. Features frequently have a score
            associated with them. The meaning of that score is dependent
            upon its implementation.
 Arguments: optional for assignment.
 Returns:   score of the feature. If the feature is complex, by default the
			sum of the subfeature scores is returned, but you can retrieve
            various other values if you like (see examples below).
 Examples:

 print $exon->score;            # prints "10"
 print $exon->score('maximum'); # prints 10, any operator would give the same
                                # value since $exon is a simple feature

 print $gene->score('minimum'); # prints "20"
 print $gene->score('maximum'); # prints "30"
 print $gene->score('sum');     # prints "50"
 print $gene->score('product'); # prints "600"
 print $gene->score('average'); # prints "25"
 print $gene->score;            # undefined value, some operation needed for
                                # complex feature

=cut

sub score {
	my ($this, $value) = @_;
	if (defined $value) {
		#----- assignment
		if ($value =~ /^\d+$/) {
			if ($this->isComplex) {
				warn $complex_warning;
				return;
			}
			$this->{'score'} = $value;
		}
		
		#----- complex accession
		if ($this->isSimple) {
			return $this->{'score'}
		}
		else {
			if    ($value eq 'minimum') {return $this->minimum('score')}
			elsif ($value eq 'maximum') {return $this->maximum('score')}
			elsif ($value eq 'sum')     {return $this->sum('score')}
			elsif ($value eq 'product') {return $this->product('score')}
			elsif ($value eq 'average') {return $this->average('score')}
			else                        {die "unknown argument\n"}
		}
	}
	else {
		#----- simple accession
		if ($this->isComplex) {return $this->sum('score')}
		else                  {return $this->{'score'}}
	}
}

#-----------------------------------------------------------------------------

=head2 P

 Function:  access/set the P attribute. P is frequently associated with
            features and refers to some probabalistic measurement against
            some kind of model.
 Arguments: optional for assignment.
 Returns:   P of the feature. If the feature is complex, by default the
			minimum P of the subfeature is returned, but you can retrieve
            various other values if you like (see examples in score).
 Examples:  see the score attribute above.

=cut

sub P {
	my ($this, $value) = @_;
	if (defined $value) {
		#----- assignment
		if ($value =~ /^\d+$/) {
			if ($this->isComplex) {
				warn $complex_warning;
				return;
			}
			$this->{'P'} = $value;
		}
		
		#----- complex accession
		if ($this->isSimple) {
			return $this->{'P'}
		}
		else {
			if    ($value eq 'minimum') {return $this->minimum('P')}
			elsif ($value eq 'maximum') {return $this->maximum('P')}
			elsif ($value eq 'sum')     {return $this->sum('P')}
			elsif ($value eq 'product') {return $this->product('P')}
			elsif ($value eq 'average') {return $this->average('P')}
			else                        {die "unknown argument\n"}
		}
	}
	else {
		#----- simple accession
		if ($this->isComplex) {return $this->sum('P')}
		else                  {return $this->{'P'}}
	}
}

#-----------------------------------------------------------------------------

=head2 percent

 Function:  access/set the percent attribute. A percent attribute usually
            indicates percent similarity or percent identity in a pairwise
            alignment. If you use percent to mean percent similar, beware
            that it is dependent upon what scoring matrix you use.
 Arguments: optional for assignment.
 Returns:   percent of the feature. If the feature is complex, by default the
			average of the subfeature percents is returned, but you can
            retrieve various other values if you like (see examples in score).
 Examples:  see the score attribute above.

=cut

sub percent {
	my ($this, $value) = @_;
	if (defined $value) {
		#----- assignment
		if ($value =~ /^\d+$/) {
			if ($this->isComplex) {
				warn $complex_warning;
				return;
			}
			$this->{'percent'} = $value;
		}
		
		#----- complex accession
		if ($this->isSimple) {
			return $this->{'percent'}
		}
		else {
			if    ($value eq 'minimum') {return $this->minimum('percent')}
			elsif ($value eq 'maximum') {return $this->maximum('percent')}
			elsif ($value eq 'sum')     {return $this->sum('percent')}
			elsif ($value eq 'product') {return $this->product('percent')}
			elsif ($value eq 'average') {return $this->average('percent')}
			else                        {die "unknown argument\n"}
		}
	}
	else {
		#----- simple accession
		if ($this->isComplex) {return $this->average('percent')}
		else                  {return $this->{'percent'}}
	}
}

#-----------------------------------------------------------------------------

=head2 match

 Function:  used for features produced by pairwise alignment algorithms.
 Arguments: optional for assignment.
 Returns:   reference to a matching feature.
 Examples:

 $feature->match->begin;

=cut

sub sbjct {
	my ($this, $value) = @_;
	if (defined $value) {
		if ($this->isComplex) {warn $complex_warning}
		else                  {$this->{match} = $value}
	}
	else {
		if ($this->isComplex) {return $this->{match}}
		else {return $this->{match}}
	}
}

#-----------------------------------------------------------------------------

=head2 sequence

 Function:  contains the sequence of the feature. This is useful for pairwise
            alignments to store the alignment sequence, which may be
            different form the actual sequence due to gaps.
 Arguments: optional for assignment.
 Returns:   reference to an alignment sequence
 Examples:

 $feature->sequence;
 $feature->match->sequence;

=cut

sub sequence {
	my ($this, $value) = @_;
	if (defined $value) {
		if ($this->isComplex) {warn $complex_warning}
		else                  {$this->{sequence} = $value}
	}
	else {
		if ($this->isComplex) {}
		else {return $this->{sequence}}
	}
}

#-----------------------------------------------------------------------------

=head2 tags

 Function:  associate the feature with some higher order concept.
 Arguments: optional for assignment.
 Returns:   reference to list of tag names.
 Examples:

 print @{$exon->tags}; # prints "exon"
 print @{$gene->tags}; # prints "mRNA"

=cut

sub tags {
	my ($this, $tags) = @_;
	if (ref $tags eq 'ARRAY') {$this->{'tags'} = $tags}
	else                      {return $this->{'tags'}}
}

#-----------------------------------------------------------------------------

=head2 members

 Function:  identifiers a list of similar features. This attribute is used in
            clustering to keep a record of those features that were "thown
            away" during the clustering. See the ClusterFeatures class method
            below.
 Arguments: optional for assignment.
 Returns:   reference to list of members. Members may be scalar values (eg. the
            name of a another feature), or Bio::SeqFeature objects.
 Examples:

 print @{$exon->members};        # accession
 $exon->members(['foo', 'bar']); # assignment

=cut

sub members {
	my ($this, $members) = @_;
	if (ref $members eq 'ARRAY') {$this->{'members'} = $members}
	else                         {return $this->{'members'}}
}

#-----------------------------------------------------------------------------
# Simple Methods
#-----------------------------------------------------------------------------

=head2 isComplex

 Function:  used to determine if a feature is complex (contains subfeatures).
 Arguments: none.
 Returns:   1 or 0 (TRUE or FALSE).
 Examples:

 print "complex" if $exon->isComplex; # doesn't print
 print "complex" if $gene->isComplex; # prints "complex"

=cut

sub isComplex {
	my ($this) = @_;
	if (exists $this->{'subfeature'}) {return 1}
	else                              {return 0}
}

#-----------------------------------------------------------------------------

=head2 isSimple

 Function:  used to determine if a feature is simple (contains no subfeatures).
 Arguments: none.
 Returns:   1 or 0 (TRUE or FALSE).
 Examples:
 
 print "simple" if $exon->isSimple; # prints "simple"
 print "simple" if $gene->isSimple; # doesn't print

=cut

sub isSimple {
	my ($this) = @_;
	if    (exists $this->{'subfeature'}) {return 0}
	elsif (exists $this->{'metric'})     {return 0}
	else                                 {return 1}
}

#-----------------------------------------------------------------------------

=head2 hasTag

 Function:  used to determine if a feature has a given tag.
 Arguments: list of tags to check.
 Returns:   1 or 0 (TRUE o FALSE)
 Examples:

 print "yup" if $exon->hasTag('exon'); # prints "yup"

=cut

sub hasTag {
	my ($this, @tags) = @_;
	return unless exists $this->{'tags'} and @{$this->{'tags'}};
	foreach my $query (@tags) {
		foreach my $tag (@{$this->{'tags'}}) {
			return 1 if $query eq $tag;
		}
	}
	return 0;
}

#-----------------------------------------------------------------------------

=head2 query

 Function:  used to ask questions about a feature. The questions can be
            quite complicated. The querries are written in pseudo perl.
            Undefined or FALSE querries return TRUE. The following rules should
            be followed:
            - string constants must use single quotes
            - pattern matches must use slashes and contain no spaces
            - method arguments may not contain spaces
 Arguments: a query string.
 Returns:   1 or 0 (TRUE or FALSE).
 Examples:

 $gene->query("length > 100");
 $gene->query("begin < 100 and end > 200");
 $gene->query("average('score') > 20");
 $gene->query("hasTags('foo', 'bar')");
 $gene->query("name =~ /sapiens/");
 $gene->query("contains(150)");

=cut

sub query {
	my ($this, $filter) = @_;
	return 1 unless defined $filter;
	return 1 unless $filter =~ /\w/;

	#----- parse filtering into eval-able code (see class variables)
	$filter =~ s/($FIELD) ($NOP) ($NUMBER)/\$this->$1 $2 $3/g;
	$filter =~ s/($FIELD) ($TOP) ($TEXT)/\$this->$1 $2 $3/g;
	$filter =~ s/($FIELD) ($POP) ($MATCH)/\$this->$1 $2 $3/g;
	$filter =~ s/($FIELD)\(($ARG)\)/\$this->$1($2)/;

	my $result = eval $filter;
	warn "ERROR in SeqFeature::query\n$@" if $@;

	# do something with undefined values?
	return $result;
}

#-----------------------------------------------------------------------------

=head2 contains

 Function:  determines if a feature contains a particular coordinate.
 Arguments: a position.
 Returns:   1 or 0 (TRUE or FALSE).
 Examples:

 $exon->contains(20);

=cut

sub contains {
	my ($this, $coor) = @_;
	my ($begin, $end) = ($this->begin, $this->end);
	if (abs($begin - $coor) + abs($end - $coor) > abs($begin - $end))
		{return FALSE}
	else
		{return TRUE}
}

#-----------------------------------------------------------------------------

=head2 cmpByCoor

 Function:  compares two features by looking at their coordinates. This
            function assumes features are simple. See cmpFeature for a more
            general method.
 Arguments: a SeqFeature object, a boundary, and an optional strand-
            sensitive flag.
 Returns:   various values.
   'identical'    features are within the boundary
   'subset'       the feature is smaller on both ends
   'superset'     the feature is longer on both ends
   'longRight'    the feature is longer on the right only
   'longLeft'     the feature is longer on the left only
   'shortRight'   the feature is shorter on the right only
   'shortLeft'    the feature is shorter on the left only
   'skewRight'    the feature is longer on the right and shorter on the left
   'skewLeft'     the feature is longer on the left and shorter on the right
 Examples:

 $exon->cmpByCoor($gene, 10);
 $exon->cmpByCoor($gene, 10, 0);

=cut

sub cmpByCoor {
	my ($f1, $f2, $fuzzy, $strand) = @_;
	$strand = TRUE unless defined $strand;
	my ($begin, $end, $result);

	# strand sensitive?
	return FALSE if $f1->strand != $f2->strand
		and $strand == TRUE;

	# get the coordinates
	my ($b1, $e1) = ($f1->begin, $f1->end);
	my ($b2, $e2) = ($f2->begin, $f2->end);

	# if strand insensitive, fake positive strand
	if ($b1 > $e1) {($b1, $e1) = ($e1, $b1)}
	if ($b2 > $e2) {($b2, $e2) = ($e2, $b2)}

	# any overlap at all? (an endpoint of 2 must fall within 1)
	return FALSE unless $f1->contains($b2) or $f1->contains($e2)
		or $f2->contains($b1) or $f2->contains($e1);

	# do the comparison
	$begin = _fcmp($b1, $b2, $fuzzy);
	$end   = _fcmp($e1, $e2, $fuzzy);

	# produce the results
	if ($begin ==  1 and $end ==  1)
		{return ($f1->strand > 0) ?   'skewRight' : 'skewLeft'}
	elsif ($begin ==  1 and $end ==  0)
		{return ($f1->strand > 0) ?   'shortLeft' : 'longLeft'}
	elsif ($begin ==  1 and $end == -1)
		{return ($f1->strand > 0) ?      'subset' : 'superset'}
	elsif ($begin ==  0 and $end ==  1)
		{return ($f1->strand > 0) ?   'longRight' : 'shortRight'}
	elsif ($begin ==  0 and $end ==  0)
		{return                              'identical'}
	elsif ($begin ==  0 and $end == -1)
		{return ($f1->strand > 0) ?  'shortRight' : 'longRight'}
	elsif ($begin == -1 and $end ==  1)
		{return ($f1->strand > 0) ?    'superset' : 'subset'}
	elsif ($begin == -1 and $end ==  0)
		{return ($f1->strand > 0) ?    'longLeft' : 'shortLeft'}
	elsif ($begin == -1 and $end == -1)
		{return ($f1->strand > 0) ?    'skewLeft' : 'skewRight'}
}

#-----------------------------------------------------------------------------

=head2 cmpFeature

 Function:  compares two features by looking at their coordinates. Simple
            and complex features can be compared to each other, but the
            result of the comparison will depend upon the mode (see below).
 Arguments: cmpFeature arguments are passed in as named parameters.
            -subject     a Bio::SeqFeature object.
            -mode        "simple", "identical" or "any".
            -boundary    some integer.
            -strand      1 or 0 [default 1].
 Returns:   1 or 0 (TRUE or FALSE).

 More Detail:

            The comparison routine looks at the coordinates of the begin
            and end attributes and determines if they are within some
            boundary. No matter how high the boundary is set, the features
            must always overlap each other.

            The return values depend upon the comparison mode. In "simple"
            mode, all features are compared as simple features. That is,
            subfeatures are not compared. In "identical mode", all features,
            including subfeatures must be identical (within some boundary
            of course). In "any" mode, the features must share some region
			in common.

            By default, the comparisons are strand sensitive. This may not
            always be appropriate (eg. ESTs), so the -strand paramater may be
            set to 0 (strand insensitive).

 Examples:

 $gene->cmpFeature(
	-feature  => $gene2,
	-boundary => 10,
	-strand   => 1,
	-mode     => 'identical');

=cut

sub cmpFeature {
	my ($f1, %p) = @_;

	my $f2     = $p{-subject};
	die "cmpFeature requires subject\n" unless defined $f2;

	my $mode   = $p{-mode};
	my $bound  = $p{-boundary};
	my $strand = $p{'-strand'} ;

	#----- defaults
	$strand = TRUE    unless defined $p{'-strand'};
	$bound =  1       unless defined $p{'-boundary'} and $p{-boundary} >= 0;
	$mode =  'simple' unless defined $p{'-mode'}     and $p{-mode};

	#----- check if they are on the same strand
	return 0 if $f1->strand != $f2->strand and $strand == TRUE;

	#----- complex features need to be expanded (only two levels now)
	my (@g1, @g2); # g is for gene - a simple way to think about it
	if ($f1->isComplex and $mode ne 'simple') {
		my $sf = $f1->subfeature;
		@g1 = @$sf;
	}
	else {push @g1, $f1};

	if ($f2->isComplex and $mode ne 'simple') {
		my $sf = $f2->subfeature;
		@g2 = @$sf;
	}
	else {push @g2, $f2};

	#----- if they each have one feature, life is easy
	if (@g1 == 1 and @g2 == 1) {
		my $result = $f1->cmpByCoor($f2, $bound, $strand);
		if ($mode eq 'any') {
			if ($result) {return 1}
			else         {return 0}
		}
		elsif ($mode eq 'identical') {
			if ($result eq 'identical') {return 1}
			else                        {return 0}
		}
		elsif ($mode eq 'simple') {
			if ($result eq 'identical') {return 1}
			else                        {return 0}
		}
		# more modes to be added later
	}

	#----- if identical mode then they must have the same number of exons
	return FALSE if @g1 != @g2 and $mode eq 'identical';

	#------------------------------------------------------------------
	# Find the first matching exons. In the 'any' mode, any TRUE value
	# will be taken and returned as TRUE (actually a 1).
	#------------------------------------------------------------------
	my $match = FALSE;
	my ($s1, $s2); # the starting offsets of the first matching exons
	FIND_EXON: for(my $i=0;$i<@g1;$i++) {
		$s2 = _exonInGene($bound, $mode, $g1[$i], @g2);
		if (defined $s2) {
			$s1 = $i;
			$match = TRUE;
			return TRUE if ($mode eq 'any');
			last FIND_EXON;
		}
	}
	return FALSE unless $match;
	return FALSE if $mode eq 'identical' and $s2 != 0;

	if ($mode eq 'identical') {
		my $id = TRUE;
		for(my $i=0;$i<@g1;$i++) {
			if ($g1[$i]->cmpByCoor($g2[$i], $bound) ne 'identical') {
				return FALSE;
			}
		}
	}
	return TRUE if $mode eq 'identical';
}

#-----------------------------------------------------------------------------
# exonInGene
# Returns the offset of the first "match" to a query exon in sbjct gene.
# Match is defined as a return value of 'identical'.
# This is a direction sensitive algorithm that assumes starting from
# the left end of a complex feature and moving right
#-----------------------------------------------------------------------------
sub _exonInGene {
	my ($bound, $mode, $e1, @gene) = @_;
	my $r; # the result of the search

	SEARCH: for(my $i=0;$i<@gene;$i++) {
		my $e2 = $gene[$i];
		my $r; # the result of the comparison

		$r = $e1->cmpByCoor($e2, $bound);

		return $r if $mode eq 'any' and $r;
		return $i if $r eq 'identical';
	}

	# If no matching exons were found, return the undefined value
	return undef;
}

sub _fcmp {
	my ($this, $that, $bound) = @_;
	if    (abs($this - $that) <= $bound) {return 0}
	elsif ($this < $that)                {return -1}
	else                                 {return 1}
}


#-----------------------------------------------------------------------------

=head2 minimum

 Function:  gets the minimum value of an attribute in a complex feature.
 Arguments: the name of the attribute.
 Returns:   minimum value.
 Examples:

 print $gene->minimum('score');

=cut

sub minimum {
	my ($this, $att) = @_;
	my @sf = @{$this->{'subfeature'}};
	my $min;
	foreach my $sub (@sf) {
		my $val = eval "\$sub->$att";
		warn "EVAL error in SeqFeature::minimum\n" if $@;
		next unless defined $val;
		$min = $val unless defined $min;
		$min = $val if $val < $min;
	}
	return $min;
}

#-----------------------------------------------------------------------------

=head2 maximum

 Function:  gets the maximum value of an attribute in a complex feature.
 Arguments: the name of the attribute.
 Returns:   maximum value.
 Examples:

 print $gene->maximum('score');

=cut

sub maximum {
	my ($this, $att) = @_;
	my @sf = @{$this->{'subfeature'}};
	my $max;
	foreach my $sub (@sf) {
		my $val = eval "\$sub->$att";
		warn "EVAL error in SeqFeature::maximum\n" if $@;
		next unless defined $val;
		$max = $val unless defined $max;
		$max = $val if $val > $max;
	}
	return $max;
}

#-----------------------------------------------------------------------------

=head2 sum

 Function:  gets the sum of values of an attribute in a complex feature.
 Arguments: the name of the attribute.
 Returns:   the sum of the attributes.
 Examples:

 print $gene->sum('score');

=cut

sub sum {
	my ($this, $att) = @_;
	my @sf = @{$this->{'subfeature'}};
	my $sum = 0;
	foreach my $sub (@sf) {
		my $val = eval "\$sub->$att";
		warn "EVAL error in SeqFeature::sum\n" if $@;
		next unless defined $val;
		$sum += $val;
	}
	return $sum;
}

#-----------------------------------------------------------------------------

=head2 product

 Function:  gets the product of the values of an attribute in a complex
			feature.
 Arguments: the name of the attribute.
 Returns:   the product of the attributes.
 Examples:

 print $gene->product('score');

=cut

sub product {
	my ($this, $att) = @_;
	my @sf = @{$this->{'subfeature'}};
	my $pro = 0;
	foreach my $sub (@sf) {
		my $val = eval "\$sub->$att";
		warn "EVAL error in SeqFeature::product\n" if $@;
		next unless defined $val;
		$pro *= $val;
	}
	return $pro;
}

#-----------------------------------------------------------------------------

=head2 average

 Function:  gets the average value of an attribute in a complex feature.
 Arguments: the name of the attribute.
 Returns:   the average of the attributes.
 Examples:

 print $gene->average('score');

=cut

sub average {
	my ($this, $att) = @_;
	my @sf = @{$this->{'subfeature'}};
	my $sum = 0;
	foreach my $sub (@sf) {
		my $val = eval "\$sub->$att";
		warn "EVAL error in SeqFeature::average\n" if $@;
		next unless defined $val;
		$sum += $val;
	}
	return $sum / @sf;
}

#-----------------------------------------------------------------------------

=head2 cat

 Function:  concatenates the values of an attribute in a complex feature.
 Arguments: the name of the attribute.
 Returns:   concatenated attributes.
 Examples:

 print $gene->cat('sequence');

=cut

sub cat {
	my ($this, $att) = @_;
	my @sf = @{$this->{'subfeature'}};
	my $cat = '';
	foreach my $sub (@sf) {
		my $val = eval "\$sub->$att";
		warn "EVAL error in SeqFeature::average\n" if $@;
		next unless defined $val;
		$cat .= $val;
	}
	return $cat;
}

#-----------------------------------------------------------------------------
# Class Methods
#-----------------------------------------------------------------------------

=head2 Class Methods

B<Bio::SeqFeature> defines several class methods for dealing with lists of
B<Bio::SeqFeature> objects. These are explained below. In the examples, $list
is a reference to a list of B<Bio::SeqFeature> objects;

=cut


#-----------------------------------------------------------------------------
# RepeatMasker
#-----------------------------------------------------------------------------

=head2 parseRepeatMasker

 Function:  parses RepeatMasker files into Bio::SeqFeature objects.
 Arguments: passed as key => value pairs (see parameters).
 Returns:   a reference to a list of SeqFeature objects.
 Parameters:
   -handle      an open object-oriented filehandle or pipe
   -filter      [optional] filtering statement for subjects
   -tags        a reference to a list of tags to assoicate with each feature
 Examples:

 $list = parseRepeatMasker(-handle => $filhandle);
 $list = parseRepeatMasker(-tags => ['repeat']);
 $list = parseRepeatMasker(-handle => $filhandle, -filter => 'length > 50');

=cut

sub parseRepeatMasker {
	my (%p) = @_;
	my ($file, $filter, $tags) = ($p{'-handle'}, $p{'-filter'}, $p{'-tags'});
	my @bag;

	# column number that contains begin, end, score, and name
	# these numbers are for RepeatMakser's *.out file
	my ($begin, $end, $score, $name) = (5,6,0,9);

	<$file>;<$file>;<$file>; # strip off header
	while(<$file>) {
		my @field = split;

		#----- create a feature
		my $feature = new Bio::SeqFeature(
			'-name'  => $field[$name],
			'-begin' => $field[$begin],
			'-end'   => $field[$end],
			'-score' => $field[$score],
			'-tags'  => $tags);

		#----- put it in the bag
		push @bag, $feature if $feature->query($filter);
	}
	return \@bag;
}


#-----------------------------------------------------------------------------
# MaskedSequence
#-----------------------------------------------------------------------------

=head2 parseMaskedSequence

 Function:  parses masked sequence into Bio::SeqFeature objects.
 Arguments: passed as key => value pairs (see parameters).
 Returns:   a reference to a list of Bio::SeqFeature objects.
 Parameters:
   -handle      an open object-oriented filehandle or pipe
   -maskchar    character used in the masking
   -filter      [optional] filtering statement for subjects
   -tags        a reference to a list of tags to assoicate with each feature
 Examples:

 $list = parseMaskedSequence(-handle => $filehandle, -maskchar => 'N');
 $list = parseMaskedSequence(-handle => $filehandle, -maskchar => 'X');
 $list = parseMaskedSequence(-handle => $filehandle, -maskchar => 'N',
 	-filter =>'length > 50', -tags => ['repeat', '>50bp']);

=cut

sub parseMaskedSequence {
	my (%p) = @_;
	my ($file, $maskchar, $filter, $tags) = ($p{'-handle'}, $p{'-maskchar'},
		$p{'-filter'}, $p{'-tags'});

	#----- read in the FASTA file
	my ($def, @seq) = <$file>;
	chomp @seq;
	my $seq = join('', @seq);

	my @bag; # holds the features

	my ($i, $start, $stop);
	SEGMENT: for($i=0;$i<length($seq)-1;$i++) {
		$start = index($seq, $maskchar, $i);
		last SEGMENT if $start == -1;
		$i = $start;
		while (substr($seq, $i, 1) eq $maskchar) {
			$stop = $i;
			$i++;
		}

        #----- create the feature
        my $feature = new Bio::SeqFeature(
			'-begin'  => $start,
			'-end'    => $stop,
			'-tags'   => $tags,
			'-strand' => 0);

		#----- put it in the bag if it's query returns true
		push @bag, $feature if $feature->query($filter);

		redo SEGMENT;
	}
	return \@bag;
}


#-----------------------------------------------------------------------------
# Genscan
#-----------------------------------------------------------------------------

=head2 parseGenscan

 Function:  parses Genscan output into Bio::SeqFeature objects.
 Arguments: passed as key => value pairs (see parameters).
 Returns:   a reference to a list of Bio::SeqFeature objects.
 Parameters:
   -handle      an open object-oriented filehandle or pipe
   -genefilter  [optional] filtering statement for genes
   -exonfilter  [optional] filtering statement for exons
   -genetags    a reference to a list of tags to assoicate with each gene
   -exontags    a reference to a list of tags to assoicate with each exon
 Examples:
 $list = parseGenscan(-handle => $filehandle, -genetags => ['genscan', 'gene']);
 $list = parseGenscan(-hanlde => $filehandle, -genefilter => "strand > 0");
 $list = parseGenscan(-handle => $filehandle,
 	-genefilter => "average('score') > 50");

=cut

sub parseGenscan {
	my (%p) = @_;
	my ($file, $geneFilter, $exonFilter, $genetags, $exontags) = 
		($p{'-handle'}, $p{'-genefilter'}, $p{'-exonfilter'}, $p{'-genetags'},
		$p{'-exontags'});
	while(<$file>) {last if $_ =~ /^-----/} # strip off header
	my $line = <$file>;                     # dump the blank line

	my @gene_list; # stores all the genes
	my $gene_id = 0;

	GENE: {

		my @exon_list; # exons for this gene

		EXON: {
			$line = <$file>;
			last GENE if $line =~ /^Predicted|^Suboptimal|^NO EXONS/;
			last EXON unless $line =~ /\w/;
			my @field = split(/\s+/, $line);

			#----- throw out non-exons
			redo EXON if $field[2] =~ /Prom|PlyA/;

			my $exon = new Bio::SeqFeature(
				'-name'  => "Genscan $field[2] exon $field[0]",
				'-begin' => $field[4],
				'-end'   => $field[5],
				'-P'     => $field[12],
				'-score' => $field[13],
				'-tags'  => $exontags);

			push @exon_list, $exon if $exon->query($exonFilter);
			redo EXON;
		}

		#----- make a gene from the exons
		$gene_id++;
		my $gene = new Bio::SeqFeature(
			'-name'       => "Genscan gene $gene_id",
			'-subfeature' => \@exon_list,
			'-tags'       => $genetags);

		push @gene_list, $gene if @exon_list and $gene->query($geneFilter);
		redo GENE;
	}

	return \@gene_list;
}


#-----------------------------------------------------------------------------
# WU-BLAST
#-----------------------------------------------------------------------------

=head2 parseBlast2

 Function:  parses Blast2 reports (WU-BLAST) into Bio::SeqFeature objects.
 Arguments: passed as key => value pairs (see parameters).
 Returns:   a reference to a list of Bio::SeqFeature objects.
 Parameters:
   -handle         an open object-oriented filehandle or pipe
   -sbjctfilter    [optional] filtering statement for subjects
   -hspfilter      [optional] filtering statement for high scoring pairs
   -clustermode    [optional] 'any' or 'identical'
   -clusterbound   [optional] an integer
   -clusterstrand  [optional] 1 or 0 (TRUE or FALSE)
   -clusterattr    [optional] attribute of redundancies to save
   -sbjcttags      a reference to a list of tags to assoicate with each SBJCT
   -hsptags        a reference to a list of tags to assoicate with each HSP
 Examples:
 $data = parseBlast2(-handle => $filehandle, -exontags => ['BLAST', 'HSP']);
 $data = parseBlast2(-handle => $filehandle, -sjbctfilter => 'P < 1e-30',
 	-clustermode => 'identical', -clusterbound => 10, -clusterstrand => 1,
	-clusterattr => 'match->name');

=cut

sub parseBlast2 {
	my (%p) = @_;
	my ($file, $sbjct_filter, $hsp_filter, $cmode, $cbound, $cstrand, $cattr,
		$sbjcttags, $hsptags) = ($p{'-handle'}, $p{'-sbjctfilter'},
		$p{'-hspfilter'}, $p{'-clustermode'}, $p{'-clusterbound'},
		$p{'-clusterstrand'}, $p{'-clusterattr'}, $p{'-sbjcttags'},
		$p{'-hsptags'});

	#----- some regular expressions used frequently in parsing
	my $Int = '\d+';
	my $Float = '[\de\.\-\+]+';
	my $Seq = '[A-Z\-\*\+\|]+';

	#----- blast file
	while(<$file>) {last if $_ =~ /^>/} # strip off header
	my $line = $_;  # $line is used somewhat globally

	#----- SBJCT parsing starts here
	my @sbjct; # holds all the SBJCTs in the report
	SBJCT: {
		last SBJCT if blastDone($line);
		#----- get the def line
		my $def = $line;
		DEF_LINE: {
			$line = <$file>;
			last SBJCT if blastDone($line);
			if ($line =~ /^\s+Length = \d+/) {last DEF_LINE}
			else {$def .= $line}
			redo DEF_LINE;
		}
		my ($len) = $line =~ /Length = ([\d,]+)/;
		$len =~ s/,//g;
		$def =~ s/\s+/ /g;

		#----- HSP parsing starts here
		my @hsp;   # holds all the HSPs for a SBJCT
		my @match; # holds all the individual matches
		HSP: {

			CHECK: {
				last CHECK if $line =~ /^\s+Score/;
				last SBJCT if blastDone($line);
				$line = <$file>;
				redo CHECK;
			}
			my ($score) = $line =~ /Score = ($Int)/;
			my ($p)     = $line =~ /[Sum ]*P[\(\d+\)]* = ($Float)/;
			$line = <$file>;
			my ($percent) = $line =~ /Identities = \d+\/\d+ \((\d+)%\)/;
			<$file>;

			my $first = 1;
			my ($q_begin, $q_seq, $q_end); # query begin, seq, end
			my $a_seq;                     # alignment between
			my ($s_begin, $s_seq, $s_end); # sbjct begin, seq, end
			$q_seq = $a_seq = $s_seq = '';
			my ($ql, $al, $sl);            # each line

			ALIGNMENT: {
				if ($first) {$ql = <$file>}
				else        {$ql = $line} 
				$al = <$file>; $sl = <$file>;
				my ($qb, $qs, $qe) = $ql =~ /($Int) ($Seq) ($Int)/;
				my $offset = index($ql, $qs);
				my $as = substr($al, $offset); chomp $as;
				my ($sb, $ss, $se) = $sl =~ /($Int) ($Seq) ($Int)/;
				($q_begin, $s_begin) = ($qb, $sb) if $first;
				($q_end, $s_end) = ($qe, $se);
				$q_seq .= $qs; $a_seq .= $as; $s_seq .= $ss; $first = 0;
				<$file>;
				$line = <$file>;
				redo ALIGNMENT if $line =~ /^Query/;
			}

			#----- save or filter HSP
			my $match = new Bio::SeqFeature('-begin' => $s_begin,
				'-end' => $s_end, '-sequence' => $s_seq);
			my $hsp = new Bio::SeqFeature('-begin' => $q_begin,
				'-end' => $q_end, '-score' => $score, '-P' => $p,
				'-percent' => $percent, '-sequence' => $q_seq,
				'-match' => $match, '-tags' => $hsptags);

			if ($hsp->query($hsp_filter)) {
				push @hsp, $hsp;
				push @match, $match;
			}

			#----- get another HSP
			redo HSP if $line =~ /^\s+Score|^\s+Minus Strand HSPs:/;

			#----- all HSP's might fall below filter
			unless (@hsp) {
				$line = <$file>;
				redo SBJCT;
			}

			#----- save or filter SBJCT
			my $hit = new Bio::SeqFeature('-subfeature' => \@match,
				'-name' => $def);
			my $sbjct = new Bio::SeqFeature('-subfeature' => \@hsp,
				'-ordered' => 0, '-match' => $hit, '-tags' => $sbjcttags);

			push @sbjct, $sbjct if $sbjct->query($sbjct_filter);
			
			#----- cluster
			if (defined $cmode and @sbjct > 1) {
				my $cluster = Bio::SeqFeature::ClusterFeatures(\@sbjct, $cmode,
					$cbound, $cstrand, $cattr);
				@sbjct = @{$cluster};
			}

			#----- get another SBJCT
			$line = <$file>;
			redo SBJCT;
		}
	}
	return \@sbjct;
}


sub blastDone {
	my ($line) = @_;
	return 1 unless defined $line;
	if ($line =~ /^Parameters/) {return 1}
	elsif ($line =~ /^FATAL/) {return 1}
	else {return 0}
}

#-----------------------------------------------------------------------------

=head2 SortFeatures

 Function:  sorts a list of features.
 Arguments: a reference to a list of Bio::SeqFeature objects and
            a list of attributes to sort. Each attribute is followed by a
            + or - indicating if the sort on this attribute should be by
            ascending or descending value. If no argument is given, the
            default sort (strand- begin-) is used.
 Returns:   a reference to a sorted list. The original list is not altered. If
            you wish to do so, you should reassign the original list (see the
            last example below).
 Examples:

 SortFeatures($list);                     # default
 SortFeatures($list, 'length+');          # ascending length
 SortFeatures($list, 'score-', 'begin+'); # descending score, ascending begin
 $list = SortFeatures($list);             # alter the original list

=cut

sub SortFeatures {
	my ($list, @att) = @_;
	my @sort;
	if (not @att) {
		@sort = sort _defaultSort @{$list->features};
		$list = \@sort;
		return;
	}
	my $code = "sort {\n";
	foreach my $att (@att) {
		my $mode = chop $att;
		if ($mode eq '+') {
			$code .= "  \$a->$att <=> \$b->$att or \n";
			$code .= "  \$a->$att cmp \$b->$att or \n";
		}
		elsif ($mode eq '-') {
			$code .= "  \$b->$att <=> \$a->$att or \n";
			$code .= "  \$b->$att cmp \$a->$att or \n";
		}
		else {die "ERROR: sorting must have + or -\n"}
	}
	$code = substr($code, 0, length($code) -4);
	$code .= "\n} \@{\$list};";
	#print $code, "\n";
	@sort = eval $code;
	print $@ if $@;
	return \@sort;
}
sub _defaultSort {$b->strand <=> $a->strand or $a->begin <=> $b->begin}

#-----------------------------------------------------------------------------

=head2 SelectFeatures

 Function:  selects features from a list of Bio::SeqFeature objects.
 Arguments: a filtering statement. See the query method above.
 Returns:   a reference to a list of Bio::SeqFeatures.
 Examples:

 $hits = SelectFeatures($list, 'score > 50');
 $yfg  = SelectFeatures($list, "hasTag('gene') and name =~ /actin/");
 foreach my $feature (SelectFeatures($list, "hasTag('exon')")) {do_something($f)}

=cut

sub SelectFeatures {
	my ($list, $filter) = @_;
	my @keep;
	foreach my $f (@$list) {push @keep, $f if $f->query($filter)}
	return \@keep;
}

#-----------------------------------------------------------------------------

=head2 ClusterFeatures

B<ClusterFeatures> is especially useful where there are many redundant features
and you only want the "best" ones (eg. BLAST searches). If you are using
B<Bio::SeqFeature> while parsing a file, you can call B<Cluster> repeatedly to
keep your object from getting too big. The algorithm is greedy and given two
similar features, it will just keep the first one in the list. You may want to
sort your list before clustering becase of the greediness. B<Cluster> uses
cmpFeature (above), so you should peruse that method to understand the behavior
of the algorithm.

 Function:  clusters similar features.
 Arguments: a reference to a list of features, a clustering mode, a boundary,
            a strand sensitive flag, and the name of an attribute to save. If
            no attribute is named, the object is stored (and memory usage will
            remain unchanged.
 Returns:   a clustered list.
 Examples:

 ClusterFeatures($list, $mode, $bound, $strand, $attribute);
 ClusterFeatures($list, 'identical', 10, 1, 'name');
 ClusterFeatures($list, 'any', 50, 0);

=cut

sub ClusterFeatures {
	my ($list, $mode, $bound, $strand, $att) = @_;
	my @list = @{$list};
	my @keep;
	QUERY: while(@list) {
		my $query = shift @list;
		my $i = -1; # splice offset
		my $unique = 1;
		SBJCT: foreach my $sbjct (@list) {
			$i++;
			if ($query->cmpFeature('-subject' => $sbjct, '-mode' => $mode,
				'-bound' => $bound, '-strand' => $strand)) {
				$unique = 0;
				if (not defined $att) {push @{$query->{'members'}}, $sbjct}
				else {push @{$query->{'members'}}, eval "\$sbjct->$att"}
				splice @list, $i, 1;   # remove the clustered feature
				unshift @list, $query; # put the query back in the list
				last SBJCT;
			}
		}
		push @keep, $query if $unique;
	}
	return \@keep;
}


1;
