package Bio::DB::SeqFeature::Store::FeatureFileLoader;


=head1 NAME

Bio::DB::SeqFeature::Store::FeatureFileLoader -- feature file loader for Bio::DB::SeqFeature::Store

=head1 SYNOPSIS

  use Bio::DB::SeqFeature::Store;
  use Bio::DB::SeqFeature::Store::FeatureFileLoader;

  # Open the sequence database
  my $db      = Bio::DB::SeqFeature::Store->new( -adaptor => 'DBI::mysql',
                                                 -dsn     => 'dbi:mysql:test',
                                                 -write   => 1 );

  my $loader = 
    Bio::DB::SeqFeature::Store::FeatureFileLoader->new(-store    => $db,
                                                       -verbose  => 1,
					               -fast     => 1);

  $loader->load('./my_genome.fff');


=head1 DESCRIPTION

The Bio::DB::SeqFeature::Store::FeatureFileLoader object parsers
FeatureFile-format sequence annotation files and loads
Bio::DB::SeqFeature::Store databases. For certain combinations of
SeqFeature classes and SeqFeature::Store databases it features a "fast
load" mode which will greatly accelerate the loading of databases by a
factor of 5-10.

FeatureFile Format (.fff) is very simple:

 mRNA B0511.1 Chr1:1..100 Type=UTR;Note="putative primase"
 mRNA B0511.1 Chr1:101..200,300..400,500..800 Type=CDS
 mRNA B0511.1 Chr1:801..1000 Type=UTR

 reference = Chr3
 Cosmid	B0511	516..619
 Cosmid	B0511	3185..3294
 Cosmid	B0511	10946..11208
 Cosmid	B0511	13126..13511
 Cosmid	B0511	11394..11539
 EST	yk260e10.5	15569..15724
 EST	yk672a12.5	537..618,3187..3294
 EST	yk595e6.5	552..618
 EST	yk595e6.5	3187..3294
 EST	yk846e07.3	11015..11208
 EST	yk53c10
 	yk53c10.3	15000..15500,15700..15800
 	yk53c10.5	18892..19154
 EST	yk53c10.5	16032..16105
 SwissProt	PECANEX	13153-13656	Note="Swedish fish"
 FGENESH	"Predicted gene 1"	1-205,518-616,661-735,3187-3365,3436-3846	"Pfam domain"
 # file ends

There are up to four columns of WHITESPACE (not necessarily tab)
delimited text. Embedded whitespace must be escaped using shell
escaping rules (quoting the column or backslashing whitespace).

  Column 1: The feature type. You may use type:subtype as a convention
            for method:source.

  Column 2: The feature name/ID.

  Column 3: The position of this feature in base pair
            coordinates. Ranges can be given as either 
            start-end or start..end. A chromosome position
            can be specified using the format "reference:start..end".
            A discontinuous feature can be specified by giving
            multiple ranges separated by commas. Minus-strand features
            are indicated by specifying a start > end.

  Column 4: Comment/attribute field. A single Note can be given, or
            a series of attribute=value pairs, separated by
            spaces or semicolons, as in "score=23;type=transmembrane"

=head2 Specifying Positions and Ranges

A feature position is specified using a sequence ID (a genbank
accession number, a chromosome name, a contig, or any other meaningful
reference system, followed by a colon and a position range. Ranges are
two integers separated by double dots or the hyphen. Examples:
"Chr1:516..11208", "ctgA:1-5000". Negative coordinates are allowed, as
in "Chr1:-187..1000".

A discontinuous range ("split location") uses commas to separate the
ranges.  For example:

 Gene B0511.1  Chr1:516..619,3185..3294,10946..11208

In the case of a split location, the sequence id only has to appear in
front of the first range.

Alternatively, a split location can be indicated by repeating the
features type and name on multiple adjacent lines:

 Gene	B0511.1	Chr1:516..619
 Gene	B0511.1	Chr1:3185..3294
 Gene	B0511.1	Chr1:10946..11208

If all the locations are on the same reference sequence, you can
specify a default chromosome using a "reference=E<lt>seqidE<gt>":

 reference=Chr1
 Gene	B0511.1	516..619
 Gene	B0511.1	3185..3294
 Gene	B0511.1	10946..11208

The default seqid is in effect until the next "reference" line
appears.

=head2 Feature Tags

Tags can be added to features by adding a fourth column consisting of
"tag=value" pairs:

 Gene  B0511.1  Chr1:516..619,3185..3294 Note="Putative primase"

Tags and their values take any form you want, and multiple tags can be
separated by semicolons. You can also repeat tags multiple times:

 Gene  B0511.1  Chr1:516..619,3185..3294 GO_Term=GO:100;GO_Term=GO:2087

Several tags have special meanings:

 Tag     Meaning
 ---     -------

 Type    The primary tag for a subfeature.
 Score   The score of a feature or subfeature.
 Phase   The phase of a feature or subfeature.
 URL     A URL to link to (via the Bio::Graphics library).
 Note    A note to attach to the feature for display by the Bio::Graphics library.

For example, in the common case of an mRNA, you can use the "Type" tag
to distinguish the parts of the mRNA into UTR and CDS:

 mRNA B0511.1 Chr1:1..100 Type=UTR
 mRNA B0511.1 Chr1:101..200,300..400,500..800 Type=CDS
 mRNA B0511.1 Chr1:801..1000 Type=UTR

The top level feature's primary tag will be "mRNA", and its subparts
will have types UTR and CDS as indicated. Additional tags that are
placed in the first line of the feature will be applied to the top
level. In this example, the note "Putative primase" will be applied to
the mRNA at the top level of the feature:

 mRNA B0511.1 Chr1:1..100 Type=UTR;Note="Putative primase"
 mRNA B0511.1 Chr1:101..200,300..400,500..800 Type=CDS
 mRNA B0511.1 Chr1:801..1000 Type=UTR

=head2 Feature Groups

Features can be grouped so that they are rendered by the "group"
glyph.  To start a group, create a two-column feature entry showing
the group type and a name for the group.  Follow this with a list of
feature entries with a blank type.  For example:

 EST	yk53c10
 	yk53c10.3	15000-15500,15700-15800
 	yk53c10.5	18892-19154

This example is declaring that the ESTs named yk53c10.3 and yk53c10.5
belong to the same group named yk53c10.

=head2 Comments and the #include Directive

Lines that begin with the # sign are treated as comments and
ignored. When a # sign appears within a line, everything to the right
of the symbol is also ignored, unless it looks like an HTML fragment or
an HTML color, e.g.:

 # this is ignored
 [Example]
 glyph   = generic   # this comment is ignored
 bgcolor = #FF0000
 link    = http://www.google.com/search?q=$name#results

Be careful, because the processing of # signs uses a regexp heuristic. To be safe, 
always put a space after the # sign to make sure it is treated as a comment.

The special comment "#include 'filename'" acts like the C preprocessor
directive and will insert the comments of a named file into the
position at which it occurs. Relative paths will be treated relative
to the file in which the #include occurs. Nested #include directives
are allowed:

 #include "/usr/local/share/my_directives.txt"
 #include 'my_directives.txt'
 #include chromosome3_features.gff3

You can enclose the file path in single or double quotes as shown
above. If there are no spaces in the filename the quotes are optional.

Include file processing is not very smart. Avoid creating circular
#include references. You have been warned!

=head2 Caveats

Note that this loader always creates denormalized features such that
subfeatures and their parents are stored as one big database
object. The GFF3 format and its loader is usually preferred for both
space and execution efficiency.

=head1 METHODS

=cut


use strict;
use Carp 'croak';
use File::Spec;
use Text::ParseWords 'shellwords','quotewords';

use base 'Bio::DB::SeqFeature::Store::Loader';

=head2 new

 Title   : new
 Usage   : $loader = Bio::DB::SeqFeature::Store::FeatureFileLoader->new(@options)
 Function: create a new parser
 Returns : a Bio::DB::SeqFeature::Store::FeatureFileLoader parser and loader
 Args    : several - see below
 Status  : public

This method creates a new FeatureFile loader and establishes its connection
with a Bio::DB::SeqFeature::Store database. Arguments are -name=E<gt>$value
pairs as described in this table:

 Name               Value
 ----               -----

 -store             A writable Bio::DB::SeqFeature::Store database handle.

 -seqfeature_class  The name of the type of Bio::SeqFeatureI object to create
                      and store in the database (Bio::DB::SeqFeature by default)

 -sf_class          A shorter alias for -seqfeature_class

 -verbose           Send progress information to standard error.

 -fast              If true, activate fast loading (see below)

 -chunk_size        Set the storage chunk size for nucleotide/protein sequences
                       (default 2000 bytes)

 -tmp               Indicate a temporary directory to use when loading non-normalized
                       features.

When you call new(), a connection to a Bio::DB::SeqFeature::Store
database should already have been established and the database
initialized (if appropriate).

Some combinations of Bio::SeqFeatures and Bio::DB::SeqFeature::Store
databases support a fast loading mode. Currently the only reliable
implementation of fast loading is the combination of DBI::mysql with
Bio::DB::SeqFeature. The other important restriction on fast loading
is the requirement that a feature that contains subfeatures must occur
in the FeatureFile file before any of its subfeatures. Otherwise the
subfeatures that occurred before the parent feature will not be
attached to the parent correctly. This restriction does not apply to
normal (slow) loading.

If you use an unnormalized feature class, such as
Bio::SeqFeature::Generic, then the loader needs to create a temporary
database in which to cache features until all their parts and subparts
have been seen. This temporary databases uses the "bdb" adaptor. The
-tmp option specifies the directory in which that database will be
created. If not present, it defaults to the system default tmp
directory specified by File::Spec-E<gt>tmpdir().

The -chunk_size option allows you to tune the representation of
DNA/Protein sequence in the Store database. By default, sequences are
split into 2000 base/residue chunks and then reassembled as
needed. This avoids the problem of pulling a whole chromosome into
memory in order to fetch a short subsequence from somewhere in the
middle. Depending on your usage patterns, you may wish to tune this
parameter using a chunk size that is larger or smaller than the
default.

=cut

# sub new {} inherited

=head2 load

 Title   : load
 Usage   : $count = $loader->load(@ARGV)
 Function: load the indicated files or filehandles
 Returns : number of feature lines loaded
 Args    : list of files or filehandles
 Status  : public

Once the loader is created, invoke its load() method with a list of
FeatureFile or FASTA file paths or previously-opened filehandles in order to
load them into the database. Compressed files ending with .gz, .Z and
.bz2 are automatically recognized and uncompressed on the fly. Paths
beginning with http: or ftp: are treated as URLs and opened using the
LWP GET program (which must be on your path).

FASTA files are recognized by their initial "E<gt>" character. Do not feed
the loader a file that is neither FeatureFile nor FASTA; I don't know what
will happen, but it will probably not be what you expect.

=cut

# sub load {} inherited

=head2 accessors

The following read-only accessors return values passed or created during new():

 store()          the long-term Bio::DB::SeqFeature::Store object

 tmp_store()      the temporary Bio::DB::SeqFeature::Store object used
                    during loading

 sfclass()        the Bio::SeqFeatureI class

 fast()           whether fast loading is active

 seq_chunk_size() the sequence chunk size

 verbose()        verbose progress messages

=cut

# sub store          {} inherited
# sub tmp_store      {} inherited
# sub sfclass        {} inherited
# sub fast           {} inherited
# sub seq_chunk_size {} inherited
# sub verbose        {} inherited

=head2 default_seqfeature_class

  $class = $loader->default_seqfeature_class

Return the default SeqFeatureI class (Bio::Graphics::Feature).

=cut

sub default_seqfeature_class { #override
  my $self = shift;
  return 'Bio::Graphics::Feature';
}


=head2 load_fh

  $count = $loader->load_fh($filehandle)

Load the FeatureFile data at the other end of the filehandle and return true
if successful. Internally, load_fh() invokes:

  start_load();
  do_load($filehandle);
  finish_load();

=cut

# sub load_fh { } inherited

=head2 start_load, finish_load

These methods are called at the start and end of a filehandle load.

=cut

sub create_load_data {
    my $self = shift;
    $self->SUPER::create_load_data();
    $self->{load_data}{mode}          = 'fff';
    $self->{load_data}{CurrentGroup}  = undef;
}

sub finish_load {
    my $self = shift;
    $self->_store_group;
    $self->SUPER::finish_load;
}

=head2 load_line

    $loader->load_line($data);

Load a line of a FeatureFile file. You must bracket this with calls to
start_load() and finish_load()!

    $loader->start_load();
    $loader->load_line($_) while <FH>;
    $loader->finish_load();

=cut

sub load_line {
    my $self = shift;
    my $line = shift;

    chomp($line);
    return unless $line =~ /\S/;     # blank line
    my $load_data = $self->{load_data};

    $load_data->{mode} = 'fff' if /\s/;  # if it has any whitespace in
                                         # it, then back to fff mode

    if ($line =~ /^\#\s?\#\s*([\#]+)/) {    ## meta instruction
      $load_data->{mode} = 'fff';
      $self->handle_meta($1);

    } elsif ($line =~ /^\#/) {
      $load_data->{mode} = 'fff';  # just to be safe
      return; # comment
    }

    elsif ($line =~ /^>\s*(\S+)/) { # FASTA lines are coming
      $load_data->{mode} = 'fasta';
      $self->start_or_finish_sequence($1);
    }

    elsif ($load_data->{mode} eq 'fasta') {
      $self->load_sequence($line);
    }

    elsif ($load_data->{mode} eq 'fff') {
      $self->handle_feature($line);
      if (++$load_data->{count} % 1000 == 0) {
	my $now = $self->time();
	my $nl = -t STDOUT && !$ENV{EMACS} ? "\r" : "\n";
	$self->msg(sprintf("%d features loaded in %5.2fs...$nl",
			   $load_data->{count},$now - $load_data->{start_time}));
	$load_data->{start_time} = $now;
      }
    }

    else {
      $self->throw("I don't know what to do with this line:\n$line");
    }
}


=head2 handle_meta

  $loader->handle_meta($meta_directive)

This method is called to handle meta-directives such as
##sequence-region. The method will receive the directive with the
initial ## stripped off.

=cut

# sub handle_meta { } inherited

=head2 handle_feature

  $loader->handle_feature($gff3_line)

This method is called to process a single FeatureFile line. It manipulates
information stored a data structure called $self-E<gt>{load_data}.

=cut

sub handle_feature {
  my $self     = shift;
  local $_     = shift;

  my $ld       = $self->{load_data};

  # handle reference line
  if (/^reference\s*=\s*(.+)/) {
      $ld->{reference} = $1;
      return;
  }

  # parse data lines
  my @tokens = quotewords('\s+',1,$_);
  for (0..2) { # remove quotes from everything but last column
      next unless defined $tokens[$_];
      $tokens[$_] =~ s/^"//;
      $tokens[$_] =~ s/"$//;
  }

  if (@tokens < 3) {      # short line; assume a group identifier
      $self->store_current_feature();
      my $type               = shift @tokens;
      my $name               = shift @tokens;
      $ld->{CurrentGroup}    = $self->_make_indexed_feature($name,$type,'',{_ff_group=>1});
      $self->_indexit($name => 1);
      return;
  }

  my($type,$name,$strand,$bounds,$attributes);

  if ($tokens[2] =~ /^([+-.]|[+-]?[01])$/) { # old version
      ($type,$name,$strand,$bounds,$attributes) = @tokens;
  } else {                                   # new version
      ($type,$name,$bounds,$attributes) = @tokens;
  }

  # handle case of there only being one value in the last column,
  # in which case we treat it the same as Note="value"
  my $attr = $self->parse_attributes($attributes);

  # @parts is an array of ([ref,start,end],[ref,start,end],...)
  my @parts = map { [/(?:(\w+):)?(-?\d+)(?:-|\.\.)(-?\d+)/]} split /(?:,| )\s*/,$bounds;

  # deal with groups -- a group is ending if $type is defined
  # and CurrentGroup is set
  if ($type && $ld->{CurrentGroup}) {
      $self->_store_group();
  }

  $type   = '' unless defined $type;
  $name   = '' unless defined $name;
  $type ||= $ld->{CurrentGroup}->primary_tag if $ld->{CurrentGroup};

  my $reference = $ld->{reference} || 'ChrUN';
  foreach (@parts) {
      if (defined $_ && ref($_) eq 'ARRAY' 
	  && defined $_->[1] 
	  && defined $_->[2]) 
      {
	  $strand     ||= $_->[1] <= $_->[2] ? '+' : '-';
	  ($_->[1],$_->[2])   = ($_->[2],$_->[1]) if $_->[1] > $_->[2];
      }
      $reference = $_->[0] if defined $_->[0];
      $_ = [@{$_}[1,2]]; # strip off the reference.
  }

  # now @parts is an array of [start,end] and $reference contains the seqid

  # apply coordinate mapper
  if ($self->{coordinate_mapper} && $reference) {
      my @remapped = $self->{coordinate_mapper}->($reference,@parts);
      ($reference,@parts) = @remapped if @remapped;
  }

  # either create a new feature or add a segment to it
  my $feature = $ld->{CurrentFeature};
  
  $ld->{OldPartType} = $ld->{PartType};
  if (exists $attr->{Type} || exists $attr->{type})  {
      $ld->{PartType}   = $attr->{Type}[0] || $attr->{type}[0];
  } else {
      $ld->{PartType}   = $type;
  }

  if ($feature) {
      local $^W = 0;  # avoid uninit warning when display_name() is called

      # if this is a different feature from what we have now, then we
      # store the current one, and create a new one
      if ($feature->display_name ne $name ||
	  $feature->method       ne $type) {
	  $self->store_current_feature;  # new feature, store old one
	  undef $feature;
      } else { # create a new multipart feature
	  $self->_multilevel_feature($feature,$ld->{OldPartType})
	      unless $feature->get_SeqFeatures;
	  my $part = $self->_make_feature($name,
					  $ld->{PartType},
					  $strand,
					  $attr,
					  $reference,
					  @{$parts[0]});
	  $feature->add_SeqFeature($part);
      }
  }

  $feature ||= $self->_make_indexed_feature($name,
					    $type,   # side effect is to set CurrentFeature
					    $strand,
					    $attr,
					    $reference,
					    @{$parts[0]});

  # add more segments to the current feature
  if (@parts > 1) {
      for my $part (@parts) {
	  $type ||= $feature->primary_tag;
	  my $sp  = $self->_make_feature($name,
					 $ld->{PartType},
					 $strand,
					 $attr,
					 $reference,
					 @{$part});
      $feature->add_SeqFeature($sp);
      }
  }
}

sub _multilevel_feature { # turn a single-level feature into a multilevel one
    my $self = shift;
    my $f    = shift;
    my $type = shift;
    my %attributes     = $f->attributes;
    $attributes{Score} = [$f->score] if defined $f->score;
    $attributes{Phase} = [$f->phase] if defined $f->phase;
    my @args = ($f->display_name,
		$type||$f->type,
		$f->strand,
		\%attributes,
		$f->seq_id,
		$f->start,
		$f->end);
    my $subpart = $self->_make_feature(@args);
    $f->add_SeqFeature($subpart);
}

sub _make_indexed_feature {
    my $self = shift;
    my $f    = $self->_make_feature(@_);
    my $name = $f->display_name;
    $self->{load_data}{CurrentFeature} = $f;
    $self->{load_data}{CurrentID}      = $name;
    $self->_indexit($name => 1);
    return $f;
}

sub _make_feature {
    my $self = shift;
    my ($name,$type,$strand,$attributes,$ref,$start,$end) = @_;

    # some basic error checking
     $self->throw("syntax error at line $.: '$_'")
 	if ($ref   && !defined $start)
 	or ($ref   && !defined $end)
 	or ($start && $start   !~  /^[-\d]+$/)
 	or ($end   && $end     !~  /^[-\d]+$/)
 	or !defined $type
 	or !defined $name;

    $strand ||= '';

    my @args = (-name        => $name,
		-strand      => $strand eq '+' ?  1 
                               :$strand eq '-' ? -1
		               :$strand eq ''  ?  0
                               :$strand eq '.' ?  0
		               :$strand == 1   ?  1
		               :$strand == -1  ? -1
                               :0,
		-attributes  => $attributes,
	);

    if (my ($method,$source) = $type =~ /(\S+):(\S+)/) {
	push @args,(-primary_tag => $method,
		    -source      => $source);
    } else {
	push @args,(-primary_tag => $type);
    }

    push @args,(-seq_id       => $ref)   if defined $ref;
    push @args,(-start        => $start) if defined $start;
    push @args,(-end          => $end)   if defined $end;

    # pull out special attributes
    if (my $score = $attributes->{Score} || $attributes->{score}) {
	push @args,(-score => $score->[0]);
	delete $attributes->{$_} foreach qw(Score score);
    }

    if (my $note  = $attributes->{Note}  || $attributes->{note}) {
	push @args,(-desc => join '; ',@$note);
	delete $attributes->{$_} foreach qw(Note note);
    }

    if (my $url = $attributes->{url} || $attributes->{Url}) {
	push @args,(-url => $url->[0]);
	delete $attributes->{$_} foreach qw (Url url);
    }

    if (my $phase = $attributes->{phase} || $attributes->{Phase}) {
	push @args,(-phase => $phase->[0]);
	delete $attributes->{$_} foreach qw (Phase phase);
    }

    $self->_indexit($name=>1)
	if $self->index_subfeatures && $name;

    return $self->sfclass->new(@args);
}

=head2 store_current_feature

  $loader->store_current_feature()

This method is called to store the currently active feature in the
database. It uses a data structure stored in $self-E<gt>{load_data}.

=cut

sub store_current_feature {  # overridden
    my $self = shift;

    # handle open groups
    # if there is an open group, then we simply add the current
    # feature to the group.
    my $ld = $self->{load_data};
    if ($ld->{CurrentGroup} && $ld->{CurrentFeature}) {
	$ld->{CurrentGroup}->add_SeqFeature($ld->{CurrentFeature})
	    unless $ld->{CurrentGroup} eq $ld->{CurrentFeature};   # paranoia - shouldn't happen
	return;
    }
    else {
	$self->SUPER::store_current_feature();
    }
}

sub _store_group {
    my $self  = shift;
    my $ld    = $self->{load_data};
    my $group = $ld->{CurrentGroup} or return;
    # if there is an unattached feature, then add it
    $self->store_current_feature() if $ld->{CurrentFeature};
    $ld->{CurrentFeature} = $group;
    $ld->{CurrentID}      = $group->display_name;
    $self->_indexit($ld->{CurrentID} => 1);
    undef $ld->{CurrentGroup};
    $self->store_current_feature();
}

=head2 build_object_tree

 $loader->build_object_tree()

This method gathers together features and subfeatures and builds the
graph that connects them.

=cut

###
# put objects together
#
sub build_object_tree {
    croak "We shouldn't be building an object tree in the FeatureFileLoader";
}

=head2 build_object_tree_in_tables

 $loader->build_object_tree_in_tables()

This method gathers together features and subfeatures and builds the
graph that connects them, assuming that parent/child relationships
will be stored in a database table.

=cut

sub build_object_tree_in_tables {
    croak "We shouldn't be building an object tree in the FeatureFileLoader";
}

=head2 build_object_tree_in_features

 $loader->build_object_tree_in_features()

This method gathers together features and subfeatures and builds the
graph that connects them, assuming that parent/child relationships are
stored in the seqfeature objects themselves.

=cut

sub build_object_tree_in_features {
  croak "We shouldn't be building an object tree in the FeatureFileLoader";
}

=head2 attach_children

 $loader->attach_children($store,$load_data,$load_id,$feature)

This recursively adds children to features and their subfeatures. It
is called when subfeatures are directly contained within other
features, rather than stored in a relational table.

=cut

sub attach_children {
    croak "We shouldn't be attaching children in the
    FeatureFileLoader!";
}

=head2 parse_attributes

 @attributes = $loader->parse_attributes($attribute_line)

This method parses the information contained in the $attribute_line
into a flattened hash (array). It may return one element, in which case it is
an implicit

=cut

sub parse_attributes {
  my $self  = shift;
  my $att   = shift;

  $att     ||= ''; # to prevent uninit variable warnings from quotewords()

  my @pairs =  quotewords('[;\s]',1,$att);
  my %attributes;
  for my $pair (@pairs) {
      unless ($pair =~ /=/) {
	  push @{$attributes{Note}},(quotewords('',0,$pair))[0] || $pair;
      } else {
	  my ($tag,$value) = quotewords('\s*=\s*',0,$pair);
	  $tag = 'Note' if $tag eq 'description';
	  push @{$attributes{$tag}},$value;
      }
  }
  return \%attributes;
}

=head2 start_or_finish_sequence

  $loader->start_or_finish_sequence('Chr9')

This method is called at the beginning and end of a fasta section.

=cut


1;

__END__


=head1 BUGS

This is an early version, so there are certainly some bugs. Please
use the BioPerl bug tracking system to report bugs.

=head1 SEE ALSO

L<Bio::DB::SeqFeature::Store>,
L<Bio::DB::SeqFeature::Segment>,
L<Bio::DB::SeqFeature::NormalizedFeature>,
L<Bio::DB::SeqFeature::GFF3Loader>,
L<Bio::DB::SeqFeature::Store::DBI::mysql>,
L<Bio::DB::SeqFeature::Store::bdb>

=head1 AUTHOR

Lincoln Stein E<lt>lstein@cshl.orgE<gt>.

Copyright (c) 2006 Cold Spring Harbor Laboratory.

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.

=cut


