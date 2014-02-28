#!/usr/bin/perl

# minimal annotation server

use strict;
use warnings;
use Apache::DBI;
use Bio::DB::GFF;
use CGI qw/header path_info param url request_method/;
use Digest::MD5 'md5_hex';
use Carp;

my $VERSION = 'DAS/1.00';
(my $BASENAME = url(-absolute=>1)) =~ s!http://[^/]+/!!;

use vars qw($DB %ERRCODES %CATEGORIES $HEADER
	    %DSN %TYPE2CATEGORY %TYPEOBJECTS
	    %EXCLUDE
	   );

#         dsn          description              db server    map master
%DSN = (
	'chr22_transcripts'     => [q(EST-predicted transcripts on chr22 from Jean Thierry-Mieg),
				    'dbi:mysql:database=tm_chr22;host=brie3.cshl.org',
				    'http://servlet.sanger.ac.uk:8080/das/ensembl110']
	);
########################################################################################

%ERRCODES = (
	     200 => 'OK',
	     400 => 'Bad command',
	     401 => 'Bad data source',
	     402 => 'Bad command arguments',
	     403 => 'Bad reference object',
	     404 => 'Bad stylesheet',
	     405 => 'Coordinate error',
	     500 => 'Internal server error (oops)',
	     501 => 'Unimplemented feature',
	     );

%CATEGORIES = (
	       component     => [qw(Sequence:Contig Sequence:Link Sequence:Genomic_canonical 
				    static_golden_path_contig:ensembl ensembl_contig:ensembl)],
	       transcription => [qw(Sequence:curated polyA_site stop CpG misc_signal intron exon transcript CDS)],
	       homology      => [qw(similarity)],
	       repeat        => [qw(Alu repeat repeat_region repeat_unit misc_feature)],
	       structural    => [qw(Clone cosmid OLIGO PCR_product structural compression Comment Conflict)],
	       experimental  => [qw(experimental RNAi)],
);

%EXCLUDE = (
	    'static_golden_path_contig:ensembl' => 1,
	    'ensembl_contig:ensembl' => 1,
	    'Sequence:Contig' => 1,
	    );

while (my($c,$v) = each %CATEGORIES) { # invert nicely
  for my $typename (@$v) {
    my $typeobj      = Bio::DB::GFF::Typename->new($typename);
    $TYPE2CATEGORY{$typeobj} = $c;
    $TYPEOBJECTS{$typeobj}   = $typeobj;
  }
}

$HEADER = 0;
my ($junk,$DSN,$OPERATION) = split '/',path_info();

do { error_header('invalid request',400); exit 0 } unless $DSN;
do { list_dsns();   exit 0 } if $DSN eq 'dsn' or $OPERATION eq 'dsn';
do { error_header('invalid data source, use the dsn command to get list',401); exit 0 } unless $DSN{$DSN};

do { error_header('Could not open database',500); exit 0 }
     unless $DB = openDB($DSN);

do { entry_points(); exit 0 } if $OPERATION eq 'entry_points';
do { types();        exit 0 } if $OPERATION eq 'types';
do { features();     exit 0 } if $OPERATION eq 'features';
do { stylesheet();   exit 0 } if $OPERATION eq 'stylesheet';

error_header('invalid request',400);
exit 0;

# -----------------------------------------------------------------
sub openDB {
  my $name = shift;
  my $db = Bio::DB::GFF->new(-adaptor=>'dbi::mysqlopt',-dsn=>$DSN{$name}[1]);
  $db->automerge(0);
  $db->debug(0);
  return $db;
}

# -----------------------------------------------------------------
sub list_dsns {
  my $j = ' 'x3;
  ok_header();
  print qq(<?xml version="1.0" standalone="yes"?>\n<!DOCTYPE DASDSN SYSTEM "http://www.biodas.org/dtd/dasdsn.dtd">\n);
  print "<DASDSN>\n";

  for my $dsn (sort keys %DSN) {
    print "$j<DSN>\n";
    print qq($j$j<SOURCE id="$dsn">$DSN{$dsn}[0]</SOURCE>\n);
    print qq($j$j<MAPMASTER>$DSN{$dsn}[2]/</MAPMASTER>\n);
    print qq($j$j<DESCRIPTION>This is the $DSN{$dsn}[0] database</DESCRIPTION>\n);
    print "$j</DSN>\n";
  }
  print "</DASDSN>\n";
}

# -----------------------------------------------------------------
sub entry_points {
  my $segments = get_segments();

  my @parts;
  if ($segments) {
    @parts = map { get_segment_obj(@$_) } @$segments;
    @parts = map { $_->contained_features(-types=>['Sequence:Link','Sequence:Contig','Sequence:Genomic_canonical'],-merge=>0) } @parts;
  } else {
    @parts = grep {$_->name =~ /^CHR/i} $DB->features(-types=>['Sequence:Link','Sequence:Contig','Sequence:Genomic_canonical'],-merge=>0);
  }

  my $url = get_url();

  ok_header();
  print <<END;
<?xml version="1.0" standalone="no"?>
<!DOCTYPE DASEP SYSTEM "http://www.biodas.org/dtd/dasep.dtd">
<DASEP>
<ENTRY_POINTS href="$url" version="1.0">
END
;

  for my $part (@parts) {
    $part->absolute(1);
    my $name  = $part->name;
    my $st    = $part->start;
    my $en    = $part->stop;
    my $class = $part->class;
    my $length = $part->length;
    my $orientation = $part->strand > 0 ? '+' : '-';
    my $subparts = $part->source =~ /Link|Chromosome|Contig/ ? 'yes' : 'no';
    print qq(<SEGMENT id="$name" size="$length" start="$st" stop="$en" class="$class" orientation="$orientation" subparts="$subparts">$name</SEGMENT>\n);
  }
  print "</ENTRY_POINTS>\n</DASEP>\n";
}

# -----------------------------------------------------------------
# get the features for the segment indicated
sub features {
  my @segments = get_segments() or return;

  my $summary = param('summary');
  my $url = get_url();
  my @filter = param('type');
  my @category = param('category');
  push @filter,make_categories(@category);


  ok_header();
  print <<END
<?xml version="1.0" standalone="yes"?>
<!DOCTYPE DASGFF SYSTEM "http://www.biodas.org/dtd/dasgff.dtd">
<DASGFF>
<GFF version="1.0" href="$url">
END
;

  foreach (@segments) {
    my ($reference,$refclass,$start,$stop) = @$_;
    my $seq = get_segment_obj($reference,$refclass,$start,$stop);
    unless ($seq) {
      print qq(<SEGMENT id="$reference" start="$start" stop="$stop" version="1.0" />\n);
      next;
    }

    if (lc(param('category')) eq 'component') {
      dump_framework($seq);
      next;
    }

    my $r = $seq->refseq;
    my $s = $seq->start;
    my $e = $seq->stop;
    ($s,$e) = ($e,$s) if $s > $e;

    print qq(<SEGMENT id="$r" start="$s" stop="$e" version="1.0">\n);

    my $iterator = $seq->features(-types=>\@filter,-merge=>0,-iterator=>1);

    while (my $f = $iterator->next_seq) {
      my $type        = $f->type;
      next if $EXCLUDE{$type};

      my $flabel      = $f->info || $f->type;
      my $source      = $f->source;
      my $method      = $f->method;
      my $start       = $f->start;
      my $end         = $f->stop;
      my $score       = $f->score;
      my $orientation = $f->strand;
      my $phase       = $f->phase;
      my $group       = $f->group;
      my $id          = $f->id;

      $phase       ||= 0;
      $orientation ||= 0;
      $score       ||= '-';
      $orientation = $orientation >= 0 ? '+' : '-';

      # hack hack hack
      my $category = transmute($type);
      ($start,$end) = ($end,$start) if $start > $end;

      # group stuff
      my $hash       = $group;
#      my @notes      = $f->notes;
      my @notes;
      my $info       = $f->info;
      my $group_info;

      if (ref($info)) {
	my $class = $info->class;
	$id = "$class:$info";
	if ($DSN eq 'elegans') {
	  $group_info = qq(<LINK href="http://www.wormbase.org/db/get?name=$info;class=$class">$info</LINK>);
	}
      } else {
	$hash       = md5_hex($group);
	$group_info = join "\n",map {qq(<NOTE>$_</NOTE>)} @notes;
      }

      my ($target,$target_info);
      if (($target = $f->target) && $target->can('start')) {
	my $start = $target->start;
	my $stop  = $target->stop;
	$target_info = qq(<TARGET id="$target" start="$start" stop="$stop" />);
      }

      if ($category eq 'component') {
	my $strt = 1;
	my $stp  = $stop - $start + 1;
	$target_info = qq(<TARGET id="$id" start="$strt" stop="$stp" />);
      }

    my $map;
    if ($type =~ /Segment|Link|Genomic_canonical|Contig/i) { $map = qq( reference="yes") } else { $map = qq() }
    $map .= qq( subparts="yes") if $type =~ /Segment|Link/i;

    ## Need not print feature for map in annotation services
    ## The next 2 lines are ucommented at Wash U:
    # if (($DSN ne "welegans") && ($c eq "structural")) {
    # } else {

      print <<END;
   <FEATURE id="$id" label="$flabel">
      <TYPE id="$type" category="$category"$map>$type</TYPE>
      <METHOD id="$method">$method</METHOD>
      <START>$start</START>
      <END>$end</END>
      <SCORE>$score</SCORE>
      <ORIENTATION>$orientation</ORIENTATION>
      <PHASE>$phase</PHASE>
END
;
      if ($hash) {
	print qq(      <GROUP id="$hash">\n);
	print qq(        $group_info\n)  if $group_info;
	print qq(        $target_info\n) if $target_info;
	print qq(      </GROUP>\n);
      }
      print <<END;
    </FEATURE>
END
      ;
      # }  # End Wash U if statement
    }

    print qq(</SEGMENT>\n);
  }

print <<END;
</GFF>
</DASGFF>
END
}

sub dump_framework {
  my $seq = shift;
  my $start = $seq->start;
  my $stop  = $seq->stop;

  my @parts = $seq->contained_features(-type=>['Sequence:Link','Sequence:Genomic_canonical','Sequence:Contig'],-merge=>0);

  print qq(<SEGMENT id="$seq" start="$start" stop="$stop" version="1.0">\n);

  for my $part (@parts) {
    my ($st,$en)    = ($part->start,$part->stop);
    my $orientation = $part->strand >= 0 ? '+1' : '-1';
    my $length = $part->length;
    my $type   = $part->type;
    my $method = $type->method;
    my $description = qq(category="component" reference="yes");
    $description .= qq( subparts="yes") unless $part->source eq 'Genomic_canonical';

    print <<END
<FEATURE id="Sequence:$part" label="$part">
   <TYPE id="$type" $description>$part</TYPE>
   <METHOD id="$method">$method</METHOD>
   <START>$st</START>
   <END>$en</END>
   <SCORE>-</SCORE>
   <ORIENTATION>$orientation</ORIENTATION>
   <PHASE>-</PHASE>
   <GROUP id="Sequence:$part">
      <TARGET id="$part" start="1" stop="$length">$part</TARGET>
   </GROUP>
</FEATURE>
END
  ;
  }
  print qq(</SEGMENT>\n);
}

sub types {
  return all_types() unless param('ref') or param('segment');

  my $type    = param('entry_type') || 'Sequence';
  my $summary = param('summary');
  my $url     = get_url();
  my @filter  = param('type');

  my @segments = get_segments() or return;

  ok_header();

  print <<END;
<?xml version="1.0" standalone="yes"?>
<!DOCTYPE DASTYPES SYSTEM "http://www.biodas.org/dtd/dastypes.dtd">
<DASTYPES>
<GFF version="1.2" summary="yes" href="$url">
END
;

  foreach (@segments) {
    my ($reference,$class,$start,$stop) = @$_;
    next unless $reference;
    my $seq = get_segment_obj($reference,$class,$start,$stop) or next;
    unless ($seq) {  #empty section
      print qq(<SEGMENT id="$reference" start="$start" stop="$stop" version="1.0">\n);
      print qq(</SEGMENT>\n);
      next;
    }

    my $s = $seq->start;
    my $e = $seq->stop;

    # use absolute coordinates -- people expect it
    my $name = $seq->refseq;

    print qq(<SEGMENT id="$name" start="$s" stop="$e" version="1.0">\n);

    my @args = (-enumerate=>1);
    push @args,(-types=>\@filter) if @filter;
    my %histogram = $seq->types(@args);
    foreach (keys %histogram) {
      my ($method,$source) = split ':';
      my $count = $histogram{$_};
      my $category  = transmute($_);
      print qq(\t<TYPE id="$_" category="$category" method="$method" source="$source">$count</TYPE>\n)
	unless $EXCLUDE{$_};
    }
    print qq(</SEGMENT>\n);
  }
print <<END;
</GFF>
</DASTYPES>
END
}

# list of all the types
sub all_types {
  my @methods = $DB->types;

  ok_header();
  my $url = get_url();
  print <<END;
<?xml version="1.0" standalone="yes"?>
<!DOCTYPE DASTYPES SYSTEM "http://www.biodas.org/dtd/dastypes.dtd">
<DASTYPES>
<GFF version="1.2" summary="yes" href="$url">
<SEGMENT>
END
    ;

  for my $id (@methods) {
    next if $EXCLUDE{$id};
    my $category = transmute($id);
    my $method   = $id->method;
    my $source   = $id->source;
    print qq(\t<TYPE id="$id" category="$category" method="$method" source="$source" />\n);
  }

  print <<END
</SEGMENT>
</GFF>
</DASTYPES>
END
    ;

}

# Big time kludge -- just outputs the prebuilt stylesheet in this
# directory.  Used primarily for testing.
sub stylesheet {
  ok_header();
  open my $STYLE, '<', "style.xml" or die "Could not read file 'style.xml': $!\n";
  while(<$STYLE>) {
    print $_;
  }
  close $STYLE;
}

# really, really bad shit
# calculate type and category from acedb type and method
sub transmute {
    my $type = shift;

    # look in $TYPE2CATEGORY first to see if we have an exact match
    my $category = $TYPE2CATEGORY{$type};
    return $category if $category;

    # otherwise do a fuzzy match using the values of %TYPEOBJECTS
    for my $typeobj (values %TYPEOBJECTS) {
      warn "comparing $typeobj to $type";

      if ($typeobj->match($type)) {
	$category = $TYPE2CATEGORY{$typeobj};  # fetch category for this object
	$TYPE2CATEGORY{$type} = $category;     # remember this match for later
	return $category;
      }
    }
    return 'miscellaneous';  # no success
}

# -----------------------------------------------------------------
sub get_url {
  my $url = url(-path=>1, -query=>1);
  $url =~ tr/&/\;/;
  return $url;
}

# -----------------------------------------------------------------
sub error_header {
  my ($message,$code) = @_;
  $code ||= 500;
#  $code = "$code $ERRCODES{$code}";
  print header(-type          =>'text/plain',
	       -X_DAS_Version => $VERSION,
	       -X_DAS_Status  => $code,
	      ) unless $HEADER++;
  return if request_method() eq 'HEAD';
  print $message;
}

sub ok_header {
  print header(-type          =>'text/plain',
	       -X_DAS_Version => $VERSION,
	       -X_DAS_Status  => "200 OK", 
	      ) unless $HEADER++;
}

# phony dtd
sub dtd {
    ok_header();
    print <<DTD;
<!-- phony dtd for debugging parsers -->
DTD
}

# -----------------------------------------------------------------
sub get_segments {
  # extended segment argument
  my @segments;
  foreach (param('segment')) {
    my ($ref,$start,$stop) = /^(\S+?)(?::(\d+),(\d+))?$/;
    push @segments,[$ref,$start,$stop];
  }
  push @segments,[scalar param('ref'),scalar param('start'),scalar param('stop')] if param('ref');
  return unless @segments;

  foreach (@segments){
    my ($reference,$start,$stop) = @$_;
    my $class = param('entry_type') || 'Sequence';
    my $name  = $reference;

    if ($reference =~ /^(\w+):(\S+)$/) {
      $class = $1;
      $name  = $2;
    }
    my @values = ($name,$class,$start,$stop);
    $_ = \@values;
  }

  return wantarray ? @segments : \@segments;
}

# -----------------------------------------------------------------
sub get_segment_obj {
  my ($reference,$class,$start,$stop) = @_;
  my @args = (-name=>$reference);
  push @args,(-class=>$class) if defined $class;
  push @args,(-start=>$start) if defined $start;
  push @args,(-stop=>$stop)   if defined $stop;

  my $segment = $DB->segment(@args) or return;
  return $segment;
}


# -----------------------------------------------------------------
sub make_categories {
  my @filter;
  for my $category (@_) {
    my $c = lc $category;
    push @filter,@{$CATEGORIES{$c}} if $CATEGORIES{$c};
    push @filter,$category         unless  $CATEGORIES{$c};
  }
  return @filter;
}
