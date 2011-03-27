#!perl -w
#
#   A client showing how to use Bio::Biblio module, a module for
#   accessing and querying a bibliographic repository.
#   It also shows how to use modules Bio::Biblio::IO::medlinexml
#   Bio::Biblio::IO::medline2ref which converts XML MEDLINE
#   citations into a simple hash table and into full Perl objects.
#
#   It has many options in order to cover as many methods as
#   possible.  Because of that, it can be also used as a fully
#   functional command-line client for querying repository and
#   retrieving citations from it.
#
#   Usage: biblio.PLS -h
#
#   martin.senger@gmail.com
#   February 2002
#
#-----------------------------------------------------------------------------

use strict;

sub get_usage {
    exec('perldoc',$0);
}

BEGIN {
    # add path to the directory with this script
    my $mylib;
    ($mylib = $0) =~ s|/[^/]+$||;
    unshift @INC, $mylib;

    # be prepare for command-line options/arguments
    use Getopt::Std;

    # general options
    use vars qw/ $opt_h $opt_v $opt_q /;
    # specialized options
    use vars qw/ $opt_a $opt_b $opt_c $opt_d $opt_D $opt_e $opt_k $opt_n $opt_p $opt_r $opt_s /;
    # options with a value
    use vars qw/ $opt_f $opt_F $opt_g $opt_i $opt_l $opt_m $opt_O $opt_V /;
    my $switches = 'fFgilmOV';   # these are switches taking an argument (a value)
    getopt ($switches);

    # help wanted?
    if ($opt_h) {
	print get_usage;
	exit 0;
    }
}

use Bio::Biblio;       # to read data via SOAP
use Bio::Biblio::IO;   # to convert resulting XML to Biblio objects
use Data::Dumper;      # to print resulting data in a raw form

# --- print version and exit
if ($opt_v) {
    print "$Bio::Biblio::VERSION\n";
    print "$Bio::Biblio::Revision\n";
    exit 0;
}

# --- deal with a local file
&convert_and_print ($opt_f) if $opt_f;

# --- create a Biblio object;
#     the new() method understands the following parameters
#     (but none of them is mandatory - unless the default service location
#      is not where you want to go today):
#
#      -location           (taken from '-l' option if given)
#      -collection_id      (taken from '-i' option if given)
#      -destroy_on_exit    (set to false if '-k' or '-p' or '-i' are given)
#
#      And just for information (these can be used from your script but
#      they are not set-able by this script):
#
#      -access => 'soap'   (not set-able here, a default value will be used)
#      -namespace => '...' (not set-able here, a default value will be used)
#      -soap               (not set-able here)
#
#      Additionally, it uses env. variable HTTPPROXY to create parameter
#      '-httpproxy'.
#
my @location   = ('-location', $opt_l) if defined $opt_l;
my @collection = ('-collection_id', $opt_i) if defined $opt_i;
my @destroy    = ('-destroy_on_exit', 0) if $opt_k or $opt_p or $opt_i;
my @httpproxy  = ('-httpproxy', $ENV{'HTTPPROXY'}) if defined $ENV{'HTTPPROXY'};
my $biblio = new Bio::Biblio (@location, @collection, @destroy, @httpproxy);

die "Stopped. No success in accessing the bibliographic repository.\n" unless $biblio;

#
# all remaining command-line arguments (if any remains after getopts) are:
#     -find <keywords> [-attrs <attributes>]
# and these (up-to-)pairs can be repeated...
#
# ...and it creates a query collection (perhaps more than one) and
# assigns it (or the last one in case on 'chained' finds) to $biblio
#
my ($keywords, $attrs, $next);
while ($next = shift) {
    if ($next eq '-find') {
	$biblio = &_find ($biblio, $keywords, $attrs) if $keywords;
	$keywords = shift;
	undef $attrs;
    } elsif ($next eq '-attrs') {
	$attrs = shift;
    }
}
$biblio = &_find ($biblio, $keywords, $attrs) if $keywords;

#
# now we have either the top-level collection (if there were no -find
# arguments), or a resulting collection from the -find queries above
# ...let's do with it what was asked by options
#

# ...print the number of citations
print $biblio->get_count . "\n" if $opt_c;

# ...get one particular citation (this method does not use any -finds above)
&convert_and_print ($biblio->get_by_id ($opt_g)) if $opt_g;

# ...print all citation IDs
print join ("\n", @{ $biblio->get_all_ids }) . "\n" if $opt_d;

# ...print all citations - returned as one big string from the server
&convert_and_print ($biblio->get_all) if $opt_s;

# ... reset iteration in the collection again to the first citation
if ($opt_r) {
    $biblio->reset_retrieval;
    print "Reset OK.\n" unless $opt_q;
}

# ...print more citations (perhaps all) - returned as an array of citations
$opt_m = 100000000 if $opt_a;
if (defined $opt_m) {
    foreach my $cit (@{ $biblio->get_more ($opt_m) }) {
	&convert_and_print ($cit);
    }
}

# ...print next citation from the current collection
&convert_and_print ($biblio->get_next) if $opt_n;

# ...check existence of a collection and completeness of its iterator
if ($opt_e) {
    my $exists = $biblio->exists;
    my $has_next = $biblio->has_next if $exists;
    $exists = '0' unless $exists;
    $has_next = '0' unless $has_next;

    if ($opt_q) {
	print "$exists\n$has_next\n";
    } else {
	print "Exists: $exists\tHas next: $has_next\n";
    }
}

# ...destroy collection
if ($opt_D) {
    $biblio->destroy;
    print "Destroyed OK.\n" unless $opt_q;
}

# ...print the collection ID
if ($opt_p) {
    my $id = $biblio->get_collection_id;
    print "$id\n" if $id;
}

# ...controlled vocabularies
if ($opt_V) {

    # ...print all vocabulary names (-Vn)
    if ($opt_V =~ /^n/) {
	print join ("\n", @{ $biblio->get_vocabulary_names }) . "\n";

    } else {
	my ($arg, $name, $value) = split (/\:\:/, $opt_V, 3);

	# ...print all values from a given vocabulary (-Vv::<name>)
	if ($opt_V =~ /^v/) {
	    print join ("\n", @{ $biblio->get_all_values ($name) }) . "\n";

	# ...print all entries from a given vocabulary (-Va::<name>)
	} elsif ($opt_V =~ /^a/) {
	    print Data::Dumper->Dump ( [$biblio->get_all_entries ($name)], ['All entries']);

	# ...print description of a given vocabulary entry (-Vd::<name>::<value>)
	} elsif ($opt_V =~ /^d/) {
	    print $biblio->get_entry_description ($name, $value) . "\n";

	# ...check existence of a vocabulary value (-Ve::<name>::<value>)
	} elsif ($opt_V =~ /^e/) {
            my $contains = $biblio->contains ($name, $value);
	    $contains = '0' unless $contains;
	    print "Value '$value' in vocabulary '$name': $contains\n" unless $opt_q;
	    print "$contains\n" if $opt_q;
	}
    }
}

sub _find {
    my ($biblio, $keywords, $attrs) = @_;
    $| = 1;
    print "Looking for '$keywords'" . ($attrs ? " in attributes '$attrs'..." : "...")
	unless $opt_q;
    my ($new_biblio) = $biblio->find ($keywords, $attrs);
    print "\tFound " . $new_biblio->get_count . "\n"
	unless $opt_q;
    print "\tReturned collection is '" . $new_biblio->get_collection_id . "'.\n"
	if $opt_k and not $opt_q;
    return $new_biblio;
}

sub convert_and_print {
    my ($citation) = @_;

    # if no -O option given or if it is -Ox we are happy returning XML string
    unless (defined $opt_O and $opt_O !~ /^x/) {
	return if $opt_f;   # we do not do a simple file reading
	&print_one ($citation);
	return;
    }

    my @args;

    # -Or means to return a raw hash, everything else means to return
    # Biblio objects - but there may be more types of them depending
    # also on -F (which format the citation is in)
    if ($opt_O =~ /^r/) {
	push (@args, ('-result' => 'raw'));
    } elsif ($opt_F and $opt_F =~ /^p/) {
	push (@args, ('-result' => 'pubmed2ref'));
    }
    # default: -result => 'medline2ref'

    # an argument to specify that we want parse XML (which we always want
    # but there can be various XML formats)
    if ($opt_F and $opt_F =~ /^p/) {
	push (@args, ('-format' => 'pubmedxml'));
    } else {
	push (@args, ('-format' => 'medlinexml'));
    }

    # where to take the citation from
    if ($opt_f) {
	push (@args, ('-file' => $citation));
    } else {
	push (@args, ('-data' => $citation));
    }

    # make an instance of a converter
    my $io = new Bio::Biblio::IO (@args);

    # and finally make the conversion
    while (my $bibref = $io->next_bibref) {
	&print_one ($bibref);
    }
#    return $io->next_bibref;
}

sub print_one {
    my ($citation) = @_;
    return unless defined $citation;
    if (ref (\$citation) eq 'SCALAR') {
	print $citation;
    } elsif (ref ($citation) =~ /^HASH|ARRAY|SCALAR$/o) {
	print Data::Dumper->Dump ( [$citation], ['Citation']);
    } else {
	print $citation->print_me;
    }
}


__END__


=head1 NAME

biblio.PLS - bioperl client for accessing and querying a bibliographic
repository.

=head1 SYNOPSIS

Usage:

  biblio.PLS [vh]
  biblio.PLS [bcFgOpq]         [-l <URL>]
  biblio.PLS [abcdDeFknmOpqrs] [-l <URL>] -i <collection-ID>
  biblio.PLS [abcdDeFknmOpqrs] [-l <URL>] - -find <keywords> \
         [-attrs <attrs>]...
  biblio.PLS [Vq]              [-l <URL>]
  biblio.PLS [FOq]             [-f <filename>]

=head1 DESCRIPTION

A client showing how to use Bio::Biblio module, a module for accessing
and querying a bibliographic repository.  It also shows how to use
modules Bio::Biblio::IO::medlinexml Bio::Biblio::IO::medline2ref which
converts XML MEDLINE citations into a simple hash table and into full
Perl objects.

It has many options in order to cover as many methods as possible.
Because of that, it can be also used as a fully functional
command-line client for querying repository and retrieving citations
from it.

=head1 OPTIONS

=head2 What service to contact:

  -l <URL> ... a location where a Bibliographic Query service is
               provided as a WebService
               (default: http://www.ebi.ac.uk/openbqs/services/MedlineSRS)

=head2 What query collection to use:

Some options do not need to specify a collection, some do.

  -i <collection_id>  ... the collection ID can be obtained in a
                          previous invocation by specifying argument
                          '-p' (print ID)
  -find <keywords> [-attrs <attrs>]
                      ... create a collection from citations
                          containing given keywords - either in all
                          default attributes, or only in the given
                          attributes;

                          it is possible to repeat it, for example:
                             -find brazma -attrs authors -find -study
                          (the repetitions refine previous results)
                          both <keywords> and <attrs> may be
                          comma-delimited multi-values;
                          note that '-find' must be separated from
                          the rest of options by '-';

                          note that this script is a bit stupid
                          regarding quoted keywords, or keywords
                          containing commans... TBD better

  what XML format is used for citations:
  -Fm     ... MEDLINE (default)
  -Fp     ... PubMed

=head2 What to do (with the query collection):

  -g <id>    ... get citation <id>
  -c         ... get count (a number of citations)
  -p         ... print collection ID (which may be used in the next
                 invocation as an '-i' argument); it implies also '-k'
  -b         ... print citations in a non-XML format (TBD)

Other options can be used only on a sub-collection - which can be
obtained directly by specifying '-i' argument, or indirectly by
specifying one or more queries by '-find' arguments:

  -d         ... get all citation IDs
  -n         ... get next citation
  -m [<how_many>] ... get 'how_many' more
  -r         ... reset iteration to the first citation in the collection
                 (now you can use '-n' or '-m' again)
  -a         ... get all citations - as an array
  -s         ... as '-a' but get it as one string
  -e         ... check if given collection exists and has more citations
  -k         ... keep resulting collection persistent (makes sense only
                 when collection IDs are being printed otherwise you
                 would not know how to contact the persistent collection
                 next time)
  -D         ... destroy given collection (makes sense together with '-i')

Options specifying output format of the results:

  -Ox        ... output in XML format (default)
  -Oo        ... output as Biblio objects
  -Or        ... output as a raw hashtable

  The options above can be used also for converting an XML MEDLINE
  local file without using any SOAP connection at all;

  -f <filename> ... an XML file to be read and converted

Options dealing with controlled vocabularies:

  -Vn                  ... get all vocabulary names
  -Vv::<name>          ... get all values from vocabulary <name>
  -Va::<name>          ... get everything from vocabulary <name>
  -Vd::<name>::<value> ... get description of <value>
                           from vocabulary <name>
  -Ve::<name>::<value> ... return 1 if <value> exists
                           in vocabulary <name>

And the remaining options:

  -h  ... get help
  -v  ... get version
  -q  ... be quiet (less verbose)

=head1 EXAMPLES

  biblio.PLS - -find Java -attrs abstract -find perl

Several separate invocations sharing the same query collection:

  biblio.PLS -p -q - -find Brazma,Robinson > b.tmp
  biblio.PLS -i `cat b.tmp` -d
  MEDLINE2005/10693778
  MEDLINE2005/10977099
  MEDLINE2005/11726920
  MEDLINE2005/12225585
  MEDLINE2005/12227734
  biblio.PLS -i `cat b.tmp` -g 10693778
  <MedlineCitation Status="Completed">
   ...
  </MedlineCitation>

  biblio.PLS -i `cat b.tmp` -e
  Exists: 1       Has next: 1

  biblio.PLS -i `cat b.tmp` -D
  Destroyed OK.

  biblio.PLS -i `cat b.tmp` -e
  Exists: 0       Has next: 0

Access to controlled vocabularies:

  biblio.PLS -Vn
  MEDLINE2005/JournalArticle/properties
  MEDLINENEW/resource_types
  MEDLINE2005/resource_types
  MEDLINE2005/Person/properties
  MEDLINE2005/*/publication_type
  MEDLINENEW/JournalArticle/properties
  repository_subsets
  MEDLINE2005/*/citation_subset

  biblio.PLS -Vv::MEDLINE2005/JournalArticle/properties
  AllText
  ID
  PMID
  ISSN
  ...

Converting local XML MEDLINE file:

  biblio.PLS -g 10693778 > a_file.xml
  biblio.PLS -f a_file.xml -Oo   ... to Perl objects
  biblio.PLS -f a_file.xml -Or   ... as a raw hash

=head1 ENVIRONMENT VARIABLES

  HTTPPROXY = <HTTP proxy server>

Use this if you use this script on a machine which needs to access
remote HTTP targets via a proxy server.  For example:

  export HTTPPROXY=http://128.243.220.41:3128
  biblio.PLS -c

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

  bioperl-l@bioperl.org                  - General discussion
  http://bioperl.org/wiki/Mailing_lists  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug trackingi system to help us keep track
of the bugs and their resolution. Bug reports can be submitted via the
web:

  https://redmine.open-bio.org/projects/bioperl/

=head1 AUTHOR - Martin Senger

Email martin.senger@gmail.com

=head1 HISTORY

Written February 2002
Updated July 2005

=cut
