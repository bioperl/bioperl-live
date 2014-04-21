#!perl
# Author:  Jason Stajich <jason@bioperl.org>
# Purpose: Retrieve the NCBI Taxa ID for organism(s)

# TODO: add rest of POD
#

use LWP::UserAgent;
use XML::Twig;
use strict;
use warnings;
use Getopt::Long;
my $verbose = 0;
my $plain   = 0;
my $help    = 0;
my $USAGE = "taxid4species: [-v] [-p] \"Genus1 species1\" \"Genus2 species2\"";

GetOptions('v|verbose' => \$verbose,
	   'p|plain'   => \$plain,
	   'h|help'    => \$help);
die("$USAGE\n") if $help;

my $ua = new LWP::UserAgent();

my $urlbase = 'http://www.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=taxonomy&term=';

my (@organisms) = @ARGV;
die("must provide valid organism") unless @organisms;
my $organismstr = join(" OR ", @organisms);
$organismstr =~ s/\s/\+/g;

my $response = $ua->get($urlbase.$organismstr);
my $t = XML::Twig->new();
print $response->content,"\n"if($verbose);
$t->parse($response->content);
my $root = $t->root;
my $list = $root->first_child('IdList');
my @data;
foreach my $child ($list->children('Id') ) {
    push @data, $child->text;
    if( $plain ) { print $child->text, "\n" }
}
unless( $plain  ) {
    $list = $root->first_child('TranslationStack');
    foreach my $set ($list->children('TermSet') ) {
	foreach my $term ( $set->children('Term') ) {
	    print "\"",$term->text(), "\", ", shift @data, "\n";
	}
    }
}

=head1 NAME

bp_taxid4species - simple script which returns the NCBI Taxonomic id for a requested species

=head1 SYNOPSIS

bp_taxid4species [-v] [-p] [-h] "Genus1 species1" "Genus2 species2"

Options:
  -v   verbose
  -p   plain
  -h   help

=head1 DESCRIPTION

This simple script shows how to get the taxa id from NCBI Entrez and
will return a list of taxa ids for requested organisms.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

  bioperl-l@bioperl.org                  - General discussion
  http://bioperl.org/wiki/Mailing_lists  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution. Bug reports can be submitted via the
web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR

 Jason Stajich jason-at-bioperl-dot-org

=cut

