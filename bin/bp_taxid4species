#!/usr/bin/perl
# Author:  Jason Stajich <jason@bioperl.org>
# Purpose: Retrieve the NCBI Taxa ID for organism(s)

# TODO: add rest of POD
#

use LWP::UserAgent;
use XML::Twig;
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
my $verbose = 0;
my $plain   = 0;
my $help    = 0;
my $USAGE = "taxid4species: [-v] [-p] \"Genus1 species1\" \"Genus2 species2\"";

GetOptions('v|verbose' => \$verbose,
	   'p|plain'   => \$plain,
	   'h|help'    => \$help);
die("$USAGE\n") if $help;

my $ua = new LWP::UserAgent();

my $urlbase = 'https://www.ncbi.nlm.nih.gov/entrez/eutils/';
my $esearch = 'esearch.fcgi?db=taxonomy&usehistory=y&term=';
my $esummary = 'esummary.fcgi?db=taxonomy&query_key=QUERYKEY&WebEnv=WEBENV';

my (@organisms) = @ARGV;
die("must provide valid organism") unless @organisms;
my $organismstr = join(" OR ", @organisms);
$organismstr =~ s/\s/\+/g;

# Esearch
my $response = $ua->get($urlbase . $esearch . $organismstr);
my $t = XML::Twig->new();
print $response->content,"\n"if($verbose);
$t->parse($response->content);
my $root = $t->root;
my $querykey = $root->first_child('QueryKey')->text;
my $webenv = $root->first_child('WebEnv')->text;

# Esummary
$esummary =~ s/QUERYKEY/$querykey/;
$esummary =~ s/WEBENV/$webenv/;
$response = $ua->get($urlbase . $esummary);
$t = XML::Twig->new();
print $response->content,"\n"if($verbose);
$t->parse($response->content);
$root = $t->root;

# Parse XML
my %taxinfo;
foreach my $docsum ($root->children) {
    foreach my $item ($docsum->children('Item')) {
        if ($item->{att}{Name} eq 'ScientificName') {
            my $sciname = $item->text;
            $taxinfo{lc $sciname}{sciname} = $sciname;
            $taxinfo{lc $sciname}{tid} = $docsum->first_child_text('Id');
            last;
        }
    }
}

# Output in same order as given on command line
foreach my $orgn (@organisms) {
    if (exists $taxinfo{lc $orgn}) {
        my $tid = $taxinfo{lc $orgn}{tid};
        
        if ($plain) { print $tid, "\n"; }
        else { print join(", ", "'$orgn'", $tid), "\n"; }
    }
    else { print "'$orgn' not found\n"; }
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

