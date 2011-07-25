#!perl
#
# OBDA Registry compliant sequence retrieval script
#
# Copyright Heikki Lehvaslaiho <heikki-at-bioperl-dot-org>
# You may distribute this program under the same terms as perl itself


use Bio::DB::Registry;
use Bio::SeqIO;
use Getopt::Long;
use strict;
use warnings;

my ($help, $format, $namespace, $dbname) = ('', 'embl', 'acc', 'embl');
GetOptions ("help" => \$help, "format=s" => \$format,
            "namespace=s" => \$namespace, "dbname=s" => \$dbname );
if ($help || !@ARGV) {
  system("perldoc $0");
  exit 0;
}

my $get_function = 'get_Seq_by_'. $namespace;

my $registry = new Bio::DB::Registry();
while (my $id = shift) {
    my $db = $registry->get_database($dbname);
    my $seq = $db->$get_function($id);

    if ($seq) {
        my $out = new Bio::SeqIO('-format' => $format);
        $out->write_seq($seq);
    } else {
        print STDERR "Could not find sequence with identifier [$id]\n";
    }
}

=head1 NAME

bp_biogetseq - sequence retrieval using OBDA registry

=head1 DESCRIPTION

This script retrives sequences from the source defined by users
registry setup.  The current alternatives are from a local indexed
file, sql database or over the web.

=head1 USAGE

  Usage: bp_biogetseq --dbname embl --format embl
                   --namespace acc id [ ids... ]
         * dbname defaults to embl
         * format defaults to embl
         * namespace defaults to 'acc' ['id', 'acc', 'version']
         * unnamed arguments are ids in the given namespace

=cut

