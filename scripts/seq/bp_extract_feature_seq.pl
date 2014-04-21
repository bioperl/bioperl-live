#!perl
use strict;
use warnings;
use Bio::SeqIO;
use Getopt::Long;
# Author Jason Stajich <jason@bioperl.org>

=head1 NAME

bp_extract_feature_seq - extract the corresponding sequence for a specified feature type

=head1 SYNOPSIS

bp_extract_feature_seq [--format FORMAT] [--feature CDS] [--output FILE] [--input] FILE

=head1 DESCRIPTION

This script will extract the sequence for all the features you specify.

=head1 OPTIONS

=over

=item B<-i>, B<--input>

Specifies the sequence file to be read.

=item B<--format>

Format of the file specifed by B<--input>. If not given, it will try to guess the
correct format from the file extension.

=item B<--feature>

Feature to be extracted. By default, it extracts the CDS feature.

=item B<-o>, B<--output>

File where the extracted features will be saved. If not specified, STDOUT is used.

=back

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list. Your participation is much appreciated.

  L<bioperl-l@bioperl.org>                  - General discussion
  L<http://bioperl.org/wiki/Mailing_lists>  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution. Bug reports can be submitted via
email or the web:

  L<https://github.com/bioperl/bioperl-live/issues>

=head1 AUTHOR

 Jason Stajich <jason-at-bioperl-dot-org>

=cut

my ($input,$format,$featuretype,$output);
$featuretype ='CDS';
GetOptions(
           'i|input:s' => \$input,
           'format:s'  => \$format,
           'feature:s' => \$featuretype,
           'o|output:s'=> \$output);

$input || shift if @ARGV;

my $in = new Bio::SeqIO(-file => $input,
                        -format => $format);
my $out;
if ($output ) {
    $out = new Bio::SeqIO(-file => ">$output", -format => 'fasta');
} else { 
    $out = new Bio::SeqIO(-format => 'fasta'); # use STDOUT for output
}

my $count = 1;
while( my $seq = $in->next_seq ) {
    foreach my $f ( grep { $_->primary_tag =~ /$featuretype/i }
                    $seq->get_SeqFeatures ) {
        my $s = $f->spliced_seq;
        if( $featuretype =~ /gene|CDS/ ) {
            $s->display_id($f->has_tag('gene') ? join(',',sort $f->each_tag_value('gene')) :
                           $f->has_tag('label') ? join(',',$f->each_tag_value('label')): 
                           $s->display_id);
        } else {
            $s->display_id(sprintf("%s_%s_%d",
                                   $seq->display_id,
                                   $f->primary_tag,
                                   $count++));
        }
        $out->write_seq($s);
    }
}
