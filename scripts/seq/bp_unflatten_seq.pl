#!perl
use strict;
use warnings;
# Author Chris Mungall <cjm-at-bioperl.org>

=head1 NAME

bp_unflatten_seq - unflatten a genbank or genbank-style feature file into
a nested SeqFeature hierarchy

=head1 SYNOPSIS

  bp_unflatten_seq.PLS -e 3 -gff ~/cvs/bioperl-live/t/data/AE003644_Adh-genomic.gb

  bp_unflatten_seq.PLS --detail ~/cvs/bioperl-live/t/data/AE003644_Adh-genomic.gb

  bp_unflatten_seq.PLS -i foo.embl --from embl --to chadoxml -o out.chado.xml

  bp_unflatten_seq.PLS --notypemap --detail --to asciitree -ethresh 2 AE003644_Adh-genomic.gb

=head1 DESCRIPTION 

This script will B<unflatten> a genbank or genbank-style file of
SeqFeatures into a nested hierarchy.

See L<Bio::SeqFeature::Tools::Unflattener>

In a GenBank/EMBL representation, features are 'flat' - for example,
there is no link between an mRNA and a CDS, other than implicit links
(eg via tags or via splice site coordinates) which may be hard to code
for.

This is most easily illustrated with the default output format,
B<asciitree>

An unflattened genbank feature set may look like this (AB077698)

  Seq: AB077698
    databank_entry                                   1..2701[+]
    gene                                             
      mRNA                                           
        CDS hCHCR-G                                  80..1144[+]
        exon                                         80..1144[+]
      five_prime_UTR                                 1..79[+]
      located_sequence_feature                       137..196[+]
      located_sequence_feature                       239..292[+]
      located_sequence_feature                       617..676[+]
      located_sequence_feature                       725..778[+]
      three_prime_UTR                                1145..2659[+]
      polyA_site                                     1606..1606[+]
      polyA_site                                     2660..2660[+]

Or like this (portion of AE003734)


  gene                                             
    mRNA CG3320-RA                                 
      CDS CG3320-PA                                53126..54971[-]
      exon                                         52204..53323[-]
      exon                                         53404..53631[-]
      exon                                         53688..53735[-]
      exon                                         53798..53918[-]
      exon                                         54949..55287[-]
    mRNA CG3320-RB                                 
      CDS CG3320-PB                                53383..54971[-]
      exon                                         52204..53631[-]
      exon                                         53688..53735[-]
      exon                                         53798..53918[-]
      exon                                         54949..55287[-]

The unflattening will also 'normalize' the containment hierarchy (in
the sense of standardising it - e.g. making sure there is always a
transcript record, even if genbank just specifies CDS and gene)

By default, the GenBank types will be mapped to SO types

See L<Bio::SeqFeature::Tools::TypeMapper>

=head1 COMMAND LINE ARGUMENTS

=over

=item -i|input FILE

input file (can also be specified as last argument)

=item -from FORMAT

input format (defaults to genbank)

probably doesnt make so much sense to use this for non-flat formats;
ie other than embl/genbank

=item -to FORMAT

output format (defaults to asciitree)

should really be a format that is nested SeqFeature aware; I think
this is only asciitree, chadoxml and gff3

=item -gff

with export to GFF3 format (pre-3 GFFs make no sense with unflattened
sequences, as they have no set way of representing feature graphs)

=item -o|output FILE

outfile defaults to STDOUT

=item -detail

show extra detail on features (asciitree mode only)

=item -e|ethresh INT

sets the error threshold on unflattening

by default this script will throw a wobbly if it encounters weird
stuff in the genbank file - raise the error threshold to signal these
to be ignored (and reported on STDERR)

=item -nomagic

suppress use_magic in unflattening (see
L<Bio::SeqFeature::Tools::Unflattener>

=item -notypemap

suppress type mapping (see
L<Bio::SeqFeature::Tools::TypeMapper>


=back

=head1 TODO

L<Bio::SeqFeature::Tools::Unflattener> allows fine-grained control
over the unflattening process - need to add more options to allow this
control at the command line


=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

  bioperl-l@bioperl.org                  - General discussion
  http://bioperl.org/wiki/Mailing_lists  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution. Bug reports can be submitted via
email or the web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR

 Chris Mungall E<lt>cjm-at-bioperl.orgE<gt>

=cut

use Bio::SeqIO;
use Bio::SeqFeature::Tools::Unflattener;
use Bio::SeqFeature::Tools::TypeMapper;
use Bio::SeqFeature::Tools::IDHandler;
use Bio::Tools::GFF;

use Getopt::Long;

my ($input,$from,$to,$output,$verbosity,$ethresh,$nomagic,$group_tag,$detail,
    $notypemap);
$from = 'genbank';
$to = 'asciitree';
$ethresh = 3;
my $gff;
my @remove_types = ();

GetOptions(
	   'i|input:s' => \$input,
	   'from:s'  => \$from,
	   'to:s'  => \$to,
	   'o|output:s'=> \$output,
	   "verbosity|v=s"=>\$verbosity,
	   "ethresh|e=s"=>\$ethresh,
	   "remove_type=s@"=>\@remove_types,
	   "nomagic"=>\$nomagic,
	   "notypemap"=>\$notypemap,
	   "group_tag"=>\$group_tag,
	   "detail"=>\$detail,
           "gff"=>\$gff,
	   "h|help"=>sub {
	       system("perldoc $0");
	       exit 0;
	   },
	  );
                       
                       
if ($to =~ /^gff/i) {
    $gff = 1;
}

$input = $input || shift if @ARGV;

my $in = new Bio::SeqIO(-file => $input,
			-format => $from);
my $out;
my @out_opt = $output ? (-file => ">$output") : ();
unless ($gff) {
    $out = new Bio::SeqIO(-format=>$to, @out_opt);
    $out->show_detail($detail) if $out->can("show_detail") && $detail;
}

my $unflattener = Bio::SeqFeature::Tools::Unflattener->new;
$unflattener->verbose($verbosity);
$unflattener->error_threshold($ethresh);
my $tm = Bio::SeqFeature::Tools::TypeMapper->new;
my $idhandler = Bio::SeqFeature::Tools::IDHandler->new;

while( my $seq = $in->next_seq ) {    
    $unflattener->remove_types(-seq=>$seq,
                               -types=>\@remove_types)
      if @remove_types;

    $unflattener->unflatten_seq(-seq=>$seq,
				-use_magic=>!$nomagic,
				-group_tag=>$group_tag,
			       );
    $unflattener->report_problems(\*STDERR);
    $tm->map_types_to_SO(-seq=>$seq) unless $notypemap;

    my @seq_args = ($seq);
    if ($to eq 'chadoxml') {
	@seq_args = (-seq=>$seq, -nounflatten=>1)
    }
    if ($gff) {
        my $gffio = Bio::Tools::GFF->new(@out_opt, -noparse=>1, -gff_version => 3);
        $idhandler->set_ParentIDs_from_hierarchy($seq);
        foreach my $feature ($seq->get_all_SeqFeatures) {
            $gffio->write_feature($feature);
        }
        $gffio->close();
    }
    else {
        $out->write_seq(@seq_args);
    }

}

__END__
