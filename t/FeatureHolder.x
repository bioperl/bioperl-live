# -*-Perl-*- mode (to keep my emacs happy)
# $Id: FeatureHolder.x,v 1.1 2004-03-06 02:02:28 cjm Exp $

use strict;
use vars qw($DEBUG $TESTCOUNT);
BEGIN {     
    eval { require Test; };
    if( $@ ) {
	use lib 't';
    }
    use Test;
    $TESTCOUNT = 6;
    plan tests => $TESTCOUNT;
}

use Bio::Seq;
use Bio::SeqIO;
use Bio::Root::IO;
use Bio::SeqFeature::Tools::Unflattener;
use Bio::SeqFeature::Tools::TypeMapper;
use Bio::Tools::GFF;

ok(1);

my $verbosity = -1;   # Set to -1 for release version, so warnings aren't printed

my ($seq, @sfs);
my $unflattener = Bio::SeqFeature::Tools::Unflattener->new;
my $tm = Bio::SeqFeature::Tools::TypeMapper->new;


if (1) {
    my @path = ("t","data","AE003644_Adh-genomic.gb");
    # allow cmd line override
    if (@ARGV) {
	@path = (shift @ARGV);
    }
    $seq = getseq(@path);
    
    ok ($seq->accession_number, 'AE003644');
    my @topsfs = $seq->get_SeqFeatures;
    
    # UNFLATTEN
    $unflattener->verbose(1);
    $unflattener->unflatten_seq(-seq=>$seq,
                                -use_magic=>1);
    $tm->map_types_to_SO(-seq=>$seq);
    $path[-1] .= ".gff3";
    my $o =
      Bio::Root::IO->catfile(
                             @path
                            );
    my $gffio = Bio::Tools::GFF->new(-file=>">$o" , -noparse=>1, -gff_version => 3);
    $seq->set_ParentIDs_from_hierarchy();
    foreach my $feature ($seq->get_all_SeqFeatures) {
        $gffio->write_feature($feature);
    }
    $gffio->close();

    $gffio = Bio::Tools::GFF->new(-file=>"$o", -gff_version => 3);
    my $seq = Bio::Seq->new;
    my $feature;
    # loop over the input stream
    while($feature = $gffio->next_feature()) {
        $seq->add_SeqFeature($feature);
    }
    $gffio->close();
    $seq->create_hierarchy_from_ParentIDs;

    $o =~ s/\.gff3$/chado\-xml/;
    my $outio = new Bio::SeqIO(-format=>'chadoxml', -file=>">$o");
    $outio->write_seq($seq);

    # no way to check chado output for now

    $o =~ s/\.chado\-xml$/chaos\-xml/;
    $outio = new Bio::SeqIO(-format=>'chaosxml', -file=>">$o");
    $outio->write_seq($seq);

    # no way to check chado output for now

}


sub getseq {
    my @path = @_;
    my $seqio =
      Bio::SeqIO->new('-file'=> Bio::Root::IO->catfile(
                                                       @path
                                                      ), 
                      '-format' => 'GenBank');
    $seqio->verbose($verbosity);

    my $seq = $seqio->next_seq();
    return $seq;
}
