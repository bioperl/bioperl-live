#!/usr/bin/perl

use strict;
use Bio::LiveSeq::IO::BioPerl;
use Bio::LiveSeq::Mutator;
use Bio::LiveSeq::Mutation;
use Bio::Variation::IO;

if ($#ARGV < 1) { # one argument input
    print <<USAGE;

LiveSeq::Mutator example by Joseph Insana

Arguments: filename containing embl entry, gene_name
           It will create that Bio::LiveSeq::Gene and print out some
           basic informations about it.
           It will then issue mutations and print results

Usage:     change_gene.pl filename genename

Example:   change_gene.pl ../../t/data/ar.embl AR
USAGE
exit;
} else {

    my $filename=$ARGV[0];
    my $loader=Bio::LiveSeq::IO::BioPerl->load(-file => "$filename");

    my $gene_name=$ARGV[1];
    my $gene=$loader->gene2liveseq(-gene_name => $gene_name,
				   -getswissprotinfo => 0);

    print STDERR "Gene: ",$gene->name,"\n";
    print STDERR "    Moltype: ", $gene->get_DNA->alphabet,  "\n";
    print STDERR "    Features:\n";
    print STDERR $gene->printfeaturesnum();
    print STDERR "    Gene has boundaries ",$gene->upbound," - ",$gene->downbound,"\n";
    print STDERR "    Gene has maxtranscript with start ",$gene->maxtranscript->start,
          " end ",$gene->maxtranscript->end," strand ",$gene->maxtranscript->strand,"\n";
    print STDERR "    DNA  has boundaries ",$gene->get_DNA->start," - ",$gene->get_DNA->end,"\n";
    print STDERR "\n";

    print STDERR "Now issuing mutations to the gene....\n";

    my $mutation = new Bio::LiveSeq::Mutation (-seq =>'A',
				       -pos => 64
				       );
    my $mutate = Bio::LiveSeq::Mutator->new(-gene => $gene,
					 -numbering => "coding"
					 );
    $mutate->add_Mutation($mutation);
    my $results=$mutate->change_gene();
    print "\n";
    if ($results) {
	my $out = Bio::Variation::IO->new( '-format' => 'flat');
	$out->write($results);
    }
}
