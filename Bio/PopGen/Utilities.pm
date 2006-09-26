# $Id$
#
# BioPerl module for Bio::PopGen::Utilities
#
# Cared for by Jason Stajich <jason-at-open-bio-dot-org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::PopGen::Utilities - Utilities for working with PopGen data and objects

=head1 SYNOPSIS

  use Bio::PopGen::Utilities;
  use Bio::AlignIO;

  my $in = new Bio::AlignIO(-file   => 't/data/t7.aln',
                            -format => 'clustalw');
  my $aln = $in->next_aln;
  # get a population, each sequence is an individual and 
  # for the default case, every site which is not monomorphic
  # is a 'marker'.  Each individual will have a 'genotype' for the
  # site which will be the specific base in the alignment at that
  # site
  my $pop = Bio::PopGen::Utilities->aln_to_population(-alignment => $aln);

  # get the synonymous sites from the alignemt only as the 'genotypes'
  # for the population
  my $synpop = Bio::PopGen::Utilities->aln_to_population(-site_model => 'syn',
                                                         -alignment  => $aln);


=head1 DESCRIPTION

This object provides some convience function to turn sequence
alignments into usable objects for the Population genetics modules
(Bio::PopGen).

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
the web:

  http://bugzilla.open-bio.org/

=head1 AUTHOR - Jason Stajich

Email jason-at-open-bio-dot-org

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::PopGen::Utilities;
use strict;

use Bio::Align::DNAStatistics;
use Bio::PopGen::Population;
use Bio::PopGen::Individual;

use base qw(Bio::Root::Root);


=head2 aln_to_population

 Title   : aln_to_population
 Usage   : my $pop = Bio::PopGen::Utilities->aln_to_population($aln);
 Function: Turn and alignment into a set of L<Bio::PopGen::Individual>
           objects grouped in a L<Bio::PopGen::Population> object

           Sites are treated as 'Markers' in the Bioperl PopGen object
           model in the sense that a site is a unique location for which
           an individual will have a genotype (a set of alleles). 
           In this implementation we are assuming that each individual 
           has a single entry in the alignment file.

           Specify a site model as one of those listed
           'all' -- every base in the alignment is considered a site
           'syn' -- Synonomous sites. Those where a seen substition do 
                    not change the amino acid [Assumes this is only 
                    coding sequence and the frame starts with first base 
                    in the alignment]
           'non' -- Non-Synonomous sites.  Those where a substitution changes
                    the encoded amino acid.

           The option -site_model
                for Non-synonymous: 'non' or 'non-synonomous' or 'NS' or 'Ka'
		    Synonymous	  : 'synonomous' or 'syn' or 'S' or 'Ks'
                    All           : 'all' 
          To see all sites, including those which are fixed in the population
          add -include_monomorphic => 1
          to the arguments
 Returns : 
 Args    : -include_monomorphic => 1   to specify all sites, 
                                       even those which are monomorphic
                                       in the population 
                                  (useful for HKA test mostly) 
                            [default is false]
           -site_model     => one-of 'all', 'syn', or 'non' 
                             to specify a site model you want to see data
                             for
                            [default is all]
           -alignment      => provide a L<Bio::SimpleAlign> object [required]

=cut

sub aln_to_population{
   my ($self,@args) = @_;
   my ($aln,
       $sitemodel,
       $includefixed) = $self->_rearrange([qw(ALIGNMENT
					      SITE_MODEL
					      INCLUDE_MONOMORPHIC)],
					  @args);
   if( ! defined $aln ) { 
       $self->warn("Must provide a valid Bio::SimpleAlign object to run aln_to_population");
       return;
   }
   if( ! $aln->is_flush ) {
       $self->warn("Must provide a Bio::SimpleAlign object with aligned sequences to aln_to_population!");
       return;
   }

   my $population = Bio::PopGen::Population->new(-source => 'alignment');
   my @seqs = map { $_->seq() } $aln->each_seq;

   if( ! defined $sitemodel ||
       $sitemodel =~ /all/i ) {
       my $ct = 0;
       my @inds;
       my @seqs;
       for my $seq ( $aln->each_seq ) {
	   my $ind = Bio::PopGen::Individual->new(-unique_id => $seq->display_id);
	   push @seqs, $seq->seq;
	   push @inds, $ind;
       }
       for( my $i = 0; $i < $aln->length; $i++ ) {
	   my $nm = "Site-$i";
	   my (@genotypes,%set);
	   # do we skip indels?
	   for my $seq ( @seqs ) {
	       my $site = substr($seq,$i,1);
	       $set{$site}++;
	       push @genotypes, $site;
	   }
	   if( keys %set > 1 || $includefixed ) {
	       for( my $i = 0; $i < scalar @genotypes; $i++ ) {
		   $inds[$i]->add_Genotype(Bio::PopGen::Genotype->new
					   (-marker_name  => $nm,
					    -individual_id=> $inds[$i]->unique_id,
					    -alleles      => [$genotypes[$i]]));
	       }
	   }
       }
       for my $ind ( @inds ) { 
	   $population->add_Individual($ind);
       }
   } else { 
       $self->throw("Can only build sites based on all the data right now!");
       my ($sitecount,@sites) = ($aln->length);
       my @sitecat;
       # ToDo: categorize site a syn, non-syn, monomorphic
       #      4-fold degenerate?
       my (@codons,@codons_v, $codon_ct);
       
       for( my $i = 0; $i < $sitecount; $i++ ) {
	   if( $i && $i % 3 == 0 ) {
	       # A A T  T T G  T C G
	       # A A A  T A G  T A G
	       # A A T  T A G  T C T	       
	       for my $cod ( @{$codons[$codon_ct]} ) {
		   $codons_v[$codon_ct]->{$cod}++;
	       }
	       $codon_ct++;	       
	   }
	   my $seqct = 0;
	   foreach my $seq ( @seqs ) {   
	       my $char = substr($seq,$i,1);
	       $sites[$i]->{'alleles'}->{$char}++;
	       $sites[$i]->{'seq'}->[$seqct] = $char;
	       $codons[$codon_ct]->[$seqct] .= $char;	       
	       $seqct++;
	   }
       }

       # at the end @sites will be full, each entry is a column and it
       # will have a hashref with 2 values, 'alleles' which will have
       # a frequency for each base as an allele. 
       # 'seq' will have the
       # participating residue for each sequence


       my ($i,$seqctr) = (0,0); 
       for my $site ( @sites ) { 
	   my %alleles = %{$site->{'alleles'}};
	   my %codons = $codons_v[$i % 3]->[$seqctr];
	   $i++;
	   $seqctr++;
       }
   }
   return $population;
}



1;
