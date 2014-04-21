#
# BioPerl module for Bio::PopGen::Utilities
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Jason Stajich <jason-at-bioperl-dot-org>
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

  my $in = Bio::AlignIO->new(-file   => 't/data/t7.aln',
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
  my $synpop = Bio::PopGen::Utilities->aln_to_population(-site_model => 'cod',
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

=head2 Support 

Please direct usage questions or support issues to the mailing list:

I<bioperl-l@bioperl.org>

rather than to the module maintainer directly. Many experienced and 
reponsive experts will be able look at the problem and quickly 
address it. Please include a thorough description of the problem 
with code and data examples if at all possible.

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution. Bug reports can be submitted via
the web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Jason Stajich

Email jason-at-bioperl-dot-org

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
use constant CodonLen => 3;


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
           'cod' -- codons 

           The option -site_model
                for All sites          : 'all' 
                    Codon sites        : 'cod' or 'codon'

          To see all sites, including those which are fixed in the population
          add -include_monomorphic => 1
          to the arguments
 Returns : 
 Args    : -include_monomorphic => 1   to specify all sites, 
                                       even those which are monomorphic
                                       in the population 
                                  (useful for HKA test mostly) 
                            [default is false]
           -phase          => specify a phase for the data, this is only
                              used if the site_mode is codon
                            [default is 0]
           -site_model     => one-of 'all', 'codon'
                             to specify a site model for the data extraction
                             from the alignment
                            [default is all]
           -alignment      => provide a L<Bio::SimpleAlign> object [required]

=cut

sub aln_to_population{
   my ($self,@args) = @_;
   my ($aln,
       $sitemodel,$phase,
       $includefixed,$checkisa) = $self->_rearrange([qw(ALIGNMENT
				                     SITE_MODEL
					             PHASE
					             INCLUDE_MONOMORPHIC
					             CHECKISA)],
					  @args);

   my %ambig_code = ('?' => ['?','?'],
		     'N' => ['?','?'],
		     '-' => ['?','?'],
                     'G' => ['G','G'],
                     'A' => ['A','A'],
                     'T' => ['T','T'],
                     'C' => ['C','C'],
                     'R' => ['A','G'],
                     'Y' => ['C','T'],
                     'W' => ['T','A'],
                     'M' => ['C','A'],
                     'S' => ['C','G'],
                     'K' => ['G','T']);
					  
   if( ! defined $aln ) { 
       $self->warn("Must provide a valid Bio::SimpleAlign object to run aln_to_population");
       return;
   }

   if( ! $aln->is_flush ) {
       $self->warn("Must provide a Bio::SimpleAlign object with aligned sequences to aln_to_population!");
       return;
   }
   $phase = 0 unless defined $phase;
   if( $phase != 0 && $phase != 1 && $phase != 2 ) { 
       warn("phase must be 0,1, or 2");
       return;
   }
   my $alength = $aln->length;
   my @inds;
   if( ! defined $sitemodel || $sitemodel =~ /all/i ) {
       my $ct = 0;       
       my @seqs;
       for my $seq ( $aln->each_seq ) {
	   push @seqs, $seq->seq;
	   push @inds, Bio::PopGen::Individual->new(-unique_id => $seq->display_id);
       }

       for( my $i = 0; $i < $alength; $i++ ) {
	   my (@genotypes,%set);
	   
           # do we skip indels?
	   # slicing vertically
	   for my $seq ( @seqs ) {
	       my $site = uc(substr($seq,$i,1));
	       push @genotypes, $ambig_code{$site};
	       $set{$site}++;
	   }
	   if( keys %set > 1 || $includefixed ) {
	       my $genoct = scalar @genotypes;
	       for( my $j = 0; $j < $genoct; $j++ ) {
		   $inds[$j]->add_Genotype(Bio::PopGen::Genotype->new
					   (-marker_name  => ($i+1),
					    -individual_id=> $inds[$j]->unique_id,
					    -alleles      => $genotypes[$j]));
	       }
	   }
       }
   } elsif( $sitemodel =~ /cod(on)?/i ) {
       my $ct = 0;
       my @seqs;
       for my $seq ( $aln->each_seq ) {
	   push @seqs, $seq->seq;
	   push @inds, Bio::PopGen::Individual->new(-unique_id => $seq->display_id);
       }
       my $codonct = 0;
       for( my $i = $phase; $i < $alength; $i += CodonLen ) {
	   my (@genotypes,%set,$genoct);
	   
	   for my $seq ( @seqs ) {
	       my @unambig_site;
	       my $site = uc(substr($seq,$i,CodonLen));
	       if( length($site) < CodonLen ) {
		   # at end of alignment and this is not in phase
		   $self->debug("phase was $phase, but got to end of alignment with overhang of $site");
		   next;
	       }
	       # do we check for gaps/indels here?
	       for (my $pos=0; $pos<CodonLen; $pos++)
               {
                    $unambig_site[0] .= $ambig_code{substr($site, $pos, 1)}[0];
                    $unambig_site[1] .= $ambig_code{substr($site, $pos, 1)}[1];
               }
               push @genotypes, [@unambig_site];
               $set{$site}++;
	   }
	   $genoct = scalar @genotypes;
	   
	   # do we include fixed sites? I think we should leave it to the user.
	   if( keys %set > 1 || $includefixed ) {
	       for( my $j = 0; $j < $genoct; $j++ ) {
		   $inds[$j]->add_Genotype(Bio::PopGen::Genotype->new
					   (-marker_name  => ($i/CodonLen),
					    -individual_id=> $inds[$j]->unique_id,
					    -alleles      => $genotypes[$j]));
	       }
	       $codonct++;
	   }
       }
   } else { 
       $self->throw("Can only build sites based on all the data right now!");
   }
   return Bio::PopGen::Population->new(-checkisa => 0,
				       -source => 'alignment',
				       -individuals=> \@inds);
}

1;
