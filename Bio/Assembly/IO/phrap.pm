# $Id$
#
# BioPerl driver for phrap.out file
#
# Copyright by Robson F. de Souza
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=head1 NAME

Bio::Assembly::IO::phrap - driver to load phrap.out files.

=head1 SYNOPSIS

    # Building an input stream
    use Bio::Assembly::IO;

    # Assembly loading methods
    $io = Bio::Assembly::IO->new(-file=>"SGC0-424.phrap.out",
                                -format=>"phrap");

    $assembly = $io->next_assembly;

=head1 DESCRIPTION

This package was developed to load the phrap.out files from the
(phred/phrap/consed) package by Phill Green. This files contain just
the messages printed to standard out by phrap when building an
assembly.  This output is redirected by phredPhrap perl-script to a
file in the project's directory and hold some bit of information
regarding assembly quality, connections between contigs and clone's
position inside contigs.  It should be noted that such files have no
data about the sequence. neither for contig consensus nor for any
aligned sequence. Anyway, such information may be loaded from Fasta
files in the projects directory and added to the assembly object
later.

Note that, because no sequence is loaded for the contig consensus and
locations for aligned sequences are only given in "ungapped consensus"
coordinates in a phrap.out file, you can't make coordinate changes in
assemblies loaded by pharp.pm, unless you add an aligned
coordinates for each sequence to each contig's features collection
yourself. See L<Bio::Assembly::Contig::Coordinate_Systems> and
L<Bio::Assembly::Contig::Feature_collection>..

This driver also loads singlets into the assembly contigs as
Bio::Assembly::Singlet objects, altough without their sequence strings.
It also adds a feature for the entire sequence, thus storing the singlet
length in its end position, and adds a tag '_nof_trimmed_nonX' to the
feature, which stores the number of non-vector bases in the singlet.

=head2 Implementation

Assemblies are loaded into Bio::Assembly::Scaffold objects composed by
Bio::Assembly::Contig objects. No features are added to Bio::Assembly::Contig
"_aligned_coord:$seqID" feature class, therefore you can't make
coordinate changes in contigs loaded by this module. Contig objects
created by this module will have the following special feature
classes, identified by their primary tags, in their features
collection:

"_main_contig_feature:$ID" : main feature for contig $ID.  This
                              feature is used to store information
                              about the entire consensus
                              sequence. This feature always start at
                              base 1 and its end position is the
                              consensus sequence length. A tag,
                              'trimmed_length' holds the length of the
                              trimmed good quality region inside the
                              consensus sequence.

"_covered_region:$index" : coordinates for valid clones inside the
                              contig. $index is the covered region
                              number, starting at 1 for the covered
                              region closest to the consensus sequence
                              first base.

"_unalign_coord:$seqID" : location of a sequence in "ungapped
                              consensus" coordinates (consensus
                              sequence without gaps).  Primary and
                              secondary scores, indel and
                              substitutions statistics are stored as
                              feature tags.

"_internal_clones:$cloneID" : clones inside contigs $cloneID should be
                              used as the unique id for each
                              clone. These features have six tags:
                              '_1st_name', which is the id of the
                              upstream (5') aligned sequence
                              delimiting the clone; '_1st_strand', the
                              upstream sequence strand in the
                              alignment; '_2nd_name', downstream (3')
                              sequence id; '_2nd_strand', the
                              downstream sequence strand in the
                              alignment; '_length', unaligned clone
                              length; '_rejected', a boolean flag,
                              which is false if the clone is valid and
                              true if it was rejected.

All coordinates for the features above are expressed as "ungapped
consensus" coordinates (See L<Bio::Assembly::Contig::Coordinate_Systems>..

=head2 Feature collection

#

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to the
Bioperl mailing lists  Your participation is much appreciated.

  bioperl-l@bioperl.org                  - General discussion
  http://bioperl.org/wiki/Mailing_lists  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via the
web:

  http://bugzilla.open-bio.org/


=head1 AUTHOR - Robson Francisco de Souza

Email rfsouza@citri.iq.usp.br

head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

package Bio::Assembly::IO::phrap;

use strict;

use Bio::Assembly::Scaffold;
use Bio::Assembly::Singlet;
use Bio::Assembly::Contig;
use Bio::LocatableSeq;
use Bio::Seq;
use Bio::SeqFeature::Generic;

use base qw(Bio::Assembly::IO);

=head2 next_assembly

 Title   : next_assembly
 Usage   : $unigene = $stream->next_assembly()
 Function: returns the next assembly in the stream
 Returns : Bio::Assembly::Scaffold object
 Args    : NONE

=cut

sub next_assembly {
  my $self = shift; # Package reference

  # Resetting assembly data structure
  my $Assembly = Bio::Assembly::Scaffold->new(-source=>'phrap');

  # Looping over all phrap out file lines
  my ($contigOBJ);
  while ($_ = $self->_readline) {
    chomp;

    # Loading exact dupicated reads list
    # /Exact duplicate reads:/ && do {
    #   my @exact_dupl;
    #   while (<FILE>) {
    #     last if (/^\s*$/);
    #     /(\S+)\s+(\S+)/ && do {
    #       push(@exact_dupl,[$1,$2]);
    #     };
    #     $self->{'assembly'}{'exact_dupl_reads'} =
    #       new Data::Table(\@exact_dupl,['included','excluded'],0);
    #   }
    # };

    # Loading singlets reads data
    /^(\d+) isolated singlet/ && do { # should it match 'singlets' and 'singletons'?
      while ($_ = $self->_readline) {
        chomp;
        last if (/^$/);
        if (/^\s+(\S+)\s+(\d+)\s+\((\d+)\)/) {
          my ($singletID, $length, $nof_trimmed_nonX) = ($1, $2, $3);
          # Create singlet object, and add it to scaffold
          my $seq = Bio::LocatableSeq->new(
            -start      => 1,
            -end        => $length,
            -strand     => 1,
            -nowarnonempty => 1,
            -id         => $singletID,
            -primary_id => $singletID,
            -alphabet   => 'dna');
          my $singletOBJ = Bio::Assembly::Singlet->new(-seqref=>$seq);
          my $feat = Bio::SeqFeature::Generic->new(
            -start   => 1,
            -end     => $length,
            -primary => "_main_contig_feature:".$singletOBJ->id(),
            -tag     => { '_nof_trimmed_nonX' => $nof_trimmed_nonX }
          );                         
          $singletOBJ->add_features([ $feat ],1);
          $Assembly->add_singlet($singletOBJ);
        }
      }
    };
  
    # Loading contig information
    /^Contig (\d+)\.\s+(\d+) reads?; (\d+) bp \(untrimmed\), (\d+) \(trimmed\)\./ && do {
      my ($contigID, $nof_reads, $length, $trimmed_length) = ($1, $2, $3, $4);
      $contigOBJ = Bio::Assembly::Contig->new( -id     => $contigID,
                                               -source => 'phrap'   );
      my $feat   = Bio::SeqFeature::Generic->new(
        -start   => 1,
        -end     => $length,
        -primary => "_main_contig_feature:".$contigOBJ->id(),
        -tag     => { '_trimmed_length' => $trimmed_length }
      );
      $contigOBJ->add_features([ $feat ],1);
      $Assembly->add_contig($contigOBJ);
    };
  
    # Loading read information
    /^(C?)\s+(-?\d+)\s+(\d+)\s+(\S+)\s+(\d+)\s+\(\s*(\d+)\)\s+(\d+\.\d*)\s+(\d+\.\d*)\s+(\d+\.\d*)/ && do {
      my ($strand, $start, $end, $readID, $primary_score, $secondary_score,
        $substitutions, $deletions, $insertions) = ($1, $2, $3, $4, $5, $6, $7,
        $8, $9);
      $strand = ($strand eq 'C' ? -1 : 1);
      my $seq = Bio::LocatableSeq->new(
        -start      => $start,
        -end        => $end,
        -nowarnonempty => 1,
        -strand     => $strand,
        -id         => $readID,
        -primary_id => $readID,
        -alphabet   => 'dna');
      my $unalign_coord = Bio::SeqFeature::Generic->new(
        -start   => $start,
        -end     => $end,
        -primary => "_unalign_coord:$readID",
        -tag     => {'_primary_score'=>$primary_score,
                     '_secondary_score'=>$secondary_score,
                     '_substitutions'=>$substitutions,
                     '_insertions'=>,$insertions,
                     '_deletions'=>$deletions }
      );
      $unalign_coord->attach_seq($seq);
      $contigOBJ->add_seq($seq);
      $contigOBJ->add_features([ $unalign_coord ]);
    };
  
    # Loading INTERNAL clones description
    /INTERNAL\s+Contig\s+(\d+)\s+opp\s+sense/ && do {
      my $contigID = $1;
      my $contig = $Assembly->get_contig_by_id($contigID);
      while ($_ = $self->_readline) {
        my (@data,$rejected,$c1_strand,$c2_strand);
  
        (@data = /\s+(\*?)\s+(C?)\s+(\S+)\s+(C?)\s+(\S+)\s+(-?\d+)\s+(-?\d+)\s+(-?\d+)/) && do {
           if ($data[0] eq '*') { $rejected = 1 } else { $rejected = 0 }
           $c1_strand = ($data[1] eq 'C' ? -1 : 1);
           $c2_strand = ($data[3] eq 'C' ? -1 : 1);
           (my $clone_name = $data[2]) =~ s/^(\S+)\.\w.*/$1/;
           my $clone = Bio::SeqFeature::Generic->new(
             -start   => $data[6],
             -end     => $data[7],
             -strand  => 0,
             -primary => "_internal_clone:$clone_name",
             -tag     => {'_1st_strand'=>,$c1_strand,
                          '_2nd_strand'=>,$c2_strand,
                          '_1st_name'=>$data[2],
                          '_2nd_name'=>$data[4],
                          '_length'=>$data[5],
                          '_rejected'=>$rejected}
          );
          $contig->add_features([ $clone ]);
        };
  
        /Covered regions:/ && do {
          my %coord  = /(\d+)/g; my $i = 0;
          foreach my $start (sort { $a <=> $b } keys %coord) {
            my $cov = Bio::SeqFeature::Generic->new(
              -start   => $start,
              -end     => $coord{$start},
              -primary => '_covered_region:'.++$i
            );
            # 1: attach feature to contig consensus, if any
            $contig->add_features([ $cov ],1);
          }
          last; # exit while loop
        }; # /Covered regions:/

      } # while ($_ = $self->_readline)
    }; # /INTERNAL\s+Contig\s+(\d+)\s+opp\s+sense/
  } # while ($_ = $self->_readline)

  return $Assembly;
}

=head2 write_assembly (NOT IMPLEMENTED)

    Title   : write_assembly
    Usage   : $ass_io->write_assembly($assembly)
    Function: Write the assembly object in Phrap compatible ACE format
    Returns : 1 on success, 0 for error
    Args    : A Bio::Assembly::Scaffold object

=cut

sub write_assembly {
    my $self = shift;
    $self->throw_not_implemented();   
}

1;

__END__
