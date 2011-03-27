#
# BioPerl module for Bio::Align::Utilities
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Jason Stajich <jason-at-bioperl.org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Align::Utilities - A collection of utilities regarding converting
and manipulating alignment objects

=head1 SYNOPSIS

  use Bio::Align::Utilities qw(:all);
  # %dnaseqs is a hash of CDS sequences (spliced)


  # Even if the protein alignments are local make sure the start/end
  # stored in the LocatableSeq objects are to the full length protein.
  # The CoDing Sequence that is passed in should still be the full 
  # length CDS as the nt alignment will be generated.
  #
  my $dna_aln = &aa_to_dna_aln($aa_aln,\%dnaseqs);


  # generate bootstraps
  my $replicates = &bootstrap_replicates($aln,$count);


=head1 DESCRIPTION

This module contains utility methods for manipulating sequence
alignments ( L<Bio::Align::AlignI>) objects.

The B<aa_to_dna_aln> utility is essentially the same as the B<mrtrans>
program by Bill Pearson available at
ftp://ftp.virginia.edu/pub/fasta/other/mrtrans.shar.  Of course this
is a pure-perl implementation, but just to mention that if anything
seems odd you can check the alignments generated against Bill's
program.

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
of the bugs and their resolution. Bug reports can be submitted via the
web:

  https://redmine.open-bio.org/projects/bioperl/

=head1 AUTHOR - Jason Stajich

Email jason-at-bioperl-dot-org

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

#' keep my emacs happy
# Let the code begin...


package Bio::Align::Utilities;
use vars qw(@EXPORT @EXPORT_OK $GAP $CODONGAP %EXPORT_TAGS);
use strict;
use Carp;
use Bio::Root::Version;
require Exporter;

use base qw(Exporter);

@EXPORT = qw();
@EXPORT_OK = qw(aa_to_dna_aln bootstrap_replicates cat);
%EXPORT_TAGS = (all =>[@EXPORT, @EXPORT_OK]);
BEGIN {
    use constant CODONSIZE => 3;
    $GAP = '-';
    $CODONGAP = $GAP x CODONSIZE;
}

=head2 aa_to_dna_aln

 Title   : aa_to_dna_aln
 Usage   : my $dnaaln = aa_to_dna_aln($aa_aln, \%seqs);
 Function: Will convert an AA alignment to DNA space given the 
           corresponding DNA sequences.  Note that this method expects 
           the DNA sequences to be in frame +1 (GFF frame 0) as it will
           start to project into coordinates starting at the first base of 
           the DNA sequence, if this alignment represents a different 
           frame for the cDNA you will need to edit the DNA sequences
           to remove the 1st or 2nd bases (and revcom if things should be).
 Returns : Bio::Align::AlignI object 
 Args    : 2 arguments, the alignment and a hashref.
           Alignment is a Bio::Align::AlignI of amino acid sequences. 
           The hash reference should have keys which are 
           the display_ids for the aa 
           sequences in the alignment and the values are a 
           Bio::PrimarySeqI object for the corresponding 
           spliced cDNA sequence. 

See also: L<Bio::Align::AlignI>, L<Bio::SimpleAlign>, L<Bio::PrimarySeq>

=cut

sub aa_to_dna_aln {
    my ($aln,$dnaseqs) = @_;
    unless( defined $aln && 
	    ref($aln) &&
	    $aln->isa('Bio::Align::AlignI') ) { 
	croak('Must provide a valid Bio::Align::AlignI object as the first argument to aa_to_dna_aln, see the documentation for proper usage and the method signature');
    }
    my $alnlen = $aln->length;
    my $dnaalign = Bio::SimpleAlign->new();
    $aln->map_chars('\.',$GAP);

    foreach my $seq ( $aln->each_seq ) {    
	my $aa_seqstr = $seq->seq();
	my $id = $seq->display_id;
	my $dnaseq = $dnaseqs->{$id} || $aln->throw("cannot find ".
						     $seq->display_id);
	my $start_offset = ($seq->start - 1) * CODONSIZE;

	$dnaseq = $dnaseq->seq();
	my $dnalen = $dnaseqs->{$id}->length;
	my $nt_seqstr;
	my $j = 0;
	for( my $i = 0; $i < $alnlen; $i++ ) {
	    my $char = substr($aa_seqstr,$i + $start_offset,1);	    
	    if ( $char eq $GAP || $j >= $dnalen )  { 
		$nt_seqstr .= $CODONGAP;
	    } else {
		$nt_seqstr .= substr($dnaseq,$j,CODONSIZE);
		$j += CODONSIZE;
	    }
	}
	$nt_seqstr .= $GAP x (($alnlen * 3) - length($nt_seqstr));

	my $newdna = Bio::LocatableSeq->new(-display_id  => $id,
					   -alphabet    => 'dna',
					   -start       => $start_offset+1,
					   -end         => ($seq->end * 
							    CODONSIZE),
					   -strand      => 1,
					   -seq         => $nt_seqstr);    
	$dnaalign->add_seq($newdna);
    }
    return $dnaalign;
}

=head2 bootstrap_replicates

 Title   : bootstrap_replicates
 Usage   : my $alns = &bootstrap_replicates($aln,100);
 Function: Generate a pseudo-replicate of the data by randomly
           sampling, with replacement, the columns from an alignment for
           the non-parametric bootstrap.
 Returns : Arrayref of L<Bio::SimpleAlign> objects
 Args    : L<Bio::SimpleAlign> object
           Number of replicates to generate

=cut

sub bootstrap_replicates {
   my ($aln,$count) = @_;
   $count ||= 1;
   my $alen = $aln->length;
   my (@seqs,@nm);
   $aln->set_displayname_flat(1);
   for my $s ( $aln->each_seq ) {
       push @seqs, $s->seq();
       push @nm, $s->id;
   }
   my (@alns,$i);
   while( $count-- > 0 ) {
       my @newseqs;
       for($i =0; $i < $alen; $i++ ) {
	   my $index = int(rand($alen));
	   my $c = 0;
	   for ( @seqs ) {
	       $newseqs[$c++] .= substr($_,$index,1);
	   }
       }
       my $newaln = Bio::SimpleAlign->new();
       my $i = 0;
       for my $s ( @newseqs ) {
       (my $tmp = $s) =~ s{[$Bio::LocatableSeq::GAP_SYMBOLS]+}{}g;
	   $newaln->add_seq( Bio::LocatableSeq->new
			     (-start         => 1,
			      -end           => length($tmp),
			      -display_id    => $nm[$i++],
			      -seq           => $s));
       }
       push @alns, $newaln;
   }
   return \@alns;
}

=head2 cat

 Title     : cat
 Usage     : $aln123 = cat($aln1, $aln2, $aln3)
 Function  : Concatenates alignment objects. Sequences are identified by id.
             An error will be thrown if the sequence ids are not unique in the
             first alignment. If any ids are not present or not unique in any
             of the additional alignments then those sequences are omitted from
             the concatenated alignment, and a warning is issued. An error will
             be thrown if any of the alignments are not flush, since
             concatenating such alignments is unlikely to make biological
             sense.
 Returns   : A new Bio::SimpleAlign object
 Args      : A list of Bio::SimpleAlign objects

=cut

sub cat {
    my ($self, @aln) = @_;
    $self->throw("cat method called with no arguments") unless $self;
    for ($self,@aln) {
	$self->throw($_->id. " not a Bio::Align::AlignI object") unless $_->isa('Bio::Align::AlignI');
	$self->throw($_->id. " is not flush") unless $_->is_flush;
    }
    my $aln = $self->new;
    $aln->id($self->id);
    $aln->annotation($self->annotation);
    my %unique;
  SEQ: foreach my $seq ( $self->each_seq() ) {
        throw("ID: ", $seq->id, " is not unique in initial alignment.") if exists $unique{$seq->id};
        $unique{$seq->id}=1;

        # Can be Bio::LocatableSeq, Bio::Seq::Meta or Bio::Seq::Meta::Array
 	my $new_seq = $seq->new(-id=> $seq->id,
                                -strand  => $seq->strand,
                                -verbose => $self->verbose);
	$new_seq->seq($seq->seq);
	$new_seq->start($seq->start);
	$new_seq->end($seq->end);
        if ($new_seq->isa('Bio::Seq::MetaI')) {
            for my $meta_name ($seq->meta_names) {
                $new_seq->named_submeta($meta_name, $new_seq->start, $new_seq->end, $seq->named_meta($meta_name));
            }
        }
 	for my $cat_aln (@aln) {
            my @cat_seq=$cat_aln->each_seq_with_id($seq->id);
 	    if (@cat_seq==0) {
                $self->warn($seq->id. " not found in alignment ". $cat_aln->id. ", skipping this sequence.");
                next SEQ;
            }
 	    if (@cat_seq>1) {
                $self->warn($seq->id. " found multiple times in alignment ".  $cat_aln->id. ", skipping this sequence.");
                next SEQ;
            }
            my $cat_seq=$cat_seq[0];
            my $old_end=$new_seq->end;
            $new_seq->seq($new_seq->seq.$cat_seq->seq);
            
            # Not sure if this is a sensible way to deal with end coordinates
            $new_seq->end($new_seq->end+$cat_seq->end+1-$cat_seq->start);

            if ($cat_seq->isa('Bio::Seq::Meta::Array')) {
                unless ($new_seq->isa('Bio::Seq::Meta::Array')) {
                    my $meta_seq=Bio::Seq::Meta::Array->new;
                    $meta_seq->seq($new_seq->seq);
                    $meta_seq->start($new_seq->start);
                    $meta_seq->end($new_seq->end);
                    if ($new_seq->isa('Bio::Seq::Meta')) {
                        for my $meta_name ($new_seq->meta_names) {
                            $meta_seq->named_submeta($meta_name,
                                                     $new_seq->start,
                                                     $old_end,
                                                     [split(//, $new_seq->named_meta($meta_name))]
                                                    );
                        }
                    }
                    $new_seq=$meta_seq;
                }
                for my $meta_name ($cat_seq->meta_names) {
                    $new_seq->named_submeta($meta_name,
                                            $old_end+1,
                                            $new_seq->end,
                                            $cat_seq->named_meta($meta_name)
                                           );
                }
            } elsif ($cat_seq->isa('Bio::Seq::Meta')) {
                if ($new_seq->isa('Bio::Seq::Meta::Array')) {
                    for my $meta_name ($cat_seq->meta_names) {
                        $new_seq->named_submeta($meta_name,
                                                $old_end+1,
                                                $new_seq->end,
                                                [split(//,$cat_seq->named_meta($meta_name))]
                                               );
                    }
                } else {
                    unless ($new_seq->isa('Bio::Seq::Meta')) {
                        my $meta_seq=Bio::Seq::Meta::Array->new;
                        $meta_seq->seq($new_seq->seq);
                        $meta_seq->start($new_seq->start);
                        $meta_seq->end($new_seq->end);
                        $new_seq=$meta_seq;
                    }
                    for my $meta_name ($cat_seq->meta_names) {
                        $new_seq->named_submeta($meta_name,
                                                $old_end+1,
                                                $new_seq->end,
                                                $cat_seq->named_meta($meta_name)
                                               );
                    }
                }
            }
        }
        $aln->add_seq($new_seq);
    }
    my $cons_meta = $self->consensus_meta;
    my $new_cons_meta;
    if ($cons_meta) {
        $new_cons_meta = Bio::Seq::Meta->new();
        for my $meta_name ($cons_meta->meta_names) {
            $new_cons_meta->named_submeta($meta_name, 1, $self->length, $cons_meta->$meta_name);
        }
    }
    my $end=$self->length;
    for my $cat_aln (@aln) {
        my $cat_cons_meta=$cat_aln->consensus_meta;
        if ($cat_cons_meta) {
            $new_cons_meta = Bio::Seq::Meta->new() if !$new_cons_meta;
            for my $meta_name ($cat_cons_meta->meta_names) {
                $new_cons_meta->named_submeta($meta_name, $end+1, $end+$cat_aln->length, $cat_cons_meta->$meta_name);
            }
        }
        $end+=$cat_aln->length;
    } 
    $aln->consensus_meta($new_cons_meta) if $new_cons_meta;
    return $aln;
}


1;
