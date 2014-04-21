
#
# BioPerl module for Bio::Tools::pSW
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Ewan Birney <birney@sanger.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Tools::pSW - pairwise Smith Waterman object

=head1 SYNOPSIS

    use Bio::Tools::pSW;
    use Bio::AlignIO;
    my $factory = Bio::Tools::pSW->new( '-matrix' => 'blosum62.bla',
				       '-gap' => 12,
				       '-ext' => 2,
				       );

    #use the factory to make some output

    $factory->align_and_show($seq1,$seq2,STDOUT);

    # make a Bio::SimpleAlign and do something with it

    my $aln = $factory->pairwise_alignment($seq1,$seq2);
    my $alnout = Bio::AlignIO->new(-format => 'msf',
				  -fh     => \*STDOUT);

    $alnout->write_aln($aln);

=head1 INSTALLATION

This module is included with the central Bioperl distribution:

   http://bio.perl.org/Core/Latest
   ftp://bio.perl.org/pub/DIST

Follow the installation instructions included in the INSTALL file.

=head1 DESCRIPTION

pSW is an Alignment Factory for protein sequences. It builds pairwise
alignments using the Smith-Waterman algorithm. The alignment algorithm is
implemented in C and added in using an XS extension. The XS extension basically
comes from the Wise2 package, but has been slimmed down to only be the
alignment part of that (this is a good thing!). The XS extension comes
from the bioperl-ext package which is distributed along with bioperl.
I<Warning:> This package will not work if you have not compiled the
bioperl-ext package.

The mixture of C and Perl is ideal for this sort of 
problem. Here are some plus points for this strategy: 

=over 2

=item Speed and Memory 

The algorithm is actually implemented in C, which means it is faster than
a pure perl implementation (I have never done one, so I have no idea
how faster) and will use considerably less memory, as it efficiently
assigns memory for the calculation.

=item Algorithm efficiency

The algorithm was written using Dynamite, and so contains an automatic
switch to the linear space divide-and-conquer method. This means you
could effectively align very large sequences without killing your machine
(it could take a while though!).

=back

=head1 FEEDBACK

=head2 Mailing Lists 

User feedback is an integral part of the evolution of this and other
Bioperl modules.  Send your comments and suggestions preferably to one
of the Bioperl mailing lists.  Your participation is much appreciated.

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
the bugs and their resolution. Bug reports can be submitted via the
web:

  https://github.com/bioperl/bioperl-live/issues           

=head1 AUTHOR

Ewan Birney, birney-at-sanger.ac.uk or birney-at-ebi.ac.uk

=head1 CONTRIBUTORS

Jason Stajich, jason-at-bioperl.org

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with an underscore "_".

=cut

# Let the code begin...

package Bio::Tools::pSW;
use strict;
no strict ( 'refs');

BEGIN {
    eval {
	require Bio::Ext::Align;
    };
    if ( $@ ) {
	die("\nThe C-compiled engine for Smith Waterman alignments (Bio::Ext::Align) has not been installed.\n Please read the install the bioperl-ext package\n\n");
	exit(1);
    }
}

use Bio::SimpleAlign;


use base qw(Bio::Tools::AlignFactory);



sub new {
  my($class,@args) = @_;
  
  my $self = $class->SUPER::new(@args);

  my($matrix,$gap,$ext) = $self->_rearrange([qw(MATRIX
						GAP
						EXT
						)],@args);
  
  #default values - we have to load matrix into memory, so 
  # we need to check it out now
  if( ! defined $matrix || !($matrix =~ /\w/) ) {
      $matrix = 'blosum62.bla';
  }

  $self->matrix($matrix); # will throw exception if it can't load it
  $self->gap(12) unless defined $gap;
  $self->ext(2) unless defined $ext;

  # I'm pretty sure I am not doing this right... ho hum...
  # This was not roght ($gap and $ext could not be 0) It is fixed now /AE
  if(  defined $gap ) {
      if( $gap =~ /^\d+$/ ) {
	  $self->gap($gap);
      } else {
	  $self->throw("Gap penalty must be a number, not [$gap]");
      }
  }
  if( defined $ext ) {
      if( $ext =~ /^\d+$/ )  {
	  $self->ext($ext);
      } else {
	  $self->throw("Extension penalty must be a number, not [$ext]");
      }
  }
 
  return $self; 
}


=head2 pairwise_alignment

 Title   : pairwise_alignment
 Usage   : $aln = $factory->pairwise_alignment($seq1,$seq2)
 Function: Makes a SimpleAlign object from two sequences
 Returns : A SimpleAlign object
 Args    :


=cut

sub pairwise_alignment{
    my ($self,$seq1,$seq2) = @_;
    my($t1,$t2,$aln,$out,@str1,@str2,@ostr1,@ostr2,$alc,$tstr,$tid,$start1,$end1,$start2,$end2,$alctemp);
    
    if( ! defined $seq1 || ! $seq1->isa('Bio::PrimarySeqI') ||
	! defined $seq2 || ! $seq2->isa('Bio::PrimarySeqI') ) {
	$self->warn("Cannot call pairwise_alignment without specifing 2 sequences (Bio::PrimarySeqI objects)");
	return;
    }
    # fix Jitterbug #1044
    if( $seq1->length() < 2 || 
	$seq2->length() < 2 ) {
	$self->warn("cannot align sequences with length less than 2");
	return;
    }
    $self->set_memory_and_report();
    # create engine objects 
    $seq1->display_id('seq1') unless ( defined $seq1->id() );
    $seq2->display_id('seq2') unless ( defined $seq2->id() );

    $t1  = &Bio::Ext::Align::new_Sequence_from_strings($seq1->id(),
						       $seq1->seq());
    $t2  = &Bio::Ext::Align::new_Sequence_from_strings($seq2->id(),
						       $seq2->seq());
    $aln = &Bio::Ext::Align::Align_Sequences_ProteinSmithWaterman($t1,$t2,$self->{'matrix'},-$self->gap,-$self->ext);
    if( ! defined $aln || $aln == 0 ) {
	$self->throw("Unable to build an alignment");
    }

    # free sequence engine objects

    $t1 = $t2 = 0;

    # now we have to get into the AlnBlock structure and
    # figure out what is aligned to what...

    # we are going to need the sequences as arrays for convience

    @str1 = split(//, $seq1->seq());
    @str2 = split(//, $seq2->seq());

    # get out start points

    # The alignment is in alignment coordinates - ie the first
    # residues starts at -1 and ends at 0. (weird I know).
    # bio-coordinates are +2 from this...

    $start1 = $aln->start()->alu(0)->start +2;
    $start2 = $aln->start()->alu(1)->start +2;

    # step along the linked list of alc units...

    for($alc = $aln->start();$alc->at_end() != 1;$alc = $alc->next()) {
	if( $alc->alu(0)->text_label eq 'SEQUENCE' ) {
	    push(@ostr1,$str1[$alc->alu(0)->start+1]);
	} else {
	    # assumme it is in insert!
	    push(@ostr1,'-');
	}

	if( $alc->alu(1)->text_label eq 'SEQUENCE' ) {
	    push(@ostr2,$str2[$alc->alu(1)->start+1]);
	} else {
	    # assumme it is in insert!
	    push(@ostr2,'-');
	}
	$alctemp = $alc;
    }

    #
    # get out end points
    #

    # end points = real residue end in 'C' coordinates = residue
    # end in biocoordinates. Oh... the wonder of coordinate systems!

    $end1 = $alctemp->alu(0)->end+1;
    $end2 = $alctemp->alu(1)->end+1;

    # get rid of the alnblock 
    $alc = 0;
    $aln = 0;

    # new SimpleAlignment
    $out = Bio::SimpleAlign->new(); # new SimpleAlignment

    $tstr = join('',@ostr1);
    $tid = $seq1->id();
    $out->add_seq(Bio::LocatableSeq->new( -seq=> $tstr,
					 -start => $start1,
					 -end   => $end1,
					 -id=>$tid ));

    $tstr = join('',@ostr2);
    $tid = $seq2->id();
    $out->add_seq(Bio::LocatableSeq->new( -seq=> $tstr,
					 -start => $start2,
					 -end => $end2,
					 -id=> $tid ));

    # give'm back the alignment

    return $out;
}

=head2 align_and_show

 Title   : align_and_show
 Usage   : $factory->align_and_show($seq1,$seq2,STDOUT)

=cut

sub align_and_show {
    my($self,$seq1,$seq2,$fh) = @_;
    my($t1,$t2,$aln,$id,$str);

if( ! defined $seq1 || ! $seq1->isa('Bio::PrimarySeqI') ||
	! defined $seq2 || ! $seq2->isa('Bio::PrimarySeqI') ) {
	$self->warn("Cannot call align_and_show without specifing 2 sequences (Bio::PrimarySeqI objects)");
	return;
    }
    # fix Jitterbug #1044
    if( $seq1->length() < 2 || 
	$seq2->length() < 2 ) {
	$self->warn("cannot align sequences with length less than 2");
	return;
    }
    if( ! defined $fh ) { 
	$fh = \*STDOUT;
    }
    $self->set_memory_and_report();
    $seq1->display_id('seq1') unless ( defined $seq1->id() );
    $seq2->display_id('seq2') unless ( defined $seq2->id() );

    $t1  = &Bio::Ext::Align::new_Sequence_from_strings($seq1->id(),$seq1->seq());

    $t2  = &Bio::Ext::Align::new_Sequence_from_strings($seq2->id(),$seq2->seq());
    $aln = &Bio::Ext::Align::Align_Sequences_ProteinSmithWaterman($t1,$t2,$self->{'matrix'},-$self->gap,-$self->ext);
    if( ! defined $aln || $aln == 0 ) {
	$self->throw("Unable to build an alignment");
    }

    &Bio::Ext::Align::write_pretty_seq_align($aln,$t1,$t2,12,50,$fh);

}

=head2 matrix

 Title     : matrix()
 Usage     : $factory->matrix('blosum62.bla');
 Function  : Reads in comparison matrix based on name
           :
 Returns   : 
 Argument  : comparison matrix

=cut

sub matrix {
    my($self,$comp) = @_;
    my $temp;

    if( !defined $comp ) {
	$self->throw("You must have a comparison matrix to set!");
    }

    # talking to the engine here...

    $temp = &Bio::Ext::Align::CompMat::read_Blast_file_CompMat($comp);

    if( !(defined $temp) || $temp == 0 ) {
	$self->throw("$comp cannot be read as a BLAST comparison matrix file");
    }

    $self->{'matrix'} = $temp;
}



=head2 gap

 Title     : gap
 Usage     : $gap = $factory->gap() #get
           : $factory->gap($value) #set
 Function  : the set get for the gap penalty
 Example   :
 Returns   : gap value 
 Arguments : new value

=cut

sub gap {
    my ($self,$val) = @_;
    

    if( defined $val ) {
	if( $val < 0 ) {    # Fixed so that gap==0 is allowed /AE
	    $self->throw("Can't have a gap penalty less than 0");
	}
	$self->{'gap'} = $val;
    }
    return $self->{'gap'};
}


=head2 ext

 Title     : ext
 Usage     : $ext = $factory->ext() #get
           : $factory->ext($value) #set
 Function  : the set get for the ext penalty
 Example   :
 Returns   : ext value 
 Arguments : new value

=cut

sub ext {
    my ($self,$val) = @_;
    
    if( defined $val ) {
	if( $val < 0 ) {    # Fixed so that gap==0 is allowed /AE
	    $self->throw("Can't have a gap penalty less than 0");
	}
	$self->{'ext'} = $val;
    }
    return $self->{'ext'};
}

1;
