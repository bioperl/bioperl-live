
#
# BioPerl module for Bio::Tools::Sim4::Results
#
# Cared for by Ewan Birney <birney@sanger.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Tools::Sim4::Results - Results of one Sim4 run

=head1 SYNOPSIS

   $sim4 = Bio::Tools::Sim4->new( -file => 'sim4.results' );
   $sim4 = Bio::Tools::Sim4->new( -fh   => \*INPUT );
   
   # parse the results.

   foreach $exonset ( $sim4->each_ExonSet() ) {
      # $exonset is-a SeqFeature::Generic with Bio::Tools::Sim4::Exons as sub features
      print "Delimited on sequence from ", $exonset->start(), " to ", $exonset->end() "\n";
      foreach $exon ( $exonset->each_Exon() ) {
	  # $exon is-a SeqFeature::Generic 
	  print "Exon from ", $exon->start, " to ", $exon->end, "\n";
          
          # you can get out what it matched using the homol_SeqFeature attribute
          $homol = $exon->homol_SeqFeature;
          print "Matched to sequence", $homol->seqname, " at ", $homol->start, " to ", $homol->end, "\n";
      }
   }


=head1 DESCRIPTION

The sim4 module provides a parser and results object for sim4 output. The
sim4 results are specialised types of SeqFeatures, meaning you can add them
to AnnSeq objects fine, and manipulate them in the "normal" seqfeature manner.

The sim4 Exon objects are Bio::SeqFeature::Homol inherieted objects. The 
$h = $exon->homol_SeqFeature() is the seqfeature on the matching object, in
which the start/end points are where the hit lies.

To make this module work sensibly you need to run

  sim4 genomic.fasta est.database.fasta

one fiddle here is that there are only two real possibilities to the matching
criteria: either one sequence needs reversing or not. Because of this, it
is impossible to tell whether the match is in the forward or reverse strand
of the genomic dna. The exon objects have their strand attribute set to 0.
However, the homol objects have their strand attribute set to 1 or -1, depending
on whether it needs reverse complementing for this to work.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists.  Your participation is much appreciated.

  vsns-bcd-perl@lists.uni-bielefeld.de          - General discussion
  vsns-bcd-perl-guts@lists.uni-bielefeld.de     - Technically-oriented discussion
  http://bio.perl.org/MailList.html             - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via email
or the web:

  bioperl-bugs@bio.perl.org
  http://bio.perl.org/bioperl-bugs/

=head1 AUTHOR - Ewan Birney

Email birney@sanger.ac.uk

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Tools::Sim4::Results;
use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::Root::Object

use Bio::Root::Object;
use Bio::Tools::Sim4::ExonSet;
use Bio::Tools::Sim4::Exon;

@ISA = qw(Bio::Root::Object);
# new() is inherited from Bio::Root::Object

# _initialize is where the heavy stuff will happen when new is called

sub _initialize {
  my($self,@args) = @_;

  my $make = $self->SUPER::_initialize(@args);

  my($fh,$file) =
      $self->_rearrange([qw(FH
			    FILE
			    )],
			@args);

  $self->{'exon_set_array'} = [];

  if( defined $fh && defined $file ) {
      $self->throw("You have defined both a filehandle and file to read from. Not good news!");
  }

  $fh && $self->_parse_fh($fh);

  $file && $self->_parse_file($file);

# set stuff in self from @args
  return $make; # success - we hope!
}



=head2 each_ExonSet

 Title   : each_ExonSet
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub each_ExonSet{
   my ($self) = @_;

   return @{$self->{'exon_set_array'}};
}

=head2 add_ExonSet

 Title   : add_ExonSet
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub add_ExonSet{
   my ($self,$es) = @_;

   $es->isa('Bio::Tools::Sim4::ExonSet') || $self->throw("Attempting to add non ExonSet to Sim4 results. Yuk!");
   
   push(@{$self->{'exon_set_array'}},$es);
}

=head2 _parse_file

 Title   : _parse_file
 Usage   : 
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub _parse_file{
   my ($self,$file) = @_;

   open(FH,$file) || $self->throw("Could not open [$file] as a filehandle!");

   $self->_parse_fh(\*FH);

   close(FH) || $self->throw("Could not close [$file] as filehandle. Probably a pipe close, in which case sim4 might not be correctly set up on your system. Try 'which sim4' and check it works");
}


=head2 _parse_fh

 Title   : _parse_fh
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub _parse_fh{
   my ($self,$fh) = @_;

   while( <$fh> ) {
       chomp;
       # bascially, each sim4 'hit' starts with seq1...

       /^seq1/ && do {
	   # main parsing code of one hit

	   my $es = Bio::Tools::Sim4::ExonSet->new();
	   $self->add_ExonSet($es);

	   /^seq1\s+=\s+(\S+)\,\s+(\d+)/ || $self->throw("Sim4 parsing error on seq1 [$_] line. Sorry!");
	   my $filename1 = $1;
	   my $length1 = $2;

	   # get the next line
	   $_ = <$fh>; 

	   # the second hit has also the database name in the >name syntex (in brackets).
	   /^seq2\s+=\s+(\S+)\s+\(>?(\S+)\s*\)\,\s+(\d+)/ || $self->throw("Sim4 parsing error on seq2 [$_] line. Sorry!");
	   my $filename2 = $1;
	   my $name2 = $2;
	   my $length2 = $3;
	   
	   # get the next line

	   $_ = <$fh>;
	   
	   # major loop over start/end points.

	   my $hit_direction = 1;

	   while ( <$fh> ) {
	       chomp;
	       /^\(complement\)/ && do { $hit_direction = -1; next; };

	       # print STDERR "Major loop with $_\n";
	       # this matches
	       # start-end (start-end) perc%
	       /(\d+)-(\d+)\s+\((\d+)-(\d+)\)\s+(\d+)%/ && do {
		   my $f1start = $1;
		   my $f1end   = $2;
		   my $f2start = $3;
		   my $f2end   = $4;
		   my $perc = $5;
		   
		   my $homol;
		   if( $hit_direction == 1 ) {
		       $homol = Bio::SeqFeature::Generic->new(
								 -start => $f2start,
								 -end   => $f2end,
								 -strand => $hit_direction,
								 );
		   } else {
		       $homol = Bio::SeqFeature::Generic->new(
								 -end => $length2 - $f2start,
								 -start => $length2 - $f2end +1,
								 -strand => $hit_direction,
								 );
		   }
							    

		   my $exon = Bio::Tools::Sim4::Exon->new();
		   $exon->start($f1start);
		   $exon->end($f1end);
		   $exon->percentage_id($5);
		   $exon->score($5);
		   $exon->homol_SeqFeature($homol);
		   $es->add_Exon($exon);
	       };
	       if( ! /\-\>/ && ! /\<\-/ && ! /\=\=/) {
		   last;
	       }
	   }
	   next; # back to while loop
       }; # end of do on seq line

       # if a blank line, move on, otherwise issue a warning/

       /^\s*$/ || do { chomp; $self->warn("Could not understand non blank line [$_] in sim4 output. Skipping"); };
   }
}

1;
