##############################################################################
# Bioperl module Bio::Tools::BPlite
##############################################################################
#
# The original BPlite.pm module has been written by Ian Korf !
# see http://sapiens.wustl.edu/~ikorf
#
# You may distribute this module under the same terms as perl itself

package Bio::Tools::BPlite;

use strict;
use vars qw(@ISA);

use Bio::Root::Object; # root object to inherit from
use Bio::Tools::BPlite::Sbjct; # we want to use Sbjct

@ISA = qw(Bio::Root::Object);

# _initialize is where the heavy stuff will happen when new is called

sub _initialize {
  my ($self, @args) = @_; 
  my $make = $self->SUPER::_initialize;

  my ($fh) = $self->_rearrange([qw(FH)],@args);

  if (ref $fh !~ /GLOB/)
    { $self->throw("Expecting a GLOB reference, not $fh!"); }

  $self->{FH} = $fh;
  $self->{LASTLINE} = "";
  
  if ($self->_parseHeader) {$self->{REPORT_DONE} = 0} # there are alignments
  else                     {$self->{REPORT_DONE} = 1} # empty report
  
  return $make; # success - we hope!
}

sub query    {shift->{QUERY}}

sub database {shift->{DATABASE}}

sub nextSbjct {
  my ($self) = @_;
  $self->_fastForward or return 0;
  
  #######################
  # get all sbjct lines #
  #######################
  my $def = $self->{LASTLINE};
  my $FH = $self->{FH};
  while(<$FH>) {
    if    ($_ !~ /\w/)            {next}
    elsif ($_ =~ /Strand HSP/)    {next} # WU-BLAST non-data
    elsif ($_ =~ /^\s{0,2}Score/) {$self->{LASTLINE} = $_; last}
    else                          {$def .= $_}
  }
  $def =~ s/\s+/ /g;
  $def =~ s/\s+$//g;
  $def =~ s/Length = [\d,]+$//g;
  return 0 unless $def =~ /^>/;
  $def =~ s/^>//;

  ####################
  # the Sbjct object #
  ####################
  my $sbjct = new Bio::Tools::BPlite::Sbjct(-name=>$def, 
                                            -fh=>$self->{FH}, 
					    -lastline=>$self->{LASTLINE}, 
					    -parent=>$self);
  return $sbjct;
}

sub _parseHeader {
  my ($self) = @_;
  my $FH = $self->{FH};
  
  while(<$FH>) {
    if ($_ =~ /^Query=\s+(.+)/)    {
      my $query = $1;
      while(<$FH>) {
        last if $_ !~ /\S/;
	$query .= $_;
      }
      $query =~ s/\s+/ /g;
      $query =~ s/^>//;
      $self->{QUERY} = $query;
    }
    elsif ($_ =~ /^Database:\s+(.+)/) {$self->{DATABASE} = $1}
    elsif ($_ =~ /^>/)                {$self->{LASTLINE} = $_; return 1}
    elsif ($_ =~ /^Parameters|^\s+Database:/) {
      $self->{LASTLINE} = $_;
      return 0; # there's nothing in the report
    }
  }
}

sub _fastForward {
  my ($self) = @_;
  return 0 if $self->{REPORT_DONE}; # empty report
  return 1 if $self->{LASTLINE} =~ /^>/;

  my $FH = $self->{FH};
  while(<$FH>) {
    if ($_ =~ /^>|^Parameters|^\s+Database:/) {
      $self->{LASTLINE} = $_;
      return 1;
    }
  }
  $self->warn("Possible error while parsing BLAST report!");
}

1;
__END__

=head1 NAME

Bio::Tools::BPlite - Lightweight BLAST parser

=head1 SYNOPSIS

 use Bio::Tools::BPlite;
 my $report = new BPlite(-fh=>\*STDIN);
 $report->query;
 $report->database;
 while(my $sbjct = $report->nextSbjct) {
     $sbjct->name;
     while (my $hsp = $sbjct->nextHSP) {
         $hsp->score;
         $hsp->bits;
         $hsp->percent;
         $hsp->P;
         $hsp->match;
         $hsp->positive;
         $hsp->length;
	 $hsp->querySeq;
	 $hsp->sbjctSeq;
	 $hsp->homologySeq;
	 $hsp->query->start;
	 $hsp->query->end;
	 $hsp->sbjct->start;
	 $hsp->sbjct->end;
	 $hsp->sbjct->seqname;
	 $hsp->sbjct->overlaps($exon);
     }
 }

=head1 DESCRIPTION

BPlite is a package for parsing BLAST reports. The BLAST programs are a family
of widely used algorithms for sequence database searches. The reports are
non-trivial to parse, and there are differences in the formats of the various
flavors of BLAST. BPlite parses BLASTN, BLASTP, BLASTX, TBLASTN, and TBLASTX
reports from both the high performance WU-BLAST, and the more generic
NCBI-BLAST.

Many people have developed BLAST parsers (I myself have made at least three).
BPlite is for those people who would rather not have a giant object
specification, but rather a simple handle to a BLAST report that works well
in pipes.

=head2 Object

BPlite has three kinds of objects, the report, the subject, and the HSP. To
create a new report, you pass a filehandle reference to the BPlite constructor.

 my $report = new BPlite(-fh=>\*STDIN); # or any other filehandle

The report has two attributes (query and database), and one method (nextSbjct).

 $report->query;     # access to the query name
 $report->database;  # access to the database name
 $report->nextSbjct; # gets the next subject
 while(my $sbjct = $report->nextSbjct) {
     # canonical form of use is in a while loop
 }

A subject is a BLAST hit, which should not be confused with an HSP (below). A
BLAST hit may have several alignments associated with it. A useful way of
thinking about it is that a subject is a gene and HSPs are the exons. Subjects
have one attribute (name) and one method (nextHSP).

 $sbjct->name;    # access to the subject name
 "$sbjct";        # overloaded to return name
 $sbjct->nextHSP; # gets the next HSP from the sbjct
 while(my $hsp = $sbjct->nextHSP) {
     # canonical form is again a while loop
 }

An HSP is a high scoring pair, or simply an alignment. 
HSP objects inherit all the useful methods from RangeI/SeqFeatureI/FeaturePair,
but provide an additional set of attributes (score, bits, percent, P, match, 
positive, length, querySeq, sbjctSeq, homologySeq) that should be familiar to
anyone who has seen a blast report. 

For lazy/efficient coders, two-letter abbreviations are available for the 
attributes with long names (qs, ss, hs). Ranges of the aligned sequences in
query/subject and other information (like seqname) are stored
in SeqFeature objects (i.e.: $hsp->query, $hsp->sbjct which is equal to
$hsp->feature1, $hsp->feature2). querySeq, sbjctSeq and homologySeq do only
contain the alignment sequences from the blast report.

 $hsp->score;
 $hsp->bits;
 $hsp->percent;
 $hsp->P;
 $hsp->match;
 $hsp->positive;
 $hsp->length;
 $hsp->querySeq;      $hsp->qs;
 $hsp->sbjctSeq;      $hsp->ss;
 $hsp->homologySeq;   $hsp->hs;
 $hsp->query->start;
 $hsp->query->end;
 $hsp->query->seqname;
 $hsp->sbjct->primary_tag; # "similarity"
 $hsp->sbjct->source_tag;  # "BLAST"
 $hsp->sbjct->start;
 $hsp->sbjct->end;
 ...
 "$hsp"; # overloaded for query->start..query->end bits

I've included a little bit of overloading for double quote variable
interpolation convenience. A subject will return its name and an HSP will
return its query->start, query->end, and bits in the alignment. Feel free 
to modify this to whatever is most frequently used by you.

So a very simple look into a BLAST report might look like this.

 my $report = new BPlite(-fh=>\*STDIN);
 while(my $sbjct = $report->nextSbjct) {
     print "$sbjct\n";
     while(my $hsp = $sbjct->nextHSP) {
	 	print "\t$hsp\n";
     }
 }

The output of such code might look like this:

 >foo
     100..155 29.5
     268..300 20.1
 >bar
     100..153 28.5
     265..290 22.1


=head1 AUTHORS

Ian Korf (ikorf@sapiens.wustl.edu, http://sapiens.wustl.edu/~ikorf), 
Lorenz Pollak (lorenz@ist.org, bioperl port)

=head1 ACKNOWLEDGEMENTS

This software was developed at the Genome Sequencing Center at Washington
Univeristy, St. Louis, MO.

=head1 COPYRIGHT

Copyright (C) 1999 Ian Korf. All Rights Reserved.

=head1 DISCLAIMER

This software is provided "as is" without warranty of any kind.

=cut










