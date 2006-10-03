# $Id$
##############################################################################
# Bioperl module Bio::Tools::BPlite
##############################################################################
#
# The original BPlite.pm module has been written by Ian Korf !
# see http://sapiens.wustl.edu/~ikorf
#
# You may distribute this module under the same terms as perl itself

=head1 NAME

Bio::Tools::BPlite - Lightweight BLAST parser

=head1 SYNOPSIS

 use Bio::Tools::BPlite;
 my $report = new Bio::Tools::BPlite(-fh=>\*STDIN);

  {
    $report->query;
    $report->database;
    while(my $sbjct = $report->nextSbjct) {
	$sbjct->name;
	while (my $hsp = $sbjct->nextHSP) {
	    $hsp->score;
	    $hsp->bits;
	    $hsp->percent;
	    $hsp->P;
            $hsp->EXP;
	    $hsp->match;
	    $hsp->positive;
	    $hsp->length;
	    $hsp->querySeq;
	    $hsp->sbjctSeq;
	    $hsp->homologySeq;
	    $hsp->query->start;
	    $hsp->query->end;
	    $hsp->hit->start;
	    $hsp->hit->end;
	    $hsp->hit->seq_id;
	    $hsp->hit->overlaps($exon);
	}
    }

    # the following line takes you to the next report in the stream/file
    # it will return 0 if that report is empty,
    # but that is valid for an empty blast report.
    # Returns -1 for EOF.

    last if ($report->_parseHeader == -1);
    redo;
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

 my $report = new Bio::Tools::BPlite(-fh=>\*STDIN); # or any other filehandle

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
 $sbjct->nextHSP; # gets the next HSP from the sbjct
 while(my $hsp = $sbjct->nextHSP) {
     # canonical form is again a while loop
 }

An HSP is a high scoring pair, or simply an alignment.  HSP objects
inherit all the useful methods from RangeI/SeqFeatureI/FeaturePair,
but provide an additional set of attributes (score, bits, percent, P,
match, EXP, positive, length, querySeq, sbjctSeq, homologySeq) that
should be familiar to anyone who has seen a blast report.

For lazy/efficient coders, two-letter abbreviations are available for the 
attributes with long names (qs, ss, hs). Ranges of the aligned sequences in
query/subject and other information (like seqname) are stored
in SeqFeature objects (i.e.: $hsp-E<gt>query, $hsp-E<gt>subject which is equal to
$hsp-E<gt>feature1, $hsp-E<gt>feature2). querySeq, sbjctSeq and homologySeq do only
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
 $hsp->query->seq_id;
 $hsp->hit->primary_tag; # "similarity"
 $hsp->hit->source_tag;  # "BLAST"
 $hsp->hit->start;
 $hsp->hit->end;
 ...

So a very simple look into a BLAST report might look like this.

 my $report = new Bio::Tools::BPlite(-fh=>\*STDIN);
 while(my $sbjct = $report->nextSbjct) {
     print ">",$sbjct->name,"\n";
     while(my $hsp = $sbjct->nextHSP) {
	 	print "\t",$hsp->start,"..",$hsp->end," ",$hsp->bits,"\n";
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

=head1 CONTRIBUTORS

Jason Stajich, jason@cgt.mc.duke.edu

=head1 COPYRIGHT

Copyright (C) 1999 Ian Korf. All Rights Reserved.

=head1 DISCLAIMER

This software is provided "as is" without warranty of any kind.

=cut

package Bio::Tools::BPlite;

use strict;

use Bio::Tools::BPlite::Sbjct; # we want to use Sbjct
use Symbol;

use base qw(Bio::Root::Root Bio::SeqAnalysisParserI Bio::Root::IO);

# new comes from a RootI now

=head2 new

 Title   : new
 Function: Create a new Bio::Tools::BPlite object
 Returns : Bio::Tools::BPlite
 Args    : -file     input file (alternative to -fh)
           -fh       input stream (alternative to -file)

=cut

sub new {
  my ($class, @args) = @_; 
  my $self = $class->SUPER::new(@args);
    $self->warn("Use of Bio::Tools::BPlite is deprecated. ".
                   "Use Bio::SearchIO classes instead");
  # initialize IO
  $self->_initialize_io(@args);

  $self->{'QPATLOCATION'} = [];  # Anonymous array of query pattern locations for PHIBLAST

  if ($self->_parseHeader) {$self->{'REPORT_DONE'} = 0} # there are alignments
  else                     {$self->{'REPORT_DONE'} = 1} # empty report
  
  return $self; # success - we hope!
}

# for SeqAnalysisParserI compliance

=head2 next_feature

 Title   : next_feature
 Usage   : while( my $feat = $res->next_feature ) { # do something }
 Function: SeqAnalysisParserI implementing function. This implementation
           iterates over all HSPs. If the HSPs of the current subject match
           are exhausted, it will automatically call nextSbjct().
 Example :
 Returns : A Bio::SeqFeatureI compliant object, in this case a
           Bio::Tools::BPlite::HSP object, and FALSE if there are no more
           HSPs.
 Args    : None

=cut

sub next_feature{
   my ($self) = @_;
   my ($sbjct, $hsp);
   $sbjct = $self->{'_current_sbjct'};
   unless( defined $sbjct ) {
       $sbjct = $self->{'_current_sbjct'} = $self->nextSbjct;
       return  unless defined $sbjct;
   }   
   $hsp = $sbjct->nextHSP;
   unless( defined $hsp ) {
       $self->{'_current_sbjct'} = undef;
       return $self->next_feature;
   }
   return $hsp || undef;
}

=head2 query

 Title    : query
 Usage    : $query = $obj->query();
 Function : returns the query object
 Example  :
 Returns  : query object
 Args     :

=cut

sub query    {shift->{'QUERY'}}

=head2 qlength

 Title    : qlength
 Usage    : $len = $obj->qlength();
 Function : returns the length of the query 
 Example  :
 Returns  : length of query
 Args     :

=cut

sub qlength  {shift->{'LENGTH'}}

=head2 pattern

 Title    : pattern
 Usage    : $pattern = $obj->pattern();
 Function : returns the pattern used in a PHIBLAST search

=cut

sub pattern {shift->{'PATTERN'}}

=head2 query_pattern_location

 Title    : query_pattern_location
 Usage    : $qpl = $obj->query_pattern_location();
 Function : returns reference to array of locations in the query sequence
            of pattern used in a PHIBLAST search

=cut

sub query_pattern_location {shift->{'QPATLOCATION'}}

=head2 database

 Title    : database
 Usage    : $db = $obj->database();
 Function : returns the database used in this search
 Example  :
 Returns  : database used for search
 Args     :

=cut

sub database {shift->{'DATABASE'}}

=head2 nextSbjct

 Title    : nextSbjct
 Usage    : $sbjct = $obj->nextSbjct();
 Function : Method of iterating through all the Sbjct retrieved 
            from parsing the report 
 Example  : while ( my $sbjct = $obj->nextSbjct ) {}
 Returns  : next Sbjct object or null if finished
 Args     :

=cut

sub nextSbjct {
  my ($self) = @_;
  
  $self->_fastForward or return;
  local $_;
  #######################
  # get all sbjct lines #
  #######################
  my $def = $self->_readline();  
  while(defined ($_ = $self->_readline() ) ) {
    if    (! /\w/)           {next}
    elsif (/Strand HSP/o)    {next} # WU-BLAST non-data
    elsif (/^\s{0,2}Score/o) {$self->_pushback($_); last}
    elsif (/^Histogram|^Searching|^Parameters|
            ^\s+Database:|
            ^\s+Posted date:/ox) {
	$self->_pushback($_); 
	last;
    } else {
	$def .= $_;
    }
  }
  if( ! $def ) { 
      return;
  }
  $def =~ s/\s+/ /g;
  $def =~ s/\s+$//g;
  
  my $length;
  if( $def =~ s/Length = ([\d,]+)$//g ) {
      $length = $1;
  }
  return unless $def =~ /^>/;
  $def =~ s/^>//;

  ####################
  # the Sbjct object #
  ####################
  my $sbjct = new Bio::Tools::BPlite::Sbjct('-name'=>$def,
					    '-length'=>$length,
                                            '-parent'=>$self);
  return $sbjct;
}

# begin private routines

sub _parseHeader {
  my ($self) = @_;

  # normally, _parseHeader will break out of the parse as soon as it
  # reaches a new Subject (i.e. the first one after the header) if you
  # call _parseHeader twice in a row, with nothing in between, all you
  # accomplish is a ->nextSubject call..  so we need a flag to
  # indicate that we have *entered* a header, before we are allowed to
  # leave it!

  my $header_flag = 0; # here is the flag/ It is "false" at first, and
                       # is set to "true" when any valid header element
                       # is encountered
  local $_;
  $self->{'REPORT_DONE'} = 0;  # reset this bit for a new report
  while(defined($_ = $self->_readline() ) ) {
      s/\(\s*\)//;      
      if (/^Query=(?:\s+(.+))?/s) {
	  $header_flag = 1;	# valid header element found
	  my $query = $1;
	  while( defined($_ = $self->_readline() ) ) {
	      # Continue reading query name until encountering either
	      # a line that starts with "Database" or a blank line.
	      # The latter condition is needed in order to be able to
	      # parse megablast output correctly, since Database comes
	      # before (not after) the query.
	      if( ($_ =~ /^Database/) || ($_ =~ /^$/) ) {
		  $self->_pushback($_); last;
	      }	      
	      $query .= $_;
	  }
	  $query =~ s/\s+/ /g;
	  $query =~ s/\s+$//;
	  $query =~ s/^>//;

	  my $length = 0;
	  if( $query =~ /\(([\d,]+)\s+\S+\)\s*$/ ) {      
	      $length = $1;
	      $length =~ s/,//g;
	  } else { 
	      $self->debug("length is 0 for '$query'\n");
	  }
	  $self->{'QUERY'} = $query;
	  $self->{'LENGTH'} = $length;
      }
      elsif (/^(<b>)?(T?BLAST[NPX])\s+([\w\.-]+)\s+(\[[\w-]*\])/o) { 
	  $self->{'BLAST_TYPE'} = $2; 
	  $self->{'BLAST_VERSION'} = $3;
      }				# BLAST report type - not a valid header element # JB949
      
      # Support Paracel BTK output
      elsif ( $_ =~ /(^[A-Z0-9_]+)\s+BTK\s+/ ) { 
	  $self->{'BLAST_TYPE'} = $1;
	  $self->{'BTK'} = 1;
      } 
      elsif ($_ =~ /^Database:\s+(.+)/) {$header_flag = 1;$self->{'DATABASE'} = $1} # valid header element found
      elsif ($_ =~ /^\s*pattern\s+(\S+).*position\s+(\d+)\D/) {   
	  # For PHIBLAST reports
	  $header_flag = 1;	# valid header element found
	  $self->{'PATTERN'} = $1;
	  push (@{$self->{'QPATLOCATION'}}, $2);
      } 
      elsif (($_ =~ /^>/) && ($header_flag==1)) {$self->_pushback($_); return 1} # only leave if we have actually parsed a valid header!
      elsif (($_ =~ /^Parameters|^\s+Database:/) && ($header_flag==1)) { 
      # if we entered a header, and saw nothing before the stats at the end, 
      # then it was empty
	  $self->_pushback($_);
	  return 0;		# there's nothing in the report
      } elsif( /Reference:\s+Aaron E\. Darling/ ) {
	  $self->{'BTK'} = 1;
      }  
      # bug fix suggested by MI Sadowski via Martin Lomas
      # see bug report #1118
      if( ref($self->_fh()) !~ /GLOB/ && 
	  $self->_fh()->can('EOF') && eof($self->_fh()) ) {
	  $self->warn("unexpected EOF in file\n");
	  return -1;
      }
  }
  return -1; # EOF
}

sub _fastForward {
    my ($self) = @_;
    return 0 if $self->{'REPORT_DONE'}; # empty report
    local $_;
    while(defined( $_ = $self->_readline() ) ) {
	if (/^Histogram|^Searching|^Parameters|^\s+Database:|
             ^\s+Posted date:/xo) {
	    return 0;
	} elsif( $self->{'BTK'} && /^BLAST/o ) {
	    return 0;
	} elsif( /^>/ ) {
	    $self->_pushback($_);	
	    return 1;
	}
    }
    unless( $self->{'BTK'} ) { # Paracel BTK reports have no footer
	$self->warn("Possible error (1) while parsing BLAST report!");
    }
}

1;
__END__
