# $Id$
#
# Bioperl module Bio::Tools::BPbl2seq
#	based closely on the Bio::Tools::BPlite modules
#	Ian Korf (ikorf@sapiens.wustl.edu, http://sapiens.wustl.edu/~ikorf),
#	Lorenz Pollak (lorenz@ist.org, bioperl port)
#
#
# Copyright Peter Schattner
#
# You may distribute this module under the same terms as perl itself
# _history
# October 20, 2000
# May 29, 2001
#	Fixed bug which prevented reading of more than one HSP / hit.
#	This fix required changing calling syntax as described below. (PS)
# POD documentation - main docs before the code

=head1 NAME

Bio::Tools::BPbl2seq - Lightweight BLAST parser for pair-wise sequence
alignment using the BLAST algorithm.

=head1 SYNOPSIS

  use Bio::Tools::BPbl2seq;
  my $report = Bio::Tools::BPbl2seq->new(-file => 't/bl2seq.out');
  $report->sbjctName;
  $report->sbjctLength;
  while(my $hsp = $report->next_feature) {
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
	 $hsp->sbjct->seq_id;
	 $hsp->sbjct->overlaps($exon);
 }

=head1 DESCRIPTION

B<NOTE:> This module's functionality has been implemented in
L<Bio::SearchIO::blast> and therefore is not actively maintained.

BPbl2seq is a package for parsing BLAST bl2seq reports. BLAST bl2seq
is a program for comparing and aligning two sequences using BLAST.
Although the report format is similar to that of a conventional BLAST,
there are a few differences so that BPlite is unable to read bl2seq
reports directly.

From the user's perspective, one difference between bl2seq and other
blast reports is that the bl2seq report does not print out the name of
the first of the two aligned sequences.  (The second sequence name is
given in the report as the name of the "hit").  Consequently, BPbl2seq
has no way of identifying the name of the initial sequence unless it
is passed to constructor as a second argument as in:

	my $report = Bio::Tools::BPbl2seq->new(\*FH, "ALEU_HORVU");

If the name of the first sequence (the "query") is not passed to
BPbl2seq.pm in this manner, the name of the first sequence will be
left as "unknown".  (Note that to preserve a common interface with the
other BLAST programs the two sequences being compared are referred to
in bl2seq as "query" and "subject" although this is perhaps a bit
misleading when simply comparing 2 sequences as opposed to querying a
database.)

In addition, since there will only be (at most) one "subject" (hit) in
a bl2seq report, one should use the method $report-E<gt>next_feature,
rather than $report-E<gt>nextSbjct-E<gt>nextHSP to obtain the next
high scoring pair.

One should note that the previous (0.7) version of BPbl2seq used
slightly different syntax. That version had a bug and consequently the
old syntax has been eliminated.  Attempts to use the old syntax will
return error messages explaining the (minor) recoding required to use
the current syntax.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists.  Your participation is much appreciated.

  bioperl-l@bioperl.org                  - General discussion
  http://bioperl.org/wiki/Mailing_lists  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via the
web:

  http://bugzilla.open-bio.org/

=head1 AUTHOR - Peter Schattner

Email: schattner@alum.mit.edu

=head1 ACKNOWLEDGEMENTS

Based on work of:
Ian Korf (ikorf@sapiens.wustl.edu, http://sapiens.wustl.edu/~ikorf),
Lorenz Pollak (lorenz@ist.org, bioperl port)

=head1 CONTRIBUTORS

Jason Stajich, jason-at-bioperl.org

=cut

#'
package Bio::Tools::BPbl2seq;

use strict;
use Bio::Tools::BPlite;
use Bio::Tools::BPlite::Sbjct; # we want to use Sbjct
use Symbol;

use base qw(Bio::Root::Root Bio::SeqAnalysisParserI Bio::Root::IO);

#@ISA = qw(Bio::Tools::BPlite);

=head2 new

 Title   : new
 Function: Create a new Bio::Tools::BPbl2seq object
 Returns : Bio::Tools::BPbl2seq
 Args    : -file     input file (alternative to -fh)
           -fh       input stream (alternative to -file)
           -queryname    name of query sequence
           -report_type What type of BLAST was run (blastn,blastp,tblastn...)

=cut

sub new {
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(@args);
    $self->warning("Use of Bio::Tools::BPbl2seq is deprecated".
                   "Use Bio::SearchIO classes instead");
    # initialize IO
    $self->_initialize_io(@args);

    my ($queryname,$rt) = $self->_rearrange([qw(QUERYNAME 
						REPORT_TYPE)], @args);
    $queryname = 'unknown' if( ! defined $queryname );
    if( $rt && $rt =~ /BLAST/i ) {
	$self->{'BLAST_TYPE'} = uc($rt);
    } else { 
	$self->warn("Must provide which type of BLAST was run (blastp,blastn, tblastn, tblastx, blastx) if you want strand information to get set properly for DNA query or subjects");
    }
    my $sbjct = $self->getSbjct();
    $self->{'_current_sbjct'} = $sbjct;

    $self->{'_query'}->{'NAME'} = $queryname;
    return $self;
}


=head2 getSbjct

 Title    :
 Usage    : $sbjct = $obj->getSbjct();
 Function : Method of obtaining single "subject" of a bl2seq report
 Example  : my $sbjct = $obj->getSbjct ) {}
 Returns  : Sbjct object or undef if finished
 Args     :

=cut

sub getSbjct {
  my ($self) = @_;
#  $self->_fastForward or return;

  #######################
  # get bl2seq "sbjct" name and length #
  #######################
  my $length;
  my $def;
 READLOOP: while(defined ($_ = $self->_readline) ) {
     if ($_ =~ /^>(.+)$/) {
	$def = $1;
	next READLOOP;
     }
    elsif ($_ =~ /^\s*Length\s.+\D(\d+)/i) {
	$length = $1;	
	next READLOOP;
     }
    elsif ($_ =~ /^\s{0,2}Score/) {
	$self->_pushback($_); 	
	last READLOOP;
     }
  }
  return if ! defined $def;
  $def =~ s/\s+/ /g;
  $def =~ s/\s+$//g;
  

  ####################
  # the Sbjct object #
  ####################
  my $sbjct = new Bio::Tools::BPlite::Sbjct('-name'=>$def,
					    '-length'=>$length,
					    '-parent'=>$self);
  return $sbjct;
}




=head2 next_feature

 Title   : next_feature
 Usage   : while( my $feat = $res->next_feature ) { # do something }
 Function: calls next_feature function from BPlite.
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
       $self->debug(" No hit object found for bl2seq report \n ");
       return;
   }
   $hsp = $sbjct->nextHSP;
   return $hsp || undef;
}

=head2  queryName

 Title    :
 Usage    : $name = $report->queryName();
 Function : get /set the name of the query
 Example  :
 Returns  : name of the query
 Args     :

=cut

sub  queryName {
    my ($self, $queryname) = @_;
    if( $queryname ) {
	$self->{'_query'}->{'NAME'} = $queryname;
    }
    $self->{'_query'}->{'NAME'};
}

=head2  sbjctName

 Title    :
 Usage    : $name = $report->sbjctName();
 Function : returns the name of the Sbjct
 Example  :
 Returns  : name of the Sbjct
 Args     :

=cut

sub  sbjctName {
	my $self = shift;
#	unless( defined  $self->{'_current_sbjct'} ) {
#       		my $sbjct = $self->{'_current_sbjct'} = $self->nextSbjct;
#       		return unless defined $sbjct;
#   	}
	$self->{'_current_sbjct'}->{'NAME'} || '';
}

=head2 sbjctLength

 Title    :  sbjctLength
 Usage    : $length = $report->sbjctLength();
 Function : returns the length of the Sbjct
 Example  :
 Returns  : name of the Sbjct
 Args     :

=cut

sub sbjctLength {
	my $self = shift;
#	unless( defined  $self->{'_current_sbjct'} ) {
#       		my $sbjct = $self->{'_current_sbjct'} = $self->nextSbjct;
#       		return unless defined $sbjct;
#   	}
	$self->{'_current_sbjct'}->{'LENGTH'};
}

=head2 P

 Title    : P
 Usage    :
 Function : Syntax no longer supported, error message only

=cut

sub P     {
	my $self = shift;
	$self->throw("Syntax used is no longer supported.\n  See BPbl2seq.pm documentation for current syntax.\n ");
}

=head2 percent

 Title    : percent
 Usage    : $hsp->percent();
 Function : Syntax no longer supported, error message only

=cut

sub percent  {
	my $self = shift;
	$self->throw("Syntax used is no longer supported.\n  See BPbl2seq.pm documentation for current syntax.\n ");
}

=head2 match

 Title    : match
 Usage    : $hsp->match();
 Function : Syntax no longer supported, error message only

=cut

sub match  {
	my $self = shift;
	$self->throw("Syntax used is no longer supported.\n  See BPbl2seq.pm documentation for current syntax.\n ");
}

=head2 positive

 Title    : positive
 Usage    : $hsp->positive();
 Function : Syntax no longer supported, error message only

=cut

sub positive  {
	my $self = shift;
	$self->throw("Syntax used is no longer supported.\n  See BPbl2seq.pm documentation for current syntax.\n ") ;
}

=head2 querySeq

 Title    : querySeq
 Usage    : $hsp->querySeq();
 Function : Syntax no longer supported, error message only

=cut

sub querySeq  {
	my $self = shift;
	$self->throw("Syntax used is no longer supported.\n  See BPbl2seq.pm documentation for current syntax.\n ") ;
}

=head2 sbjctSeq

 Title    : sbjctSeq
 Usage    : $hsp->sbjctSeq();
 Function : Syntax no longer supported, error message only

=cut

sub sbjctSeq  {
	my $self = shift;
	$self->throw("Syntax used is no longer supported.\n  See BPbl2seq.pm documentation for current syntax.\n ") ;
}

=head2 homologySeq

 Title    : homologySeq
 Usage    : $hsp->homologySeq();
 Function : Syntax no longer supported, error message only

=cut

sub homologySeq  {
	my $self = shift;
	$self->throw("Syntax used is no longer supported.\n  See BPbl2seq.pm documentation for current syntax.\n ") ;
}

=head2 qs

 Title    : qs
 Usage    : $hsp->qs();
 Function : Syntax no longer supported, error message only

=cut

sub qs        {
	my $self = shift;
	$self->throw("Syntax used is no longer supported.\n  See BPbl2seq.pm documentation for current syntax.\n ") ;
}

=head2 ss

 Title    : ss
 Usage    : $hsp->ss();
 Function : Syntax no longer supported, error message only

=cut

sub ss     {
	my $self = shift;
	$self->throw("Syntax used is no longer supported.\n  See BPbl2seq.pm documentation for current syntax.\n ") ;
}

=head2 hs

 Title    : hs
 Usage    : $hsp->hs();
 Function : Syntax no longer supported, error message only

=cut

sub hs   {
	my $self = shift;
	$self->throw("Syntax used is no longer supported.\n  See BPbl2seq.pm documentation for current syntax.\n ") ;
}

sub _fastForward {
    my ($self) = @_;
    return 0 if $self->{'REPORT_DONE'}; # empty report
    while(defined( $_ = $self->_readline() ) ) {
	if ($_ =~ /^>|^Parameters|^\s+Database:|^\s+Posted date:|^\s*Lambda/) {
	    $self->_pushback($_);	
	    return 1;
	}
    }
    $self->warn("Possible error (1) while parsing BLAST report!");
}

sub DESTROY { 
    my $self = shift; 
    if( defined  $self->{'_current_sbjct'} ) { 
	$self->{'_current_sbjct'}->{'PARENT'} = undef;
	$self->{'_current_sbjct'} = undef;
    }
    $self->_io_cleanup(); 
}

1;
__END__
