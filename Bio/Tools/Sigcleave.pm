#-----------------------------------------------------------------------------
# PACKAGE : Bio::Tools::Sigcleave.pm
# AUTHOR  : Chris Dagdigian, dag@sonsorol.org
# CREATED : Jan 28 1999
# REVISION: $Id$
#            
# Copyright (c) 1997-9 bioperl, Chris Dagdigian and others. All Rights Reserved.
#           This module is free software; you can redistribute it and/or 
#           modify it under the same terms as Perl itself.
#
# _History_
#
# Object framework ripped from Steve Chervits's SeqPattern.pm
# 
# Core EGCG Sigcleave emulation from perl code developed by
# Danh Nguyen & Kamalakar Gulukota which itself was based 
# loosely on Colgrove's signal.c program.
#
# The overall idea is to replicate the output of the sigcleave
# program which was distributed with the EGCG extension to the GCG sequence
# analysis package. There is also an accessor method for just getting at
# the raw results.
#
#-----------------------------------------------------------------------------

package Bio::Tools::Sigcleave;

use Bio::Root::Global qw(:devel);
use Bio::Seq ();

@ISA = qw(Bio::Seq);
use strict;
use vars qw ($ID $VERSION %WeightTable);
$ID  = 'Bio::Tools::Sigcleave';
$VERSION = 0.01;

=head1 NAME

Bio::Tools::Sigcleave.pm - Bioperl object for sigcleave analysis

=head1 SYNOPSIS

=head2 Object Creation

    use Bio::Tools::Sigcleave ();

    $sigcleave_object = new Bio::Tools::Sigcleave(-file=>'sigtest.aa',
                                                  -desc=>'test sigcleave protein seq',
                                                  -type=>'AMINO',
                                                  -threshold=>'3.5',
                                                 );

Sigcleave objects can be created via the same methods as Bio::Seq objects. The
one additional parameter is "-threshold" which sets the score reporting limit
for the algorithim. The above exmple shows a sigcleave object being created
from a protein sequence file. See the Bio::Seq documention to see the other ways
that objects can be created.

=head2 Object Methods & Accessors

     %raw_results      = $sigcleave_object->signals;

     $formatted_output = $sigcleave_object->pretty_print;

=head1 INSTALLATION

This module is included with the central Bioperl distribution:

   http://bio.perl.org/Core/Latest
   ftp://bio.perl.org/pub/DIST

Follow the installation instructions included in the README file.

=head1 DESCRIPTION

"Sigcleave" was a program distributed as part of the free EGCG add-on to
earlier versions of the GCG Sequence Analysis package. 

From the EGCG documentation:
  SigCleave uses the von Heijne method to locate signal sequences, and to identify 
  the cleavage site. The method is 95% accurate in resolving signal sequences from 
  non-signal sequences with a cutoff score of 3.5, and 75-80% accurate in identifying 
  the cleavage site. The program reports all hits above a minimum value. 

The EGCG Sigcleave program was written by Peter Rice 
(E-mail: pmr@sanger.ac.uk Post: Informatics Division, The Sanger Centre,
Wellcome Trust Genome Campus, Hinxton, Cambs, CB10 1SA, UK). 

Since EGCG is no longer distributed for the latest versions of GCG, this code
was developed to emulate the output of the original program as much as possible for
those who lost access to sigcleave when upgrading to newer versions of GCG.

The EGCG website can be found at: http://www.sanger.ac.uk/Software/EGCG/

There are 2 accessor methods for this object. "signals" will return a perl
associative array containing the sigcleave scores keyed by amino acid position.
"pretty_print" returns a formatted string similar to the output of the original
sigcleave utility.

In both cases, the "threshold" setting controls the score reporting level. If no
value for threshold is passed in by the user, the code defaults to a reporting value
of 3.5. 

In this implemntation the accessor will never return any score/position pair which does not
meet the threshold limit. This is the slightly different from the behaviour of
the 8.1 EGCG sigcleave program which will report the highest of the under-threshold
results if nothing else is found.

Example of pretty_print output:

	SIGCLEAVE of sigtest from: 1 to 146

	Report scores over 3.5
	Maximum score 4.9 at residue 131

	 Sequence:  FVILAAMSIQGSA-NLQTQWKSTASLALET
        	    | (signal)    | (mature peptide)
          	118            131

	 Other entries above 3.5

	Maximum score 3.7 at residue 112

	 Sequence:  CSRQLFGWLFCKV-HPGAIVFVILAAMSIQGSANLQTQWKSTASLALET
         	   | (signal)    | (mature peptide)
           	99            112

=head1 USAGE

No warranty implied or expressed. Use at your own risk :) Users unfamiliar
with the original Sigcleave application should read the von Heijne papers. 

The emphasis here is on correctly replicating the calls that 8.1 EGCG sigcleave
would make. This code has been tested against a non-redundant curated set
of 405 Swissprot proteins representing secreted, non-secreted, membrane and
transit proteins. Except for the EGCG sigcleave habit of reporting an
under-threshold score if nothing better is found the output was identical.

The weight matrix in this code is for eukaryote signal sequences.

Please see the example script located in the bioperl distribution
to see how this code can be used.

=head1 FEEDBACK

When updating and maintaining a module, it helps to know that people
are actually using it. Let us know if you find a bug, think this code
is useful or have any improvements/features to suggest.

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track the bugs and 
their resolution. Bug reports can be submitted via email or the web:

    bioperl-bugs@bio.perl.org                   
    http://bio.perl.org/bioperl-bugs/           

=head1 AUTHOR

Chris Dagdigian, dag@sonsorol.org  & others

=head1 VERSION

Bio::Tools::Sigcleave.pm, $Id$

=head1 COPYRIGHT

Copyright (c) 1999 Chris Dagdigian & others. All Rights Reserved.
This module is free software; you can redistribute it and/or 
modify it under the same terms as Perl itself.

=head1 REFERENCES / SEE ALSO

von Heijne G. (1986) "A new method for predicting signal sequences cleavage sites." 
Nucleic Acids Res. 14, 4683-4690. 

von Heijne G. (1987) in "Sequence Analysis in Molecular Biology: Treasure Trove or Trivial Pursuit" 
(Acad. Press, (1987), 113-117). 

EGCG website: http://www.sanger.ac.uk/Software/EGCG/

=head1 APPENDIX

The following documentation describes the various functions
contained in this module. Some functions are for internal 
use and are not meant to be called by the user; they are 
preceded by an underscore ("_").

=cut

#
##
###
#### END of main POD documentation.
###
##
#

  %Bio::Tools::Sigcleave::WeightTable = (
   'A' => [16 ,13 ,14 ,15 ,20 ,18 ,18 ,17 ,25 ,15 ,47 , 6 ,80 ,18 , 6 ,14.5],
   'C' => [3 , 6 , 9 , 7 , 9 ,14 , 6 , 8 , 5 , 6 ,19 , 3 , 9 , 8 , 3 , 4.5],
   'D' => [0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 5 , 3 , 0 , 5 , 0 ,10 ,11 , 8.9],
   'E' => [0 , 0 , 0 , 1 , 0 , 0 , 0 , 0 , 3 , 7 , 0 , 7 , 0 ,13 ,14 ,10.0],
   'F' => [13 , 9 ,11 ,11 , 6 , 7 ,18 ,13 , 4 , 5 , 0 ,13 , 0 , 6 , 4 , 5.6],
   'G' => [4 , 4 , 3 , 6 , 3 ,13 , 3 , 2 ,19 ,34 , 5 , 7 ,39 ,10 , 7 ,12.1],
   'H' => [0 , 0 , 0 , 0 , 0 , 1 , 1 , 0 , 5 , 0 , 0 , 6 , 0 , 4 , 2 , 3.4],
   'I' => [15 ,15 , 8 , 6 ,11 , 5 , 4 , 8 , 5 , 1 ,10 , 5 , 0 , 8 , 7 , 7.4],
   'K' => [0 , 0 , 0 , 1 , 0 , 0 , 1 , 0 , 0 , 4 , 0 , 2 , 0 ,11 , 9 ,11.3],
   'L' => [71 ,68 ,72 ,79 ,78 ,45 ,64 ,49 ,10 ,23 , 8 ,20 , 1 , 8 , 4 ,12.1],
   'M' => [0 , 3 , 7 , 4 , 1 , 6 , 2 , 2 , 0 , 0 , 0 , 1 , 0 , 1 , 2 , 2.7],
   'N' => [0 , 1 , 0 , 1 , 1 , 0 , 0 , 0 , 3 , 3 , 0 ,10 , 0 , 4 , 7 , 7.1],
   'P' => [2 , 0 , 2 , 0 , 0 , 4 , 1 , 8 ,20 ,14 , 0 , 1 , 3 , 0 ,22 , 7.4],
   'Q' => [0 , 0 , 0 , 1 , 0 , 6 , 1 , 0 ,10 , 8 , 0 ,18 , 3 ,19 ,10 , 6.3],
   'R' => [2 , 0 , 0 , 0 , 0 , 1 , 0 , 0 , 7 , 4 , 0 ,15 , 0 ,12 , 9 , 7.6],
   'S' => [9 , 3 , 8 , 6 ,13 ,10 ,15 ,16 ,26 ,11 ,23 ,17 ,20 ,15 ,10 ,11.4],
   'T' => [2 ,10 , 5 , 4 , 5 ,13 , 7 , 7 ,12 , 6 ,17 , 8 , 6 , 3 ,10 , 9.7],
   'V' => [20 ,25 ,15 ,18 ,13 ,15 ,11 ,27 , 0 ,12 ,32 , 3 , 0 , 8 ,17 ,11.1],
   'W' => [4 , 3 , 3 , 1 , 1 , 2 , 6 , 3 , 1 , 3 , 0 , 9 , 0 , 2 , 0 , 1.8],
   'Y' => [0 , 1 , 4 , 0 , 0 , 1 , 3 , 1 , 1 , 2 , 0 , 5 , 0 , 1 , 7 , 5.6],
);

##
## Now we calculate the _real_ values for the weight table
##
##
## yeah yeah yeah there is lots of math here that gets repeated
## every single time a sigcleave object gets created. This is
## a quick hack to make sure that we get the scores as accurate as 
## possible. Need all those significant digits....
##
## suggestions for speedup aproaches welcome
##

my ($i, $j, $expected);
foreach $i (keys %WeightTable)
 {
    $expected = $WeightTable{$i}[15];
    if ($expected > 0)
     {
	for ($j=0; $j<16; $j++)
	  {
	    if ($WeightTable{$i}[$j] == 0)
	     {
		$WeightTable{$i}[$j] = 1; 
		if ($j == 10 || $j == 12) { $WeightTable{$i}[$j] = 1.e-10; }
	     }
	    $WeightTable{$i}[$j] = log($WeightTable{$i}[$j]/$expected); 
	  } 
     }
 }

#####################################################################################
##                                 CONSTRUCTOR                                     ##
#####################################################################################

=head1 _initialize

 Title     : _initialize
 Usage     : n/a; automatically called by Bio::Root::Object::new()
 Purpose   : Verifies that the type is correct for superclass (Bio::Seq.pm)
           : and calls superclass constructor last.
 Returns   : n/a
 Argument  : Parameters passed to new()
 Throws    : Exception if given type is not protein.
 Comments  : 
See Also   : B<Bio::Root::Object::new()>, B<Bio::Seq::_initialize()>

=cut

#----------------
sub _initialize {
#----------------
    my($self, %param) = @_;
    
    my($threshold,$type) = $self->_rearrange([qw(THRESHOLD
                                                 TYPE)], %param);

    my $make = $self->SUPER::_initialize(%param);

    # complain if not protein
    if ($type =~ /nuc|[dr]na/i) {
	$self->throw("Sigcleave.pm only supports protein sequences.");
    } elsif ($type =~ /amino|pep|prot/i) {
	$self->{type} = 'amino';
    }

    # set threshold if supplied, otherwise default to 3.5
    if (defined $threshold) {
	$self->{threshold} = $threshold;
    } else {
      $self->{threshold} = 3.5;
    }

    $self->_Analyze;
    $make;
}

=head1 _Analyze

 Title     : _Analyze
 Usage     : N/A This is an internal method. Not meant to be called from outside
           : the package
           :
 Purpose   : calculates sigcleave score and amino acid position for the
           : given protein sequence. The score reporting threshold can
           : be adjusted by passing in the "threshold" parameter during
           : object construction. If no threshold is passed in, the code
           : defaults to reporting any scores equal to or above 3.5
           :
 Returns   : nothing. results are added to the object
 Argument  : none.
 Throws    : nothing.
 Comments  : nothing.
See Also   : n/a

=cut

#----------------
sub _Analyze {
#----------------
my($self) = @_;

## need to shut strict() up

my %signals;
my @hitWeight = ();
my @hitsort   = ();
my @hitpos    = ();
my $maxSite   = "";
my $seqPos    = "";
my $istart    = "";
my $iend      = "";
my $icol      = "";
my $i         = "";
my $weight    = "";
my $k         = 0;
my $c         = 0;
my $seqBegin  = 0; 
my $pVal      = -13; 
my $nVal      = 2;
my $nHits     = 0;
my $seqEnd    = $self->length;
my $pep       = $self->seq;
my $minWeight = $self->threshold;

 ## The weight table is keyed by UPPERCASE letters so we uppercase
 ## the pep string because we don't want to alter the actual object
 ## sequence.

 $pep =~ tr/a-z/A-Z/;
 
for($seqPos = $seqBegin; $seqPos < $seqEnd; $seqPos++)
  {
    $istart = (0 > $seqPos + $pVal)? 0 : $seqPos + $pVal;
    $iend = ($seqPos + $nVal - 1 < $seqEnd)? $seqPos + $nVal - 1 : $seqEnd;
    $icol= $iend - $istart + 1; 
    $weight = 0.00;
    for ($k=0; $k<$icol; $k++)
     {
	$c = substr($pep, $istart + $k, 1);

        ## CD: The if(defined) stuff was put in here because Sigcleave.pm
        ## CD: kept getting warnings about undefined vals during 'make test' ... 
	if(defined $WeightTable{$c}[$k]) { $weight += $WeightTable{$c}[$k]; }
        
     }

    if ($weight >= $minWeight)
        {  $signals{$seqPos+1} = sprintf ("%.1f", $weight);   }
  }

  $self->{"signal_scores"} = { %signals };
}

=head1 threshold

 Title     : threshold
 Usage     : $value = $self->threshold 
           :
 Purpose   : Accessor method sigcleave score reporting threshold.
           : 
 Returns   : float.
           : 
 Argument  : none. 
 Throws    : none.
 Comments  : none.
See Also   : n/a

=cut

#----------------
sub threshold {
#----------------
my $self = shift;
return $self->{threshold};
}

=head1 signals

 Title     : signals
 Usage     : %sigcleave_results = $sigcleave_object->signals;
           :
 Purpose   : Accessor method for sigcleave results
           : 
 Returns   : Associative array. The key value represents the amino acid position
           : and the value represents the score. Only scores that
           : are greater than or equal to the THRESHOLD value are reported.
           : 
 Argument  : none. 
 Throws    : none.
 Comments  : none.
See Also   : THRESHOLD

=cut

#----------------
sub signals {
#----------------
my $self = shift;
my %results;
my $position;

 foreach $position ( sort keys %{ $self->{signal_scores} } ) {
     $results{$position} = $self->{signal_scores}{$position};    
 }

return %results;
}

=head1 pretty_print

 Title     : pretty_print
 Usage     : $output = $sigcleave_object->pretty_print;
           : print $sigcleave_object->pretty_print;
           :
 Purpose   : Emulates the output of the EGCG Sigcleave
           : utility.  
           : 
 Returns   : A formatted string.
 Argument  : none.
 Throws    : none.
 Comments  : none.
See Also   : n/a

=cut

#----------------
sub pretty_print {
#----------------
my $self = shift;
my $pos;
my $output;
my $cnt = 1;
my %results  = $self->signals;
my @hits     = keys %results;
my $hitcount = $#hits; $hitcount++;
my $thresh   = $self->threshold;
my $seqlen   = $self->length;
my $name     = $self->id;
my $pep      = $self->seq;
   $pep      =~ tr/a-z/A-Z/;

$output = "SIGCLEAVE of $name from: 1 to $seqlen\n\n";

if($hitcount > 0) {
   $output .= "Report scores over $thresh\n";
	foreach $pos ((sort { $results{$b} cmp $results{$a} } keys %results)) {
 	  my $start = $pos - 13;
          $output .= sprintf ("Maximum score %1.1f at residue %3d\n",$results{$pos},$pos);
    	  $output .= "\n";
    	  $output .= " Sequence:  ";
    	  $output .= substr($pep,$start,13);
    	  $output .= "-";
    	  $output .= substr($pep,$pos,50);
    	  $output .= "\n";
    	  $output .= "            | \(signal\)    | \(mature peptide\)\n";
    	  $output .= sprintf("          %3d            %3d\n\n",$start,$pos);

   		if(($hitcount > 1) && ($cnt == 1)) { 
    		$output .= " Other entries above $thresh\n\n";
   		} 
    	$cnt++;
 	}
}
$output;
}

1;
__END__

#########################################################################
#  End of class 
#########################################################################

