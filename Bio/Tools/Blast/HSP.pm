#----------------------------------------------------------------------------
# PACKAGE : Bio::Tools::Blast::HSP.pm
# AUTHOR  : Steve Chervitz (sac@bioperl.org)
# CREATED : March 1996
# STATUS  : Alpha
# REVISION: $Id$
#
# For the latest version and documentation, visit the distribution site:
#    http://genome-www.stanford.edu/perlOOP/bioperl/blast/
#
# To generate documentation, run this module through pod2html
# (preferably from Perl v5.004 or better).
#
# Copyright (c) 1996-2000 Steve Chervitz. All Rights Reserved.
#           This module is free software; you can redistribute it and/or 
#           modify it under the same terms as Perl itself.
#----------------------------------------------------------------------------

package Bio::Tools::Blast::HSP;

use Bio::Root::Global qw(:devel);
use Bio::Root::Object ();
#use Bio::Root::Err   qw(:std);

@ISA = qw( Bio::Root::Object);  

use strict;
use vars qw($ID $GAP_SYMBOL @SCORE_CUTOFFS $Revision %STRAND_SYMBOL);
$ID       = 'Bio::Tools::Blast::HSP';
$Revision = '$Id$';  #'

$GAP_SYMBOL    = '-';          # Need a more general way to handle gap symbols.
@SCORE_CUTOFFS = ( 100, 30 );  # Bit score cutoffs (see homol_score()).
%STRAND_SYMBOL = ('Plus' => 1, 'Minus' => -1);

## POD Documentation:

=head1 NAME

Bio::Tools::Blast::HSP.pm - Bioperl BLAST High-Scoring Segment Pair object

=head1 SYNOPSIS

=head2 Object Creation

The construction of HSP objects is handled by Bio::Tools::Blast:: Sbjct.pm. 
You should not need to use this package directly. See L<_initialize()|_initialize>
for a description of constructor parameters.

    require Bio::Tools::Blast::HSP;

    $hspObj = eval{ new Bio::Tools::Blast::HSP(-DATA    =>\@hspData, 
					       -PARENT  =>$sbjct_object, 
					       -NAME    =>$hspCount,
 		                               -PROGRAM =>'TBLASTN',
					       );
		};

@hspData includes the raw BLAST report data for a specific HSP,
and is prepared by Bio::Tools::Blast::Sbjct.pm. 

=head1 INSTALLATION

This module is included with the central Bioperl distribution:

   http://bio.perl.org/Core/Latest
   ftp://bio.perl.org/pub/DIST

Follow the installation instructions included in the README file.

=head1 DESCRIPTION

The Bio::Tools::Blast::HSP.pm module encapsulates data and methods for 
manipulating, parsing, and analyzing HSPs ("High-scoring Segment Pairs")
derived from BLAST sequence analysis.

This module is a utility module used by the B<Bio::Tools::Blast::Sbjct.pm>
and is not intended for separate use. Please see documentation for 
B<Bio::Tools::Blast.pm> for some basic information about using
HSP objects (L<Links:>).

=over 0

=item * Supports BLAST versions 1.x and 2.x, gapped and ungapped.

=back

Bio::Tools::Blast::HSP.pm has the ability to extract a list of all
residue indices for identical and conservative matches along both
query and sbjct sequences. Since this degree of detail is not always
needed, this behavior does not occur during construction of the HSP
object.  These data will automatically be collected as necessary as
the HSP.pm object is used.

=head1 DEPENDENCIES

Bio::Tools::Blast::HSP.pm is a concrete class that inherits from
B<Bio::Root::Object.pm> and relies on B<Bio::Tools::Sbjct.pm> as a
container for HSP.pm objects.  B<Bio::Seq.pm> and B<Bio::UnivAln.pm>
are employed for creating sequence and alignment objects,
respectively.


=head2 Relationship to UnivAln.pm & Seq.pm

HSP.pm can provide the query or sbjct sequence as a B<Bio::Seq.pm>
object via the L<seq()|seq> method. The HSP.pm object can also create a
two-sequence B<Bio::UnivAln.pm> alignment object using the the query
and sbjct sequences via the L<get_aln()|get_aln> method. Creation of alignment
objects is not automatic when constructing the HSP.pm object since
this level of functionality is not always required and would generate
a lot of extra overhead when crunching many reports.


=head1 FEEDBACK

=head2 Mailing Lists 

User feedback is an integral part of the evolution of this and other
Bioperl modules.  Send your comments and suggestions preferably to one
of the Bioperl mailing lists.  Your participation is much appreciated.

    bioperl-l@bioperl.org          - General discussion
    http://bio.perl.org/MailList.html             - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution. Bug reports can be submitted via email
or the web:

    bioperl-bugs@bio.perl.org                   
    http://bio.perl.org/bioperl-bugs/           

=head1 AUTHOR

Steve Chervitz, E<lt>sac@bioperl.orgE<gt>

=head1 SEE ALSO

 Bio::Tools::Blast::Sbjct.pm   - Blast hit object.
 Bio::Tools::Blast.pm          - Blast object.
 Bio::Seq.pm                   - Biosequence object  
 Bio::UnivAln.pm               - Biosequence alignment object.
 Bio::Root::Object.pm          - Proposed base class for all Bioperl objects.

=head2 Links:

 http://bio.perl.org/Core/POD/Tools/Blast/Sbjct.pm.html

 http://bio.perl.org/Projects/modules.html  - Online module documentation
 http://bio.perl.org/Projects/Blast/        - Bioperl Blast Project     
 http://bio.perl.org/                       - Bioperl Project Homepage

=head1 COPYRIGHT

Copyright (c) 1996-98 Steve Chervitz. All Rights Reserved.
This module is free software; you can redistribute it and/or 
modify it under the same terms as Perl itself.

=cut



#
##
###
#### END of main POD documentation.
###
##
#

=head1 APPENDIX

Methods beginning with a leading underscore are considered private
and are intended for internal use by this module. They are
B<not> considered part of the public interface and are described here
for documentation purposes only.

=cut

#####################################################################################
##                                 CONSTRUCTOR                                     ##
#####################################################################################

=head2 _initialize

 Usage     : n/a; automatically called by Bio::Root::Object::new()
           : Bio::Tools::Blast::HSP.pm objects are constructed 
           : automatically by Bio::Tools::Sbjct.pm, so there is no need
           : for direct consumption.
 Purpose   : Initializes HSP data and calls private methods to extract 
           : the data for a given HSP.
           : Calls superclass constructor first (Bio::Root::Object.pm).
 Returns   : n/a
 Argument  : Named parameters passed from new():
           : All tags must be uppercase (does not call _rearrange()).
           :      -DATA    => array ref containing raw data for one HSP.
           :      -PARENT  => Sbjct.pm object ref.
           :      -NAME    => integer (1..n).
           :      -PROGRAM => string ('TBLASTN', 'BLASTP', etc.).

See Also   : L<_set_data()|_set_data>, B<Bio::Root::Object::new()>, B<Bio::Tools::Blast::Sbjct::_set_hsps()>

=cut

#----------------
sub _initialize {
#----------------
    my( $self, %param ) = @_;
    
    $self->SUPER::_initialize( %param );

    # The gapped and program booleans may be needed after the HSP object
    # is built.
#    $self->{'_gapped'} = $param{-GAPPED} || 0;
    $self->{'_prog'} = $param{-PROGRAM} || 0;  
    $self->_set_data( @{$param{-DATA}} );
}

#####################################################################################
##                                  ACCESSORS                                      ##
#####################################################################################


=head2 _set_data

 Usage     : n/a; called automatically during object construction.
 Purpose   : Sets the query sequence, sbjct sequence, and the "match" data
           : which consists of the symbols between the query and sbjct lines
           : in the alignment.
 Argument  : Array (all lines from a single, complete HSP, one line per element)
 Throws    : Propagates any exceptions from the methods called ("See Also")

See Also   : L<_set_seq()|_set_seq>, L<_set_residues()|_set_residues>, L<_set_score_stats()|_set_score_stats>, L<_set_match_stats()|_set_match_stats>, L<_initialize()|_initialize>

=cut

#--------------
sub _set_data {
#--------------
    my $self = shift;
    my @data = @_;
    my @queryList  = ();  # 'Query' = SEQUENCE USED TO QUERY THE DATABASE.
    my @sbjctList  = ();  # 'Sbjct' = HOMOLOGOUS SEQUENCE FOUND IN THE DATABASE.
    my @matchList  = ();
    my $matchLine  = 0;   # Alternating boolean: when true, load 'match' data.
    my @linedat = ();
    
    $DEBUG and print STDERR "$ID: set_data()\n";

    my($line, $aln_row_len, $length_diff);
    $length_diff = 0;

    # Collecting data for all lines in the alignment
    # and then storing the collections for possible processing later.
    #
    # Note that "match" lines may not be properly padded with spaces.
    # This loop now properly handles such cases:
    # Query: 1141 PSLVELTIRDCPRLEVGPMIRSLPKFPMLKKLDLAVANIIEEDLDVIGSLEELVIXXXXX 1200
    #             PSLVELTIRDCPRLEVGPMIRSLPKFPMLKKLDLAVANIIEEDLDVIGSLEELVI
    # Sbjct: 1141 PSLVELTIRDCPRLEVGPMIRSLPKFPMLKKLDLAVANIIEEDLDVIGSLEELVILSLKL 1200

    foreach $line( @data ) {
	next if $line =~ /^\s*$/;

	if( $line =~ /^ ?Score/ ) {
	    $self->_set_score_stats( $line );
	} elsif( $line =~ /^ ?(Identities|Positives|Strand)/ ) {
	    $self->_set_match_stats( $line );
	} elsif( $line =~ /^ ?Frame = ([\d+-]+)/ ) {
	    # Version 2.0.8 has Frame information on a separate line.
	   $self->{'_frame'} = $1;
	} elsif( $line =~ /^(Query:?[\s\d]+)([^\s\d]+)/ ) {
	    push @queryList, $line;
	    $self->{'_match_indent'} = CORE::length $1;
	    $aln_row_len = (CORE::length $1) + (CORE::length $2);
	    $matchLine = 1;
	} elsif( $matchLine ) {
	    # Pad the match line with spaces if necessary.
	    $length_diff = $aln_row_len - CORE::length $line;
	    $length_diff and $line .= ' 'x $length_diff;
	    push @matchList, $line;
	    $matchLine = 0;
	} elsif( $line =~ /^Sbjct/ ) {
	    push @sbjctList, $line;
	}
    }

    # Storing the query and sbjct lists in case they are needed later.
    # We could make this conditional to save memory.
    $self->{'_queryList'} = \@queryList; 
    $self->{'_sbjctList'} = \@sbjctList; 

    # Storing the match list in case it is needed later.
    $self->{'_matchList'} = \@matchList;

    if(not defined ($self->{'_numIdentical'})) {
      $self->throw("Can't parse match statistics.",
		   "Possibly a new or unrecognized Blast format.");
    }

    if(!scalar @queryList or !scalar @sbjctList) {
        $self->throw("Can't find query or sbjct alignment lines.",
	  	     "Possibly unrecognized Blast format.");
      }
}



=head2 _set_score_stats

 Usage     : n/a; called automatically by _set_data()
 Purpose   : Sets various score statistics obtained from the HSP listing.
 Argument  : String with any of the following formats:
           : blast2:  Score = 30.1 bits (66), Expect = 9.2
           : blast2:  Score = 158.2 bits (544), Expect(2) = e-110
           : blast1:  Score = 410 (144.3 bits), Expect = 1.7e-40, P = 1.7e-40
           : blast1:  Score = 55 (19.4 bits), Expect = 5.3, Sum P(3) = 0.99
 Throws    : Exception if the stats cannot be parsed, probably due to a change
           : in the Blast report format.

See Also   : L<_set_data()|_set_data>

=cut

#--------------------
sub _set_score_stats {
#--------------------
    my ($self, $data) = @_;

    my ($expect, $p);

    if($data =~ /Score = +([\d.e+-]+) bits \(([\d.e+-]+)\), +Expect = +([\d.e+-]+)/) {
	# blast2 format n = 1
	$self->{'_bits'}   = $1;
	$self->{'_score'}  = $2;
	$expect            = $3;
    } elsif($data =~ /Score = +([\d.e+-]+) bits \(([\d.e+-]+)\), +Expect\((\d+)\) = +([\d.e+-]+)/) {
	# blast2 format n > 1
	$self->{'_bits'}   = $1;
	$self->{'_score'}  = $2;
	$self->{'_n'}      = $3;
	$expect            = $4;

    } elsif($data =~ /Score = +([\d.e+-]+) \(([\d.e+-]+) bits\), +Expect = +([\d.e+-]+), P = +([\d.e-]+)/) {
	# blast1 format, n = 1
	$self->{'_score'}  = $1;
	$self->{'_bits'}   = $2;
	$expect            = $3;
	$p                 = $4;

    } elsif($data =~ /Score = +([\d.e+-]+) \(([\d.e+-]+) bits\), +Expect = +([\d.e+-]+), +Sum P\((\d+)\) = +([\d.e-]+)/) {
	# blast1 format, n > 1
	$self->{'_score'}  = $1;
	$self->{'_bits'}   = $2;
	$expect            = $3;
	$self->{'_n'}      = $4;
	$p                 = $5;

    } else {
	$self->throw("Can't parse score statistics: unrecognized format.", "$data");
    }

    $expect = "1$expect" if $expect =~ /^e/i;    
    $p      = "1$p"      if defined $p and $p=~ /^e/i; 

    $self->{'_expect'} = $expect;
    $self->{'_p'}      = $p || undef;

}



=head2 _set_match_stats

 Usage     : n/a; called automatically by _set_data()
 Purpose   : Sets various matching statistics obtained from the HSP listing.
 Argument  : blast2: Identities = 23/74 (31%), Positives = 29/74 (39%), Gaps = 17/74 (22%)
           : blast2: Identities = 57/98 (58%), Positives = 74/98 (75%)
           : blast1: Identities = 87/204 (42%), Positives = 126/204 (61%)
           : blast1: Identities = 87/204 (42%), Positives = 126/204 (61%), Frame = -3
           : WU-blast: Identities = 310/553 (56%), Positives = 310/553 (56%), Strand = Minus / Plus
 Throws    : Exception if the stats cannot be parsed, probably due to a change
           : in the Blast report format.
 Comments  : The "Gaps = " data in the HSP header has a different meaning depending
           : on the type of Blast: for BLASTP, this number is the total number of
           : gaps in query+sbjct; for TBLASTN, it is the number of gaps in the
           : query sequence only. Thus, it is safer to collect the data
           : separately by examining the actual sequence strings as is done
           : in _set_seq().

See Also   : L<_set_data()|_set_data>, L<_set_seq()|_set_seq>

=cut

#--------------------
sub _set_match_stats {
#--------------------
    my ($self, $data) = @_;

    if($data =~ m!Identities = (\d+)/(\d+)!) {
      # blast1 or 2 format
      $self->{'_numIdentical'} = $1;
      $self->{'_totalLength'}  = $2;
    }
    
    if($data =~ m!Positives = (\d+)/(\d+)!) {
      # blast1 or 2 format
      $self->{'_numConserved'} = $1;
      $self->{'_totalLength'}  = $2;
    }
    
    if($data =~ m!Frame = ([\d+-]+)!) { 
      $self->{'_frame'} = $1; 
    }

    # Strand data is not always present in this line.
    # _set_seq() will also set strand information.
    if($data =~ m!Strand = (\w+) / (\w+)!) { 
	$self->{'_queryStrand'} = $1; 
	$self->{'_sbjctStrand'} = $2; 
    }

#    if($data =~ m!Gaps = (\d+)/(\d+)!) {
#	 $self->{'_totalGaps'} = $1;
#    } else {
#	 $self->{'_totalGaps'} = 0;
#    }
}



=head2 _set_seq_data

 Usage     : n/a; called automatically when sequence data is requested.
 Purpose   : Sets the HSP sequence data for both query and sbjct sequences.
           : Includes: start, stop, length, gaps, and raw sequence.
 Argument  : n/a
 Throws    : Propagates any exception thrown by _set_match_seq()
 Comments  : Uses raw data stored by _set_data() during object construction.
           : These data are not always needed, so it is conditionally
           : executed only upon demand by methods such as gaps(), _set_residues(),
           : etc. _set_seq() does the dirty work.

See Also   : L<_set_seq()|_set_seq>

=cut

sub _set_seq_data {
    my $self = shift;

    $self->_set_seq('query', @{$self->{'_queryList'}});
    $self->_set_seq('sbjct', @{$self->{'_sbjctList'}});

    # Liberate some memory.
    @{$self->{'_queryList'}} = @{$self->{'_sbjctList'}} = ();
    undef $self->{'_queryList'};
    undef $self->{'_sbjctList'};

    $self->{'_set_seq_data'} = 1;
}



=head2 _set_seq

 Usage     : n/a; called automatically by _set_seq_data()
           : $hsp_obj->($seq_type, @data);
 Purpose   : Sets sequence information for both the query and sbjct sequences.
           : Directly counts the number of gaps in each sequence (if gapped Blast).
 Argument  : $seq_type = 'query' or 'sbjct'
           : @data = all seq lines with the form:
           : Query: 61  SPHNVKDRKEQNGSINNAISPTATANTSGSQQINIDSALRDRSSNVAAQPSLSDASSGSN 120
 Throws    : Exception if data strings cannot be parsed, probably due to a change
           : in the Blast report format.
 Comments  : Uses first argument to determine which data members to set
           : making this method sensitive data member name changes.
           : Behavior is dependent on the type of BLAST analysis (TBLASTN, BLASTP, etc).
 Warning   : Sequence endpoints are normalized so that start < end. This affects HSPs
           : for TBLASTN/X hits on the minus strand. Normalization facilitates use
           : of range information by methods such as match().

See Also   : L<_set_seq_data()|_set_seq_data>, L<matches()|matches>, L<range()|range>, L<start()|start>, L<end()|end>

=cut

#-------------
sub _set_seq {
#-------------
    my $self      = shift;
    my $seqType   = shift;
    my @data      = @_;
    my @ranges    = ();
    my @sequence  = ();
    my $numGaps   = 0;

    foreach( @data ) {
	if( m/(\d+) *(\D+) *(\d+)/) {
	    push @ranges, ( $1, $3 ) ;
	    push @sequence, $2;
	} else {
	    $self->warn("Bad sequence data: $_");
	}
    }

    (scalar(@sequence) and scalar(@ranges)) || $self->throw("Can't set sequence: missing data",
							    "Possibly unrecognized Blast format.");
 
    # Sensitive to member name changes.
    $seqType = "_\L$seqType\E";
    $self->{$seqType.'Start'} = $ranges[0];
    $self->{$seqType.'Stop'}  = $ranges[ $#ranges ];
    $self->{$seqType.'Seq'}   = \@sequence;
	
    $self->{$seqType.'Length'} = abs($ranges[ $#ranges ] - $ranges[0]) + 1;

    # Adjust lengths for BLASTX, TBLASTN, TBLASTX sequences
    # Converting nucl coords to amino acid coords.

    my $prog = $self->{'_prog'};
    if($prog eq 'TBLASTN' and $seqType eq '_sbjct') {
	$self->{$seqType.'Length'} /= 3;
    } elsif($prog eq 'BLASTX' and $seqType eq '_query') {
	$self->{$seqType.'Length'} /= 3;
    } elsif($prog eq 'TBLASTX') {
	$self->{$seqType.'Length'} /= 3;
    }

    $self->{$seqType.'Strand'} = 'Plus' if $prog =~ /BLAST[NX]/;

    # Normalize sequence endpoints so that start < end.
    # Reverse complement or 'minus strand' HSPs get flipped here.
    if($self->{$seqType.'Start'} > $self->{$seqType.'Stop'}) {
	($self->{$seqType.'Start'}, $self->{$seqType.'Stop'}) = 
	    ($self->{$seqType.'Stop'}, $self->{$seqType.'Start'});
	$self->{$seqType.'Strand'} = 'Minus';
    } 

    ## Count number of gaps in each seq. Only need to do this for gapped Blasts.
#    if($self->{'_gapped'}) {
	my $seqstr = join('', @sequence);
	$seqstr =~ s/\s//g;
        my $num_gaps = CORE::length($seqstr) - $self->{$seqType.'Length'};
	$self->{$seqType.'Gaps'} = $num_gaps if $num_gaps > 0;
#    }
}


=head2 _set_residues

 Usage     : n/a; called automatically when residue data is requested.
 Purpose   : Sets the residue numbers representing the identical and
           : conserved positions. These data are obtained by analyzing the
           : symbols between query and sbjct lines of the alignments.
 Argument  : n/a
 Throws    : Propagates any exception thrown by _set_seq_data() and _set_match_seq().
 Comments  : These data are not always needed, so it is conditionally
           : executed only upon demand by methods such as seq_inds().
           : Behavior is dependent on the type of BLAST analysis (TBLASTN, BLASTP, etc).

See Also   : L<_set_seq_data()|_set_seq_data>, L<_set_match_seq()|_set_match_seq>, L<seq_inds()|seq_inds>

=cut

#------------------
sub _set_residues {
#------------------
    my $self      = shift;
    my @sequence  = ();

    $self->_set_seq_data() unless $self->{'_set_seq_data'};

    # Using hashes to avoid saving duplicate residue numbers.
    my %identicalList_query = ();
    my %identicalList_sbjct = ();
    my %conservedList_query = ();
    my %conservedList_sbjct = ();
    
    my $aref = $self->_set_match_seq() if not ref $self->{'_matchSeq'};
    $aref  ||= $self->{'_matchSeq'};
    my $seqString = join('', @$aref );

    my $qseq = join('',@{$self->{'_querySeq'}});
    my $sseq = join('',@{$self->{'_sbjctSeq'}});
    my $resCount_query = $self->{'_queryStop'} || 0;
    my $resCount_sbjct = $self->{'_sbjctStop'} || 0;

    my $prog = $self->{'_prog'};
    if($prog !~ /^BLASTP|^BLASTN/) {
	if($prog eq 'TBLASTN') {
	    $resCount_sbjct /= 3;
	} elsif($prog eq 'BLASTX') {
	    $resCount_query /= 3;
	} elsif($prog eq 'TBLASTX') {
	    $resCount_query /= 3;
	    $resCount_sbjct /= 3;
	}
    }

    my ($mchar, $schar, $qchar);
    while( $mchar = chop($seqString) ) {
	($qchar, $schar) = (chop($qseq), chop($sseq));
	if( $mchar eq '+' ) { 
	    $conservedList_query{ $resCount_query } = 1; 
	    $conservedList_sbjct{ $resCount_sbjct } = 1; 
	} elsif( $mchar ne ' ' ) { 
	    $identicalList_query{ $resCount_query } = 1; 
	    $identicalList_sbjct{ $resCount_sbjct } = 1;
	}
	$resCount_query-- if $qchar ne $GAP_SYMBOL;
	$resCount_sbjct-- if $schar ne $GAP_SYMBOL;
    }
    $self->{'_identicalRes_query'} = \%identicalList_query;
    $self->{'_conservedRes_query'} = \%conservedList_query;
    $self->{'_identicalRes_sbjct'} = \%identicalList_sbjct;
    $self->{'_conservedRes_sbjct'} = \%conservedList_sbjct;

}




=head2 _set_match_seq

 Usage     : n/a. Internal method.
           : $hsp_obj->_set_match_seq()
 Purpose   : Set the 'match' sequence for the current HSP (symbols in between
	   : the query and sbjct lines.)				
 Returns   : Array reference holding the match sequences lines.
 Argument  : n/a
 Throws    : Exception if the _matchList field is not set.
 Comments  : The match information is not always necessary. This method
           : allows it to be conditionally prepared.
           : Called by _set_residues>() and seq_str().

See Also   : L<_set_residues()|_set_residues>, L<seq_str()|seq_str>

=cut

#-------------------
sub _set_match_seq {
#-------------------
    my $self = shift;

## DEBUGGING CODE:
#    if($self->parent->name eq '1AK5_' and $self->parent->parent->name eq 'YAR073W') {
#	print "\n_set_match_seq() called for HSP ", $self->name, " of hit ${\$self->parent->name} in query ${\$self->parent->parent->name}"; <STDIN>;
 #   }

    ref($self->{'_matchList'}) || $self->throw("Can't set HSP match sequence: No data");
    
    my @data = @{$self->{'_matchList'}};

    my(@sequence);
    foreach( @data ) {
	chomp($_);
	## Remove leading spaces; (note: aln may begin with a space
	## which is why we can't use s/^ +//).
	s/^ {$self->{'_match_indent'}}//;   
	push @sequence, $_;
    }
    # Liberate some memory.
    @{$self->{'_matchList'}} = undef;
    $self->{'_matchList'} = undef;

    $self->{'_matchSeq'} = \@sequence;

## DEBUGGING CODE:
#    if($self->parent->name eq '1AK5_' and $self->parent->parent->name eq 'YAR073W') {
#	print "RETURNING: $self->{'_matchSeq'}:\n   @{$self->{'_matchSeq'}}";<STDIN>;
#    }

    $self->{'_matchSeq'};
}



=head2 score

 Usage     : $hsp_obj->score()
 Purpose   : Get the Blast score for the HSP.
 Returns   : Integer
 Argument  : n/a
 Throws    : n/a

See Also   : L<bits()|bits>

=cut

#---------
sub score { my $self = shift; $self->{'_score'}; }
#---------



=head2 bits

 Usage     : $hsp_obj->bits()
 Purpose   : Get the Blast score in bits for the HSP.
 Returns   : Float
 Argument  : n/a
 Throws    : n/a


See Also   : L<score()|score>

=cut

#--------
sub bits { my $self = shift; $self->{'_bits'}; }
#--------



=head2 n

 Usage     : $hsp_obj->n()
 Purpose   : Get the N value (num HSPs on which P/Expect is based).
           : This value is not defined with NCBI Blast2 with gapping.
 Returns   : Integer or null string if not defined.
 Argument  : n/a
 Throws    : n/a
 Comments  : The 'N' value is listed in parenthesis with P/Expect value:
           : e.g., P(3) = 1.2e-30  ---> (N = 3).
           : Not defined in NCBI Blast2 with gaps.
           : This typically is equal to the number of HSPs but not always.
           : To obtain the number of HSPs, use Bio::Tools::Blast::Sbjct::num_hsps().

See Also   : L<score()|score>

=cut

#-----
sub n { my $self = shift; $self->{'_n'} || ''; }
#-----



=head2 frame

 Usage     : $hsp_obj->frame()
 Purpose   : Get the reading frame number (-/+ 1, 2, 3) (TBLASTN/X only).
 Returns   : Integer or null string if not defined.
 Argument  : n/a
 Throws    : n/a

=cut

#---------
sub frame { my $self = shift; $self->{'_frame'} || ''; }
#---------



=head2 signif()

 Usage     : $hsp_obj->signif()
 Purpose   : Get the P-value or Expect value for the HSP.
 Returns   : Float (0.001 or 1.3e-43)
           : Returns P-value if it is defined, otherwise, Expect value.
 Argument  : n/a
 Throws    : n/a
 Comments  : Provided for consistency with Sbjct::signif()
           : Support for returning the significance data in different
           : formats (e.g., exponent only), is not provided for HSP objects.
           : This is only available for the Sbjct or Blast object.

See Also   : L<p()|p>, L<expect()|expect>, B<Bio::Tools::Blast::Sbjct::signif()>

=cut

#-----------
sub signif { 
#-----------
    my $self = shift; 
    my $val ||= defined($self->{'_p'}) ? $self->{'_p'} : $self->{'_expect'};
    $val; 
}



=head2 expect

 Usage     : $hsp_obj->expect()
 Purpose   : Get the Expect value for the HSP.
 Returns   : Float (0.001 or 1.3e-43)
 Argument  : n/a
 Throws    : n/a
 Comments  : Support for returning the expectation data in different
           : formats (e.g., exponent only), is not provided for HSP objects.
           : This is only available for the Sbjct or Blast object.

See Also   : L<p()|p>

=cut

#----------
sub expect { my $self = shift; $self->{'_expect'}; }
#----------



=head2 p

 Usage     : $hsp_obj->p()
 Purpose   : Get the P-value for the HSP.
 Returns   : Float (0.001 or 1.3e-43) or undef if not defined.
 Argument  : n/a
 Throws    : n/a
 Comments  : P-value is not defined with NCBI Blast2 reports.
           : Support for returning the expectation data in different
           : formats (e.g., exponent only) is not provided for HSP objects.
           : This is only available for the Sbjct or Blast object.

See Also   : L<expect()|expect>

=cut

#-----
sub p { my $self = shift; $self->{'_p'}; }
#-----


=head2 length

 Usage     : $hsp->length( [seq_type] )
 Purpose   : Get the length of the aligned portion of the query or sbjct.
 Example   : $hsp->length('query')
 Returns   : integer
 Argument  : seq_type: 'query' | 'sbjct' | 'total'  (default = 'total')
 Throws    : n/a
 Comments  : 'total' length is the full length of the alignment
           : as reported in the denominators in the alignment section: 
           : "Identical = 34/120 Positives = 67/120".
           : Developer note: when using the built-in length function within
           : this module, call it as CORE::length().

See Also   : L<gaps()|gaps>

=cut

#-----------
sub length {
#-----------
    my( $self, $type ) = @_;
    $type  ||= 'total';

    $type ne 'total' and $self->_set_seq_data() unless $self->{'_set_seq_data'};

    ## Sensitive to member name format.
    $type = "_\L$type\E";
    $self->{$type.'Length'};
}



=head2 gaps

 Usage     : $hsp->gaps( [seq_type] )
 Purpose   : Get the number of gaps in the query, sbjct, or total alignment.
           : Also can return query gaps and sbjct gaps as a two-element list
           : when in array context.
 Example   : $total_gaps      = $hsp->gaps();
           : ($qgaps, $sgaps) = $hsp->gaps();
           : $qgaps           = $hsp->gaps('query');
 Returns   : scalar context: integer
           : array context without args: (int, int)  =  ('queryGaps', 'sbjctGaps')
 Argument  : seq_type: 'query' | 'sbjct' | 'total'  
           : (default = 'total', scalar context)
           : Array context can be "induced" by providing an argument of 'list' or 'array'.
 Throws    : n/a

See Also   : L<length()|length>, L<matches()|matches>

=cut

#---------
sub gaps {
#---------
    my( $self, $seqType ) = @_;

    $self->_set_seq_data() unless $self->{'_set_seq_data'};

    $seqType  ||= (wantarray ? 'list' : 'total');

    if($seqType =~ /list|array/i) {
	return (($self->{'_queryGaps'} || 0), ($self->{'_sbjctGaps'} || 0));
    }

    if($seqType eq 'total') {
	return ($self->{'_queryGaps'} + $self->{'_sbjctGaps'}) || 0;
    } else {
	## Sensitive to member name format.
	$seqType = "_\L$seqType\E";
	return $self->{$seqType.'Gaps'} || 0;
    }
}



=head2 matches

 Usage     : $hsp->matches([seq_type], [start], [stop]);
 Purpose   : Get the total number of identical and conservative matches 
           : in the query or sbjct sequence for the given HSP. Optionally can
           : report data within a defined interval along the seq.
           : (Note: 'conservative' matches are called 'positives' in the
	   : Blast report.)
 Example   : ($id,$cons) = $hsp_object->matches('sbjct');
           : ($id,$cons) = $hsp_object->matches('query',300,400);
 Returns   : 2-element array of integers 
 Argument  : (1) seq_type = 'query' | 'sbjct' (default = query)
           : (2) start = Starting coordinate (optional)
           : (3) stop  = Ending coordinate (optional)
 Throws    : Exception if the supplied coordinates are out of range.
 Comments  : Relies on seq_str('match') to get the string of alignment symbols
           : between the query and sbjct lines which are used for determining
           : the number of identical and conservative matches.

See Also   : L<length()|length>, L<gaps()|gaps>, L<seq_str()|seq_str>, B<Bio::Tools::Blast::Sbjct::_adjust_contigs()>

=cut

#-----------
sub matches {
#-----------
    my( $self, %param ) = @_;
    my(@data);
    my($seqType, $beg, $end) = ($param{-SEQ}, $param{-START}, $param{-STOP});
    $seqType ||= 'query';

    if(!defined $beg && !defined $end) {
	## Get data for the whole alignment.
	push @data, ($self->{'_numIdentical'}, $self->{'_numConserved'});
    } else {
	## Get the substring representing the desired sub-section of aln.
	$beg ||= 0;
	$end ||= 0;
	my($start,$stop) = $self->range($seqType);
	if($beg == 0) { $beg = $start; $end = $beg+$end; }
	elsif($end == 0) { $end = $stop; $beg = $end-$beg; }

	if($end >= $stop) { $end = $stop; } ##ML changed from if (end >stop)
	else { $end += 1;}   ##ML moved from commented position below, makes
                             ##more sense here
#	if($end > $stop) { $end = $stop; }
	if($beg < $start) { $beg = $start; }
#	else { $end += 1;}

#	my $seq = substr($self->seq_str('match'), $beg-$start, ($end-$beg));

	## ML: START fix for substr out of range error ------------------
	my $seq = "";
	if (($self->{'_prog'} eq 'TBLASTN') and ($seqType eq 'sbjct'))
	{
	    $seq = substr($self->seq_str('match'),
			  int(($beg-$start)/3), int(($end-$beg+1)/3));

	} elsif (($self->{'_prog'} eq 'BLASTX') and ($seqType eq 'query'))
	{
	    $seq = substr($self->seq_str('match'),
			  int(($beg-$start)/3), int(($end-$beg+1)/3));
	} else {
	    $seq = substr($self->seq_str('match'), 
			  $beg-$start, ($end-$beg));
	}
	## ML: End of fix for  substr out of range error -----------------

	
	## ML: debugging code
	## This is where we get our exception.  Try printing out the values going
	## into this:
	##
#	 print STDERR 
#	     qq(*------------MY EXCEPTION --------------------\nSeq: ") , 
#	     $self->seq_str("$seqType"), qq("\n),$self->name,",(  index:";
#	 print STDERR  $beg-$start, ", len: ", $end-$beg," ), (HSPRealLen:", 
#	     CORE::length $self->seq_str("$seqType");
#	 print STDERR ", HSPCalcLen: ", $stop - $start +1 ," ), 
#	     ( beg: $beg, end: $end ), ( start: $start, stop: stop )\n";
	 ## ML: END DEBUGGING CODE----------

	if(!CORE::length $seq) {
	    $self->throw("Undefined sub-sequence ($beg,$end). Valid range = $start - $stop");
	}
	## Get data for a substring.
#	printf "Collecting HSP subsection data: beg,end = %d,%d; start,stop = %d,%d\n%s<---\n", $beg, $end, $start, $stop, $seq;
#	printf "Original match seq:\n%s\n",$self->seq_str('match');
	$seq =~ s/ //g;  # remove space (no info).
	my $len_cons = CORE::length $seq;
	$seq =~ s/\+//g;  # remove '+' characters (conservative substitutions)
	my $len_id = CORE::length $seq;
	push @data, ($len_id, $len_cons);
#	printf "  HSP = %s\n  id = %d; cons = %d\n", $self->name, $len_id, $len_cons; <STDIN>;
    }
    @data;
}



=head2 frac_identical

 Usage     : $hsp_object->frac_identical( [seq_type] );
 Purpose   : Get the fraction of identical positions within the given HSP.
 Example   : $frac_iden = $hsp_object->frac_identical('query');
 Returns   : Float (2-decimal precision, e.g., 0.75).
 Argument  : seq_type: 'query' | 'sbjct' | 'total'
           : default = 'total' (but see comments below).
 Throws    : n/a
 Comments  : Different versions of Blast report different values for the total
           : length of the alignment. This is the number reported in the
           : denominators in the stats section:
           : "Identical = 34/120 Positives = 67/120".
           : BLAST-GP uses the total length of the alignment (with gaps)
           : WU-BLAST uses the length of the query sequence (without gaps).
           : Therefore, when called without an argument or an argument of 'total',
           : this method will report different values depending on the
           : version of BLAST used.
           : 
           : To get the fraction identical among only the aligned residues,
           : ignoring the gaps, call this method with an argument of 'query'
           : or 'sbjct'.

See Also   : L<frac_conserved()|frac_conserved>, L<num_identical()|num_identical>, L<matches()|matches>

=cut

#-------------------
sub frac_identical {
#-------------------
# The value is calculated as opposed to storing it from the parsed results.
# This saves storage and also permits flexibility in determining for which
# sequence (query or sbjct) the figure is to be calculated.

    my( $self, $seqType ) = @_;
    $seqType ||= 'total';

    if($seqType ne 'total') {
      $self->_set_seq_data() unless $self->{'_set_seq_data'};
    }
    ## Sensitive to member name format.
    $seqType = "_\L$seqType\E";

    sprintf( "%.2f", $self->{'_numIdentical'}/$self->{$seqType.'Length'});
}


=head2 frac_conserved

 Usage     : $hsp_object->frac_conserved( [seq_type] );
 Purpose   : Get the fraction of conserved positions within the given HSP.
           : (Note: 'conservative' positions are called 'positives' in the
	   : Blast report.)
 Example   : $frac_cons = $hsp_object->frac_conserved('query');
 Returns   : Float (2-decimal precision, e.g., 0.75).
 Argument  : seq_type: 'query' | 'sbjct' 
           : default = 'total' (but see comments below).
 Throws    : n/a
 Comments  : Different versions of Blast report different values for the total
           : length of the alignment. This is the number reported in the
           : denominators in the stats section:
           : "Identical = 34/120 Positives = 67/120".
           : BLAST-GP uses the total length of the alignment (with gaps)
           : WU-BLAST uses the length of the query sequence (without gaps).
           : Therefore, when called without an argument or an argument of 'total',
           : this method will report different values depending on the
           : version of BLAST used.
           :
           : To get the fraction conserved among only the aligned residues,
           : ignoring the gaps, call this method with an argument of 'query'
           : or 'sbjct'.

See Also   : L<frac_conserved()|frac_conserved>, L<num_conserved()|num_conserved>, L<matches()|matches>

=cut

#--------------------
sub frac_conserved {
#--------------------
# The value is calculated as opposed to storing it from the parsed results.
# This saves storage and also permits flexibility in determining for which
# sequence (query or sbjct) the figure is to be calculated.
 
    my( $self, $seqType ) = @_;
    $seqType ||= 'total';

    if($seqType ne 'total') {
      $self->_set_seq_data() unless $self->{'_set_seq_data'};
    }

    ## Sensitive to member name format.
    $seqType = "_\L$seqType\E";

    sprintf( "%.2f", $self->{'_numConserved'}/$self->{$seqType.'Length'});
}


=head2 num_identical

 Usage     : $hsp_object->num_identical();
 Purpose   : Get the number of identical positions within the given HSP.
 Example   : $num_iden = $hsp_object->num_identical();
 Returns   : integer
 Argument  : n/a
 Throws    : n/a

See Also   : L<num_conserved()|num_conserved>, L<frac_identical()|frac_identical>

=cut

#-------------------
sub num_identical {
#-------------------
    my( $self) = shift;

    $self->{'_numIdentical'};
}


=head2 num_conserved

 Usage     : $hsp_object->num_conserved();
 Purpose   : Get the number of conserved positions within the given HSP.
 Example   : $num_iden = $hsp_object->num_conserved();
 Returns   : integer
 Argument  : n/a
 Throws    : n/a

See Also   : L<num_identical()|num_identical>, L<frac_conserved()|frac_conserved>

=cut

#-------------------
sub num_conserved {
#-------------------
    my( $self) = shift;

    $self->{'_numConserved'};
}



=head2 range

 Usage     : $hsp->range( [seq_type] );
 Purpose   : Gets the (start, end) coordinates for the query or sbjct sequence
           : in the HSP alignment.
 Example   : ($qbeg, $qend) = $hsp->range('query');
           : ($sbeg, $send) = $hsp->range('sbjct');
 Returns   : Two-element array of integers 
 Argument  : seq_type = string, 'query' or 'sbjct'  (default = 'query')
           : (case insensitive).
 Throws    : n/a

See Also   : L<start()|start>, L<end()|end>

=cut

#----------
sub range {
#----------
    my ($self, $seqType) = @_;

    $self->_set_seq_data() unless $self->{'_set_seq_data'};

    $seqType ||= 'query';
    ## Sensitive to member name changes.
    $seqType = "_\L$seqType\E";

    return ($self->{$seqType.'Start'},$self->{$seqType.'Stop'});
}

=head2 start

 Usage     : $hsp->start( [seq_type] );
 Purpose   : Gets the start coordinate for the query, sbjct, or both sequences
           : in the HSP alignment.
 Example   : $qbeg = $hsp->start('query');
           : $sbeg = $hsp->start('sbjct');
           : ($qbeg, $sbeg) = $hsp->start();
 Returns   : scalar context: integer
           : array context without args: list of two integers
 Argument  : In scalar context: seq_type = 'query' or 'sbjct'
           :  (case insensitive). If not supplied, 'query' is used.
           : Array context can be "induced" by providing an argument of 'list' or 'array'.
 Throws    : n/a

See Also   : L<end()|end>, L<range()|range>

=cut

#----------
sub start {
#----------
    my ($self, $seqType) = @_;

    $seqType ||= (wantarray ? 'list' : 'query');

    $self->_set_seq_data() unless $self->{'_set_seq_data'};

    if($seqType =~ /list|array/i) {
	return ($self->{'_queryStart'}, $self->{'_sbjctStart'});
    } else {
	## Sensitive to member name changes.
	$seqType = "_\L$seqType\E";
	return $self->{$seqType.'Start'};
    }
}

=head2 end

 Usage     : $hsp->end( [seq_type] );
 Purpose   : Gets the end coordinate for the query, sbjct, or both sequences
           : in the HSP alignment.
 Example   : $qbeg = $hsp->end('query');
           : $sbeg = $hsp->end('sbjct');
           : ($qbeg, $sbeg) = $hsp->end();
 Returns   : scalar context: integer
           : array context without args: list of two integers
 Argument  : In scalar context: seq_type = 'query' or 'sbjct'
           :  (case insensitive). If not supplied, 'query' is used.
           : Array context can be "induced" by providing an argument of 'list' or 'array'.
 Throws    : n/a

See Also   : L<start()|start>, L<range()|range>

=cut

#----------
sub end {
#----------
    my ($self, $seqType) = @_;

    $seqType ||= (wantarray ? 'list' : 'query');

    $self->_set_seq_data() unless $self->{'_set_seq_data'};

    if($seqType =~ /list|array/i) {
	return ($self->{'_queryStop'}, $self->{'_sbjctStop'});
    } else {
	## Sensitive to member name changes.
	$seqType = "_\L$seqType\E";
	return $self->{$seqType.'Stop'};
    }
}



=head2 strand

 Usage     : $hsp_object->strand( [seq_type] )
 Purpose   : Get the strand of the query or sbjct sequence.
 Example   : print $hsp->strand('query');
           : ($qstrand, $sstrand) = $hsp->strand();
 Returns   : -1, 0, or 1
           : -1 = Minus strand, +1 = Plus strand
           : Returns 0 if strand is not defined, which occurs
           : for non-TBLASTN/X reports.
           : In scalar context without arguments, returns queryStrand value.
           : In array context without arguments, returns a two-element list 
           :    of strings (queryStrand, sbjctStrand).
           : Array context can be "induced" by providing an argument of 'list' or 'array'.
 Argument  : seq_type: 'query' | 'sbjct' or undef
 Throws    : n/a

See Also   : L<_set_seq()|_set_seq>, L<_set_match_stats()|_set_match_stats>

=cut

#-----------
sub strand {
#-----------
    my( $self, $seqType ) = @_;
    $seqType  ||= (wantarray ? 'list' : 'query');

    return '' if $seqType eq 'query' and $self->{'_prog'} eq 'TBLASTN';

    ## Sensitive to member name format.
    $seqType = "_\L$seqType\E";

    # $seqType could be '_list'.
    $self->{'_queryStrand'} or $self->_set_seq_data() unless $self->{'_set_seq_data'};

    if($seqType  =~ /list|array/i) {
	return ('','') unless defined $self->{'_queryStrand'};
	return ($self->{'_queryStrand'}, $self->{'_sbjctStrand'});
    }
    local $^W = 0;
    $STRAND_SYMBOL{$self->{$seqType.'Strand'}} || 0;
}


#####################################################################################
##                            INSTANCE METHODS                                     ##
#####################################################################################


=head2 seq

 Usage     : $hsp->seq( [seq_type] );
 Purpose   : Get the query or sbjct sequence as a Bio::Seq.pm object.
 Example   : $seqObj = $hsp->seq('query');
 Returns   : Object reference for a Bio::Seq.pm object.
 Argument  : seq_type = 'query' or 'sbjct' (default = 'query').
 Throws    : Propagates any exception that occurs during construction
           : of the Bio::Seq.pm object.
 Comments  : The sequence is returned in an array of strings corresponding
           : to the strings in the original format of the Blast alignment.
           : (i.e., same spacing).

See Also   : L<seq_str()|seq_str>, L<seq_inds()|seq_inds>, B<Bio::Seq.pm>

=cut

#-------
sub seq {  
#-------
    my($self,$seqType) = @_; 
    $seqType ||= 'query';
    my $str = $self->seq_str($seqType);
    my $num =  $self->name;
    my $name = $seqType =~ /query/i 
	            ? $self->parent->parent->name 
	            : $self->parent->name; 
	
    require Bio::Seq;

    new Bio::Seq (-ID   => $name,
		  -SEQ  => $str,
		  -DESC => "Blast HSP #$num, $seqType sequence",
		  );
}



=head2 seq_str

 Usage     : $hsp->seq_str( seq_type );
 Purpose   : Get the full query, sbjct, or 'match' sequence as a string.
           : The 'match' sequence is the string of symbols in between the 
           : query and sbjct sequences.
 Example   : $str = $hsp->seq_str('query');
 Returns   : String
 Argument  : seq_Type = 'query' or 'sbjct' or 'match'
 Throws    : Exception if the argument does not match an accepted seq_type.
 Comments  : Calls _set_residues() to set the 'match' sequence if it has
           : not been set already.

See Also   : L<seq()|seq>, L<seq_inds()|seq_inds>, L<_set_match_seq()|_set_match_seq>

=cut

#------------
sub seq_str {  
#------------
    my($self,$seqType) = @_; 

    ## Sensitive to member name changes.
    $seqType = "_\L$seqType\E";

    $self->_set_seq_data() unless $self->{'_set_seq_data'};

    if($seqType =~ /sbjct|query/) {
	my $seq = join('',@{$self->{$seqType.'Seq'}}); 
	$seq =~ s/\s+//g;
	return $seq;

    } elsif( $seqType =~ /match/i) {
	# Only need to call _set_match_seq() if the match seq is requested.
	my $aref = $self->_set_match_seq() unless ref $self->{'_matchSeq'};
	$aref =  $self->{'_matchSeq'};

## DEBUGGING CODE:
#	if($self->parent->name eq '1AK5_' and $self->parent->parent->name eq 'YAR073W') {
#	    print "seq_str():\n   @$aref";<STDIN>;
#	}

	return join('',@$aref); 

    } else {
	$self->throw("Invalid or undefined sequence type: $seqType",
		     "Valid types: query, sbjct, match");
    }
}




=head2 seq_inds

 Usage     : $hsp->seq_inds( seq_type, class, collapse );
 Purpose   : Get a list of residue positions (indices) for all identical 
           : or conserved residues in the query or sbjct sequence.
 Example   : @ind = $hsp->seq_inds('query', 'identical');
           : @ind = $hsp->seq_inds('sbjct', 'conserved');
           : @ind = $hsp->seq_inds('sbjct', 'conserved', 1);
 Returns   : List of integers 
           : May include ranges if collapse is true.
 Argument  : seq_type  = 'query' or 'sbjct'  (default = query)
           : class     = 'identical' or 'conserved' (default = identical)
           :              (can be shortened to 'id' or 'cons')
           :              (actually, anything not 'id' will evaluate to 'conserved').
           : collapse  = boolean, if true, consecutive positions are merged
           :             using a range notation, e.g., "1 2 3 4 5 7 9 10 11" 
           :             collapses to "1-5 7 9-11". This is useful for 
           :             consolidating long lists. Default = no collapse.
 Throws    : n/a.
 Comments  : Calls _set_residues() to set the 'match' sequence if it has
           : not been set already.

See Also   : L<seq()|seq>, L<_set_residues()|_set_residues>, L<collapse_nums()|collapse_nums>, B<Bio::Tools::Blast::Sbjct::seq_inds()>

=cut

#---------------
sub seq_inds {
#---------------
    my ($self, $seq, $class, $collapse) = @_;

    $seq  ||= 'query';
    $class ||= 'identical';
    $collapse ||= 0;

    $self->_set_residues() unless defined $self->{'_identicalRes_query'};

    $seq  = ($seq !~ /^q/i ? 'sbjct' : 'query');
    $class = ($class !~ /^id/i ? 'conserved' : 'identical');

    ## Sensitive to member name changes.
    $seq  = "_\L$seq\E";
    $class = "_\L$class\E";

    my @ary = sort { $a <=> $b } keys %{ $self->{"${class}Res$seq"}};

    return $collapse ? &collapse_nums(@ary) : @ary;
}




=head2 get_aln

 Usage     : $hsp->get_aln()
 Purpose   : Get a Bio::UnivAln.pm object constructed from the query + sbjct 
           : sequences of the present HSP object.
 Example   : $aln_obj = $hsp->get_aln();
 Returns   : Object reference for a Bio::UnivAln.pm object.
 Argument  : n/a.
 Throws    : Propagates any exception ocurring during the construction of
           : the Bio::UnivAln object.
 Comments  : Requires Bio::UnivAln.pm.
           : The Bio::UnivAln.pm object is constructed from the query + sbjct 
           : sequence objects obtained by calling seq().
           : Gap residues are included (see $GAP_SYMBOL). It is important that 
           : Bio::UnivAln.pm recognizes the gaps correctly. A strategy for doing 
           : this is being considered. Currently it is hard-wired.

See Also   : L<seq()|seq>, B<Bio::UnivAln.pm>

=cut

#------------
sub get_aln {
#------------
    my $self = shift;
    
    require Bio::UnivAln;
    
    my $qseq = $self->seq('query');
    my $sseq = $self->seq('sbjct');

    my $desc = sprintf "HSP #%s of query %s vs. sbjct %s", 
	      $self->name, $self->parent->parent->name, $self->parent->name;

    my $type = $self->{'_prog'} =~ /P$|^T/ ? 'amino' : 'dna';

    Bio::UnivAln->new( -seqs => [$qseq, $sseq],
		       -desc => $desc,
		       -type => $type,
		       );
}


=head2 display

 Usage     : $sbjct_object->display( %named_parameters );
 Purpose   : Display information about Bio::Tools::Blast::Sbjct.pm data members
           : including:  length, gaps, score, significance value, 
           :             sequences and sequence indices.
 Example   : $object->display(-SHOW=>'stats');
 Argument  : Named parameters: (TAGS CAN BE UPPER OR LOWER CASE)
           :          -SHOW  => 'hsp',
           :          -WHERE => filehandle (default = STDOUT)
 Returns   : n/a
 Status    : Experimental
 Comments  : For more control over the display of sequence data,
           : use seq(), seq_str(), seq_inds().

See Also   : L<_display_seq()|_display_seq>, L<seq()|seq>, L<seq_str()|seq_str>, L<seq_inds()|seq_inds>, L<_display_matches()|_display_matches>, B<Bio::Root::Object::display()>

=cut

#-----------
sub display {
#-----------
    my( $self, %param ) = @_;
    
     my $sbjctName = $self->parent->name();
    my $queryName = $self->parent->parent->name();
    my $layout = $self->parent->parent->_layout();
    
    my $OUT = $self->set_display(%param);
    
    printf( $OUT "%-15s: %d\n", "LENGTH TOTAL", $self->length('total') );
    printf( $OUT "%-15s: %d\n", "LENGTH QUERY", $self->length('query') );
    printf( $OUT "%-15s: %d\n", "LENGTH SBJCT", $self->length('sbjct') );
    printf( $OUT "%-15s: %d\n", "GAPS QUERY", $self->gaps('query') );
    printf( $OUT "%-15s: %d\n", "GAPS SBJCT", $self->gaps('sbjct') );
    printf( $OUT "%-15s: %d\n", "SCORE", $self->{'_score'} );
    printf( $OUT "%-15s: %0.1f\n", "BITS", $self->{'_bits'} );
    if($layout == 1) {
	printf( $OUT "%-15s: %.1e\n", "P-VAL", $self->{'_p'} );
	printf( $OUT "%-15s: %.1e\n", "EXPECT", $self->{'_expect'} );
    } else {
	printf( $OUT "%-15s: %.1e\n", "EXPECT", $self->{'_expect'} );
    }

    my $queryLength = $self->length('query'); 

    printf( $OUT "%-15s: %d (%0.0f%%)\n", "IDENTICAL", $self->{'_numIdentical'}, 
	   $self->{'_numIdentical'}/$queryLength * 100 );
    printf( $OUT "%-15s: %d (%0.0f%%)  %s \n", "CONSERVED", $self->{'_numConserved'}, 
	   $self->{'_numConserved'}/$queryLength * 100,
	   "includes identical" );
    
    $self->_display_seq('query', $queryName, $OUT);
    $self->_display_seq('sbjct', $sbjctName, $OUT);
    $self->_display_matches($queryName, $sbjctName, $OUT);
}




=head2 _display_seq

 Usage     : n/a; called automatically by display()
 Purpose   : Display information about query and sbjct HSP sequences.
           : Prints the start, stop coordinates and the actual sequence.
 Example   : n/a
 Argument  : 
 Returns   : printf call.
 Status    : Experimental
 Comments  : For more control, use seq(), seq_str(), or seq_inds().

See Also   : L<display()|display>, L<seq()|seq>, L<seq_str()|seq_str>, L<seq_inds()|seq_inds>, L<_display_matches()|_display_matches>

=cut

#------------------
sub _display_seq {
#------------------
    my( $self, $seqType, $name, $OUT ) = @_;
    
    $self->_set_seq_data() unless $self->{'_set_seq_data'};

    # Sensitive to member name changes.
    my $mem = "_\L$seqType\E";
    printf( $OUT "\n%10s: %s\n%10s  %s\n", "\U$seqType\E", "$name", "-----",
	   ('-'x ((CORE::length $name) + 2)) );
    printf( $OUT "%13s: %d\n", "START", $self->{$mem.'Start'} );
    printf( $OUT "%13s: %d\n", "STOP", $self->{$mem.'Stop'} );
    printf( $OUT "%13s: \n", "SEQ" );
    foreach( @{ $self->{$mem.'Seq'}} ) {
	printf( $OUT "%15s%s\n", "", $_);
    }
}


=head2 _display_matches

 Usage     : n/a; called automatically by display()
 Purpose   : Display information about identical and conserved positions
           : within both the query and sbjct sequences.
 Example   : n/a
 Argument  :
 Returns   : printf call.
 Status    : Experimental
 Comments  : For more control, use seq_inds().

See Also   : L<display()|display>, L<seq_inds()|seq_inds>, L<_display_seq()|_display_seq>, 

=cut

#--------------------
sub _display_matches {
#--------------------
    my( $self, $queryName, $sbjctName, $OUT) = @_;
    my($resNum, $count);
    
    $self->_set_residues() unless defined $self->{'_identicalRes_query'};

    printf( $OUT "\n%10s:    \n%10s\n", "HITS", "-----" );
    foreach( @{ $self->{'_matchSeq'}} ) {
	printf( $OUT "%15s%s\n", "", $_ );
    }		
    
    print $OUT "\n\U$queryName\E\n------------\n";
    printf( $OUT "\n%5s%s:\n%5s%s\n\t", "", "IDENTICAL RESIDUES IN $queryName (n=$self->{'_numIdentical'})", 
	   "", "--------------------------------------------" );
    $count = 0;
    foreach $resNum ( sort keys %{ $self->{'_identicalRes_query' }} ) {
	$count++;
	print $OUT "$resNum";
	$count > 0 and print $OUT +( $count % 15 ? ", " : "\n\t");
    }
    
    print $OUT "\n";

    my $justConserved = ($self->{'_numConserved'})-($self->{'_numIdentical'});
    printf( $OUT "\n%5s%s:\n%5s%s\n\t", "","CONSERVED RESIDUES IN $queryName (n=$justConserved)", 
	   "", "--------------------------------------------" );
    $count = 0;
    foreach $resNum ( sort keys %{ $self->{'_conservedRes_query' }} ) { 
	$count++;
	print $OUT "$resNum";
	$count > 0 and print $OUT +( $count % 15 ? ", " : "\n\t");
    }
    
    
    print $OUT "\n\n\U$sbjctName\E\n------------\n";
    printf( $OUT "\n%5s%s:\n%5s%s\n\t", "", "IDENTICAL RESIDUES IN $sbjctName (n=$self->{'_numIdentical'})", 
	   "", "--------------------------------------------" );
    $count = 0;
    foreach $resNum ( sort keys %{ $self->{'_identicalRes_sbjct' }} ) {
	$count++;
	print $OUT "$resNum";
	$count > 0 and print $OUT +( $count % 15 ? ", " : "\n\t");
    }
    
    print $OUT "\n";
    $justConserved = ($self->{'_numConserved'})-($self->{'_numIdentical'});
    printf( $OUT "\n%5s%s:\n%5s%s\n\t", "","CONSERVED RESIDUES IN $sbjctName (n=$justConserved)", 
	   "", "--------------------------------------------" );
    $count       = 0;
    foreach $resNum ( sort keys %{ $self->{'_conservedRes_sbjct' }} ) { 
	$count++;
	print $OUT "$resNum";
	$count > 0 and print $OUT +( $count % 15 ? ", " : "\n\t");
    }
}




=head2 homol_data

 Usage     : $data = $hsp_object->homo_data( %named_params );
 Purpose   : Gets similarity data for a single HSP.
 Returns   : String:
           : "Homology data" for each HSP is in the format:
           :  "<integer> <start> <stop>"
           :  where integer is the value returned by homol_score().
 Argument  : Named params:  (UPPER OR LOWERCASE TAGS)
           :   currently just one param is used:
           :   -SEQ =>'query' or 'sbjct'
 Throws    : n/a
 Status    : Experimental
 Comments  : This is a very experimental method used for obtaining a
           : coarse indication of:
           :   1) how strong the similarity is between the sequences in the HSP,
           :   3) the endpoints of the alignment (sequence monomer numbers)

See Also   : L<homol_score()|homol_score>, B<Bio::Tools::Blast.::homol_data()>, B<Bio::Tools::Blast::Sbjct::homol_data()>

=cut

#---------------
sub homol_data {
#---------------
    my ($self, %param) = @_;
    my $seq = $param{-SEQ} || $param{'-seq'} || 'sbjct';  # 'query' or 'sbjct'
    my $homolScore = $self->homol_score();
    # Sensitive to member name changes.
    $seq = "_\L$seq\E";
    
    $self->_set_seq_data() unless $self->{'_set_seq_data'};
    return ( $homolScore.' '.$self->{$seq.'Start'}.' '.$self->{$seq.'Stop'});
}


=head2 homol_score

 Usage     : $self->homol_score();
 Purpose   : Get a homology score (integer 1 - 3) as a coarse representation of 
           : the strength of the similarity independent of sequence composition. 
           : Based on the Blast bit score.
 Example   : $hscore = $hsp->homol_score();
 Returns   : Integer 
 Argument  : n/a
 Throws    : n/a
 Status    : Experimental
 Comments  : See @Bio::Tools::Blast::HSP::SCORE_CUTOFFS for the specific values. 
           : Currently,  BIT_SCORE  HOMOL_SCORE
           :             ---------  -----------
           :             >=100     --> 3 
           :               30-100  --> 2
           :             < 30      --> 1

See Also   : L<homol_data()|homol_data>

=cut

#----------------
sub homol_score {
#----------------
    my $self = shift;
    
    if( $self->{'_bits'} >= $SCORE_CUTOFFS[0] ) { 1 }
    elsif($self->{'_bits'} < $SCORE_CUTOFFS[0] and 
	  $self->{'_bits'} >= $SCORE_CUTOFFS[1] ) { 2 }
    else  { 3 }
}


#####################################################################################
##                                  CLASS METHODS                                  ##
#####################################################################################

=head1 CLASS METHODS

=head2 collapse_nums

 Usage     : @cnums = collapse_nums( @numbers );
 Purpose   : Collapses a list of numbers into a set of ranges of consecutive terms:
           : Useful for condensing long lists of consecutive numbers.
           :  EXPANDED:
           :     1 2 3 4 5 6 10 12 13 14 15 17 18 20 21 22 24 26 30 31 32
           :  COLLAPSED:
           :     1-6 10 12-15 17 18 20-22 24 26 30-32
 Argument  : List of numbers and sorted numerically.
 Returns   : List of numbers mixed with ranges of numbers (see above).
 Throws    : n/a
 Comments  : Probably belongs in a more general utility class.

See Also   : L<seq_inds()|seq_inds>

=cut

#------------------
sub collapse_nums {
#------------------
# This is not the slickest connectivity algorithm, but will do for now.
    my @a = @_;
    my ($from, $to, $i, @ca, $consec);
    
    $consec = 0;
    for($i=0; $i < @a; $i++) {
	not $from and do{ $from = $a[$i]; next; };
	if($a[$i] == $a[$i-1]+1) {
	    $to = $a[$i];
	    $consec++;
	} else {
	    if($consec == 1) { $from .= ",$to"; }
	    else { $from .= $consec>1 ? "\-$to" : ""; }
	    push @ca, split(',', $from);
	    $from =  $a[$i];
	    $consec = 0;
	    $to = undef;
	}
    }
    if(defined $to) {
	if($consec == 1) { $from .= ",$to"; }
	else { $from .= $consec>1 ? "\-$to" : ""; }
    }
    push @ca, split(',', $from) if $from;

    @ca;
}


1;
__END__

#####################################################################################
#                                  END OF CLASS
#####################################################################################

=head1 FOR DEVELOPERS ONLY

=head2 Data Members

Information about the various data members of this module is provided for those 
wishing to modify or understand the code. Two things to bear in mind: 

=over 4

=item 1 Do NOT rely on these in any code outside of this module. 

All data members are prefixed with an underscore to signify that they are private.
Always use accessor methods. If the accessor doesn't exist or is inadequate, 
create or modify an accessor (and let me know, too!). 

=item 2 This documentation may be incomplete and out of date.

It is easy for these data member descriptions to become obsolete as 
this module is still evolving. Always double check this info and search 
for members not described here.

=back

An instance of Bio::Tools::Blast::HSP.pm is a blessed reference to a hash containing
all or some of the following fields:

 FIELD           VALUE
 --------------------------------------------------------------
 (member names are mostly self-explanatory)

 _score              :
 _bits               :
 _p                  :
 _n                  : Integer. The 'N' value listed in parenthesis with P/Expect value:
                     : e.g., P(3) = 1.2e-30  ---> (N = 3).
                     : Not defined in NCBI Blast2 with gaps.
                     : To obtain the number of HSPs, use Bio::Tools::Blast::Sbjct::num_hsps().
 _expect             :
 _queryLength        : 
 _queryGaps          : 
 _queryStart         :
 _queryStop          :
 _querySeq           :
 _sbjctLength        :
 _sbjctGaps          :
 _sbjctStart         :
 _sbjctStop          :
 _sbjctSeq           :
 _matchSeq           : String. Contains the symbols between the query and sbjct lines
                       which indicate identical (letter) and conserved ('+') matches
                       or a mismatch (' ').
 _numIdentical       :
 _numConserved       :
 _identicalRes_query :
 _identicalRes_sbjct :
 _conservedRes_query :
 _conservedRes_sbjct :
 _match_indent       : The number of leading space characters on each line containing
                       the match symbols. _match_indent is 13 in this example:
                         Query:   285 QNSAPWGLARISHRERLNLGSFNKYLYDDDAG
                                      Q +APWGLARIS       G+ + Y YD+ AG
                         ^^^^^^^^^^^^^ 

 INHERITED DATA MEMBERS 

 _name          : From Bio::Root::Object.pm. 
		:
 _parent        : From Bio::Root::Object.pm. This member contains a reference to the
		: Bio::Tools::Blast::Sbjct.pm object to which this hit belongs.


=cut

1;

