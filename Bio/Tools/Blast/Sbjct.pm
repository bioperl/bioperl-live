#------------------------------------------------------------------------------
# PACKAGE : Bio::Tools::Blast::Sbjct.pm
# AUTHOR  : Steve Chervitz (sac@bioperl.org) 
# CREATED : 7 October 1996
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
#------------------------------------------------------------------------------

package Bio::Tools::Blast::Sbjct;

use Bio::Root::Global  qw(:devel);
use Bio::Root::Object  ();

@ISA        = qw( Bio::Root::Object Exporter );

use strict;
use vars qw($ID %SUMMARY_OFFSET $Revision);
$ID = 'Bio::Tools::Blast::Sbjct';
$Revision = '$Id$';  #'

my $_prog       = '';
my $_signif_fmt = '';

## POD Documentation:

=head1 NAME

Bio::Tools::Blast::Sbjct - Bioperl BLAST "Hit" object

=head1 SYNOPSIS

=head2 Object Creation

The construction of HSP objects is handled by B<Bio::Tools::Blast>.
You should not need to use this package directly. See L<_initialize()|_initialize>
for a description of constructor parameters.

    require Bio::Tools::Blast::Sbjct;

    $hit = new Bio::Tools::Blast::Sbjct (-DATA    =>\@hitData, 
					 -PARENT  =>$self, 
					 -NAME    =>5,
					 -RANK    =>5,
					 -RANK_BY =>'order',
					 -MAKE    =>'query' (or 'sbjct'),
					 -OVERLAP =>2,
					 -PROGRAM =>'TBLASTN'
					 );

@hitData includes the summary line for the hit as element [0], plus 
all lines from the HSP alignment section of the BLAST report for
the present hit. 

=head1 INSTALLATION

This module is included with the central Bioperl distribution:

   http://bio.perl.org/Core/Latest
   ftp://bio.perl.org/pub/DIST

Follow the installation instructions included in the README file.

=head1 DESCRIPTION

The Bio::Tools::Blast::Sbjct.pm module encapsulates data and methods for 
parsing and manipulating "hits" from a BLAST report.
This module is a utility module used by the Bio::Tools::Blast.pm
and is not intended for separate use.

In Blast lingo, the "sbjct" sequences are all the sequences 
in a target database which were compared against a "query" sequence.
The terms "sbjct" and "hit" will be used interchangeably in this and related modules. 

This module supports BLAST versions 1.x and 2.x, gapped and ungapped.


=head2 HSP Tiling and Ambiguous Alignments

If a Blast hit has more than one HSP, the Bio::Tools::Blast::Sbjct.pm
object has the ability to merge overlapping HSPs into contiguous
blocks. This permits the Sbjct object to sum data across all HSPs
without counting data in the overlapping regions multiple times, which
would happen if data from each overlapping HSP are simply summed.  HSP
tiling is performed automatically when methods of the Sbjct object
that rely on tiled data are invoked. These include
L<frac_identical()|frac_identical>, L<frac_conserved()|frac_conserved>, L<gaps()|gaps>,
L<frac_aligned_query()|frac_aligned_query>, L<frac_aligned_hit()|frac_aligned_hit>,
L<num_unaligned_query()|num_unaligned_query>, L<num_unaligned_hit()|num_unaligned_hit>.

It also permits the assessment of an "ambiguous alignment" if the
query (or sbjct) sequences from different HSPs overlap. The existence
of an overlap could indicate a biologically interesting region in the
sequence, such as a repeated domain.  The Sbjct object uses the
-OVERLAP parameter to determine when two sequences overlap; if this is
set to 2 -- the default -- then any two sbjct or query HSP sequences
must overlap by more than two residues to get merged into the same
contig and counted as an overlap. See the L<BUGS | BUGS> section below for
"issues" with HSP tiling.


The results of the HSP tiling is reported with the following ambiguity codes:

   'q' = Query sequence contains multiple sub-sequences matching
         a single region in the sbjct sequence. 

   's' = Sbjct sequence contains multiple sub-sequences matching
         a single region in the query sequence. 

   'qs' = Both query and sbjct sequences contain more than one
          sub-sequence with similarity to the other sequence.


For addition information about ambiguous BLAST alignments, see
L<_tile_hsps()|_tile_hsps> and 

 http://www-genome.stanford.edu/Sacch3D/help/ambig_aln.html

=head1 DEPENDENCIES

Bio::Tools::Blast::Sbjct.pm is a concrete class that inherits from B<Bio::Root::Object>
and relies on two other modules:

=over 4

=item B<Bio::Tools::Blast::HSP> 

Encapsulates a single high-scoring segment pair within a hit.

=item B<Bio::Tools::Blast>

Provides a container for Sbjct.pm objects.

=back


Bio::Tools::Blast::Sbjct.pm does not currently inherit from
Bio::Root::Vector.pm since Bio::Root::Vector.pm may be re-designed to
make it usable via delegation.  Thus, a Blast.pm object would manage a
vector of Sbjct.pm objects.  Stay tuned.


=head1 BUGS

One consequence of the HSP tiling is that methods that rely on HSP
tiling such as L<frac_identical()|frac_identical>, L<frac_conserved()|frac_conserved>, L<gaps()|gaps>
etc. may report misleading numbers when C<-OVERLAP> is set to a large
number.  For example, say we have two HSPs and the query sequence tile
as follows:

            1      8             22      30        40             60 
 Full seq:  ------------------------------------------------------------
                    *  ** *   **
 HSP1:             ---------------                    (6 identical matches)
                              **   **  **
 HSP2:                        -------------           (6 identical matches)


If C<-OVERLAP> is set to some number over 4, HSP1 and HSP2 will not be
tiled into a single contig and their numbers of identical matches will
be added, giving a total of 12, not 10 if they had be combined into
one contig. This can lead to number greater than 1.0 for methods
L<frac_identical()|frac_identical> and L<frac_conserved()|frac_conserved>. This is less of an issue
with gapped Blast since it tends to combine HSPs that would be listed
separately without gapping.  (Fractions E<gt>1.0 can be viewed as a
signal for an interesting alignment that warrants further inspection,
thus turning this bug into a feature).

Using large values for C<-OVERLAP> can lead to incorrect numbers
reported by methods that rely on HSP tiling but can be useful if you
care more about detecting ambiguous alignments.  Setting C<-OVERLAP>
to zero will lead to the most accurate numbers for the
tiling-dependent methods but will be useless for detecting overlapping
HSPs since all HSPs will appear to overlap.


=head1 SEE ALSO

 Bio::Tools::Blast::HSP.pm     - Blast HSP object.
 Bio::Tools::Blast.pm          - Blast object.
 Bio::Root::Object.pm          - Proposed base class for all Bioperl objects.

Links:

 http://bio.perl.org/Core/POD/Tools/Blast/HSP.pm.html

 http://bio.perl.org/Projects/modules.html  - Online module documentation
 http://bio.perl.org/Projects/Blast/        - Bioperl Blast Project     
 http://bio.perl.org/                       - Bioperl Project Homepage


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

Steve Chervitz, sac@bioperl.org

See the L<FEEDBACK | FEEDBACK> section for where to send bug reports and comments.

=head1 COPYRIGHT

Copyright (c) 1996-2000 Steve Chervitz. All Rights Reserved.
This module is free software; you can redistribute it and/or 
modify it under the same terms as Perl itself.

=cut



#
##
###
#### END of main POD documentation.
###
##
#'


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
           : Bio::Tools::Blast::Sbjct.pm objects are constructed 
           : automatically by Bio::Tools::Blast.pm, so there is no need
           : for direct consumption.
 Purpose   : Initializes key varaiables and calls methods to parse a single Blast hit.
           : Constructs Bio::Tools::Blast::HSP.pm objects for each 
           : high-scoring segment pair (HSP).
           : Calls superclass constructor first (Bio::Root::Object.pm).
 Returns   : n/a
 Argument  : Named Parameters passed from new()
           : ALL TAGS MUST BE UPPERCASE (does not call _rearrange()).
           :     -DATA       => array reference holding all data for a single hit.
           :                    The first element should hold the description
           :                    line (from the desctiption section at the top of
           :                    the report), remaining lines should hold all lines
           :                    within the HSP alignment listing section of report.
	   :	 -PARENT     => object reference to a Bio::Tools::Blast.pm object.
	   :	 -NAME       => string (typically same as -RANK, just a temporary
           :                    name to use until the actual name of hit is parsed),
	   :	 -RANK       => integer,
	   :	 -RANK_BY    => 'order',
	   :	 -OVERLAP    => integer (maximum overlap between adjacent
           :                    HSPs when tiling)
	   :	 -PROGRAM    => string (type of Blast blastp, blastn, etc).

See Also   : L<_set_id()|_set_id>, L<_set_hsps()|_set_hsps>, L<_tile_hsps()|_tile_hsps>, B<Bio::Root::Object::new()>, B<Bio::Tools::Blast::_set_hit_db()>

=cut

#-------------------
sub _initialize {
#-------------------
    my( $self, %param ) = @_;
    
    # $make not currently used.
    my $make = $self->SUPER::_initialize( %param );
    
    # Set various class data.
    $_prog       = $param{-PROGRAM} || '';
    $_signif_fmt = $param{-SIGNIF_FMT};

    $self->{'_rank'} = $param{-RANK} || '';
    $self->_set_id( $param{-DATA}->[0]);
    $self->_set_hsps( @{$param{-DATA}} );

    $self->{'_overlap'} = $param{-OVERLAP} || 0;
}

#--------------
sub destroy {
#--------------
    my $self=shift; 
    if($self->{'_hsps'}) {
	foreach($self->hsps) { 
	    $_->destroy; 
	    undef $_; 
	}
	undef $self->{'_hsps'};
    }
    $DEBUG==2 && print STDERR "DESTROYING $self ${\$self->name}";
    $self->SUPER::destroy;
}

#####################################################################################
##                                  ACCESSORS                                      ##
#####################################################################################


=head2 rank

 Usage     : $sbjct->rank( integer or string );
 Purpose   : Sets/Gets the rank of the current Sbjct object relative to 
           : other Sbjct objects managed by a given Blast object.
 Example   : $sbjct->rank(1);
 Returns   : The current rank value.
 Argument  : Integer or string to be used for ranking the hit 
           : relative to other hits.
 Throws    : n/a
 Comments  : The rank usually corresponds to the order the listing
           : of hits in the BLAST report from lowest to highest p-value.
           : Rank need not be restricted to this value.
           : rank() may be provided by a delegated or inherited
           : iterator class in the future (such as Bio::Root::Vector.pm).

=cut

#-----------
sub rank {
#-----------
    my $self = shift;
    if(@_) {$self->{'_rank'} = shift; }
    $self->{'_rank'};
}



=head2 _set_id

 Usage     : n/a; automatically called by _initialize()
 Purpose   : Sets the name of the Sbjct sequence from the BLAST summary line.
           : The identifier is assumed to be the first
           : chunk of non-whitespace characters in the description line
           : Does not assume any semantics in the structure of the identifier
           : (Formerly, this method attempted to extract database name from
           : the seq identifiers, but this was prone to break).
 Returns   : n/a
 Argument  : String containing description line of the hit from Blast report
           : or first line of an alignment section.
 Throws    : Warning if cannot locate sequence ID.

See Also   : L<_initialize()|_initialize>, B<Bio::Tools::Blast::database()>

=cut

#---------------
sub _set_id {
#---------------
    my( $self, $desc ) = @_;
    my ($seqID1, $seqID2, $dbID, @seqDat);

    local $_ = $desc;
    my @linedat = split();
    my $data = $linedat[0];

# New strategy: Assume only that the ID is the first white space
# delimited chunk. Not attempting to extract database name.
# Clients will have to interpret it as necessary.
    if($data =~ /^(\S+)\s*/) {
        $self->name($1);
    }
    else {
        $self->warn("Can't locate sequence identifier in summary line.", "Line = $data");
        $data = 'Unknown sequence ID' if not $data;
        $self->name($data);
    }
    $self->{'_db'}  = '-';

# Old strategy: assumes semantics in the identifier
# and tries to separate out database and id components.
# Too fancy and fragile!  SAC, 2000-02-18

#    # Proceeding from more standard (NCBI-like) to less standard.
#    if($data =~ /(\S+?)[:\|]+(\S+?)[:\|]+(\S*)/) {
#       # matches: database|id1|id2 or database||id1||id2 or database:id1:id2 
#        $dbID    = $1;
#        $seqID1  = $2;
#        $seqID2  = $3;
#        if($seqID2 eq $seqID1) { undef($seqID2); }
#
#    } elsif($data =~ /(\S+?)[:\|]+(\S+)/) {
#       # matches: database|id1  or database:id1
#        $dbID    = $1;
#        $seqID1  = $2;
#
#    } elsif($data =~ /^(\S+)\s+([gb|emb|dbj|sp|pir])\s+(\S+)*/) {
#       # matches: id1 database id2 
#        $seqID1  = $1;
#        $dbID    = $2;
#        $seqID2  = $3;
#
#    } elsif($data =~ /^(\S+)\s*/) {
#        $seqID1 = $1;
#    }
#
#    ## Combine the multiple IDs.
#    $seqID2 = scalar($seqID2) ? "/$seqID2" : '';
#
#    if( !scalar $seqID1) {
#        $self->warn("Can't locate sequence identifier in summary line.", "Line = $data");
#        $self->name('Unknown sequence ID');
#    } else {
#        $self->name($seqID1.$seqID2);
#    }
#    $self->{'_db'}  = $dbID || '-';
}


=head2 _set_hsps

 Usage     : n/a; called automatically during object construction.
 Purpose   : Creates HSP.pm objects for each HSP in a BLAST hit alignment.
           : Also collects the full description of the hit from the
           : HSP alignment section.
 Returns   : n/a
 Argument  : List of strings containing raw BLAST report data for 
           : a single hit's HSP alignment data.
 Throws    : Warnings for each HSP.pm object that fails to be constructed.
           : Exception if no HSP.pm objects can be constructed.
           : Exception if can't parse length data for hit sequence.
 Comments  : Requires Bio::Tools::Blast::HSP.pm.
           : Sets the description using the full string present in 
           : the alignment data.
           : Also sets Expect and P-values for the Sbjct object by
           : copying from the HSP object. 
           : While this sacrifices some memory efficiency, it
           : improves access speed for these critical data.

See Also   : L<_initialize()|_initialize>, L<_set_desc()|_set_desc>

=cut

#--------------
sub _set_hsps { 
#--------------
    
    my( $self, @data ) = @_;
    my $start     = 0;
    my $hspCount  = 0;

    require Bio::Tools::Blast::HSP;

#    printf STDERR "$ID _set_hsps(). DATA (%d lines) =\n@data\n", scalar(@data); <STDIN>;

    my( @hspData, @hspList, @errs, @bad_names );
    my($line, $set_desc, @desc);
    $set_desc = 0;

    hit_loop:
   foreach $line( @data ) {

       if( $line =~ /^\s*Length = ([\d,]+)/ ) {
	   $self->_set_desc(@desc);
	   $set_desc = 1;
	   ($self->{'_length'} = $1) =~ s/,//g; # get rid of commas
	   next hit_loop;
       } elsif( !$set_desc) {
	   $line =~ s/^\s+|\s+$//g;
	   push @desc, $line;
	   next hit_loop;
       } elsif( $line =~ /^\s*Score/ ) {
	   ## This block is for setting multiple HSPs.

	   if( not scalar @hspData ) {
	       $start = 1; 
	       push @hspData, $line; 
	       next hit_loop;

	    } elsif( scalar @hspData) {  
		$hspCount++;
		$DEBUG and do{ print STDERR +( $hspCount % 10 ? "+" : "+\n" ); };

#		print STDERR "\n$ID: setting HSP: ${\$self->name}\n";
		my $hspObj = eval { new Bio::Tools::Blast::HSP(-DATA    =>\@hspData, 
							       -PARENT  =>$self, 
							       -NAME    =>$hspCount,
							       -PROGRAM =>$_prog,
							       ); 
				};
		if($@) {
#		   print "$ID: ERROR:\n$@";<STDIN>;
		  push @errs, $@;
		  push @bad_names, "#$hspCount";
		    $hspObj->destroy if ref $hspObj;
   		    undef $hspObj;
		} else {
		    push @hspList, $hspObj;
                    if (!defined($self->{'_expect'}) || $hspObj->expect() < $self->{'_expect'}) {
                        $self->{'_expect'} = $hspObj->expect();
                    }
                    if (!defined($self->{'_p'}) || $hspObj->p() < $self->{'_p'}) {
                        $self->{'_p'} = $hspObj->p();
                    }
		}
		@hspData = ();
		push @hspData, $line;
		next;
	   } else {
	       push @hspData, $line;
	   }
       } elsif( $start ) {
	   ## This block is for setting the last HSP (which may be the first as well!).
	   if( $line =~ /^(end|>|Parameters|CPU|Database:)/ ) {
	       $hspCount++;
	       $DEBUG and do{ print STDERR +( $hspCount % 10 ? "+" : "+\n" ); };

#	       print STDERR "\n$ID: setting HSP: ${\$self->name}\n"; 

	       my $hspObj = eval { new Bio::Tools::Blast::HSP(-DATA    =>\@hspData, 
							      -PARENT  =>$self, 
							      -NAME    =>$hspCount,
							      -PROGRAM =>$_prog,
							      );
			       };
	       if($@) {
#		   print "$ID: ERROR:\n$@";<STDIN>;
		  push @errs, $@;
		  push @bad_names, "#$hspCount";
		  $hspObj->destroy if ref $hspObj;
   		  undef $hspObj;
	       } else {
		   push @hspList, $hspObj;
                   if (!defined($self->{'_expect'}) || $hspObj->expect() < $self->{'_expect'}) {
                       $self->{'_expect'} = $hspObj->expect();
                   }
                   if (!defined($self->{'_p'}) || $hspObj->p() < $self->{'_p'}) {
                       $self->{'_p'} = $hspObj->p();
                   }
	       }
	   } else {
	       push @hspData, $line;
	   }
       }
   }		

#    print STDERR "\n--------> Done building HSPs for ${\$self->name}\n";
    
    $self->{'_length'} or $self->throw( "Can't determine hit sequence length.");

    # Adjust logical length based on BLAST flavor.
    if($_prog =~ /TBLAST[NX]/) {
	$self->{'_logical_length'} = $self->{'_length'} / 3;
    }
    
    # Handling errors as done in Blast.pm. (as of version 0.073)

    if(@errs) {
	my ($str);
	# When there are many errors, in most of the cases, they are
	# caused by the same problem. Only need to see full data for
	# the first one.
	if(@errs > 2) {
	    $str = "SHOWING FIRST EXCEPTION ONLY:\n$errs[0]";
	    $self->clear_err();  # clearing the existing set of errors.
	                         # Not necessary, unless the -RECORD_ERR =>1
	                         # constructor option was used for Blast object.
	} else {
	    $str = join("\n",@errs);
	}

    if( not scalar @hspList) {
      $self->throw("Failed to create any HSP objects for $hspCount potential HSP(s).",
		   "\n\nTRAPPED EXCEPTION(S):\n$str\nEND TRAPPED EXCEPTION(S)\n"
			 );
      } else {
	    $self->warn(sprintf("Could not create HSP objects for %d HSP(s): %s", scalar(@errs), join(', ',@bad_names)), 
			"\n\nTRAPPED EXCEPTION(S):\n$str\nEND TRAPPED EXCEPTION(S)\n"
		       );
	  }

    } else {
	$self->{'_hsps'} = \@hspList;
    }
}

=head2 _set_desc

 Usage     : n/a; called automatically by _set_hsps()
 Purpose   : Sets the description of the hit sequence.
           : For sequence without descriptions, sets description to "-".
 Argument  : Array containing description (multiple lines).
 Comments  : _set_hsps() calls this method with the data from the 
           : HSP alignment listing, which contains the complete description.
           : (Formerly, this was called from the _set_desc_data() method initially.)

See Also   : _set_hsps()

=cut

#--------------
sub _set_desc {
#--------------
    my( $self, @desc ) = @_;
    my( $desc);
    
#    print "$ID: RAW DESC:\n@desc";<STDIN>;
    
    $desc = join(" ", @desc);
    
    if($desc) {
	$desc =~ s/^\s*\S+\s+//; # remove the sequence ID(s)
	$desc =~ s/^[\s!]+//;
	$desc =~ s/ \d+$//;
	$desc =~ s/\.+$//;
	$self->{'_desc'} = $desc || '-';
    } else {
	$self->{'_desc'} = '-';
    }

#    print "$ID: _set_desc =  $desc";<STDIN>;
}


=head2 _tile_hsps

 Usage     : n/a; called automatically during object construction or
           : as needed by methods that rely on having tiled data.
 Purpose   : Collect statistics about the aligned sequences in a set of HSPs.
           : Calculates the following data across all HSPs: 
           :    -- total alignment length 
           :    -- total identical residues 
           :    -- total conserved residues
 Returns   : n/a
 Argument  : n/a
 Throws    : n/a
 Status    : Experimental
 Comments  :
           : This method performs more careful summing of data across
           : all HSPs in the Sbjct object. Simply summing the data from all HSPs
           : will overestimate the actual length of the alignment if there is 
           : overlap between different HSPs (often the case).
           : The strategy is to tile the HSPs and sum over the
           : contigs, collecting data separately from overlapping and
           : non-overlapping regions of each HSP. To facilitate this, the
           : HSP.pm object now permits extraction of data from sub-sections
           : of an HSP.
           : 
           : Additional useful information is collected from the results
           : of the tiling. It is possible that sub-sequences in
           : different HSPs will overlap significantly. In this case, it
           : is impossible to create a single unambiguous alignment by
           : concatenating the HSPs. The ambiguity may indicate the
           : presence of multiple, similar domains in one or both of the
           : aligned sequences. This ambiguity is recorded using the
           : ambiguous_aln() method.
           : 
           : This method does not attempt to discern biologically
           : significant vs. insignificant overlaps. The allowable amount of 
           : overlap can be set with the overlap() method or with the -OVERLAP
           : parameter used when constructing the Blast & Sbjct objects. 
           : 
           : For a given hit, both the query and the sbjct sequences are
           : tiled independently.
           : 
           :    -- If only query sequence HSPs overlap, 
           :          this may suggest multiple domains in the sbjct.
           :    -- If only sbjct sequence HSPs overlap, 
           :          this may suggest multiple domains in the query.
           :    -- If both query & sbjct sequence HSPs overlap, 
           :          this suggests multiple domains in both.
           :    -- If neither query & sbjct sequence HSPs overlap, 
           :          this suggests either no multiple domains in either
           :          sequence OR that both sequences have the same
           :          distribution of multiple similar domains.
           : 
           : This method can deal with the special case of when multiple
           : HSPs exactly overlap.
           : 
           : Efficiency concerns:
           :  Speed will be an issue for sequences with numerous HSPs.
           : 
 Bugs      : Currently, _tile_hsps() does not properly account for
           : the number of non-tiled but overlapping HSPs, which becomes a problem
           : as overlap() grows. Large values overlap() may thus lead to 
           : incorrect statistics for some hits. For best results, keep overlap()
           : below 5 (DEFAULT IS 2). For more about this, see the "HSP Tiling and
           : Ambiguous Alignments" section.

See Also   : L<_adjust_contigs()|_adjust_contigs>, L<ambiguous_aln()|ambiguous_aln>, L<overlap()|overlap>, L<frac_identical()|frac_identical>, L<frac_conserved()|frac_conserved>, L<frac_aligned_query()|frac_aligned_query>, L<frac_aligned_hit()|frac_aligned_hit>, L<num_unaligned_query()|num_unaligned_query>, L<num_unaligned_hit()|num_unaligned_hit>, L<HSP Tiling and Ambiguous Alignments>

=cut

#--------------
sub _tile_hsps {
#--------------
    my $self = shift;
#    my $gapped = $self->parent->gapped || 0;   # no special treatment

    $self->{'_tile_hsps'} = 1;
    $self->{'_gaps_query'} = 0;
    $self->{'_gaps_sbjct'} = 0;

    ## Simple summation scheme. Valid if there is only one HSP.
    if((defined($self->{'_n'}) and $self->{'_n'} == 1) or $self->num_hsps == 1) {
	my $hsp = $self->hsp;
	$self->{'_length_aln_query'} = $hsp->length('query');
	$self->{'_length_aln_sbjct'} = $hsp->length('sbjct');
	$self->{'_length_aln_total'} = $hsp->length('total');
	($self->{'_totalIdentical'},$self->{'_totalConserved'}) = $hsp->matches();
	$self->{'_gaps_query'} = $hsp->gaps('query');
	$self->{'_gaps_sbjct'} = $hsp->gaps('sbjct');

#	print "_tile_hsps(): single HSP, easy stats.\n";
	return;
    } else {
#	print STDERR "$ID: _tile_hsps: summing multiple HSPs\n";
	$self->{'_length_aln_query'} = 0;
	$self->{'_length_aln_sbjct'} = 0;
	$self->{'_length_aln_total'} = 0;
	$self->{'_totalIdentical'}   = 0;
	$self->{'_totalConserved'}   = 0;
    }

    ## More than one HSP. Must tile HSPs.
#    printf "\nTiling HSPs for %s (BLAST: %s)\n",$self->name, $self->parent->name;
    my($hsp, $qstart, $qstop, $sstart, $sstop);
    my(@qcontigs, @scontigs);
    my $qoverlap = 0;
    my $soverlap = 0;
    my $max_overlap = $self->{'_overlap'};

    foreach $hsp ($self->hsps()) {
#	printf "  HSP: %s\n%s\n",$hsp->name, $hsp->str('query');
#	printf "  Length = %d; Identical = %d; Conserved = %d; Conserved(1-10): %d",$hsp->length, $hsp->length(-TYPE=>'iden'), $hsp->length(-TYPE=>'cons'), $hsp->length(-TYPE=>'cons',-START=>0,-STOP=>10); <STDIN>;
	($qstart, $qstop) = $hsp->range('query');
	($sstart, $sstop) = $hsp->range('sbjct');

	my ($qgaps, $sgaps)  = $hsp->gaps();
	$self->{'_gaps_query'} += $qgaps;
	$self->{'_gaps_sbjct'} += $sgaps;

	$self->{'_length_aln_total'} += $hsp->length;
	## Collect contigs in the query sequence.
	$qoverlap = &_adjust_contigs('query', $hsp, $qstart, $qstop, \@qcontigs, $max_overlap);

	## Collect contigs in the sbjct sequence (needed for domain data and gapped Blast).
	$soverlap = &_adjust_contigs('sbjct', $hsp, $sstart, $sstop, \@scontigs, $max_overlap);

	## Collect overall start and stop data for query and sbjct over all HSPs.
	if(not defined $self->{'_queryStart'}) {
	    $self->{'_queryStart'} = $qstart;
	    $self->{'_queryStop'}  = $qstop;
	    $self->{'_sbjctStart'} = $sstart;
	    $self->{'_sbjctStop'}  = $sstop;
	} else {
	    $self->{'_queryStart'} = ($qstart < $self->{'_queryStart'} ? $qstart : $self->{'_queryStart'});
	    $self->{'_queryStop'}  = ($qstop  > $self->{'_queryStop'}  ? $qstop  : $self->{'_queryStop'});
	    $self->{'_sbjctStart'} = ($sstart < $self->{'_sbjctStart'} ? $sstart : $self->{'_sbjctStart'});
	    $self->{'_sbjctStop'}  = ($sstop  > $self->{'_sbjctStop'}  ? $sstop  : $self->{'_sbjctStop'});
	}	    
    }

    ## Collect data across the collected contigs.

#    print "\nQUERY CONTIGS:\n";
#    print "  gaps = $self->{'_gaps_query'}\n";

    foreach(@qcontigs) {
#	print "  query contig: $_->{'start'} - $_->{'stop'}\n";
#	print "         iden = $_->{'iden'}; cons = $_->{'cons'}\n";
	$self->{'_length_aln_query'} += $_->{'stop'} - $_->{'start'} + 1;
	$self->{'_totalIdentical'}   += $_->{'iden'};
	$self->{'_totalConserved'}   += $_->{'cons'};
    }

    ## Collect data for sbjct contigs. Important for gapped Blast.
    ## The totalIdentical and totalConserved numbers will be the same
    ## as determined for the query contigs.

#    print "\nSBJCT CONTIGS:\n";
#    print "  gaps = $self->{'_gaps_sbjct'}\n";

    foreach(@scontigs) {
#	print "  sbjct contig: $_->{'start'} - $_->{'stop'}\n";
#	print "         iden = $_->{'iden'}; cons = $_->{'cons'}\n";
	$self->{'_length_aln_sbjct'} += $_->{'stop'} - $_->{'start'} + 1;
    }
#   <STDIN>;
    
    if($qoverlap) {
	if($soverlap) { $self->ambiguous_aln('qs'); 
#			print "\n*** AMBIGUOUS ALIGNMENT: Query and Sbjct\n\n";
		      }
	else { $self->ambiguous_aln('q');
#	       print "\n*** AMBIGUOUS ALIGNMENT: Query\n\n";
	   }
    } elsif($soverlap) { 
	$self->ambiguous_aln('s'); 
#	print "\n*** AMBIGUOUS ALIGNMENT: Sbjct\n\n";
    }

    # Adjust length based on BLAST flavor.
    my $prog = $self->parent->program;
    if($prog eq 'TBLASTN') {
	$self->{'_length_aln_sbjct'} /= 3;
    } elsif($prog eq 'BLASTX' ) {
	$self->{'_length_aln_query'} /= 3;
    } elsif($prog eq 'TBLASTX') {
	$self->{'_length_aln_query'} /= 3;
	$self->{'_length_aln_sbjct'} /= 3;
    }
}



=head2 _adjust_contigs

 Usage     : n/a; called automatically during object construction.
 Purpose   : Builds HSP contigs for a given BLAST hit.
           : Utility method called by _tile_hsps()
 Returns   : 
 Argument  : 
 Throws    : Exceptions propagated from Bio::Tools::Blast::HSP::matches()
           : for invalid sub-sequence ranges.
 Status    : Experimental
 Comments  : This method does not currently support gapped alignments.
           : Also, it does not keep track of the number of HSPs that
           : overlap within the amount specified by overlap().
           : This will lead to significant tracking errors for large
           : overlap values.

See Also   : L<overlap()|overlap>, L<_tile_hsps()|_tile_hsps>, B<Bio::Tools::Blast::HSP>::matches

=cut

#-------------------
sub _adjust_contigs {
#-------------------
    my ($seqType, $hsp, $start, $stop, $contigs_ref, $max_overlap) = @_;

    my $overlap = 0;
    my ($numID, $numCons);

#    print "Testing $seqType data: HSP (${\$hsp->name});  $start, $stop\n"; 
    foreach(@$contigs_ref) {
#	print "  Contig: $_->{'start'} - $_->{'stop'}, iden= $_->{'iden'}, cons= $_->{'cons'}\n";

	## Test special case of a nested HSP. Skip it.
	if($start >= $_->{'start'} and $stop <= $_->{'stop'}) { 
#	    print "----> Nested HSP. Skipping.\n";
	    $overlap = 1; 
	    next;
	}

	## Test for overlap at beginning of contig.
	if($start < $_->{'start'} and $stop > ($_->{'start'} + $max_overlap)) { 
#	    print "----> Overlaps beg: existing beg,end: $_->{'start'},$_->{'stop'}, new beg,end: $start,$stop\n";<STDIN>;
	    # Collect stats over the non-overlapping region.
	    eval {
		($numID, $numCons) = $hsp->matches(-SEQ   =>$seqType, 
						   -START =>$start, 
						   -STOP  =>$_->{'start'}-1); 
	    };
	    if($@) { warn "\a\n$@\n"; }
	    else {
		$_->{'start'} = $start; # Assign a new start coordinate to the contig
		$_->{'iden'} += $numID; # and add new data to #identical, #conserved.
		$_->{'cons'} += $numCons;
		$overlap     = 1; 
	    }
	}

	## Test for overlap at end of contig.
	if($stop > $_->{'stop'} and $start < ($_->{'stop'} - $max_overlap)) { 
#	    print "----> Overlaps end: existing beg,end: $_->{'start'},$_->{'stop'}, new beg,end: $start,$stop\n";<STDIN>;
	    # Collect stats over the non-overlapping region.
	    eval {
		($numID,$numCons) = $hsp->matches(-SEQ   =>$seqType, 
						  -START =>$_->{'stop'}, 
						  -STOP  =>$stop); 
	    };
	    if($@) { warn "\a\n$@\n"; }
	    else {
		$_->{'stop'}  = $stop;  # Assign a new stop coordinate to the contig
		$_->{'iden'} += $numID; # and add new data to #identical, #conserved.
		$_->{'cons'} += $numCons;
		$overlap    = 1; 
	    }
	}
	$overlap && do {
#		print " New Contig data:\n";
#		print "  Contig: $_->{'start'} - $_->{'stop'}, iden= $_->{'iden'}, cons= $_->{'cons'}\n";
		last;
	    };
    }
    ## If there is no overlap, add the complete HSP data.
    !$overlap && do {
#	print "No overlap. Adding new contig.\n";
	($numID,$numCons) = $hsp->matches(-SEQ=>$seqType); 
	push @$contigs_ref, {'start'=>$start, 'stop'=>$stop,
			     'iden'=>$numID,  'cons'=>$numCons };
    };
#    <STDIN>;
    $overlap;
}


=head2 ambiguous_aln

 Usage     : $ambig_code = $sbjct_object->ambiguous_aln();
 Purpose   : Sets/Gets ambiguity code data member.
 Example   : (see usage)
 Returns   : String = 'q', 's', 'qs', '-'
           :   'q'  = query sequence contains overlapping sub-sequences 
           :          while sbjct does not.
           :   's'  = sbjct sequence contains overlapping sub-sequences 
           :          while query does not.
           :   'qs' = query and sbjct sequence contains overlapping sub-sequences
           :          relative to each other.
           :   '-'  = query and sbjct sequence do not contains multiple domains 
           :          relative to each other OR both contain the same distribution
           :          of similar domains.
 Argument  : n/a
 Throws    : n/a
 Status    : Experimental

See Also   : L<_tile_hsps()|_tile_hsps>,  L<HSP Tiling and Ambiguous Alignments>

=cut

#--------------------
sub ambiguous_aln { 
#--------------------
    my $self = shift;
    if(@_) { $self->{'_ambiguous_aln'} = shift; }
    $self->{'_ambiguous_aln'} || '-';
}



=head2 overlap

 Usage     : $blast_object->overlap( [integer] );
 Purpose   : Gets/Sets the allowable amount overlap between different HSP sequences.
 Example   : $blast_object->overlap(5);
           : $overlap = $blast_object->overlap();
 Returns   : Integer.
 Argument  : integer.
 Throws    : n/a
 Status    : Experimental
 Comments  : Any two HSPs whose sequences overlap by less than or equal
           : to the overlap() number of resides will be considered separate HSPs
           : and will not get tiled by _adjust_contigs().

See Also   : L<_adjust_contigs()|_adjust_contigs>, L<BUGS | BUGS>

=cut

#-------------
sub overlap { 
#-------------
    my $self = shift; 
    if(@_) { $self->{'_overlap'} = shift; }
    defined $self->{'_overlap'} ? $self->{'_overlap'} : 0;
}




=head2 score

 Usage     : $sbjct_object->score();
 Purpose   : Gets the BLAST score of the best HSP for the current Blast hit.
 Example   : $score = $sbjct_object->score();
 Returns   : Integer
 Argument  : n/a
 Throws    : n/a

See Also   : L<bits()|bits>

=cut

#----------
sub score { 
#----------
    my $self = shift;  

    # The check for $self->{'_score'} is a remnant from the 'query' mode days
    # in which the sbjct object would collect data from the description line only.

    my ($score);
    if(not defined($self->{'_score'})) {
	$score = $self->hsp->score;
    } else {
	$score = $self->{'_score'}; 
    } 
    return $score;
}



=head2 bits

 Usage     : $sbjct_object->bits();
 Purpose   : Gets the BLAST bit score of the best HSP for the current Blast hit.
 Example   : $bits = $sbjct_object->bits();
 Returns   : Integer
 Argument  : n/a
 Throws    : Exception if bit score is not set.
 Comments  : For BLAST1, the non-bit score is listed in the summary line.

See Also   : L<score()|score>

=cut

#---------
sub bits { 
#---------
    my $self = shift; 

    # The check for $self->{'_bits'} is a remnant from the 'query' mode days
    # in which the sbjct object would collect data from the description line only.

    my ($bits);
    if(not defined($self->{'_bits'})) {
	$bits = $self->hsp->bits;
    } else {
	$bits = $self->{'_bits'}; 
    } 
    return $bits;
}



=head2 n

 Usage     : $sbjct_object->n();
 Purpose   : Gets the N number for the current Blast hit.
           : This is the number of HSPs in the set which was ascribed
           : the lowest P-value (listed on the description line).
           : This number is not the same as the total number of HSPs.
           : To get the total number of HSPs, use num_hsps().
 Example   : $n = $sbjct_object->n();
 Returns   : Integer
 Argument  : n/a
 Throws    : Exception if HSPs have not been set (BLAST2 reports).
 Comments  : Note that the N parameter is not reported in gapped BLAST2.
           : Calling n() on such reports will result in a call to num_hsps().
           : The num_hsps() method will count the actual number of
           : HSPs in the alignment listing, which may exceed N in
           : some cases.

See Also   : L<num_hsps()|num_hsps>

=cut

#-----
sub n { 
#-----
    my $self = shift; 

    # The check for $self->{'_n'} is a remnant from the 'query' mode days
    # in which the sbjct object would collect data from the description line only.

    my ($n);
    if(not defined($self->{'_n'})) {
	$n = $self->hsp->n;
    } else {
	$n = $self->{'_n'}; 
    } 
    $n ||= $self->num_hsps;

    return $n;
}



=head2 frame

 Usage     : $sbjct_object->frame();
 Purpose   : Gets the reading frame for the hit sequence (TBLASTN/X only).
 Example   : $frame = $sbjct_object->frame();
 Returns   : Integer (-3 .. +3).
 Argument  : n/a
 Throws    : Exception if HSPs have not been set (BLAST2 reports).

See Also   : L<hsps()|hsps>

=cut

#----------
sub frame { 
#----------
    my $self = shift; 

    # The check for $self->{'_frame'} is a remnant from the 'query' mode days
    # in which the sbjct object would collect data from the description line only.

    my ($frame);
    if(not defined($self->{'_frame'})) {
	$frame = $self->hsp->frame;
    } else {
	$frame = $self->{'_frame'}; 
    } 
    return $frame;
}





=head2 p

 Usage     : $sbjct_object->p( [format] );
 Purpose   : Get the P-value for the given BLAST hit.
           : (Note that P-values are not provided with NCBI Blast2 reports).
 Example   : $p =  $sbjct->p;
           : $p =  $sbjct->p('exp');  # get exponent only.
           : ($num, $exp) =  $sbjct->p('parts');  # split sci notation into parts
 Returns   : Float or scientific notation number (the raw P-value, DEFAULT).
           : Integer if format == 'exp' (the magnitude of the base 10 exponent).
           : 2-element list (float, int) if format == 'parts' and P-value
           :                is in scientific notation (See Comments).
 Argument  : format: string of 'raw' | 'exp' | 'parts'
           :    'raw' returns value given in report. Default. (1.2e-34)
           :    'exp' returns exponent value only (34)
           :    'parts' returns the decimal and exponent as a 
           :            2-element list (1.2, -34) (See Comments).
 Throws    : Exception if the P-value is not defined, which will occur
           : with any NCBI Blast2 report.
 Comments  : Using the 'parts' argument is not recommended since it will not
           : work as expected if the P-value is not in scientific notation.
           : That is, floats are not converted into sci notation before
           : splitting into parts.

See Also   : L<expect()|expect>, L<signif()|signif>, L<get_exponent()|get_exponent>

=cut

#--------
sub p { 
#--------
# Some duplication of logic for p(), expect() and signif() for the sake of performance.
    my ($self, $fmt) = @_;

    my ($val);
    $fmt ||= $_signif_fmt;

    # $val can be zero.
    if(not defined($val = $self->{'_p'})) {
	## P-value not defined, must be a NCBI Blast2 report.
	my $note = '';
	if($self->parent->_layout() == 2) {
	    $note = "Blast2 does not report P-values. Use expect() instead.";
	}
	$self->throw("Can't get P-value: undefined.", $note); 
    }

    return $val if not $fmt or $fmt =~ /^raw/i;
    ## Special formats: exponent-only or as list.
    return &get_exponent($val) if $fmt =~ /^exp/i;
    return (split (/eE/, $val)) if $fmt =~ /^parts/i;

    ## Default: return the raw P-value.
    return $val;
}



=head2 expect

 Usage     : $sbjct_object->expect( [format] );
 Purpose   : Get the Expect value for the given BLAST hit.
 Example   : $e =  $sbjct->expect;
           : $e =  $sbjct->expect('exp');  # get exponent only.
           : ($num, $exp) = $sbjct->expect('parts');  # split sci notation into parts
 Returns   : Float or scientific notation number (the raw expect value, DEFAULT).
           : Integer if format == 'exp' (the magnitude of the base 10 exponent).
           : 2-element list (float, int) if format == 'parts' and Expect 
           :                is in scientific notation (see Comments).
 Argument  : format: string of 'raw' | 'exp' | 'parts'
           :    'raw' returns value given in report. Default. (1.2e-34)
           :    'exp' returns exponent value only (34)
           :    'parts' returns the decimal and exponent as a 
           :            2-element list (1.2, -34)  (see Comments).
 Throws    : Exception if the Expect value is not defined.
 Comments  : Using the 'parts' argument is not recommended since it will not
           : work as expected if the expect value is not in scientific notation.
           : That is, floats are not converted into sci notation before
           : splitting into parts.

See Also   : L<p()|p>, L<signif()|signif>, L<get_exponent()|get_exponent>

=cut

#-----------
sub expect { 
#-----------
# Some duplication of logic for p(), expect() and signif() for the sake of performance.
    my ($self, $fmt) = @_;

    my $val = $self->{'_expect'};
    $fmt ||= $_signif_fmt;

    # $val can be zero.
    defined($val) or $self->throw("Can't get Expect value: HSPs may not have been set.");

    return $val if not $fmt or $fmt =~ /^raw/i;
    ## Special formats: exponent-only or as list.
    return &get_exponent($val) if $fmt =~ /^exp/i;
    return (split (/eE/, $val)) if $fmt =~ /^parts/i;

    ## Default: return the raw Expect-value.
    return $val;
}



=head2 signif

 Usage     : $sbjct_object->signif( [format] );
 Purpose   : Get the P or Expect value for the given BLAST hit.
           : The value returned is the one which is reported in the description
           : section of the Blast report. For Blast1 and WU-Blast2, this
           : is a P-value, for Blast2, it is an Expect value.
 Example   : $obj->signif()        # returns 1.3e-34
           : $obj->signif('exp')   # returns -34
           : $obj->signif('parts') # returns (1.3, -34)
 Returns   : Float or scientific notation number (the raw P/Expect value, DEFAULT).
           : Integer if format == 'exp' (the magnitude of the base 10 exponent).
           : 2-element list (float, int) if format == 'parts' and P/Expect value
           :                is in scientific notation (see Comments).
 Argument  : format: string of 'raw' | 'exp' | 'parts'
           :    'raw' returns value given in report. Default. (1.2e-34)
           :    'exp' returns exponent value only (34)
           :    'parts' returns the decimal and exponent as a 
           :            2-element list (1.2, -34)  (see Comments).
 Throws    : n/a
 Status    : Deprecated. Use p() or expect().
 Comments  : The signif() method provides a way to deal with the fact that
           : Blast1 and Blast2 formats differ in what is reported in the
           : description lines of each hit in the Blast report. The signif()
           : method frees any client code from having to know if this is a P-value
           : or an Expect value, making it easier to write code that can process 
           : both Blast1 and Blast2 reports. This is not necessarily a good thing, since
           : one should always know when one is working with P-values or
           : Expect values (hence the deprecated status).
           : Use of expect() is recommended since all hits will have an Expect value.
           :
           : Using the 'parts' argument is not recommended since it will not
           : work as expected if the expect value is not in scientific notation.
           : That is, floats are not converted into sci notation before
           : splitting into parts.

See Also   : L<p()|p>, L<expect()|expect>, L<get_exponent()|get_exponent>

=cut

#-------------
sub signif {
#-------------
# Some duplication of logic for p(), expect() and signif() for the sake of performance.
    my ($self, $fmt) = @_;

    my $val = defined($self->{'_p'}) ? $self->{'_p'} : $self->{'_expect'};
    $fmt ||= $_signif_fmt;

    # $val can be zero.
    defined($val) or $self->throw("Can't get P- or Expect value: HSPs may not have been set.");

    return $val if not $fmt or $fmt =~ /^raw/i;
    ## Special formats: exponent-only or as list.
    return &get_exponent($val) if $fmt =~ /^exp/i;
    return (split (/eE/, $val)) if $fmt =~ /^parts/i;

    ## Default: return the raw P/Expect-value.
    return $val;
}



=head2 desc

 Usage     : $sbjct_object->desc( [integer] );
 Purpose   : Get the description for the given BLAST hit.
 Example   : (see usage)
 Returns   : String
 Argument  : Integer (optional) indicating the desired length of the
           : description string to be returned.
 Throws    : n/a

See Also   : L<_set_desc()|_set_desc>

=cut

#---------
sub desc { 
#---------
    my( $self, $len ) = @_; 
    $len = (defined $len) ? $len : (CORE::length $self->{'_desc'});
    substr( $self->{'_desc'}, 0 ,$len ); 
}



=head2 database

 Usage     : $sbjct_object->database();
 Purpose   : Get the name of the database for the hit sequence.
 Example   : (see usage)
 Returns   : String
 Argument  : n/a
 Throws    : n/a
 Status    : Deprecated. Use Bio::Tools::Blast::database()
             Extracting database name from the seq identifier is error prone.
 Comments  : Database id should be the same for all hits in a given 
           : BLAST report, however, they do not always have the same
           : name as the database name extraced by the Blast.pm object.
           : The Sbjct.pm database id is obtained from the summary line.

=cut

#--------------
sub database { 
    my $self = shift; 
    $self->warn("Bio::Tools::Sbjct::database() is deprecated.\nNo useful information is provided by this method.\nUse Bio::Tools::Blast::database().\n");
    return $self->{'_db'};
}
#--------------




=head2 hsps

 Usage     : $sbjct_object->hsps();
 Purpose   : Get a list containing all HSP objects.
           : Get the numbers of HSPs for the current hit.
 Example   : @hsps = $sbjct_object->hsps();
           : $num  = $sbjct_object->hsps();  # alternatively, use num_hsps()
 Returns   : Array context : list of Bio::Tools::Blast::HSP.pm objects.
           : Scalar context: integer (number of HSPs).
           :                 (Equivalent to num_hsps()).
 Argument  : n/a. Relies on wantarray
 Throws    : Exception if the HSPs have not been collected.

See Also   : L<hsp()|hsp>, L<num_hsps()|num_hsps>, L<_set_hsps()|_set_hsps>

=cut

#---------
sub hsps {
#---------
    my $self = shift;

    if (not ref $self->{'_hsps'}) {
	$self->throw("Can't get HSPs: data not collected.");
    }

    return wantarray 
        #  returning list containing all HSPs.
	? @{$self->{'_hsps'}}
        #  returning number of HSPs.
        : scalar(@{$self->{'_hsps'}});
}



=head2 hsp

 Usage     : $sbjct_object->hsp( [string] );
 Purpose   : Get a single HSP.pm object for the present Sbjct.pm object.
 Example   : $hspObj  = $sbjct_object->hsp;  # same as 'best'
           : $hspObj  = $sbjct_object->hsp('best');
           : $hspObj  = $sbjct_object->hsp('worst');
 Returns   : Object reference for a Bio::Tools::Blast::HSP.pm object.
 Argument  : String (or no argument).
           :   No argument (default) = highest scoring HSP (same as 'best').
           :   'best' or 'first' = highest scoring HSP.
           :   'worst' or 'last' = lowest scoring HSP.
 Throws    : Exception if the HSPs have not been collected.
           : Exception if an unrecognized argument is used.

See Also   : L<hsps()|hsps>, L<num_hsps()|num_hsps>, L<_set_hsps()|_set_hsps>

=cut

#----------
sub hsp {
#----------
    my( $self, $option ) = @_;
    $option ||= 'best';
    
    if (not ref $self->{'_hsps'}) {
	$self->throw("Can't get HSPs: data not collected.");
    }

    my @hsps = @{$self->{'_hsps'}};
    
    return $hsps[0]      if $option =~ /best|first|1/i;
    return $hsps[$#hsps] if $option =~ /worst|last/i;

    $self->throw("Can't get HSP for: $option", 
		 "Valid arguments: 'best', 'worst'");
}



=head2 num_hsps

 Usage     : $sbjct_object->num_hsps();
 Purpose   : Get the number of HSPs for the present Blast hit.
 Example   : $nhsps = $sbjct_object->num_hsps();
 Returns   : Integer
 Argument  : n/a
 Throws    : Exception if the HSPs have not been collected.

See Also   : L<hsps()|hsps>

=cut

#-------------
sub num_hsps {
#-------------
    my $self = shift;
    
    if (not defined $self->{'_hsps'}) {
	$self->throw("Can't get HSPs: data not collected.");
    }

    return scalar(@{$self->{'_hsps'}});
}



=head2 length

 Usage     : $sbjct_object->length();
 Purpose   : Get the total length of the hit sequence.
 Example   : $len    = $sbjct_object->length();
 Returns   : Integer 
 Argument  : n/a
 Throws    : n/a
 Comments  : Developer note: when using the built-in length function within
           : this module, call it as CORE::length().

See Also   : L<logical_length()|logical_length>,  L<length_aln()|length_aln>

=cut

#-----------
sub length {
#-----------
    my $self = shift;
    $self->{'_length'}; 
}    


=head2 logical_length

 Usage     : $sbjct_object->logical_length( [seq_type] );
           : (mostly intended for internal use).
 Purpose   : Get the logical length of the hit sequence.
           : If the Blast is a TBLASTN or TBLASTX, the returned length 
           : is the length of the would-be amino acid sequence (length/3).
           : For all other BLAST flavors, this function is the same as length().
 Example   : $len    = $sbjct_object->logical_length();
 Returns   : Integer 
 Argument  : seq_type = 'query' or 'sbjct' (default = 'query')
 Throws    : n/a
 Comments  : This is important for functions like frac_aligned_query()
           : which need to operate in amino acid coordinate space when dealing
           : with [T]BLAST[NX] type reports.

See Also   : L<length()|length>, L<frac_aligned_query()|frac_aligned_query>, L<frac_aligned_hit()|frac_aligned_hit>

=cut

#--------------------
sub logical_length {
#--------------------
    my $self = shift;
    my $seqType = shift || 'query';

    # Return logical sbjct length
    $seqType eq 'sbjct' and return 
	$self->{'_logical_length'} || $self->{'_length'}; 

    # Otherwise, return logical query length
    my $qlen = $self->parent->length;    

    # Adjust length based on BLAST flavor.
    my $prog = $self->parent->program;
    if($prog =~ /T?BLASTX/ ) {
	$qlen /= 3;
    }
    return $qlen;
}    



=head2 length_aln

 Usage     : $sbjct_object->length_aln( [seq_type] );
 Purpose   : Get the total length of the aligned region for query or sbjct seq.
           : This number will include all HSPs
 Example   : $len    = $sbjct_object->length_aln(); # default = query
           : $lenAln = $sbjct_object->length_aln('query');
 Returns   : Integer 
 Argument  : seq_Type = 'query' | 'sbjct'  (Default = 'query')
 Throws    : Exception if the argument is not recognized.
 Comments  : This method will report the logical length of the alignment,
           : meaning that for TBLAST[NX] reports, the length is reported
           : using amino acid coordinate space (i.e., nucleotides / 3).
           : 
           : This method requires that all HSPs be tiled. If they have not
           : already been tiled, they will be tiled first.
           : If you don't want the tiled data, iterate through each HSP
           : calling length() on each (use hsps() to get the HSPs).

See Also   : L<length()|length>, L<frac_aligned_query()|frac_aligned_query>, L<frac_aligned_hit()|frac_aligned_hit>, L<gaps()|gaps>, L<_tile_hsps()|_tile_hsps>, B<Bio::Tools::Blast::HSP::length()>

=cut

#---------------'
sub length_aln {
#---------------
    my( $self, $type ) = @_;
    
    $type ||= 'query';

    $self->_tile_hsps() if not $self->{'_tile_hsps'};

    my $data = $self->{'_length_aln_'.$type};
    
    ## If we don't have data, figure out what went wrong.
    if(!$data) {
	$self->throw("Can't get length aln for sequence type \"$type\"",
		     "Valid types are 'query', 'sbjct'");
    }		
    $data;
}    


=head2 gaps

 Usage     : $sbjct_object->gaps( [seq_type] );
 Purpose   : Get the number of gaps in the aligned query, sbjct, or both sequences.
           : Data is summed across all HSPs.
 Example   : $qgaps = $sbjct_object->gaps('query');
           : $sgaps = $sbjct_object->gaps('sbjct');
           : $tgaps = $sbjct_object->gaps();    # default = total (query + sbjct)
 Returns   : scalar context: integer
           : array context without args: two-element list of integers  
           :    (queryGaps, sbjctGaps)
           : Array context can be "induced" by providing an argument of 'list' or 'array'.
 Argument  : seq_type: 'query' | 'sbjct' | 'total' | 'list'  (default = 'total')
 Throws    : n/a
 Comments  : If you need data for each HSP, use hsps() and then interate
           : through each HSP object.
           : This method requires that all HSPs be tiled. If they have not
           : already been tiled, they will be tiled first.
           : Not relying on wantarray since that will fail in situations 
           : such as printf "%d", $hit->gaps() in which you might expect to 
           : be printing the total gaps, but evaluates to array context.

See Also   : L<length_aln()|length_aln>

=cut

#----------
sub gaps {
#----------
    my( $self, $seqType ) = @_;

    $seqType ||= (wantarray ? 'list' : 'total');

    $self->_tile_hsps() if not $self->{'_tile_hsps'};

    $seqType = lc($seqType);

    if($seqType =~ /list|array/i) {
	return ($self->{'_gaps_query'}, $self->{'_gaps_sbjct'});
    }

    if($seqType eq 'total') {
	return ($self->{'_gaps_query'} + $self->{'_gaps_sbjct'}) || 0;
    } else {
	return $self->{'_gaps_'.$seqType} || 0;
    }
}    



=head2 matches

 Usage     : $sbjct_object->matches( [class] );
 Purpose   : Get the total number of identical or conserved matches 
           : (or both) across all HSPs.
           : (Note: 'conservative' matches are indicated as 'positives' 
	   :         in the Blast report.)
 Example   : ($id,$cons) = $sbjct_object->matches(); # no argument
           : $id = $sbjct_object->matches('id');
           : $cons = $sbjct_object->matches('cons'); 
 Returns   : Integer or a 2-element array of integers 
 Argument  : class = 'id' | 'cons' OR none. 
           : If no argument is provided, both identical and conservative 
           : numbers are returned in a two element list.
           : (Other terms can be used to refer to the conservative
           :  matches, e.g., 'positive'. All that is checked is whether or
           :  not the supplied string starts with 'id'. If not, the 
           : conservative matches are returned.)
 Throws    : Exception if the requested data cannot be obtained.
 Comments  : If you need data for each HSP, use hsps() and then interate
           : through the HSP objects.
           : Does not rely on wantarray to return a list. Only checks for
           : the presence of an argument (no arg = return list).

See Also   : B<Bio::Tools::Blast::HSP::matches()>, L<hsps()|hsps>

=cut

#---------------
sub matches {
#---------------
    my( $self, $arg) = @_;
    my(@data,$data);

    if(!$arg) {
	@data = ($self->{'_totalIdentical'}, $self->{'_totalConserved'});

	return @data if @data;

    } else {

	if($arg =~ /^id/i) { 
	    $data = $self->{'_totalIdentical'};
	} else {
	    $data = $self->{'_totalConserved'};
	}
	return $data if $data;
    }
    
    ## Something went wrong if we make it to here.
    $self->throw("Can't get identical or conserved data: no data.");
}


=head2 start

 Usage     : $sbjct->start( [seq_type] );
 Purpose   : Gets the start coordinate for the query, sbjct, or both sequences
           : in the Sbjct object. If there is more than one HSP, the lowest start
           : value of all HSPs is returned.
 Example   : $qbeg = $sbjct->start('query');
           : $sbeg = $sbjct->start('sbjct');
           : ($qbeg, $sbeg) = $sbjct->start();
 Returns   : scalar context: integer 
           : array context without args: list of two integers (queryStart, sbjctStart)
           : Array context can be "induced" by providing an argument of 'list' or 'array'.
 Argument  : In scalar context: seq_type = 'query' or 'sbjct'
           :  (case insensitive). If not supplied, 'query' is used.
 Throws    : n/a
 Comments  : This method requires that all HSPs be tiled. If there is more than one
           : HSP and they have not already been tiled, they will be tiled first.
           : Remember that the start and end coordinates of all HSPs are 
           : normalized so that start < end. Strand information can only be
           : obtained on an HSP-by-HSP basis by calling $hsp->strand().

See Also   : L<end()|end>, L<range()|range>, L<HSP Tiling and Ambiguous Alignments>, B<Bio::Tools::Blast::HSP::start>()

=cut

#----------
sub start {
#----------
    my ($self, $seqType) = @_;

    $seqType ||= (wantarray ? 'list' : 'query');

    # If there is only one HSP, defer this call to the solitary HSP.
    if($self->num_hsps == 1) {
	return $self->hsp->start($seqType);
    } else {
	$self->_tile_hsps() if not $self->{'_tile_hsps'};
	if($seqType =~ /list|array/i) {
	    return ($self->{'_queryStart'}, $self->{'_sbjctStart'});
	} else {
	    ## Sensitive to member name changes.
	    $seqType = "_\L$seqType\E";
	    return $self->{$seqType.'Start'};
	}
    }
}


=head2 end

 Usage     : $sbjct->end( [seq_type] );
 Purpose   : Gets the end coordinate for the query, sbjct, or both sequences
           : in the Sbjct object. If there is more than one HSP, the largest end
           : value of all HSPs is returned.
 Example   : $qend = $sbjct->end('query');
           : $send = $sbjct->end('sbjct');
           : ($qend, $send) = $sbjct->end();
 Returns   : scalar context: integer
           : array context without args: list of two integers (queryEnd, sbjctEnd)
           : Array context can be "induced" by providing an argument of 'list' or 'array'.
 Argument  : In scalar context: seq_type = 'query' or 'sbjct'
           :  (case insensitive). If not supplied, 'query' is used.
 Throws    : n/a
 Comments  : This method requires that all HSPs be tiled. If there is more than one
           : HSP and they have not already been tiled, they will be tiled first.
           : Remember that the start and end coordinates of all HSPs are 
           : normalized so that start < end. Strand information can only be
           : obtained on an HSP-by-HSP basis by calling $hsp->strand().

See Also   : L<start()|start>, L<range()|range>, L<HSP Tiling and Ambiguous Alignments>, B<Bio::Tools::Blast::HSP::end>()

=cut

#----------
sub end {
#----------
    my ($self, $seqType) = @_;

    $seqType ||= (wantarray ? 'list' : 'query');

    # If there is only one HSP, defer this call to the solitary HSP.
    if($self->num_hsps == 1) {
	return $self->hsp->end($seqType);
    } else {
	$self->_tile_hsps() if not $self->{'_tile_hsps'};
	if($seqType =~ /list|array/i) {
	    return ($self->{'_queryStop'}, $self->{'_sbjctStop'});
	} else {
	    ## Sensitive to member name changes.
	    $seqType = "_\L$seqType\E";
	    return $self->{$seqType.'Stop'};
	}
    }
}

=head2 range

 Usage     : $sbjct->range( [seq_type] );
 Purpose   : Gets the (start, end) coordinates for the query or sbjct sequence
           : in the HSP alignment.
 Example   : ($qbeg, $qend) = $sbjct->range('query');
           : ($sbeg, $send) = $sbjct->range('sbjct');
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
    $seqType ||= 'query';
    return ($self->start($seqType), $self->end($seqType));
}


=head2 frac_identical

 Usage     : $sbjct_object->frac_identical( [seq_type] );
 Purpose   : Get the overall fraction of identical positions across all HSPs.
           : The number refers to only the aligned regions and does not
           : account for unaligned regions in between the HSPs, if any.
 Example   : $frac_iden = $sbjct_object->frac_identical('query');
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
           :
           : If you need data for each HSP, use hsps() and then iterate
           : through the HSP objects.
           : This method requires that all HSPs be tiled. If they have not
           : already been tiled, they will be tiled first.

See Also   : L<frac_conserved()|frac_conserved>, L<frac_aligned_query()|frac_aligned_query>, L<matches()|matches>, L<_tile_hsps()|_tile_hsps>

=cut

#------------------
sub frac_identical {
#------------------
    my ($self, $seqType) = @_;
    $seqType ||= 'total';

    ## Sensitive to member name format.
    $seqType = lc($seqType);

    $self->_tile_hsps() if not $self->{'_tile_hsps'};

    sprintf( "%.2f", $self->{'_totalIdentical'}/$self->{'_length_aln_'.$seqType});
}



=head2 frac_conserved

 Usage     : $sbjct_object->frac_conserved( [seq_type] );
 Purpose   : Get the overall fraction of conserved positions across all HSPs.
           : The number refers to only the aligned regions and does not
           : account for unaligned regions in between the HSPs, if any.
 Example   : $frac_cons = $sbjct_object->frac_conserved('sbjct');
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
           : To get the fraction conserved among only the aligned residues,
           : ignoring the gaps, call this method with an argument of 'query'
           : or 'sbjct'.
           :
           : If you need data for each HSP, use hsps() and then interate
           : through the HSP objects.
           : This method requires that all HSPs be tiled. If they have not
           : already been tiled, they will be tiled first.

See Also   : L<frac_identical()|frac_identical>, L<matches()|matches>, L<_tile_hsps()|_tile_hsps>

=cut

#--------------------
sub frac_conserved {
#--------------------
    my ($self, $seqType) = @_;
    $seqType ||= 'total';

    ## Sensitive to member name format.
    $seqType = lc($seqType);

    $self->_tile_hsps() if not $self->{'_tile_hsps'};

    sprintf( "%.2f", $self->{'_totalConserved'}/$self->{'_length_aln_'.$seqType});
}




=head2 frac_aligned_query

 Usage     : $sbjct_object->frac_aligned_query();
 Purpose   : Get the fraction of the query sequence which has been aligned
           : across all HSPs (not including intervals between non-overlapping
           : HSPs).
 Example   : $frac_alnq = $sbjct_object->frac_aligned_query();
 Returns   : Float (2-decimal precision, e.g., 0.75).
 Argument  : n/a
 Throws    : n/a
 Comments  : If you need data for each HSP, use hsps() and then interate
           : through the HSP objects.
           : To compute the fraction aligned, the logical length of the query
           : sequence is used, meaning that for [T]BLASTX reports, the 
           : full length of the query sequence is converted into amino acids
           : by dividing by 3. This is necessary because of the way 
           : the lengths of aligned sequences are computed.
           : This method requires that all HSPs be tiled. If they have not
           : already been tiled, they will be tiled first.

See Also   : L<frac_aligned_hit()|frac_aligned_hit>, L<_tile_hsps()|_tile_hsps>, L<logical_length()|logical_length>, L<length_aln()|length_aln>

=cut

#----------------------
sub frac_aligned_query {
#----------------------
    my $self = shift;

    $self->_tile_hsps() if not $self->{'_tile_hsps'};

    sprintf( "%.2f", $self->{'_length_aln_query'}/$self->logical_length('query'));
}



=head2 frac_aligned_hit

 Usage     : $sbjct_object->frac_aligned_hit();
 Purpose   : Get the fraction of the hit (sbjct) sequence which has been aligned
           : across all HSPs (not including intervals between non-overlapping
           : HSPs).
 Example   : $frac_alnq = $sbjct_object->frac_aligned_hit();
 Returns   : Float (2-decimal precision, e.g., 0.75).
 Argument  : n/a
 Throws    : n/a
 Comments  : If you need data for each HSP, use hsps() and then interate
           : through the HSP objects.
           : To compute the fraction aligned, the logical length of the sbjct
           : sequence is used, meaning that for TBLAST[NX] reports, the 
           : full length of the sbjct sequence is converted into amino acids
           : by dividing by 3. This is necessary because of the way 
           : the lengths of aligned sequences are computed.
           : This method requires that all HSPs be tiled. If they have not
           : already been tiled, they will be tiled first.

See Also   : L<frac_aligned_query()|frac_aligned_query>, L<matches()|matches>, L<_tile_hsps()|_tile_hsps>, L<logical_length()|logical_length>, L<length_aln()|length_aln>

=cut

#--------------------
sub frac_aligned_hit {
#--------------------
    my $self = shift;

    $self->_tile_hsps() if not $self->{'_tile_hsps'};

    sprintf( "%.2f", $self->{'_length_aln_sbjct'}/$self->logical_length('sbjct'));
}

# Safety-net methods for those who try don't read or remember the API.
# Redirecting to the proper method. These 'sbjct' versions may be more
# consistent with the API of this module since there are numerous other
# instances of using 'sbjct' in arguments. However, 'sbjct' is a bit tech-ee.

#-----------------------
sub frac_aligned_sbjct {  my $self=shift; $self->frac_aligned_hit(@_); }
#-----------------------  
sub num_unaligned_sbjct {  my $self=shift; $self->num_unaligned_hit(@_); }
#-----------------------  



=head2 num_unaligned_hit

 Usage     : $sbjct_object->num_unaligned_hit();
 Purpose   : Get the number of the unaligned residues in the hit sequence.
           : Sums across all all HSPs.
 Example   : $num_unaln = $sbjct_object->num_unaligned_hit();
 Returns   : Integer
 Argument  : n/a
 Throws    : n/a
 Comments  : See notes regarding logical lengths in the comments for frac_aligned_hit().
           : They apply here as well.
           : If you need data for each HSP, use hsps() and then interate
           : through the HSP objects.
           : This method requires that all HSPs be tiled. If they have not
           : already been tiled, they will be tiled first.

See Also   : L<num_unaligned_query()|num_unaligned_query>, L<_tile_hsps()|_tile_hsps>, L<frac_aligned_hit()|frac_aligned_hit>

=cut

#---------------------
sub num_unaligned_hit {
#---------------------
    my $self = shift;

    $self->_tile_hsps() if not $self->{'_tile_hsps'};

    my $num = $self->logical_length('sbjct') - $self->{'_length_aln_sbjct'};
    ($num < 0 ? 0 : $num );
}


=head2 num_unaligned_query

 Usage     : $sbjct_object->num_unaligned_query();
 Purpose   : Get the number of the unaligned residues in the query sequence.
           : Sums across all all HSPs.
 Example   : $num_unaln = $sbjct_object->num_unaligned_query();
 Returns   : Integer
 Argument  : n/a
 Throws    : n/a
 Comments  : See notes regarding logical lengths in the comments for frac_aligned_query().
           : They apply here as well.
           : If you need data for each HSP, use hsps() and then interate
           : through the HSP objects.
           : This method requires that all HSPs be tiled. If they have not
           : already been tiled, they will be tiled first.

See Also   : L<num_unaligned_hit()|num_unaligned_hit>, L<_tile_hsps()|_tile_hsps>, L<frac_aligned_query()|frac_aligned_query>

=cut

#-----------------------
sub num_unaligned_query {
#-----------------------
    my $self = shift;

    $self->_tile_hsps() if not $self->{'_tile_hsps'};

    my $num = $self->logical_length('query') - $self->{'_length_aln_query'};
    ($num < 0 ? 0 : $num );
}



=head2 seq_inds

 Usage     : $hit->seq_inds( seq_type, class, collapse );
 Purpose   : Get a list of residue positions (indices) across all HSPs
           : for identical or conserved residues in the query or sbjct sequence.
 Example   : @ind = $hit->seq_inds('query', 'identical');
           : @ind = $hit->seq_inds('sbjct', 'conserved');
           : @ind = $hit->seq_inds('sbjct', 'conserved', 1);
 Returns   : Array of integers 
           : May include ranges if collapse is non-zero.
 Argument  : seq_type  = 'query' or 'sbjct'  (default = query)
           : class     = 'identical' or 'conserved' (default = identical)
           :              (can be shortened to 'id' or 'cons')
           :              (actually, anything not 'id' will evaluate to 'conserved').
           : collapse  = boolean, if non-zero, consecutive positions are merged
           :             using a range notation, e.g., "1 2 3 4 5 7 9 10 11" 
           :             collapses to "1-5 7 9-11". This is useful for 
           :             consolidating long lists. Default = no collapse.
 Throws    : n/a.

See Also   : B<Bio::Tools::Blast::HSP::seq_inds()>

=cut

#-------------
sub seq_inds {
#-------------
    my ($self, $seq, $class, $collapse) = @_;

    $seq  ||= 'query';
    $class ||= 'identical';
    $collapse ||= 0;

    my (@inds, $hsp);
    foreach $hsp ($self->hsps) {
	# This will merge data for all HSPs together.
	push @inds, $hsp->seq_inds($seq, $class);
    }
    
    # Need to remove duplicates and sort the merged positions.
    if(@inds) {
	my %tmp = map { $_, 1 } @inds;
	@inds = sort {$a <=> $b} keys %tmp;
    }

    require Bio::Tools::Blast::HSP;

    $collapse ?  &Bio::Tools::Blast::HSP::collapse_nums(@inds) : @inds; 
}


#####################################################################################
##                                  INSTANCE METHODS                               ##
#####################################################################################


=head2 display

 Usage     : $sbjct_object->display( %named_parameters );
 Purpose   : Display information about Bio::Tools::Blast::Sbjct.pm data members
 Example   : $object->display(-SHOW=>'stats');
 Argument  : Named parameters: -SHOW  => 'hsp',
           :                   -WHERE => filehandle (default = STDOUT)
 Returns   : n/a
 Status    : Deprecated, Buggy.
           : Use Blast::table() or Blast::table_tiled() instead.

See Also   : L<_display_stats()|_display_stats>, L<_display_hsps()|_display_hsps>, B<Bio::Root::Object>::display

=cut

#------------
sub display {
#------------
    my( $self, %param) = @_;
     
    $param{-HEADER} = 0;
    $self->SUPER::display(%param);
    
    $self->show =~ /hsp/i and $self->_display_hsps( %param);
}


=head2 _display_stats

 Usage     : n/a; called automatically by display()
 Purpose   : Display information about Bio::Tools::Blast.pm data members.
           : Not tab-delimited.
           : Prints the rank, name, database, score, p, n, length
           : of the hit sequence, length of the aligned region,
           : fraction identical, fraction conserved, and the fraction aligned 
           : for both the query and hit sequences.
 Example   : n/a
 Argument  : one argument = filehandle object.
 Returns   : printf call.
 Status    : Deprecated, Buggy.
           : Use Blast::table() or Blast::table_tiled() instead.

See Also   : L<display()|display>  

=cut

#-------------------
sub _display_stats {
#-------------------
    my( $self, $OUT) = @_;
    my $layout = $self->parent->_layout();
    
    if($layout == 1) {
	printf( $OUT "%-3d %-20s %-11s %-5d %-5d %-9.1e %-9.1e %-4d %-3d %-5d %-5d %-5s %-6.2f %-6.2f  %-4d(%.2f)  %-4d(%.2f)\n",
		$self->rank(), $self->name(),
		($self->database() || 'UNKNOWN DB') ,
		$self->score(),$self->bits(),$self->p(),$self->expect(),
		$self->gaps(), $self->n(), 
		$self->length(), $self->length_aln('query'),
		$self->ambiguous_aln(),
		$self->frac_aligned_query, $self->frac_aligned_hit, 
		$self->matches('iden'), $self->frac_identical('query'), 
		$self->matches('cons'), $self->frac_conserved('query'));
    } else {
	printf( $OUT "%-3d %-20s %-11s %-5d %-5d %-9.1e %-4d %-3d %-5d %-5d %-5s %-6.2f %-6.2f  %-4d(%.2f)  %-4d(%.2f)\n",
		$self->rank(), $self->name(),
		($self->database()  || 'UNKNOWN DB'),
		$self->score(),$self->bits(),$self->expect(),
		$self->gaps(), $self->num_hsps, 
		$self->length(), $self->length_aln('query'),
		$self->ambiguous_aln(),
		$self->frac_aligned_query, $self->frac_aligned_hit, 
		$self->matches('iden'), $self->frac_identical('query'), 
		$self->matches('cons'), $self->frac_conserved('query') );
    }

}


=head2 _display_hsps

 Usage     : n/a; called automatically by display()
 Purpose   : Display information about each HSP in the current BLAST hit.
 Example   : n/a
 Argument  : one argument = filehandle object.
 Returns   : printf call.
 Status    : Experimental

See Also   : L<display()|display>, B<Bio::Tools::Blast::HSP::display()> 

=cut

#----------------
sub _display_hsps { 
#----------------
    my( $self, %param) = @_; 
    my $OUT = $self->fh();
    my $hspCount = 0;
    my $reply = undef;

    not defined $self->{'_hsps'} and do{ print $OUT "\nHSP data not loaded.\n\n"; return; };

#    print $OUT "\n",$self->num_hsps, " HSPs\n\n";
    
    my($hsp);
    foreach $hsp ( $self->hsps() ) {
	$hspCount++;
	print $OUT "\n   ", '-'x25, "< HSP #$hspCount >", '-'x25, "\n";
	    
	    $hsp->display( %param );
	
	if( $hspCount < $self->num_hsps ) {
	    print "\n\n\t--------> <RET> FOR NEXT HSP, q TO QUIT <--------\n";
	    chop( $reply = <STDIN>);
	    $reply =~ /^q/i and return;
	}
    }
}



=head2 homol_data

 Usage     : $data = $sbjct_object->homo_data( %named_params );
 Purpose   : Gets specific similarity data about all HSPs.
 Returns   : String
 Argument  : named parameters forwarded to Bio::Tools::Blast::HSP::homol_data().
 Throws    : n/a
 Status    : Experimental
 Comments  : This is an experimental method used for obtaining an 
           : indication of:
           :   1) how many HSPs are in a Blast alignment
           :   2) how strong the similarity is between sequences in the HSP
           :   3) the endpoints of the alignment (sequence monomer numbers)
           : "Homology data" for each HSP is in the format:
           :  "<integer> <start> <stop>"
           : Data for different HSPs are tab-delimited.

See Also   : B<Bio::Tools::Blast::homol_data()>, B<Bio::Tools::Blast::HSP::homol_data()>

=cut

#---------------
sub homol_data {
#---------------
    my ($self,%param) = @_;
    my $data = $self->name();
    
    foreach ($self->hsps) {
	$data .="\t".$_->homol_data(%param);
    }
    ## Record ambiguous alignment status.
    $data .= "\t".$self->ambiguous_aln();
    $data;
}
	

=head2 is_signif

 Usage     : $sbjct_object->is_signif();
 Purpose   : Determine if the given BLAST hit is significant.
 Example   : 
 Returns   : Boolean
 Argument  : n/a
 Throws    : n/a
 Comments  : Uses criteria defined in the parent Blast.pm object
           : to assess significance. Currently, only relies on
           : P-value and length criteria.
           : This mehtod is largely obsolete since are hits are now by
           : definition significant.

=cut

#---------------
sub is_signif { 
#---------------
    my $self = shift; 
    return ($self->{'_significance'} <= $self->parent->signif and 
	    $self->length > $self->parent->signif_len);
}

#####################################################################################
##                                  CLASS METHODS                                  ##
#####################################################################################

=head1 CLASS METHODS

=head2 get_exponent

 Usage     : &get_exponent( number );
 Purpose   : Determines the power of 10 exponent of an integer, float, 
           : or scientific notation number.
 Example   : &get_exponent("4.0e-206");
           : &get_exponent("0.00032");
           : &get_exponent("10.");
           : &get_exponent("1000.0");
           : &get_exponent("e+83");
 Argument  : Float, Integer, or scientific notation number
 Returns   : Integer representing the exponent part of the number (+ or -).
           : If argument == 0 (zero), return value is "-999".
 Comments  : Exponents are rounded up (less negative) if the mantissa is >= 5.
           : Exponents are rounded down (more negative) if the mantissa is <= -5.
           : This method probably belongs in a more general utility class.

=cut

#------------------
sub get_exponent {
#------------------
    my $data = shift;

    my($num, $exp) = split /[eE]/, $data;

    if( defined $exp) { 
	$num = 1 if not $num;
	$num >= 5 and $exp++;
	$num <= -5 and $exp--;
    } elsif( $num == 0) {
	$exp = -999;
    } elsif( not $num =~ /\./) {
	$exp = CORE::length($num) -1;
    } else {
	$exp = 0;
	$num .= '0' if $num =~ /\.$/;
	my ($c);
	my $rev = 0;
	if($num !~ /^0/) {
	    $num = reverse($num);
	    $rev = 1;
	}
	do { $c = chop($num);
	     $c == 0 && $exp++; 
	 } while( $c ne '.');

	$exp = -$exp if $num == 0 and not $rev;
	$exp -= 1 if $rev;
    }
    return $exp;
}

1;
__END__

#####################################################################################
#                                END OF CLASS                                       #
#####################################################################################


=head1 FOR DEVELOPERS ONLY

=head2 Data Members

Information about the various data members of this module is provided for those 
wishing to modify or understand the code. Two things to bear in mind: 

=over 4

=item 1 Do NOT rely on these in any code outside of this module. 

All data members are prefixed with an underscore to signify that they are private.
Always use accessor methods. If the accessor doesn't exist or is inadequate, 
create or modify an accessor (and let me know, too!). (An exception to this might
be for HSP.pm which is more tightly coupled to Sbjct.pm and
may access Sbjct data members directly for efficiency purposes, but probably 
should not).

=item 2 This documentation may be incomplete and out of date.

It is easy for these data member descriptions to become obsolete as 
this module is still evolving. Always double check this info and search 
for members not described here.

=back

An instance of Bio::Tools::Blast::Sbjct.pm is a blessed reference to a hash containing
all or some of the following fields:

 FIELD           VALUE
 --------------------------------------------------------------
 _hsps          : Array ref for a list of Bio::Tools::Blast::HSP.pm objects.
		:
 _db            : Database identifier from the summary line.
		:
 _desc          : Description data for the hit from the summary line.
		:
 _length        : Total length of the hit sequence. 
		:
 _score         : BLAST score.
		:
 _bits          : BLAST score (in bits). Matrix-independent.
		:
 _p             : BLAST P value. Obtained from summary section. (Blast1/WU-Blast only)
		:
 _expect        : BLAST Expect value. Obtained from summary section.
		:
 _n             : BLAST N value (number of HSPs) (Blast1/WU-Blast2 only)
		:
 _frame         : Reading frame for TBLASTN and TBLASTX analyses.
		:
 _totalIdentical: Total number of identical aligned monomers.
		:
 _totalConserved: Total number of conserved aligned monomers (a.k.a. "positives").
		:
 _overlap       : Maximum number of overlapping residues between adjacent HSPs
		: before considering the alignment to be ambiguous. 
		:
 _ambiguous_aln : Boolean. True if the alignment of all HSPs is ambiguous.
		:
 _length_aln_query : Length of the aligned region of the query sequence.
		   :
 _length_aln_sbjct : Length of the aligned region of the sbjct sequence.


 INHERITED DATA MEMBERS 
 ----------------------
 _name          : From Bio::Root::Object.pm. String representing the name of the 
		: sbjct sequence obtained from the BLAST report.
		:
 _parent        : From Bio::Root::Object.pm. This member contains a reference to the
		: Bio::Tools::Blast.pm object to which this hit belongs.


=cut

1;
