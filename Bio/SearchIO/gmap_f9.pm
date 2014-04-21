# $Id: gmap_f9.pm 15987 2009-08-18 21:08:55Z lstein $
#
# BioPerl module for Bio::SearchIO::gmap_f9
#
# Cared for by George Hartzell <hartzell@alerce.com>
#
# Copyright George Hartzell
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::SearchIO::gmap_f9 - Event generator for parsing gmap reports (Z format)

=head1 SYNOPSIS

   # Do not use this object directly - it is used as part of the
   # Bio::SearchIO system.

    use Bio::SearchIO;
    my $searchio = Bio::SearchIO->new(-format => 'gmap',
                                      -file   => 't/data/her2.gmapz');
    while( my $result = $searchio->next_result ) {
        while( my $hit = $result->next_hit ) {
            while( my $hsp = $hit->next_hsp ) {
                # ...
            }
        }
    }


=head1 DESCRIPTION

This object encapsulated the necessary methods for generating events
suitable for building Bio::Search objects from a GMAP "compressed"
report (from gmap run with -Z flag) Read the L<Bio::SearchIO> for more
information about how to use this.

=head2 REVERSE STRAND AND BIOPERL COORDINATES

I believe that I'm doing the correct thing when reporting hits on the
negative strand of the genome.  In particular, I've compared the
"exons" this code generates with the set returned by ncbi's megablast
web service.  NCBI's hsp's are ordered differently and have a
different genomic location (off by ~18,000,000 bases, padding?) but
the starts, ends, and lengths were similar and my strand handling
matches theirs.  E.g.

   CDNA                            GENOME
 start  end    strand   start           end             strand

blast
  1913	2989	1	86236731	86237808	-1
  1	475	1	86260509	86260983	-1
  1510	1727	1	86240259	86240476	-1
  841	989	1	86243034	86243182	-1
  1381	1514	1	86240630	86240763	-1
  989	1122	1	86242457	86242590	-1
  599	729	1	86247470	86247600	-1
  473	608	1	86259972	86260107	-1
  1255	1382	1	86240837	86240964	-1
  730	842	1	86244040	86244152	-1
  1813	1921	1	86238123	86238231	-1
  1725	1814	1	86239747	86239836	-1
  1167	1256	1	86241294	86241383	-1
  1120	1188	1	86242319	86242387	-1

gmap
  1	475	1	104330509	104330983	-1
  476	600	1	104329980	104330104	-1
  601	729	1	104317470	104317598	-1
  730	841	1	104314041	104314152	-1
  842	989	1	104313034	104313181	-1
  990	1121	1	104312458	104312589	-1
  1122	1187	1	104312320	104312385	-1
  1188	1256	1	104311294	104311362	-1
  1257	1382	1	104310837	104310962	-1
  1383	1511	1	104310633	104310761	-1
  1512	1726	1	104310260	104310474	-1
  1727	1814	1	104309747	104309834	-1
  1815	1917	1	104308127	104308229	-1
  1918	2989	1	104306731	104307802	-1

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

  bioperl-l@bioperl.org                  - General discussion
  http://bioperl.org/wiki/Mailing_lists  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution. Bug reports can be submitted via
the web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - George Hartzell

Email hartzell@alerce.com

=head1 CONTRIBUTORS

Additional contributors names and emails here

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with an underscore (_).

=cut


# Let the code begin...


package Bio::SearchIO::gmap_f9;
use strict;
use warnings;

use Bio::Search::Hit::GenericHit;
use Bio::Search::HSP::GenericHSP;

use base qw(Bio::SearchIO );

use Data::Dumper;

=head2 next_result

 Title   : next_result
 Usage   : $result = stream->next_result
 Function: Reads the next ResultI object from the stream and returns it.
 Returns : A Bio::Search::Result::ResultI object
 Args    : n/a

=cut

sub next_result {
  my $self = shift;

  my $info = [];
  my $result;
  my $hit;
  my @hsp_info;
  my $previous_hit_pos;

  while ( $_ = $self->_readline ) {
      if ( $_ =~ /^>/ ) {	# looking at the start of a result
	  if ($result) {	#    and done if there's one in progress
	      $self->_pushback($_);
	      goto DONE;
	  }
	  else {		#    otherwise start a new one.
	      my ($id, $desc, $md5) = m|>([^ ]*)\s*(.*)\s*(?:md5:(.*))?|;

	      $result = Bio::Search::Result::GenericResult->new();
	      $result->algorithm('gmap');
	      $result->query_name($id);
	      $result->query_accession($id);
	      $result->query_description($desc);
	      #$self->warn("Take care of MD5!\n");

	      $hit ||= Bio::Search::Hit::GenericHit->new( -name => 
							  "NONE_SPECIFIED");
	  }
      }
      else {			# add another position to the hit/hsp
	  # 468 H	1956 C	-14:104307764 2298317517 C	H
	  # 468 	1957 A	-14:104307763 2298317516 A
	  my $c;		# info about a column
	  ($c->{query_aa_pos}, $c->{query_aa}, $c->{query_pos},
	   $c->{query_base},
           $c->{hit_strand}, $c->{hit_chromo}, $c->{hit_pos},
	   $c->{hit_concat_pos}, $c->{hit_base}, $c->{hit_aa})
	      = ($_ =~
		 m|
                   (\d+)[ ]?(.)?[\t]
                   (\d+)[ ]?(.)?[\t]
                   # TODO chromosome isn't a number... X, Y, MT....
                   (\+\|\-)([\dxXyY]+\|MT):(\d+)[ ](\d+)[ ](.)
                   [\t]?(.)?
                  |xo
		);

	  if ($previous_hit_pos &&
	      (abs($c->{hit_pos} - $previous_hit_pos) > 1)) {
	      $hit ||= Bio::Search::Hit::GenericHit->new( -name =>
							  "NONE_SPECIFIED",
							);
	      $hit->add_hsp( $self->_hsp_from_info(\@hsp_info) );
	      @hsp_info = ();
	  }
	  push @hsp_info, $c;
	  $previous_hit_pos = $c->{hit_pos};
      }
  }

 DONE:
  if ($result) {
      $hit->add_hsp( $self->_hsp_from_info(\@hsp_info) ) if (@hsp_info);

      my ($hit_length,$query_length);
      for my $hsp ($hit->hsps) {
	  $hit_length += $hsp->length();
      $query_length += $hsp->length('query');
      }
      $hit->length($hit_length);
      $hit->query_length($query_length);
      # update this now that we actually know something useful.q
      $hit->name($hsp_info[0]->{hit_chromo}); 

      $result->add_hit($hit) if ($hit);
  }

  return($result);
}



sub _hsp_from_info {
    my $self = shift;
    my $info = shift;
    my $a = {};			# args w/ which we'll create hsp
    my $hsp;
    my $identical;

    $a->{-algorithm} = 'GMAP';

    for my $c (@{$info}) {
	$a->{-query_seq} .= $c->{query_base};
	$a->{-hit_seq} .= $c->{hit_base};
    $a->{-homology_seq} .= $c->{query_base} eq $c->{hit_base} ? $c->{hit_base} : ' ';
	$identical++ if ( $c->{query_base} eq $c->{hit_base} );
    }

    $a->{-query_seq} =~ s| |\-|g;		# switch to bioperl gaps.
    $a->{-hit_seq} =~ s| |\-|g;

    $a->{-conserved} = $a->{-identical} = $identical;

    # use the coordinates from from gmap's -f 9 output to
    # determine whether gmap revcomped the query sequence
    # to generate the alignment.  Note that this is not
    # the same as the cDNA's sense/anti-sense-ness.
    $a->{-stranded} = 'both';

    $a->{-query_start} = $info->[0]->{query_pos};
    $a->{-query_end} = $info->[-1]->{query_pos};
    $a->{-query_length} = $a->{-query_end} - $a->{-query_start} + 1;

    # hit can be either strand, -f 9 output tells us which.
    # we don't have to worry about it here, but telling the generichsp code
    # that this hit is 'stranded', it compares the start and end positions
    # sets it for us.
    $a->{-hit_start} = $info->[0]->{hit_pos};
    $a->{-hit_end} = $info->[-1]->{hit_pos};

    $a->{-hit_length} = abs($a->{-hit_end} - $a->{-hit_start}) + 1;

    $a->{-hsp_length} =
	$a->{-query_length} > $a->{-hit_length} ?
	    $a->{-query_length} : $a->{-hit_length};

    $hsp = Bio::Search::HSP::GenericHSP->new( %$a );

    return $hsp;
}

# TODO (adjust regexp to swallow lines w/out md5 sig's.
sub _parse_path_header {
    my $self = shift;
    my $path_line = shift;
    my $path = {};

    (
     $path->{query},
     $path->{db},
     $path->{path_num},
     $path->{path_total_num},
     $path->{query_length},
     $path->{exon_count},
     $path->{trimmed_coverage},
     $path->{percent_identity},
     $path->{query_start},
     $path->{query_end},
     $path->{whole_genome_start},
     $path->{whole_genome_end},
     $path->{chromosome},
     $path->{chromo_start},
     $path->{chromo_end},
     $path->{strand},
     $path->{sense},
     $path->{md5},
    ) =
	($_ =~ qr|
                >
                ([^ ]*)[ ]	# the query id}, followed by a space
	        ([^ ]*)[ ]      # the genome database, followed by a space
                (\d+)/(\d+)[ ]	# path_num/path_total_num (e.g. 3/12)
	        (\d+)[ ]	# query length, followed by a space
	        (\d+)[ ]	# hsp/exon count, followed by a space
                (\d+\.\d*)[ ]	# trimmed coverage
                (\d+\.\d*)[ ]	# percent identity
                (\d+)\.\.(\d+)[ ] # query start .. query end, followed by space
                (\d+)\.\.(\d+)[ ] # whole genome s..e, followed by space
                (\d+):		# chromosome number
                (\d+)\.\.(\d+)[ ] # chromo s..e, followed by a space
                ([+-])[ ]	# strand, followed by a space
                dir:(.*) # dir:sense or dir:antisense
                [ ]md5:([\dabcdefg]+) # md5 signature
	 |x
	);

    $path->{query} or $self->throw("query was not found in path line.");
    $path->{db} or $self->throw("db was not found in path line.");
    $path->{path_num} or $self->throw("path_num was not found in path line.");
    $path->{path_total_num} or
	$self->throw("path_total_num was not found in path line.");
    $path->{query_length} or
	$self->throw("query_length was not found in path line.");
    $path->{exon_count} or
	$self->throw("exon_count was not found in path line.");
    $path->{trimmed_coverage} or
	$self->throw("trimmed_coverage was not found in path line.");
    $path->{percent_identity} or
	$self->throw("percent_identity was not found in path line.");
    $path->{query_start} or
	$self->throw("query_start was not found in path line.");
    $path->{query_end} or
	$self->throw("query_end was not found in path line.");
    $path->{whole_genome_start} or
	$self->throw("whole_genome_start was not found in path line.");
    $path->{whole_genome_end} or
	$self->throw("whole_genome_end was not found in path line.");
    $path->{chromosome} or
	$self->throw("chromosome was not found in path line.");
    $path->{chromo_start} or
	$self->throw("chromo_start was not found in path line.");
    $path->{chromo_end} or
	$self->throw("chromo_end was not found in path line.");
    $path->{strand} or $self->throw("strand was not found in path line.");
    $path->{sense} or $self->throw("sense was not found in path line.");

    return $path;
}

sub _parse_alignment_line {
    my $self = shift;
    my $a_line = shift;
    my $align = {};

    (
     $align->{chromo_start},
     $align->{chromo_end},
     $align->{query_start},
     $align->{query_end},
     $align->{percent_identity},
     $align->{align_length},
     $align->{intron_length},
    ) =
	($_ =~ qr|
                [\t]
                ([\d]+)[ ]	# start in chromosome coord.
                ([\d]+)[ ]	# end in chromosome coord.
                ([\d]+)[ ]	# start in query coord.
                ([\d]+)[ ]	# end in query coord.
                ([\d]+) 	# percent identity (as integer)
                [\t].*[\t]	# skip the edit script
                ([\d]+) 	# length of alignment block.
                [\t]*([\d]+)* 	# length of following intron.
	 |x
	);

     $align->{chromo_start}
 	or $self->throw("chromo_start missing in alignment line.");
     $align->{chromo_end},
 	or $self->throw("chromo_end was missing in alignment line.");
     $align->{query_start},
 	or $self->throw("query_start was missing in alignment line.");
     $align->{query_end},
 	or $self->throw("query_end was missing in alignment line.");
     $align->{percent_identity},
 	or $self->throw("percent_identity was missing in alignment line.");
     $align->{align_length},
 	or $self->throw("align_length was missing in alignment line.");

    return $align;
}

1;
