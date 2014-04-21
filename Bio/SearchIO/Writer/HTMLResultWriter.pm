#
# BioPerl module for Bio::SearchIO::Writer::HTMLResultWriter
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Jason Stajich <jason@bioperl.org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# Changes 2003-07-31 (jason)
# Gary has cleaned up the code a lot to produce better looking HTML

# POD documentation - main docs before the code

=head1 NAME

Bio::SearchIO::Writer::HTMLResultWriter - write a Bio::Search::ResultI in HTML

=head1 SYNOPSIS

  use Bio::SearchIO;
  use Bio::SearchIO::Writer::HTMLResultWriter;

  my $in = Bio::SearchIO->new(-format => 'blast',
			     -file   => shift @ARGV);

  my $writer = Bio::SearchIO::Writer::HTMLResultWriter->new();
  my $out = Bio::SearchIO->new(-writer => $writer);
  $out->write_result($in->next_result);


  # to filter your output
  my $MinLength = 100; # need a variable with scope outside the method
  sub hsp_filter { 
      my $hsp = shift;
      return 1 if $hsp->length('total') > $MinLength;
  }
  sub result_filter { 
      my $result = shift;
      return $hsp->num_hits > 0;
  }

  my $writer = Bio::SearchIO::Writer::HTMLResultWriter->new
                     (-filters => { 'HSP' => \&hsp_filter} );
  my $out = Bio::SearchIO->new(-writer => $writer);
  $out->write_result($in->next_result);

  # can also set the filter via the writer object
  $writer->filter('RESULT', \&result_filter);

=head1 DESCRIPTION

This object implements the SearchWriterI interface which will produce
a set of HTML for a specific L<Bio::Search::Report::ReportI> interface.

See L<Bio::SearchIO::SearchWriterI> for more info on the filter method.

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

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Jason Stajich

Email jason-at-bioperl-dot-org

=head1 CONTRIBUTORS

Gary Williams G.Williams@hgmp.mrc.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


package Bio::SearchIO::Writer::HTMLResultWriter;
use strict;
use vars qw(%RemoteURLDefault
            $MaxDescLen $DATE $AlignmentLineWidth $Revision);

# Object preamble - inherits from Bio::Root::RootI

BEGIN {
    $Revision = '$Id$';
    $DATE = localtime(time);
    %RemoteURLDefault = ( 
      'PROTEIN' => 'http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?db=protein&cmd=search&term=%s',			  
      'NUCLEOTIDE' => 'http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?db=nucleotide&cmd=search&term=%s'
    );
    $MaxDescLen = 60;
    $AlignmentLineWidth = 60;
}


use base qw(Bio::Root::Root Bio::SearchIO::SearchWriterI);

=head2 new

 Title   : new
 Usage   : my $obj = Bio::SearchIO::Writer::HTMLResultWriter->new();
 Function: Builds a new Bio::SearchIO::Writer::HTMLResultWriter object 
 Returns : Bio::SearchIO::Writer::HTMLResultWriter
 Args    : -filters => hashref with any or all of the keys (HSP HIT RESULT)
           which have values pointing to a subroutine reference
           which will expect to get a 
           -nucleotide_url => URL sprintf string base for the nt sequences
           -protein_url => URL sprintf string base for the aa sequences
           -no_wublastlinks => boolean. Do not display WU-BLAST lines 
                               even if they are parsed out.
                               Links = (1) 

=cut

sub new {
  my($class,@args) = @_;

  my $self = $class->SUPER::new(@args);
  my ($p,$n,$filters,
      $nowublastlinks) = $self->_rearrange([qw(PROTEIN_URL 
					       NUCLEOTIDE_URL 
					       FILTERS
					       NO_WUBLASTLINKS)],@args);
  $self->remote_database_url('p',$p || $RemoteURLDefault{'PROTEIN'});
  $self->remote_database_url('n',$n || $RemoteURLDefault{'NUCLEOTIDE'});
  $self->no_wublastlinks(! $nowublastlinks);
  if( defined $filters ) {
      if( !ref($filters) =~ /HASH/i ) { 
	  $self->warn("Did not provide a hashref for the FILTERS option, ignoring.");
      } else { 
	  while( my ($type,$code) = each %{$filters} ) {
	      $self->filter($type,$code);
	  }
      }
  }

  return $self;
}

=head2 remote_database_url

 Title   : remote_database_url
 Usage   : $obj->remote_database_url($type,$newval)
 Function: This should return or set a string that contains a %s which can be
           filled in with sprintf.
 Returns : value of remote_database_url
 Args    : $type - 'PROTEIN' or 'P' for protein URLS
                   'NUCLEOTIDE' or 'N' for nucleotide URLS
           $value - new value to set [optional]


=cut

sub remote_database_url{
   my ($self,$type,$value) = @_;
   if( ! defined $type || $type !~ /^(P|N)/i ) { 
       $self->warn("Must provide a type (PROTEIN or NUCLEOTIDE)");
       return '';
   }
   $type = uc $1;
   if( defined $value) {
      $self->{'remote_database_url'}->{$type} = $value;
    }
   return $self->{'remote_database_url'}->{$type};
}

=head2 to_string

 Purpose   : Produces data for each Search::Result::ResultI in a string.
           : This is an abstract method. For some useful implementations,
           : see ResultTableWriter.pm, HitTableWriter.pm, 
           : and HSPTableWriter.pm.
 Usage     : print $writer->to_string( $result_obj, @args );
 Argument  : $result_obj = A Bio::Search::Result::ResultI object
           : @args = any additional arguments used by your implementation.
 Returns   : String containing data for each search Result or any of its
           : sub-objects (Hits and HSPs).
 Throws    : n/a

=cut

sub to_string {
    my ($self,$result,$num) = @_; 
    $num ||= 0;
    return unless defined $result;
    my $links = $self->no_wublastlinks;
    my ($resultfilter,$hitfilter, $hspfilter) = ( $self->filter('RESULT'),
						  $self->filter('HIT'),
						  $self->filter('HSP') );
    return '' if( defined $resultfilter && ! &{$resultfilter}($result) );    

    my ($qtype,$dbtype,$dbseqtype,$type);
    my $alg = $result->algorithm;
    # This is actually wrong for the FASTAs I think
    if(  $alg =~ /T(FAST|BLAST)([XY])/i ) {
	$qtype      = $dbtype = 'translated';
	$dbseqtype = $type       = 'PROTEIN';
    } elsif( $alg =~ /T(FAST|BLAST)N/i ) {
	$qtype      = '';
	$dbtype     = 'translated';
	$type       = 'PROTEIN';
	$dbseqtype  = 'NUCLEOTIDE';
    } elsif( $alg =~ /(FAST|BLAST)N/i || 
	     $alg =~ /(WABA|EXONERATE)/i ) {
	$qtype      = $dbtype = '';
	$type = $dbseqtype  = 'NUCLEOTIDE';
    } elsif( $alg =~ /(FAST|BLAST)P/  || 
	     $alg =~ /SSEARCH|HMM(PFAM|SEARCH)/i ) {
	$qtype      = $dbtype = '';
	$type = $dbseqtype  = 'PROTEIN';
    } elsif( $alg =~ /(FAST|BLAST)[XY]/i ) {
	$qtype      = 'translated';
        $dbtype     = 'PROTEIN';
	$dbseqtype  = $type      = 'PROTEIN';
    } else { 
	$self->warn("algorithm was ", $result->algorithm, " couldn't match\n");
    }
    
    
    my %baselens = ( 'Sbjct:'   => ( $dbtype eq 'translated' )  ? 3 : 1,
		     'Query:'   => ( $qtype  eq 'translated' )  ? 3 : 1);

    my $str;
    if( $num <= 1 ) { 
	$str = &{$self->start_report}($result);
    }

    $str .= &{$self->title}($result);

    $str .= $result->algorithm_reference || $self->algorithm_reference($result);
    $str .= &{$self->introduction}($result);

    $str .= "<table border=0>
            <tr><th>Sequences producing significant alignments:</th>
            <th>Score<br>(bits)</th><th>E<br>value</th></tr>";

    my $hspstr = '<p><p>';
    if( $result->can('rewind')) {
        $result->rewind(); # support stream based parsing routines
    }

    while( my $hit = $result->next_hit ) {
	next if( $hitfilter && ! &{$hitfilter}($hit) );
	my $nm = $hit->name();
	
	$self->debug( "no $nm for name (".$hit->description(). "\n") 
	    unless $nm;
	my ($gi,$acc) = &{$self->id_parser}($nm);
	my $p = "%-$MaxDescLen". "s";
	my $descsub;
	if( length($hit->description) > ($MaxDescLen - 3) ) {
	    $descsub = sprintf($p,
		substr($hit->description,0,$MaxDescLen-3) . "...");
	} else { 
	    $descsub = sprintf($p,$hit->description);
	}

	my $url_desc  = &{$self->hit_link_desc()}($self,$hit, $result);
	my $url_align = &{$self->hit_link_align()}($self,$hit, $result);

	my @hsps = $hit->hsps;
	
	if( ! @hsps ) {
	    # no HSPs so no link 
	    $str .= sprintf('<tr><td>%s %s</td><td>%s</td><td>%.2g</td></tr>'."\n",
			    $url_desc, $descsub, 
			    ($hit->bits ? $hit->bits : 
			     (defined $hsps[0] ? $hsps[0]->bits : ' ')),
			    ( $hit->significance ? $hit->significance :
			      (defined $hsps[0] ? $hsps[0]->evalue : ' ')) 
			    );
	} else { 
	    # failover to first HSP if the data does not contain a 
	    # bitscore/significance value for the Hit (NCBI XML data for one)

	    $str .= sprintf('<tr><td>%s %s</td><td>%s</td><td><a href="#%s">%.2g</a></td></tr>'."\n",
			    $url_desc, $descsub, 
			    ($hit->bits ? $hit->bits : 
			     (defined $hsps[0] ? $hsps[0]->bits : ' ')),
			    $acc,
			    ( $hit->significance ? $hit->significance :
			      (defined $hsps[0] ? $hsps[0]->evalue : ' ')) 
			    );
        my $dline = &{$self->hit_desc_line}($self, $hit, $result);
	    $hspstr .= "<a name=\"$acc\">\n".
		sprintf("><b>%s</b> %s</br><dd>Length = %s</dd><p>\n\n", $url_align, 
			$dline , &_numwithcommas($hit->length));
	    my $ct = 0;
	    foreach my $hsp (@hsps ) {
		next if( $hspfilter && ! &{$hspfilter}($hsp) );
		$hspstr .= sprintf(" Score = %s bits (%s), Expect = %s",
				   $hsp->bits || $hsp->score, 
				   $hsp->score || $hsp->bits, 
				   $hsp->evalue || '');
		if( defined $hsp->pvalue ) {
		    $hspstr .= ", P = ".$hsp->pvalue;
		}
		$hspstr .= "<br>\n";
		$hspstr .= sprintf(" Identities = %d/%d (%d%%)",
				   ( $hsp->frac_identical('total') * 
				     $hsp->length('total')),
				   $hsp->length('total'),
				   $hsp->frac_identical('total') * 100);

		if( $type eq 'PROTEIN' ) {
		    $hspstr .= sprintf(", Positives = %d/%d (%d%%)",
				       ( $hsp->frac_conserved('total') * 
					 $hsp->length('total')),
				       $hsp->length('total'),
				       $hsp->frac_conserved('total') * 100);
		}
		if( $hsp->gaps ) {
		    $hspstr .= sprintf(", Gaps = %d/%d (%d%%)",
				       $hsp->gaps('total'),
				       $hsp->length('total'),
				       (100 * $hsp->gaps('total') / 
					$hsp->length('total')));
		}

		my ($hframe,$qframe)   = ( $hsp->hit->frame, $hsp->query->frame);
		my ($hstrand,$qstrand) = ($hsp->hit->strand,$hsp->query->strand);
		# so TBLASTX will have Query/Hit frames
		#    BLASTX  will have Query frame
		#    TBLASTN will have Hit frame
		if( $hstrand || $qstrand ) {
		    $hspstr .= ", Frame = ";
		    my ($signq, $signh);
		    unless( $hstrand ) {
			$hframe = undef;
			# if strand is null or 0 then it is protein
			# and this no frame
		    } else { 
			$signh = $hstrand < 0 ? '-' : '+';
		    }
		    unless( $qstrand  ) {
			$qframe = undef;
			# if strand is null or 0 then it is protein
		    } else { 
			$signq =$qstrand < 0 ? '-' : '+';
		    }
		    # remember bioperl stores frames as 0,1,2 (GFF way)
		    # BLAST reports reports as 1,2,3 so
		    # we have to add 1 to the frame values
		    if( defined $hframe && ! defined $qframe) {  
			$hspstr .= "$signh".($hframe+1);
		    } elsif( defined $qframe && ! defined $hframe) {  
			$hspstr .= "$signq".($qframe+1);
		    } else { 
			$hspstr .= sprintf(" %s%d / %s%d",
					   $signq,$qframe+1,
					   $signh, $hframe+1);
		    }
		}
		if($links && 
		   $hsp->can('links') && defined(my $lnks = $hsp->links) ) {
		    $hspstr .= sprintf("<br>\nLinks = %s\n",$lnks);
		}

		$hspstr .= "</a><p>\n<pre>";

		my @hspvals = ( {'name' => 'Query:',
				 'seq'  => $hsp->query_string,
				 'start' => ($qstrand >= 0 ? 
					     $hsp->query->start : 
					     $hsp->query->end),
					     'end'   => ($qstrand >= 0 ? 
							 $hsp->query->end : 
							 $hsp->query->start),
							 'index' => 0,
							 'direction' => $qstrand || 1
						     },
				{ 'name' => ' 'x6,
				  'seq'  => $hsp->homology_string,
				  'start' => undef,
				  'end'   => undef,
				  'index' => 0,
				  'direction' => 1
				  },
				{ 'name'  => 'Sbjct:',
				  'seq'   => $hsp->hit_string,
				  'start' => ($hstrand >= 0 ? 
					      $hsp->hit->start : 
					      $hsp->hit->end),
					      'end'   => ($hstrand >= 0 ? 
							  $hsp->hit->end : 
							  $hsp->hit->start),
							  'index' => 0, 
							  'direction' => $hstrand || 1
						      }
				);	    


		# let's set the expected length (in chars) of the starting number
		# in an alignment block so we can have things line up
		# Just going to try and set to the largest

		my ($numwidth) = sort { $b <=> $a }(length($hspvals[0]->{'start'}),
						    length($hspvals[0]->{'end'}),
						    length($hspvals[2]->{'start'}),
						    length($hspvals[2]->{'end'}));
		my $count = 0;
		while ( $count < $hsp->length('total') ) {
		    foreach my $v ( @hspvals ) {
			my $piece = substr($v->{'seq'}, $v->{'index'} + $count,
					   $AlignmentLineWidth);
			my $cp = $piece;
			my $plen = scalar ( $cp =~ tr/\-//);
			my ($start,$end) = ('','');
			if( defined $v->{'start'} ) { 
			    $start = $v->{'start'};
			    # since strand can be + or - use the direction
			    # to signify which whether to add or substract from end
			    my $d = $v->{'direction'} * ( $AlignmentLineWidth - $plen )*
				$baselens{$v->{'name'}};
			    if( length($piece) < $AlignmentLineWidth ) {
				$d = (length($piece) - $plen) * $v->{'direction'} * 
				    $baselens{$v->{'name'}};
			    }
			    $end   = $v->{'start'} + $d - $v->{'direction'};
			    $v->{'start'} += $d;
			}
			$hspstr .= sprintf("%s %-".$numwidth."s %s %s\n",
					   $v->{'name'},
					   $start,
					   $piece,
					   $end
					   );
		    }
		    $count += $AlignmentLineWidth;
		    $hspstr .= "\n\n";
		}
		$hspstr .= "</pre>\n";
	    }
	}
#	$hspstr .= "</pre>\n";
    }

    $str .= "</table><p>\n".$hspstr;
    my ($pav, $sav) = ($result->available_parameters, $result->available_statistics);
    if ($pav || $sav) {
        # make table of search statistics and end the web page
        $str .= "<p><p><hr><h2>Search Parameters</h2>";
        if ($pav) {
        $str .= "<table border=1><tr><th>Parameter</th><th>Value</th>\n";
        foreach my $param ( sort $result->available_parameters ) {
            $str .= "<tr><td>$param</td><td>". $result->get_parameter($param) ."</td></tr>\n";
        }
        $str .= "</table>";
        }
        
        if ($sav) {
        $str .= "<p><h2>Search Statistics</h2><table border=1><tr><th>Statistic</th><th>Value</th></tr>\n";
        foreach my $stat ( sort $result->available_statistics ) {
            $str .= "<tr><td>$stat</td><td>". $result->get_statistic($stat). "</td>\n";
        }
        $str .=  "</tr></table>";
        }
    }
    $str .= $self->footer() . "<P>\n";
    return $str;
}

=head2 hit_link_desc

 Title   : hit_link_desc
 Usage   : $self->hit_link_desc(\&link_function);
 Function: Get/Set the function which provides an HTML 
           link(s) for the given hit to be used
           within the description section at the top of the BLAST report.
           This allows a person reading the report within
           a web browser to go to one or more database entries for
           the given hit from the description section.
 Returns : Function reference
 Args    : Function reference
 See Also: L<default_hit_link_desc()>

=cut

sub hit_link_desc{
    my( $self, $code ) = @_; 
    if ($code) {
        $self->{'_hit_link_desc'} = $code;
    }
    return $self->{'_hit_link_desc'} || \&default_hit_link_desc;
}

=head2 default_hit_link_desc

 Title   : default_hit_link_desc
 Usage   : $self->default_hit_link_desc($hit, $result)
 Function: Provides an HTML link(s) for the given hit to be used
           within the description section at the top of the BLAST report.
           This allows a person reading the report within
           a web browser to go to one or more database entries for
           the given hit from the description section.
 Returns : string containing HTML markup "<a href...")

           The default implementation returns an HTML link to the
           URL supplied by the remote_database_url() method
           and using the identifier supplied by the id_parser() method.
           It will use the NCBI GI if present, and the accession if not.

 Args    : First argument is a Bio::Search::Hit::HitI
           Second argument is a Bio::Search::Result::ResultI

See Also: L<hit_link_align>, L<remote_database>, L<id_parser>

=cut

sub default_hit_link_desc {
    my($self, $hit, $result) = @_;
    my $type = ( $result->algorithm =~ /(P|X|Y)$/i ) ? 'PROTEIN' : 'NUCLEOTIDE';
    my ($gi,$acc) = &{$self->id_parser}($hit->name);

    my $url = length($self->remote_database_url($type)) > 0 ? 
              sprintf('<a href="%s">%s</a>',
                      sprintf($self->remote_database_url($type),$gi || $acc), 
                      $hit->name()) :  $hit->name();

    return $url;
}


=head2 hit_link_align

 Title   : hit_link_align
 Usage   : $self->hit_link_align(\&link_function);
 Function: Get/Set the function which provides an HTML link(s) 
           for the given hit to be used
           within the HSP alignment section of the BLAST report.
           This allows a person reading the report within
           a web browser to go to one or more database entries for
           the given hit from the alignment section.
 Returns : string containing HTML markup "<a href...")

           The default implementation delegates to hit_link_desc().

 Args    : First argument is a Bio::Search::Hit::HitI
           Second argument is a Bio::Search::Result::ResultI

See Also: L<hit_link_desc>, L<remote_database>, L<id_parser>

=cut

sub hit_link_align {
    my ($self,$code) = @_;
    if ($code) {
        $self->{'_hit_link_align'} = $code;
    }
    return $self->{'_hit_link_align'} || \&default_hit_link_desc;
}

=head2 hit_desc_line

 Title   : hit_desc_line
 Usage   : $self->hit_desc_line(\&link_function);
 Function: Get/Set the function which provides HTML for the description
           information from a hit. This allows one to parse
           the rest of the description and split up lines, add links, etc.
 Returns : Function reference
 Args    : Function reference
 See Also: L<default_hit_link_desc()>

=cut

sub hit_desc_line{
    my( $self, $code ) = @_; 
    if ($code) {
        $self->{'_hit_desc_line'} = $code;
    }
    return $self->{'_hit_desc_line'} || \&default_hit_desc_line;
}

=head2 default_hit_desc_line

 Title   : default_hit_desc_line
 Usage   : $self->default_hit_desc_line($hit, $result)
 Function: Parses the description line information, splits based on the
           hidden \x01 between independent descriptions, checks the lines for
           possible web links, and adds HTML link(s) for the given hit to be
           used.

 Returns : string containing HTML markup "<a href...")
           The default implementation returns an HTML link to the
           URL supplied by the remote_database_url() method
           and using the identifier supplied by the id_parser() method.
           It will use the NCBI GI if present, and the accession if not.

 Args    : First argument is a Bio::Search::Hit::HitI
           Second argument is a Bio::Search::Result::ResultI

See Also: L<hit_link_align>, L<remote_database>, L<id_parser>

=cut

sub default_hit_desc_line {
    my($self, $hit, $result) = @_;
    my $type = ( $result->algorithm =~ /(P|X|Y)$/i ) ? 'PROTEIN' : 'NUCLEOTIDE';
    my @descs = split /\x01/, $hit->description;
    #my $descline = join("</br>",@descs)."</br>";
    my $descline = '';
    #return $descline;
    for my $sec (@descs) {
        my $url = '';
        if ($sec =~ s/((?:gi\|(\d+)\|)?        # optional GI
                     (\w+)\|([A-Z\d\.\_]+) # main 
                     (\|[A-Z\d\_]+)?) # optional secondary ID//xms) {
            my ($name, $gi, $db, $acc) = ($1, $2, $3, $4);
            #$acc ||= ($rest) ? $rest : $gi;
            $acc =~ s/^\s+(\S+)/$1/;
            $acc =~ s/(\S+)\s+$/$1/;
            $url =
            length($self->remote_database_url($type)) > 0 ? 
              sprintf('<a href="%s">%s</a> %s',
                      sprintf($self->remote_database_url($type),
                      $gi || $acc || $db), 
                      $name, $sec) :  $sec;
        } else {
            $url = $sec;
        }
        $descline .= "$url</br>\n";
    }
    return $descline;
}

=head2 start_report

  Title   : start_report
  Usage   : $index->start_report( CODE )
  Function: Stores or returns the code to
            write the start of the <HTML> block, the <TITLE> block
            and the start of the <BODY> block of HTML.   Useful
            for (for instance) specifying alternative
            HTML if you are embedding the output in
            an HTML page which you have already started.
            (For example a routine returning a null string).
            Returns \&default_start_report (see below) if not
            set. 
  Example : $index->start_report( \&my_start_report )
  Returns : ref to CODE if called without arguments
  Args    : CODE

=cut

sub start_report {
    my( $self, $code ) = @_; 
    if ($code) {
        $self->{'_start_report'} = $code;
    }
    return $self->{'_start_report'} || \&default_start_report;
}

=head2 default_start_report

 Title   : default_start_report
 Usage   : $self->default_start_report($result)
 Function: The default method to call when starting a report.
 Returns : sting
 Args    : First argument is a Bio::Search::Result::ResultI

=cut

sub default_start_report {
    my ($result) = @_;
    return sprintf(
    qq{<HTML>
      <HEAD> <CENTER><TITLE>Bioperl Reformatted HTML of %s output with Bioperl Bio::SearchIO system</TITLE></CENTER></HEAD>
      <!------------------------------------------------------------------->
      <!-- Generated by Bio::SearchIO::Writer::HTMLResultWriter          -->
      <!-- %s -->
      <!-- http://bioperl.org                                            -->
      <!------------------------------------------------------------------->
      <BODY BGCOLOR="WHITE">
    },$result->algorithm,$Revision);
    
}

=head2 title

 Title   : title
 Usage   : $self->title($CODE)

  Function: Stores or returns the code to provide HTML for the given
            BLAST report that will appear at the top of the BLAST report
            HTML output.  Useful for (for instance) specifying
            alternative routines to write your own titles.
            Returns \&default_title (see below) if not
            set. 
  Example : $index->title( \&my_title )
  Returns : ref to CODE if called without arguments
  Args    : CODE

=cut

sub title {
    my( $self, $code ) = @_; 
    if ($code) {
        $self->{'_title'} = $code;
    }
    return $self->{'_title'} || \&default_title;
}

=head2 default_title

 Title   : default_title
 Usage   : $self->default_title($result)
 Function: Provides HTML for the given BLAST report that will appear
           at the top of the BLAST report HTML output.
 Returns : string containing HTML markup
           The default implementation returns <CENTER> <H1> HTML
           containing text such as:
           "Bioperl Reformatted HTML of BLASTP Search Report
                     for gi|1786183|gb|AAC73113.1|"
 Args    : First argument is a Bio::Search::Result::ResultI

=cut

sub default_title {
    my ($result) = @_;

    return sprintf(
        qq{<CENTER><H1><a href="http://bioperl.org">Bioperl</a> Reformatted HTML of %s Search Report<br> for %s</H1></CENTER>},
		    $result->algorithm,
		    $result->query_name());
}


=head2 introduction

 Title   : introduction
 Usage   : $self->introduction($CODE)

  Function: Stores or returns the code to provide HTML for the given
            BLAST report detailing the query and the
            database information.
            Useful for (for instance) specifying
            routines returning alternative introductions.
            Returns \&default_introduction (see below) if not
            set. 
  Example : $index->introduction( \&my_introduction )
  Returns : ref to CODE if called without arguments
  Args    : CODE

=cut

sub introduction {
    my( $self, $code ) = @_; 
    if ($code) {
        $self->{'_introduction'} = $code;
    }
    return $self->{'_introduction'} || \&default_introduction;
}

=head2 default_introduction

 Title   : default_introduction
 Usage   : $self->default_introduction($result)
 Function: Outputs HTML to provide the query
           and the database information
 Returns : string containing HTML
 Args    : First argument is a Bio::Search::Result::ResultI
           Second argument is string holding literature citation

=cut

sub default_introduction {
    my ($result) = @_;

    return sprintf(
    qq{
    <b>Query=</b> %s %s<br><dd>(%s letters)</dd>
    <p>
    <b>Database:</b> %s<br><dd>%s sequences; %s total letters<p></dd>
    <p>
  }, 
		   $result->query_name, 
		   $result->query_description, 
		   &_numwithcommas($result->query_length), 
		   $result->database_name(),
		   &_numwithcommas($result->database_entries()), 
		   &_numwithcommas($result->database_letters()),
		   );
}

=head2 end_report

 Title   : end_report
 Usage   : $self->end_report()
 Function: The method to call when ending a report, this is
           mostly for cleanup for formats which require you to 
           have something at the end of the document (</BODY></HTML>)
           for HTML
 Returns : string
 Args    : none

=cut

sub end_report {
    return "</BODY>\n</HTML>\n";
}

# copied from Bio::Index::Fasta
# useful here as well

=head2 id_parser

  Title   : id_parser
  Usage   : $index->id_parser( CODE )
  Function: Stores or returns the code used by record_id to
            parse the ID for record from a string.  Useful
            for (for instance) specifying a different
            parser for different flavours of FASTA file. 
            Returns \&default_id_parser (see below) if not
            set. If you supply your own id_parser
            subroutine, then it should expect a fasta
            description line.  An entry will be added to
            the index for each string in the list returned.
  Example : $index->id_parser( \&my_id_parser )
  Returns : ref to CODE if called without arguments
  Args    : CODE

=cut

sub id_parser {
    my( $self, $code ) = @_;
    
    if ($code) {
        $self->{'_id_parser'} = $code;
    }
    return $self->{'_id_parser'} || \&default_id_parser;
}



=head2 default_id_parser

  Title   : default_id_parser
  Usage   : $id = default_id_parser( $header )
  Function: The default Fasta ID parser for Fasta.pm
            Returns $1 from applying the regexp /^>\s*(\S+)/
            to $header.
  Returns : ID string
            The default implementation checks for NCBI-style
            identifiers in the given string ('gi|12345|AA54321').
            For these IDs, it extracts the GI and accession and
            returns a two-element list of strings (GI, acc).
  Args    : a fasta header line string

=cut

sub default_id_parser {    
    my ($string) = @_;
    my ($gi,$acc);
    if( $string =~ s/gi\|(\d+)\|?// ) 
    { $gi = $1; $acc = $1;}
    
    if( $string =~ /(\w+)\|([A-Z\d\.\_]+)(\|[A-Z\d\_]+)?/ ) {
	$acc = defined $2 ? $2 : $1;
    } else {
        $acc = $string;
	$acc =~ s/^\s+(\S+)/$1/;
	$acc =~ s/(\S+)\s+$/$1/;	
    } 
    return ($gi,$acc);
}
	
sub MIN { $a <=> $b ? $a : $b; }
sub MAX { $a <=> $b ? $b : $a; }

sub footer { 
    my ($self) = @_;
    return "<hr><h5>Produced by Bioperl module ".ref($self)." on $DATE<br>Revision: $Revision</h5>\n"
    
}

=head2 algorithm_reference

 Title   : algorithm_reference
 Usage   : my $reference = $writer->algorithm_reference($result);
 Function: Returns the appropriate Bibliographic reference for the 
           algorithm format being produced
 Returns : String
 Args    : L<Bio::Search::Result::ResultI> to reference


=cut

sub algorithm_reference {
   my ($self,$result) = @_;
   return '' if( ! defined $result || !ref($result) ||
		 ! $result->isa('Bio::Search::Result::ResultI')) ;   
   if( $result->algorithm =~ /BLAST/i ) {
       my $res = $result->algorithm . ' ' . $result->algorithm_version . "<p>";
       if( $result->algorithm_version =~ /WashU/i ) {
	   return $res .
"Copyright (C) 1996-2000 Washington University, Saint Louis, Missouri USA.<br>
All Rights Reserved.<p>
<b>Reference:</b>  Gish, W. (1996-2000) <a href=\"http://blast.wustl.edu\">http://blast.wustl.edu</a><p>";	   
       } else {
	   return $res . 
"<b>Reference:</b> Altschul, Stephen F., Thomas L. Madden, Alejandro A. Schaffer,<br>
Jinghui Zhang, Zheng Zhang, Webb Miller, and David J. Lipman (1997),<br>
\"Gapped BLAST and PSI-BLAST: a new generation of protein database search<br>
programs\",  Nucleic Acids Res. 25:3389-3402.<p>";

       }       
   } elsif( $result->algorithm =~ /FAST/i ) {
       return $result->algorithm . " " . $result->algorithm_version . "<br>" .
	   "\n<b>Reference:</b> Pearson et al, Genomics (1997) 46:24-36<p>";
   } else { 
       return '';
   }
}

# from Perl Cookbook 2.17
sub _numwithcommas {
    my $num = reverse( $_[0] );
    $num =~ s/(\d{3})(?=\d)(?!\d*\.)/$1,/g;
    return scalar reverse $num;
}

=head2 Methods Bio::SearchIO::SearchWriterI

L<Bio::SearchIO::SearchWriterI> inherited methods.

=head2 filter

 Title   : filter
 Usage   : $writer->filter('hsp', \&hsp_filter);
 Function: Filter out either at HSP,Hit,or Result level
 Returns : none
 Args    : string => data type,
           CODE reference


=cut


=head2 no_wublastlinks

 Title   : no_wublastlinks
 Usage   : $obj->no_wublastlinks($newval)
 Function: Get/Set boolean value regarding whether or not to display
           Link = (1) 
           type output in the report output (WU-BLAST only)
 Returns : boolean
 Args    : on set, new boolean value (a scalar or undef, optional)


=cut

sub no_wublastlinks{
    my $self = shift;

    return $self->{'no_wublastlinks'} = shift if @_;
    return $self->{'no_wublastlinks'};
}

1;
