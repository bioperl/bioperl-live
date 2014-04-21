#
# BioPerl module for Bio::SearchIO::Writer::TextResultWriter
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Jason Stajich <jason@bioperl.org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::SearchIO::Writer::TextResultWriter - Object to implement writing
a Bio::Search::ResultI in Text.

=head1 SYNOPSIS

  use Bio::SearchIO;
  use Bio::SearchIO::Writer::TextResultWriter;

  my $in = Bio::SearchIO->new(-format => 'blast',
			     -file   => shift @ARGV);

  my $writer = Bio::SearchIO::Writer::TextResultWriter->new();
  my $out = Bio::SearchIO->new(-writer => $writer);
  $out->write_result($in->next_result);

=head1 DESCRIPTION

This object implements the SearchWriterI interface which will produce
a set of Text for a specific Bio::Search::Report::ReportI interface.

You can also provide the argument -filters =E<gt> \%hash to filter the at
the hsp, hit, or result level.  %hash is an associative array which
contains any or all of the keys (HSP, HIT, RESULT).  The values
pointed to by these keys would be references to a subroutine which
expects to be passed an object - one of Bio::Search::HSP::HSPI,
Bio::Search::Hit::HitI, and Bio::Search::Result::ResultI respectively.
Each function needs to return a boolean value as to whether or not the
passed element should be included in the output report - true if it is
to be included, false if it to be omitted.

For example to filter on sequences in the database which are too short
for your criteria you would do the following.

Define a hit filter method 

  sub hit_filter { 
      my $hit = shift;
      return $hit->length E<gt> 100; # test if length of the hit sequence
                                     # long enough    
  }
  my $writer = Bio::SearchIO::Writer::TextResultWriter->new(
       -filters => { 'HIT' =E<gt> \&hit_filter }  
      );

Another example would be to filter HSPs on percent identity, let's
only include HSPs which are 75% identical or better.

   sub hsp_filter {
       my $hsp = shift;
       return $hsp->percent_identity E<gt> 75;
   }
   my $writer = Bio::SearchIO::Writer::TextResultWriter->new(
       -filters => { 'HSP' =E<gt> \&hsp_filter }  
      );

See L<Bio::SearchIO::SearchWriterI> for more info on the filter method.


This module will use the module Text::Wrap if it is installed to wrap
the Query description line.  If you do not have Text::Wrap installed
this module will work fine but you won't have the Query line wrapped.
You will see a warning about this when you first instantiate a
TextResultWriter - to avoid these warnings from showing up, simply set
the verbosity upon initialization to -1 like this: my $writer = new
Bio::SearchIO::Writer::TextResultWriter(-verbose =E<gt> -1);

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

Email jason@bioperl.org

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::SearchIO::Writer::TextResultWriter;
use vars qw($MaxNameLen $MaxDescLen $AlignmentLineWidth 	    $DescLineLen $TextWrapLoaded);
use strict;

# Object preamble - inherits from Bio::Root::RootI

BEGIN {
    $MaxDescLen = 65;
    $AlignmentLineWidth = 60;    
    eval { require Text::Wrap; $TextWrapLoaded = 1;};
    if( $@ ) {
	$TextWrapLoaded = 0;
    }
}

use POSIX;

use base qw(Bio::Root::Root Bio::SearchIO::SearchWriterI);

=head2 new

 Title   : new
 Usage   : my $obj = Bio::SearchIO::Writer::TextResultWriter->new();
 Function: Builds a new Bio::SearchIO::Writer::TextResultWriter object 
 Returns : Bio::SearchIO::Writer::TextResultWriter
 Args    : -filters => hashref with any or all of the keys (HSP HIT RESULT)
           which have values pointing to a subroutine reference
           which will expect to get a Hit,HSP, Result object respectively
           -no_wublastlinks => boolean. Do not display WU-BLAST lines even if 
                               they are parsed out
                               Links = (1) 

=cut

sub new {
  my($class,@args) = @_;

  my $self = $class->SUPER::new(@args);
  my ($filters,$nowublastlinks) = $self->_rearrange([qw(FILTERS 
							NO_WUBLASTLINKS)],
						    @args);
  if( defined $filters ) {
      if( !ref($filters) =~ /HASH/i ) { 
	  $self->warn("Did not provide a hashref for the FILTERS option, ignoring.");
      } else { 
	  while( my ($type,$code) = each %{$filters} ) {
	      $self->filter($type,$code);
	  }
      }
  }
  $self->no_wublastlinks(! $nowublastlinks);
  unless( $TextWrapLoaded ) {
      $self->warn("Could not load Text::Wrap - the Query Description will not be line wrapped\n");
  } else { 
      $Text::Wrap::columns =  $MaxDescLen;
  }
  return $self;
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
    
    my $wublast = ($result->algorithm_version =~ /WashU/) ? 1 : 0;
    
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
	     $alg =~ /SSEARCH|(HMM|SEARCH|PFAM)/i ) {
	$qtype      = $dbtype = '';
	$type = $dbseqtype  = 'PROTEIN';
    } elsif( $alg =~ /(FAST|BLAST)[XY]/i ) {
	$qtype      = 'translated';
        $dbtype     = 'PROTEIN';
	$dbseqtype  = $type      = 'PROTEIN';
    } else { 
	print STDERR "algorithm was ", $result->algorithm, " couldn't match\n";
    }
    
    
    my %baselens = ( 'Sbjct:'   => ( $dbtype eq 'translated' )  ? 3 : 1,
		     'Query:'   => ( $qtype  eq 'translated' )  ? 3 : 1);

    my $str;
    if( ! defined $num || $num <= 1 ) { 
	$str = &{$self->start_report}($result);
    }

    $str .= &{$self->title}($result);
    $str .= $result->algorithm . " " . $result->algorithm_version . "\n\n\n";
    $str .= $result->algorithm_reference || $self->algorithm_reference($result);
    $str .= &{$self->introduction}($result);


    $str .= qq{
                                                                 Score       E
Sequences producing significant alignments:                      (bits)    value
};
    my $hspstr = '';
    if( $result->can('rewind')) {
        $result->rewind(); # support stream based parsing routines
    }
    while( my $hit = $result->next_hit ) {
	next if( defined $hitfilter && ! &{$hitfilter}($hit) );
	my $nm = $hit->name();
	$self->debug( "no $nm for name (".$hit->description(). "\n") 
	    unless $nm;
	my ($gi,$acc) = &{$self->id_parser}($nm);
	my $p = "%-$MaxDescLen". "s";
	my $descsub;
	my $desc = sprintf("%s %s",$nm,$hit->description);
	if( length($desc) - 3 > $MaxDescLen) {
	    $descsub = sprintf($p,
			       substr($desc,0,$MaxDescLen-3) . 
			       "...");
	} else { 
	    $descsub = sprintf($p,$desc);
	}
	$str .= $wublast ? sprintf("%s   %-4s  %s\n",
			$descsub,
			defined $hit->raw_score ? $hit->raw_score : ' ',
			defined $hit->significance ? $hit->significance : '?') :
                    sprintf("%s   %-4s  %s\n",
			$descsub,
			defined $hit->bits ? $hit->bits: ' ',
			defined $hit->significance ? $hit->significance : '?');        
	my @hsps = $hit->hsps;
	if( @hsps ) { 
	    $hspstr .= sprintf(">%s %s\n%9sLength = %d\n\n",
			       $hit->name, 
			       defined $hit->description ? $hit->description : '', 
			       '', # empty is for the %9s in the str formatting 
			       $hit->length);

	    foreach my $hsp ( @hsps ) { 
		next if( defined $hspfilter && ! &{$hspfilter}($hsp) );
		$hspstr .= sprintf(" Score = %4s bits (%s), Expect = %s",
				   $hsp->bits, $hsp->score, $hsp->evalue);
		if( $hsp->pvalue ) {
		    $hspstr .= ", P = ".$hsp->pvalue;
		}
		$hspstr .= "\n";
		$hspstr .= sprintf(" Identities = %d/%d (%d%%)",
				   ( $hsp->frac_identical('total') * 
				     $hsp->length('total')),
				   $hsp->length('total'),
				   POSIX::floor($hsp->frac_identical('total') 
						* 100));
		
		if( $type eq 'PROTEIN' ) {
		    $hspstr .= sprintf(", Positives = %d/%d (%d%%)",
				       ( $hsp->frac_conserved('total') * 
					 $hsp->length('total')),
				       $hsp->length('total'),
				       POSIX::floor($hsp->frac_conserved('total') * 100));

		}
		if( $hsp->gaps ) {
		    $hspstr .= sprintf(", Gaps = %d/%d (%d%%)",
				       $hsp->gaps('total'),
				       $hsp->length('total'),
				       POSIX::floor(100 * $hsp->gaps('total') / 
						    $hsp->length('total')));
		}
		$hspstr .= "\n";
		my ($hframe,$qframe)   = ( $hsp->hit->frame, 
					   $hsp->query->frame);
		my ($hstrand,$qstrand) = ($hsp->hit->strand,$hsp->query->strand);
		# so TBLASTX will have Query/Hit frames
		#    BLASTX  will have Query frame
		#    TBLASTN will have Hit frame
		if( $hstrand || $qstrand ) {
		    $hspstr .= " Frame = ";
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
		
		if( $links && 
		    $hsp->can('links') && defined(my $lnks = $hsp->links) ) {
		    $hspstr .= sprintf(" Links = %s\n",$lnks);
		}
		$hspstr .= "\n\n";

		my @hspvals = ( {'name'  => 'Query:',
				 'seq'   => $hsp->query_string,
				 'start' => ( $qstrand >= 0 ? 
					      $hsp->query->start : 
					      $hsp->query->end),
					      'end'   => ($qstrand >= 0 ? 
							  $hsp->query->end : 
							  $hsp->query->start),
							  'index' => 0,
							  'direction' => $qstrand || 1
						      },
				{ 'name' => ' 'x6, # this might need to adjust for long coordinates??
				  'seq'  => $hsp->homology_string,
				  'start' => undef,
				  'end'   => undef,
				  'index' => 0,
				  'direction' => 1
				  },
				{ 'name'  => 'Sbjct:',
				  'seq'   => $hsp->hit_string,
				  'start' => ($hstrand >= 0 ? 
					      $hsp->hit->start : $hsp->hit->end),
				      'end'   => ($hstrand >= 0 ? 
						  $hsp->hit->end : $hsp->hit->start),
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
		while ( $count <= $hsp->length('total') ) {
		    foreach my $v ( @hspvals ) {
			my $piece = substr($v->{'seq'}, $v->{'index'} +$count,
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
		    $hspstr .= "\n";
		}
	    }
	    $hspstr .= "\n";
	}
    }
    $str .= "\n\n".$hspstr;
    
    $str .= sprintf(qq{  Database: %s
    Posted date:  %s
	Number of letters in database: %s
	Number of sequences in database: %s

  Matrix: %s
  },
	        $result->database_name(),
	        $result->get_statistic('posted_date') ||
	        POSIX::strftime("%b %d, %Y %I:%M %p",localtime),
	        &_numwithcommas($result->database_letters()),
	        &_numwithcommas($result->database_entries()),
	        $result->get_parameter('matrix') || '');

    if( defined (my $open = $result->get_parameter('gapopen')) ) {
	$str .= sprintf("Gap Penalties Existence: %d, Extension: %d\n",
			$open || 0, $result->get_parameter('gapext') || 0);
    }

    # skip those params we've already output
    foreach my $param ( grep { ! /matrix|gapopen|gapext/i } 
			$result->available_parameters ) {
	$str .= "$param: ". $result->get_parameter($param) ."\n";
	
    }
    $str .= "Search Statistics\n";
    # skip posted date, we already output it
   foreach my $stat ( sort grep { ! /posted_date/ } 
		      $result->available_statistics ) {
       my $expect = $result->get_parameter('expect');
       my $v = $result->get_statistic($stat);
       if( $v =~ /^\d+$/ ) {
	   $v = &_numwithcommas($v);
       }
       if( defined $expect && 
	   $stat eq 'seqs_better_than_cutoff' ) {
	   $str .= "seqs_better_than_$expect: $v\n";
       } else { 
	   my $v = 
	   $str .= "$stat: $v\n";
       }
    }
    $str .=  "\n\n";
    return $str;
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
    return "";    
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
           at the top of the BLAST report output.
 Returns : empty for text implementation
 Args    : First argument is a Bio::Search::Result::ResultI

=cut

sub default_title {
    my ($result) = @_;
    return "";
# The HTML implementation
#    return sprintf(
#        qq{<CENTER><H1><a href="http://bioperl.org">Bioperl</a> Reformatted HTML of %s Search Report<br> for %s</H1></CENTER>},
#		    $result->algorithm,
#		    $result->query_name());
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
Query= %s
       (%s letters)

Database: %s
           %s sequences; %s total letters
}, 
		   &_linewrap($result->query_name . " " . 
			      $result->query_description), 
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
    return "";
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


=head2 algorithm_reference

 Title   : algorithm_reference
 Usage   : my $reference = $writer->algorithm_reference($result);
 Function: Returns the appropriate Bibliographic reference for the 
           algorithm format being produced
 Returns : String
 Args    : L<Bio::Search::Result::ResultI> to reference


=cut

sub algorithm_reference{
   my ($self,$result) = @_;
   return '' if( ! defined $result || !ref($result) ||
		 ! $result->isa('Bio::Search::Result::ResultI')) ;   
   if( $result->algorithm =~ /BLAST/i ) {
       my $res = $result->algorithm . ' '. $result->algorithm_version. "\n";
       if( $result->algorithm_version =~ /WashU/i ) {
	   return $res .qq{
Copyright (C) 1996-2000 Washington University, Saint Louis, Missouri USA.
All Rights Reserved.
 
Reference:  Gish, W. (1996-2000) http://blast.wustl.edu
};	   
       } else {
	   return $res . qq{
Reference: Altschul, Stephen F., Thomas L. Madden, Alejandro A. Schaffer,
Jinghui Zhang, Zheng Zhang, Webb Miller, and David J. Lipman (1997),
"Gapped BLAST and PSI-BLAST: a new generation of protein database search
programs",  Nucleic Acids Res. 25:3389-3402.
};
       }       
   } elsif( $result->algorithm =~ /FAST/i ) {
       return $result->algorithm. " ". $result->algorithm_version . "\n".
	   "\nReference: Pearson et al, Genomics (1997) 46:24-36\n";
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

sub _linewrap {
    my ($str) = @_;
    if($TextWrapLoaded) {
	return Text::Wrap::wrap("","",$str); # use Text::Wrap
    } else { return $str; }     # cannot wrap
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
