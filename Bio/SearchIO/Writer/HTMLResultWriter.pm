# $Id$
#
# BioPerl module for Bio::SearchIO::Writer::HTMLResultWriter
#
# Cared for by Jason Stajich <jason@bioperl.org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::SearchIO::Writer::HTMLResultWriter - Object to implement writing a Bio::Search::ResultI in HTML.

=head1 SYNOPSIS

  use Bio::SearchIO;
  use Bio::SearchIO::Writer::HTMLResultWriter;

  my $in = new Bio::SearchIO(-format => 'blast',
			     -file   => shift @ARGV);

  my $writer = new Bio::SearchIO::Writer::HTMLResultWriter();
  my $out = new Bio::SearchIO(-writer => $writer);
  $out->write_result($in->next_result);

=head1 DESCRIPTION

This object implements the SearchWriterI interface which will produce
a set of HTML for a specific Bio::Search::Report::ReportI interface.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

  bioperl-l@bioperl.org              - General discussion
  http://bioperl.org/MailList.shtml  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution. Bug reports can be submitted via
email or the web:

  bioperl-bugs@bioperl.org
  http://bugzilla.bioperl.org/

=head1 AUTHOR - Jason Stajich

Email jason@bioperl.org

Describe contact details here

=head1 CONTRIBUTORS

Additional contributors names and emails here

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::SearchIO::Writer::HTMLResultWriter;
use vars qw(@ISA %RemoteURLDefault $MaxDescLen $DATE $AlignmentLineWidth);
use strict;

# Object preamble - inherits from Bio::Root::RootI

BEGIN {
    $DATE = localtime(time);
    %RemoteURLDefault = ( 'PROTEIN' => 'http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?db=protein&cmd=search&term=%s',			  
			  'NUCLEOTIDE' => 'http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?db=nucleotide&cmd=search&term=%s'
			  );

    $MaxDescLen = 60;
    $AlignmentLineWidth = 60;
}

use Bio::Root::Root;
use Bio::SearchIO::SearchWriterI;

@ISA = qw(Bio::Root::Root Bio::SearchIO::SearchWriterI);

=head2 new

 Title   : new
 Usage   : my $obj = new Bio::SearchIO::Writer::HTMLResultWriter();
 Function: Builds a new Bio::SearchIO::Writer::HTMLResultWriter object 
 Returns : Bio::SearchIO::Writer::HTMLResultWriter
 Args    :


=cut

sub new {
  my($class,@args) = @_;

  my $self = $class->SUPER::new(@args);
  my ($p,$n,$filters) = $self->_rearrange([qw(PROTEIN_URL 
					     NUCLEOTIDE_URL 
					     FILTERS)],@args);
  $self->remote_database_url('p',$p || $RemoteURLDefault{'PROTEIN'});
  $self->remote_database_url('n',$n || $RemoteURLDefault{'NUCLEOTIDE'});

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
    return unless defined $result;
    my ($resultfilter,$hitfilter, $hspfilter) = ( $self->filter('RESULT'),
						  $self->filter('HIT'),
						  $self->filter('HSP') );
    return '' if( defined $resultfilter && ! &{$resultfilter}($result) );    

    my $type = ( $result->algorithm =~ /(P|X|Y)$/i ) ? 'PROTEIN' : 'NUCLEOTIDE';
    my $reference = $result->algorithm_reference || $self->algorithm_reference($result);
    $reference =~ s/\~/\n/g;
    my $str;
    if( $num <= 1 ) { 
	$str = sprintf(
qq{<HTML>
    <HEAD><CENTER><TITLE>Bioperl Reformatted HTML of %s Search output with Bioperl Bio::SearchIO system</TITLE></CENTER></HEAD>
    <BODY BGCOLOR="WHITE">
},$result->algorithm);

    }
    $str .= sprintf(	    
qq{
    <CENTER><H1>Bioperl Reformatted HTML of %s Search Report<br> for %s</H1></CENTER>
    <hr>
	<pre>%s</pre>
    <b>Query=</b>%s %s<br><dd>(%s letters)</dd>
    <p>
    <b>Database:</b> %s<br><dd>%d sequences; %d total letters<p></dd>
    <p>
    <table border=0>
    <tr><th>Sequences producing significant alignments:</th>
	<th>Score<br>(bits)</th><th>E<br>value</th></tr>
  }, 
		    $result->algorithm,		      
		    $result->query_name(), 
		    $reference,
		    $result->query_name, 
		    $result->query_description, $result->query_length, 
		    $result->database_name(),
		    $result->database_entries(),$result->database_letters(),
		    );
    my $hspstr = '<p><p>';
    if( $result->can('rewind')) {
        $result->rewind(); # support stream based parsing routines
    }
    
    while( my $hit = $result->next_hit ) {
	next if( $hitfilter && ! &{$hitfilter}($hit) );
	my $nm = $hit->name();
	my $id_parser = $self->id_parser;
	print STDERR "no $nm for name (",$hit->description(), "\n" unless $nm;
	my ($gi,$acc) = &$id_parser($nm);
	my $p = "%-$MaxDescLen". "s";
	my $descsub;
	if( length($hit->description) > ($MaxDescLen - 3) ) {
	    $descsub = sprintf($p,
		substr($hit->description,0,$MaxDescLen-3) . "...");
	} else { 
	    $descsub = sprintf($p,$hit->description);
	}
	
	my $url = length($self->remote_database_url($type)) > 0 ? 
	    sprintf('<a href="%s">%s</a>',
		    sprintf($self->remote_database_url($type),$gi || $acc), 
		    $hit->name()) :  $hit->name();

	my @hsps = $hit->hsps;
	
	# failover to first HSP if the data does not contain a 
	# bitscore/significance value for the Hit (NCBI XML data for one)
	
	$str .= sprintf('<tr><td>%s %s</td><td>%s</td><td><a href="#%s">%.2g</a></td></tr>'."\n",
			$url, $descsub, 
			($hit->raw_score ? $hit->raw_score : 
			(defined $hsps[0] ? $hsps[0]->score : ' ')),
			$acc,
			( $hit->significance ? $hit->significance :
			 (defined $hsps[0] ? $hsps[0]->evalue : ' ')) 
			);

	$hspstr .= "<a name=\"$acc\"><pre>\n".
	    sprintf("><b>%s</b> %s\n<dd>Length = %d</dd><p>\n\n", $hit->name, 
			defined $hit->description ? $hit->description : '', 
		    $hit->length);
	
	foreach my $hsp (@hsps ) {
	    next if( $hspfilter && ! &{$hspfilter}($hsp) );
	    $hspstr .= sprintf(" Score = %s bits (%s), Expect = %s",
			       $hsp->bits, $hsp->score, $hsp->evalue);
	    if( $hsp->pvalue ) {
		$hspstr .= ", P = ".$hsp->pvalue;
	    }
	    $hspstr .= "<br>\n";
	    $hspstr .= sprintf(" Identities = %d/%d (%d%%)",
			         ( $hsp->frac_identical('total') * 
				   $hsp->length('total')),
			       $hsp->length('total'),
			       $hsp->frac_identical('total') * 100);

	    if( $hit->algorithm =~ /(T)?BLAST(X|P)/ ||
		$hit->algorithm =~ /(T)?FAST(X|Y|P)/ ||
		$hit->algorithm eq 'TBLASTN' ||
		$hit->algorithm eq 'TFASTN' 
		) {
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
	    
	    my ($h,$q) = ( $hsp->hit->frame ,$hsp->query->frame);
		
	    if( $h || $q ) {
		$hspstr .= " Frame = ";
		
		if( $h && ! $q) {  $hspstr .= $h }
		elsif( $q && ! $h) {  $hspstr .= $q }
		else { 
		    $hspstr .= " $h / $q ";
		}
	    }
	    $hspstr .= "</pre></a><p>\n<pre>";
	    
	    my @hspvals = ( {'name' => 'Query:',
			     'seq'  => $hsp->query_string,
			     'start' => $hsp->query->strand >= 0 ? $hsp->query->start : $hsp->query->end,
			     'end'   => $hsp->query->strand >= 0 ? $hsp->query->end : $hsp->query->start,
			     'index' => 0,
			     'direction' => $hsp->query->strand || 1
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
			      'start' => $hsp->hit->strand >= 0 ? $hsp->hit->start : $hsp->hit->end,
			      'end'   => $hsp->hit->strand >= 0 ? $hsp->hit->end : $hsp->hit->start,
			      'index' => 0, 
			      'direction' => $hsp->hit->strand || 1
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
		    my $piece = substr($v->{'seq'}, $v->{'index'} + $count,
				       $AlignmentLineWidth);
		    my $cp = $piece;
		    my $plen = scalar ( $cp =~ tr/\-//);
		    my ($start,$end) = ('','');
		    if( defined $v->{'start'} ) { 
			$start = $v->{'start'};
			# since strand can be + or - use the direction
			# to signify which whether to add or substract from end
			my $d = ( ( $v->{'direction'} ) * 
				  ( $AlignmentLineWidth - $plen ));
			if( length($piece) < $AlignmentLineWidth ) {
			    $d = (length($piece) - $plen) * $v->{'direction'};
			}
			$end   = $v->{'start'} + $d - 1;


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
	}
	$hspstr .= "</pre>\n";
    }
    $str .= "</table><p>\n".$hspstr."<p><p><hr><h2>Search Parameters</h2><table border=1><tr><th>Parameter</th><th>Value</th>\n";
    
    
    foreach my $param ( $result->available_parameters ) {
	$str .= "<tr><td>$param</td><td>". $result->get_parameter($param) ."</td></tr>\n";
	
    }
    $str .= "</table><p><h2>Search Statistics</h2><table border=1><tr><th>Statistic</th><th>Value</th></tr>\n";
   foreach my $stat ( sort $result->available_statistics ) {
	$str .= "<tr><td>$stat</td><td>". $result->get_statistic($stat). "</td></th>\n";
    }
    $str .=  "</table>".$self->footer() . "</BODY>\n</HTML>\n";
    
    return $str;
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
    return "<hr><h5>Produced by Bioperl module ".ref($self)." on $DATE</h5>\n"
    
}

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
	   return $res .
"Copyright (C) 1996-2000 Washington University, Saint Louis, Missouri USA.\n
All Rights Reserved.\n
\n 
Reference:  Gish, W. (1996-2000) <a href=\"http://blast.wustl.edu\">http://blast.wustl.edu</a>\n";	   
       } else {
	   return $res . 
"Reference: Altschul, Stephen F., Thomas L. Madden, Alejandro A. Schaffer,\n
Jinghui Zhang, Zheng Zhang, Webb Miller, and David J. Lipman (1997),\n
\"Gapped BLAST and PSI-BLAST: a new generation of protein database search\n
programs\",  Nucleic Acids Res. 25:3389-3402.\n";

       }       
   } elsif( $result->algorithm =~ /FAST/i ) {
       return $result->algorithm. " ". $result->algorithm_version . "\n".
	   "\nReference: Pearson et al, Genomics (1997) 46:24-36\n";
   } else { 
       return '';
   }
}

=head2 filter

 Title   : filter
 Usage   : $writer->filter('hsp', \&hsp_filter);
 Function: Filter out either at HSP,Hit,or Result level
 Returns : none
 Args    : string => data type,
           CODE reference


=cut

1;
