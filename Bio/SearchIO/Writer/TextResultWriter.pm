# $Id$
#
# BioPerl module for Bio::SearchIO::Writer::TextResultWriter
#
# Cared for by Jason Stajich <jason@bioperl.org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::SearchIO::Writer::TextResultWriter - Object to implement writing a Bio::Search::ResultI in Text.

=head1 SYNOPSIS

  use Bio::SearchIO;
  use Bio::SearchIO::Writer::TextResultWriter;

  my $in = new Bio::SearchIO(-format => 'blast',
			     -file   => shift @ARGV);

  my $writer = new Bio::SearchIO::Writer::TextResultWriter();
  my $out = new Bio::SearchIO(-writer => $writer);
  $out->write_result($in->next_result);

=head1 DESCRIPTION

This object implements the SearchWriterI interface which will produce
a set of Text for a specific Bio::Search::Report::ReportI interface.

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


package Bio::SearchIO::Writer::TextResultWriter;
use vars qw(@ISA $MaxNameLen $MaxDescLen $AlignmentLineWidth);
use strict;

# Object preamble - inherits from Bio::Root::RootI

BEGIN {
    $MaxDescLen = 65;
    $AlignmentLineWidth = 60;
}

use Bio::Root::Root;
use Bio::SearchIO::SearchWriterI;
use POSIX;

@ISA = qw(Bio::Root::Root Bio::SearchIO::SearchWriterI);

=head2 new

 Title   : new
 Usage   : my $obj = new Bio::SearchIO::Writer::TextResultWriter();
 Function: Builds a new Bio::SearchIO::Writer::TextResultWriter object 
 Returns : Bio::SearchIO::Writer::TextResultWriter
 Args    :


=cut

sub new {
  my($class,@args) = @_;

  my $self = $class->SUPER::new(@args);
  my ($filters) = $self->_rearrange([qw(FILTERS)],@args);
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
    my ($self,$result) = @_; 
    return '' unless defined $result;
    my ($resultfilter,$hitfilter,
	$hspfilter) = ( $self->filter('RESULT'),
			$self->filter('HIT'),
			$self->filter('HSP'));
    return '' if( defined $resultfilter && ! &{$resultfilter}($result) );
    
    my $type = ( $result->algorithm =~ /(P|X|Y)$/i ) ? 'PROTEIN' : 'NUCLEOTIDE';
    my $str = sprintf(
qq{%s %s 

Bioperl Reformatted BLAST Format report of %s Search output 
with Bio::SearchIO Bioperl Reformatted BLAST Format Module %s

%s

Query=%s %s 
         (%d letters)

Database: %s 
          %d sequences; %d total letters

Sequences producing significant alignments:         Score       E
                                                    (bits)    value
},
		      $result->algorithm,  $result->algorithm_version,
		      $result->algorithm, ref($result),
		      $result->program_reference() || $self->algorithm_reference($result),
		      $result->query_name, $result->query_description, 
		      $result->query_length, 
		      $result->database_name(),
		      $result->database_entries(),
		      $result->database_letters(),
		      );
    my $hspstr = '';
    $result->rewind();
    while( my $hit = $result->next_hit ) {
	next if( defined $hitfilter && ! &{$hitfilter}($hit) );
	my $nm = $hit->name();
	my $id_parser = $self->id_parser;
	print STDERR "no $nm for name (",$hit->description(), ")\n" unless $nm;
	my ($gi,$acc) = &$id_parser($nm);
	my $p = "%-".$MaxDescLen. "s";
	my $descsub = $hit->name . " " . $hit->description;
	if( length($descsub) > ($MaxDescLen) ) {
	    $descsub = substr($descsub,0,$MaxDescLen-3) . "...";
	}

	$str .= sprintf("$p   %s %-4s  %s\n",
			$descsub,
			defined $hit->raw_score ? $hit->raw_score : ' ',
			defined $hit->significance ? $hit->significance : '?');

	$hspstr .= sprintf(">%s %s\n%9sLength = %d\n\n",
			   $hit->name, 
			   defined $hit->description ? $hit->description : '', 
			   '', # empty is for the %9s in the str formatting 
			   $hit->length);
	$hit->rewind(); # make sure we are at the beginning of the list
	while( my $hsp = $hit->next_hsp ) {
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
			     POSIX::floor($hsp->frac_identical('total') * 100));

	    if( $hit->algorithm =~ /(T)?BLAST(X|P)/ ||
		$hit->algorithm =~ /(T)?FAST(X|Y|P)/ ||
		$hit->algorithm eq 'TBLASTN' ||
		$hit->algorithm eq 'TFASTN' 
		) {
		$hspstr .= sprintf(", Positives = %d/%d (%d%%)",
				   ( $hsp->frac_conserved('total') * 
				     $hsp->length('total')),
				   $hsp->length('total'),
				   POSIX::floor($hsp->frac_conserved('total') * 100));
		
	    }
	    if( $hsp->gaps ) {
		$hspstr .= sprintf(", Gaps = %d/%d (%d%%)\n",
				   $hsp->gaps('total'),
				   $hsp->length('total'),
				   POSIX::floor(100 * $hsp->gaps('total') / 
					       $hsp->length('total')));
	    }
	    
	    my ($h,$q) = ( $hsp->hit->frame ,$hsp->query->frame);
		
	    if( $h || $q ) {
		$hspstr .= " Frame = ";
		
		if( $h && ! $q) {  $hspstr .= $h }
		elsif( $q && ! $h) {  $hspstr .= $q }
		else { 
		    $hspstr .= " $h / $q \n";
		}
	    }
	    $hspstr .= "\n";
	    
	    my @hspvals = ( {'name' => 'Query:',
			     'seq'  => $hsp->query_string,
			     'start' => $hsp->query->strand >= 0 ? 
				 $hsp->query->start : $hsp->query->end,
			     'end'   => $hsp->query->strand >= 0 ? 
				  $hsp->query->end : $hsp->query->start,
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
			      'direction' => $hsp->query->strand || 1
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
		$hspstr .= "\n";
	    }
	}
	$hspstr .= "\n";
    }
    $str .= "\n\n".$hspstr;

    foreach my $param ( $result->available_parameters ) {
	$str .= "$param: ". $result->get_parameter($param) ."\n";
	
    }
    $str .= "Search Statistics\n";

   foreach my $stat ( sort $result->available_statistics ) {
	$str .= "$stat: ". $result->get_statistic($stat). "\n";
    }
    $str .=  "\n\n";
    return $str;
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

=head2 filter

 Title   : filter
 Usage   : $writer->filter('hsp', \&hsp_filter);
 Function: Filter out either at HSP,Hit,or Result level
 Returns : none
 Args    : string => data type,
           CODE reference


=cut


1;
