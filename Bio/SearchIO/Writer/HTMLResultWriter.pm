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

Give standard usage here

=head1 DESCRIPTION

Describe the object here

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
  http://bioperl.org/bioperl-bugs/

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
use vars qw(@ISA $RemoteURLDefault $MaxDescLen $AlignmentLineWidth);
use strict;

# Object preamble - inherits from Bio::Root::RootI

BEGIN {
    $RemoteURLDefault = 'http://www.ebi.ac.uk/cgi-bin/dbfetch?id=%s';    
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
  my ($url) = $self->_rearrange([qw(REMOTEDBURL)],@args);
  $self->remote_database_url($url || $RemoteURLDefault);
  return $self;
}

=head2 remote_database_url

 Title   : remote_database_url
 Usage   : $obj->remote_database_url($newval)
 Function: This should return or set a string that contains a %s which can be
           filled in with sprintf.
 Returns : value of remote_database_url
 Args    : newvalue (optional)


=cut

sub remote_database_url{
   my ($self,$value) = @_;
   if( defined $value) {
      $self->{'remote_database_url'} = $value;
    }
    return $self->{'remote_database_url'};
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
    return unless defined $result;
    my $str = sprintf(
qq{<HTML>
    <HEAD><TITLE>Bioperl HTML %s Search output for Bio::Search results</TITLE></HEAD>
    <BODY BGCOLOR="WHITE">
    <H1>Bioperl HTML %s Search Report for %s</H1>
    <hr>
    <b>Query=</b>%s %s<br><dd>(%d letters)</dd>
    <p>
    <b>Database:</b> %s<br><dd>%d sequences; %d total letters<p></dd>
    <p>
    <table border=0>
    <tr><th>Sequences producing significant alignments:</th>
	<th>Score<br>(bits)</th><th>E<br>value</th></tr>
  }, 
		      $result->algorithm, $result->algorithm,
		      $result->query_name(), $result->query_name, 
		      $result->query_description, $result->query_length, 
		      $result->database_name(),
		      $result->database_entries(),$result->database_letters(),
		      );
    my $hspstr = '<p><p>';
    while( my $hit = $result->next_hit ) {
	my $nm = $hit->name();
	my $id_parser = $self->id_parser;
	my ($gi,$acc) = &$id_parser($nm);
	my $descsub = substr($hit->description,0,$MaxDescLen);
	$descsub .= " ..." if length($descsub) < length($hit->description);

	$str .= sprintf('<tr><td><a href="%s">%s</a> %s</td><td><a href="#%s">%s</a></td><td>%s</td></tr>'."\n",
			sprintf($self->remote_database_url, $gi || $acc),
			$hit->name, $descsub,$acc,
			defined $hit->raw_score ? $hit->raw_score : '?',
			defined $hit->significance ? $hit->significance : '?');
	$hspstr .= "<a name=\"$acc\"><pre>\n".
	    sprintf(">%s %s\n<dd>Length = %d</dd><p>\n\n", $hit->name, 
			$hit->description, $hit->length);
	
	while( my $hsp = $hit->next_hsp ) {
	    $hspstr .= sprintf(" Score = %s bits (%s), Expect = %s",
			       $hsp->bits, $hsp->score, $hsp->evalue);
	    if( $hsp->pvalue ) {
		$hspstr .= ", P = ".$hsp->pvalue;
	    }
	    $hspstr .= "<br>\n";
	    $hspstr .= sprintf(" Identities = %d / %d (%d%%)",
			         ( $hsp->frac_identical('total') * 
				   $hsp->length('total')),
			       $hsp->length('total'),
			       $hsp->frac_identical('total') * 100);

	    if( $hit->algorithm =~ /(T)?BLAST(X|P)/ ||
		$hit->algorithm =~ /(T)?FAST(X|Y|P)/ ||
		$hit->algorithm eq 'TBLASTN' ||
		$hit->algorithm eq 'TFASTN' 
		) {
		$hspstr .= sprintf(", Positives = %d / %d (%d%%)",
				   ( $hsp->frac_conserved('total') * 
				     $hsp->length('total')),
				   $hsp->length('total'),
				   $hsp->frac_conserved('total') * 100);
		
	    }
	    if( $hsp->gaps ) {
		$hspstr .= sprintf(", Gaps = %d / %d (%d%%)",
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
			     'start' => $hsp->query->strand > 0 ? $hsp->query->start : $hsp->query->end,
			     'end'   => $hsp->query->strand > 0 ? $hsp->query->end : $hsp->query->start,
			     'index' => 0,
			     'direction' => $hsp->query->strand
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
			      'start' => $hsp->hit->strand > 0 ? $hsp->hit->start : $hsp->hit->end,
			      'end'   => $hsp->hit->strand > 0 ? $hsp->hit->end : $hsp->hit->start,
			      'index' => 0,
			      'direction' => $hsp->query->strand
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
		    my $piece = substr($v->{'seq'}, $v->{'index'},
				       $AlignmentLineWidth);
		    my $cp = $piece;
		    my $plen = scalar ( $cp =~ tr/\-//);
		    my ($start,$end) = ('','');
		    if( defined $v->{'start'} ) { 
			$start = $v->{'start'};
			# since strand can be + or - use the direction
			# to signify which whether to add or substract from end
			$end   = $v->{'start'} + ( ( $v->{'direction'} ) * 
					   ( $AlignmentLineWidth - $plen ));
			$v->{'start'} += $AlignmentLineWidth;
		    }
		    $hspstr .= sprintf("%s %-".$numwidth."s %s %s\n",
				       $v->{'name'},
				       $start,
				       $piece,
				       $end
				       );
		}
		$count += $AlignmentLineWidth;
		$hspstr .= '<p>';
	    }
	}
	$hspstr .= "</pre>\n";
    }
    $str .= "</table><p>\n".$hspstr;

    foreach my $param ( $result->available_parameters ) {
	$str .= "$param=". $result->get_parameter($param) ."<br>\n";
	
    }
    $str .= "<p>";
   foreach my $stat ( sort $result->available_statistics ) {
	$str .= "$stat=". $result->get_statistic($stat). "<br>\n";
    }
    
    $str .=  "<hr><h6>Produced by Bioperl module".ref($self)."</h6>\n</BODY>\n</HTML>\n";
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
    { $gi = $1; }
    
    if( $string =~ /(\w+)\|([A-Z\d\.\_]+)(\|[A-Z\d\_]+)?/ ) {
	$acc = defined $2 ? $2 : $1;
    } elsif ($string =~ /^\s*(\S+)/) {
        $acc = $1;
    } 

    return ($gi,$acc);
}
	
sub MIN { $a <=> $b ? $a : $b; }
sub MAX { $a <=> $b ? $b : $a; }

sub footer { 
    
}
1;
