# $Id$
#
# BioPerl module for Bio::SearchIO::Writer::BSMLResultWriter
#
# Cared for by Jason Stajich <jason-at-bioperl-dot-org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::SearchIO::Writer::BSMLResultWriter - DESCRIPTION of Object

=head1 SYNOPSIS

  use Bio::SearchIO;
  my $in = new Bio::SearchIO(-file   => 'result.blast',
                             -format => 'blast');
  my $out = new Bio::SearchIO(-output_format  => 'BSMLResultWriter',
                              -file           => ">result.bsml");
  while( my $r = $in->next_result ) {
    $out->write_result($r);
  }

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
the web:

  http://bugzilla.bioperl.org/

=head1 AUTHOR - Jason Stajich

Email jason-at-bioperl-dot-org

Describe contact details here

=head1 CONTRIBUTORS

Additional contributors names and emails here

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::SearchIO::Writer::BSMLResultWriter;
use vars qw(@ISA);
use strict;

use Bio::Root::Root;
use Bio::SearchIO::SearchWriterI;
use XML::Writer;
use IO::String;

@ISA = qw(Bio::Root::Root Bio::SearchIO::SearchWriterI);


=head2 new

 Title   : new
 Usage   : my $obj = new Bio::SearchIO::Writer::BSMLResultWriter();
 Function: Builds a new Bio::SearchIO::Writer::BSMLResultWriter object 
 Returns : an instance of Bio::SearchIO::Writer::BSMLResultWriter
 Args    :


=cut

sub new {
  my($class,@args) = @_;

  my $self = $class->SUPER::new(@args);
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
    my $str = new IO::String();
    my $newlines = 1; # can parameterize this later
    my $writer = new XML::Writer(OUTPUT     => $str,
				 NEWLINES   => $newlines);
    $writer->xmlDecl('UTF-8');
    $writer->doctype('Bsml','-//EBI//Labbook, Inc. BSML DTD//EN',
		     'http://www.labbook.com/dtd/bsml3_1.dtd');
    $writer->startTag('Bsml');
    $writer->startTag('Definitions');
    $writer->startTag('Sequences');
    my $reporttype = $result->algorithm;
    my ($qmoltype,$hmoltype);
    my $hit = $result->next_hit;
    my $hsp = $hit->next_hsp;
    if( $hsp->query->strand == 0 ) { $qmoltype = 'aa' }
    else { $qmoltype = 'nt' }
    
    if( $hsp->hit->strand == 0 ) { $hmoltype = 'aa' }
    else { $hmoltype = 'nt' }
	
    $writer->startTag('Sequence',
		      'length' => $result->query_length,
		      'title'  => $result->query_name,
		      'molecule' => $qmoltype,
		      'id'     => $result->query_name
		      );
    $writer->endTag('Sequence');
    foreach my $hit ( $result->hits ) {
	$writer->startTag('Sequence',
			  'length' => $hit->length,
			  'title'  => $hit->name,
			  'molecule' => $hmoltype,
			  'id'     => $hit->name
		      );
	$writer->endTag('Sequence');
    }
    $writer->endTag('Sequences');

    $writer->startTag('Tables');
    $writer->startTag('Sequence-search-table',
		      'searcg-type' => $result->algorithm,
		      'query-start' => 0,
		      'query-length' => $result->query_length);
		      
    foreach my $hit ( $result->hits ) {
	$writer->startTag('Seq-pair-alignment',
			  'compxref' => sprintf("%s:%s",
						'',$result->query_name),
			  'compseq' => $result->query_name,
			  'method'  => join(' ',$result->algorithm, 
					    $result->algorithm_version),
			  'refxref'  => sprintf("%s:%s",
						$result->database_name,
						$hit->name),
			  'refseq'  => $hit->name,
			  'refstart' => 0,
			  'refend'   => $hit->length -1,
			  'reflength' => $hit->length);
	foreach my $hsp ( $hit->hsps ) {
	    $writer->startTag('Seq-pair-run',
			      'runlength' => $hit->hit->length,
			      'runprob' => '0',
			      'translated' => $result->query_algorithm =~ /[TX]/,
			      'comprunlength' => $hsp->hsp_length,
			      'comppos' => $hsp->query->start,
			      'compcomplement' => '0',
			      'complength' => $hsp->hit->length,
			      'refpos' => $hsp->hit->start, 
			      'runscore' => $hsp->score,
			      'refcomplement' => '0');
	    $writer->emptyTag('Attribute',
			      'name'    => 'percent_identity',
			      'content' => $hsp->percent_identity);
	    $writer->emptyTag('Attribute',
			      'name'    => 'percent_similarity',
			      'content' => $hsp->frac_conserved * 100);	    
	    $writer->emptyTag('Attribute',
			      'name'    => 'e_value',
			      'content' => $hsp->evalue);
	}
	$writer->endTag('Seq-pair-alignment');
    }
    $writer->endTag('Sequence-search-tables');
    $writer->endTag('Tables');
    
    $writer->endTag('Definitions');    
    $writer->endTag('Bsml');   
    $writer->end();
    return ${$str->string_ref};
}
1;
