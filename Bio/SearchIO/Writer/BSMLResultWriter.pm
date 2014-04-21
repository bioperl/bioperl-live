#
# BioPerl module for Bio::SearchIO::Writer::BSMLResultWriter
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Jason Stajich <jason-at-bioperl-dot-org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::SearchIO::Writer::BSMLResultWriter - BSML output writer

=head1 SYNOPSIS

  use Bio::SearchIO;
  my $in = Bio::SearchIO->new(-file   => 'result.blast',
                             -format => 'blast');
  my $out = Bio::SearchIO->new(-output_format  => 'BSMLResultWriter',
                              -file           => ">result.bsml");
  while( my $r = $in->next_result ) {
    $out->write_result($r);
  }

=head1 DESCRIPTION

This is a writer to produce BSML for a search result.

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
of the bugs and their resolution. Bug reports can be submitted via
the web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Jason Stajich

Email jason-at-bioperl-dot-org

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::SearchIO::Writer::BSMLResultWriter;
use strict;

use XML::Writer;
use IO::String;

use base qw(Bio::Root::Root Bio::SearchIO::SearchWriterI);


=head2 new

 Title   : new
 Usage   : my $obj = Bio::SearchIO::Writer::BSMLResultWriter->new();
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

# this implementation is largely adapted from the Incogen XSLT stylesheet
# to convert NCBI BLAST XML to BSML

sub to_string {
    my ($self,$result,$num) = @_;
    my $str = new IO::String();
    my $writer = new XML::Writer(OUTPUT     => $str,
				 DATA_INDENT => 1,
				 DATA_MODE   => 1);
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
		      'title'  => $result->query_name . " ". $result->query_description,
		      'molecule' => $qmoltype,
		      'representation' => 'virtual',
		      'id'     => $result->query_name
		      );
    # Here we're annotating the Query sequence with hits
    # hence the Feature-table
    $writer->startTag('Feature-tables');
    $writer->startTag('Feature-table',
		      'title' => "$reporttype Result", 
		      'class' => $reporttype);
    my ($hitnum,$hspnum) = (1,1);
    foreach my $hit ( $result->hits ) {	
	$hspnum = 1;
	foreach my $hsp ( $hit->hsps ) {
	    $writer->startTag('Feature',
			      'class'  => $reporttype,
			      'value-type' => 'alignment',
			      'title'  => $hit->name. " ". $hit->description,
			      );

	    $writer->emptyTag('Interval-loc',
			      'startpos' => $hsp->query->start,
			      'endpos'   => $hsp->query->end);
	    $writer->emptyTag('Qualifier',
			      'value-type' => 'score',
			      'value'      => $hsp->score,
			      );
	    
	    $writer->emptyTag('Qualifier',
			      'value-type' => 'target-start',
			      'value'      => $hsp->hit->start,
			      );
	    $writer->emptyTag('Qualifier',
			      'value-type' => 'target-end',
			      'value'      => $hsp->hit->end,
			      );
	    $writer->emptyTag('Link',
			      'title' => 'alignment',
			      'href'  => sprintf("#SPA%d.%d",$hitnum,$hspnum)
			      );
	    
	    if( $hsp->hit->strand < 0 ) {
		$writer->emptyTag('Qualifier',
				  'value-type' => 'target-on-complement',
				  'value'      => 1,
				  );
	    }
	    $hspnum++;
	    $writer->endTag('Feature');
	}
	$hitnum++;
    }
    $writer->endTag('Feature-table');
    $writer->endTag('Feature-tables');
    $writer->endTag('Sequence');
    $writer->endTag('Sequences');

    $writer->startTag('Tables');
    $writer->startTag('Sequence-search-table',
		      'search-type' => $reporttype,
		      'query-length' => $result->query_length);
    $hitnum = $hspnum = 1;
    foreach my $hit ( $result->hits ) {
	$hspnum = 1;
	foreach my $hsp ( $hit->hsps ) {
	    $writer->startTag('Seq-pair-alignment',
			      'id' => sprintf("SPA%d.%d",$hitnum,$hspnum),
			      'method'       => join(' ',$result->algorithm), 
			      'compxref'     => sprintf("%s:%s",
						'',$result->query_name),
			      'refxref'      => sprintf("%s:%s",
							$result->database_name,
							$hit->name),
			      'refseq'       => $hit->name,
			      'title'        => $result->query_name,
			      'compseq'      => $result->query_name,
			      'compcaption'  => $result->query_name . ' ' .
			                         $result->query_description,
			      'refcaption'   => $hit->name . " ". 
                                                 $hit->description,
			      'totalscore'   => $hsp->score,
			      'refstart'     => $hsp->query->start,
			      'refend'       => $hsp->query->end,
			      'compstart'    => $hsp->hit->start,
			      'compend'      => $hsp->hit->end,
			      'complength'   => $hit->length,
			      'reflength'    => $result->query_length);

	    $writer->emptyTag('Attribute',
			      'name'    => 'hit-num',
			      'content' => $hitnum);
	    $writer->emptyTag('Attribute',
			      'name'    => 'hit-id',
			      'content' => $hit->name);
	    $writer->emptyTag('Attribute',
			      'name'    => 'hsp-num',
			      'content' => $hspnum);
	    $writer->emptyTag('Attribute',
			      'name'    => 'hsp-bit-score',
			      'content' => $hsp->bits);
	    $writer->emptyTag('Attribute',
			      'name'    => 'hsp-evalue',
			      'content' => $hsp->evalue);
	    $writer->emptyTag('Attribute',
			      'name'    => 'pattern-from',
			      'content' => 0);
	    $writer->emptyTag('Attribute',
			      'name'    => 'pattern-to',
			      'content' => 0);
	    $writer->emptyTag('Attribute',
			      'name'    => 'query-frame',
			      'content' => $hsp->query->frame);
	    $writer->emptyTag('Attribute',
			      'name'    => 'hit-frame',
			      'content' => $hsp->hit->frame * $hsp->hit->strand);
	    $writer->emptyTag('Attribute',
			      'name'    => 'percent_identity',
			      'content' => sprintf("%.2f",$hsp->percent_identity));
	    $writer->emptyTag('Attribute',
			      'name'    => 'percent_similarity',
			      'content' => sprintf("%.2f",$hsp->frac_conserved('total') * 100));	    
	    my $cons = $hsp->frac_conserved('total') * $hsp->length('total');
	    my $ident = $hsp->frac_identical('total') * $hsp->length('total');
	    
	    $writer->emptyTag('Attribute',
			      'name'    => 'identity',
			      'content' => $ident);
	    $writer->emptyTag('Attribute',
			      'name'    => 'positive',
			      'content' => $cons);
	    $writer->emptyTag('Attribute',
			      'name'    => 'gaps',
			      'content' => $hsp->gaps('total'));
	    $writer->emptyTag('Attribute',
			      'name'    => 'align-len',
			      'content' => $hsp->length('total'));
	    $writer->emptyTag('Attribute',
			      'name'    => 'density',
			      'content' => 0);
	    $writer->emptyTag('Attribute',
			      'name'    => 'hit-len',
			      'content' => $hit->length);
	    my @extrafields;

	    $writer->emptyTag('Seq-pair-run',
			      'runlength'     => $hsp->hit->length,
			      'comprunlength' => $hsp->hsp_length,
			      'complength'    => $hsp->hit->length,
			      'compcomplement'=> $hsp->hit->strand < 0 ? 1 :0,
			      'refcomplement' => $hsp->query->strand < 0 ? 1 :0,
			      'refdata'       => $hsp->query_string,
			      'compdata'      => $hsp->hit_string,
			      'alignment'     => $hsp->homology_string,
			      );
	    $hspnum++;
	    $writer->endTag('Seq-pair-alignment');
	}
	$hitnum++;
    }
    $writer->endTag('Sequence-search-table');
    $writer->endTag('Tables');
    
    $writer->startTag('Research');
    $writer->startTag('Analyses');
    $writer->startTag('Analysis');
    $writer->emptyTag('Attribute',
		      'name'    => 'program',
		      'content' => $reporttype);
    $writer->emptyTag('Attribute',
		      'name'    => 'version',
		      'content' => join(' ',$reporttype, 
					$result->algorithm_version));
    $writer->emptyTag('Attribute',
		      'name'     => 'reference',
		      'content'  => $result->algorithm_reference);
    $writer->emptyTag('Attribute',
		      'name'     => 'db',
		      'content'  => $result->database_name);
    $writer->emptyTag('Attribute',
		      'name'     => 'db-size',
		      'content'  => $result->database_entries);
    $writer->emptyTag('Attribute',
		      'name'     => 'db-length',
		      'content'  => $result->database_letters);
    # $writer->emptyTag('Attribute',
    # 'name'     => 'iter-num',
    # 'content'  => $result->iteration_num);
    foreach my $attr ( $result->available_parameters ) {
	$writer->emptyTag('Attribute',
			  'name'     => $attr,
			  'content'  => $result->get_parameter($attr));
    }
    foreach my $attr ( $result->available_statistics ) {
	$writer->emptyTag('Attribute',
			  'name'     => $attr,
			  'content'  => $result->get_statistic($attr));
    }
    $writer->endTag('Analysis');    
    $writer->endTag('Analyses');    
    $writer->endTag('Research');
    
    $writer->endTag('Definitions');   
    $writer->endTag('Bsml');   
    $writer->end();
    return ${$str->string_ref};
}
1;
