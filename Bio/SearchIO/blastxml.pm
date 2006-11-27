# $Id$
#
# BioPerl module for Bio::SearchIO::blastxml
#
# Cared for by Jason Stajich <jason@bioperl.org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::SearchIO::blastxml - A SearchIO implementation of NCBI Blast XML parsing. 

=head1 SYNOPSIS

    use Bio::SearchIO;
    my $searchin = new Bio::SearchIO(-format => 'blastxml',
				     -file   => 't/data/plague_yeast.bls.xml');
    while( my $result = $searchin->next_result ) {
    }

    # one can also request that the parser NOT keep the XML data in memory
    # by using the tempfile initialization flag.
    my $searchin = new Bio::SearchIO(-tempfile => 1,
				     -format => 'blastxml',
				     -file   => 't/data/plague_yeast.bls.xml');
    while( my $result = $searchin->next_result ) {
    }

=head1 DESCRIPTION

This object implements a NCBI Blast XML parser.  It requires XML::SAX; it is
also recommended (for faster parsing) that XML::SAX::ExpatXS be installed and
set as the default parser in ParserDetails.ini.  This file is located in the
SAX subdirectory of XML in your local perl library (normally in the 'site'
directory).  Currently, XML::SAX::Expat will NOT work as expected if set as
default; you must have local copies of the NCBI DTDs if using XML::SAX::Expat.

There is one additional initialization flag from the SearchIO defaults
- that is the -tempfile flag.  If specified as true, then the parser
will write out each report to a temporary filehandle rather than
holding the entire report as a string in memory.  The reason this is
done in the first place is NCBI reports have an uncessary E<lt>?xml
version="1.0"?E<gt> at the beginning of each report and RPS-BLAST reports
have an additional unecessary RPS-BLAST tag at the top of each report.
So we currently have implemented the work around by preparsing the
file (yes it makes the process slower, but it works).

=head1 DEPENDENCIES

In addition to parts of the Bio:: hierarchy, this module uses:

 XML::SAX

It is also recommended that XML::SAX::ExpatXS be installed and made the default
XML::SAX parser using , along with the
Expat library () for faster parsing.  XML::SAX::Expat is not recommended; 
XML::SAX::ExpatXS is considered the current replacement for XML::SAX:Expat
and is actively being considered to replace XML::SAX::Expat.  XML::SAX::Expat
will work, but only if you have local copies of the NCBI BLAST DTDs. This is
due to issues with NCBI's BLAST XML format.  The DTDs and the web address to
obtain them are:

  NCBI_BlastOutput.dtd	    
  NCBI_BlastOutput.mod.dtd

  http://www.ncbi.nlm.nih.gov/data_specs/dtd/

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

  bioperl-l@bioperl.org                  - General discussion
  http://bioperl.org/wiki/Mailing_lists  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution. Bug reports can be submitted via the
web:

  http://bugzilla.open-bio.org/

=head1 AUTHOR - Jason Stajich

Email jason-at-bioperl.org

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::SearchIO::blastxml;
use vars qw($DTD %MAPPING %MODEMAP $DEBUG);
use strict;

$DTD = 'ftp://ftp.ncbi.nlm.nih.gov/blast/documents/NCBI_BlastOutput.dtd';
# Object preamble - inherits from Bio::Root::Root

use Bio::Root::Root;
use XML::SAX;
use HTML::Entities;
use IO::File;

BEGIN {
    # uncomment only for testing; trying to get XML::SAX::Expat to play nice...
    #$XML::SAX::ParserPackage = 'XML::SAX::Expat';
    # mapping of NCBI Blast terms to Bioperl hash keys
    %MODEMAP = ('BlastOutput' => 'result',
		'Hit'         => 'hit',
		'Hsp'         => 'hsp'
		);

    %MAPPING = ( 
		 # HSP specific fields
		 'Hsp_bit-score'  => 'HSP-bits',
		 'Hsp_score'      => 'HSP-score',
		 'Hsp_evalue'     => 'HSP-evalue',
		 'Hsp_query-from' => 'HSP-query_start',
		 'Hsp_query-to'   => 'HSP-query_end',
		 'Hsp_hit-from'   => 'HSP-hit_start',
		 'Hsp_hit-to'     => 'HSP-hit_end',
		 'Hsp_positive'   => 'HSP-conserved',
		 'Hsp_identity'   => 'HSP-identical',
		 'Hsp_gaps'       => 'HSP-gaps',
		 'Hsp_hitgaps'    => 'HSP-hit_gaps',
		 'Hsp_querygaps'  => 'HSP-query_gaps',
		 'Hsp_qseq'       => 'HSP-query_seq',
		 'Hsp_hseq'       => 'HSP-hit_seq',
		 'Hsp_midline'    => 'HSP-homology_seq',
		 'Hsp_align-len'  => 'HSP-hsp_length',
		 'Hsp_query-frame'=> 'HSP-query_frame',
		 'Hsp_hit-frame'  => 'HSP-hit_frame',

		 # these are ignored for now
		 'Hsp_num'          => 'HSP-order',
		 'Hsp_pattern-from' => 'patternend',
		 'Hsp_pattern-to'   => 'patternstart',
		 'Hsp_density'      => 'hspdensity',

		 # Hit specific fields
		 'Hit_id'               => 'HIT-name',
		 'Hit_len'              => 'HIT-length',
		 'Hit_accession'        => 'HIT-accession',
		 'Hit_def'              => 'HIT-description',
		 'Hit_num'              => 'HIT-order',
		 'Iteration_iter-num'   => 'HIT-iteration',
		 'Iteration_stat'       => 'HIT-iteration_statistic',
		 
		 'BlastOutput_program'   => 'RESULT-algorithm_name',
		 'BlastOutput_version'   => 'RESULT-algorithm_version',
		 'BlastOutput_query-def' => 'RESULT-query_description',
		 'BlastOutput_query-len' => 'RESULT-query_length',
		 'BlastOutput_db'        => 'RESULT-database_name',
		 'BlastOutput_reference' => 'RESULT-program_reference',
		 'BlastOutput_query-ID'  => 'runid',
		 
		 'Parameters_matrix'    => { 'RESULT-parameters' => 'matrix'},
		 'Parameters_expect'    => { 'RESULT-parameters' => 'expect'},
		 'Parameters_include'   => { 'RESULT-parameters' => 'include'},
		 'Parameters_sc-match'  => { 'RESULT-parameters' => 'match'},
		 'Parameters_sc-mismatch' => { 'RESULT-parameters' => 'mismatch'},
		 'Parameters_gap-open'  => { 'RESULT-parameters' => 'gapopen'},
		 'Parameters_gap-extend'=> { 'RESULT-parameters' => 'gapext'},
		 'Parameters_filter'    => {'RESULT-parameters' => 'filter'},
		 'Statistics_db-num'    => 'RESULT-database_entries',
		 'Statistics_db-len'    => 'RESULT-database_letters',
		 'Statistics_hsp-len'   => { 'RESULT-statistics' => 'hsplength'},
		 'Statistics_eff-space' => { 'RESULT-statistics' => 'effectivespace'},
		 'Statistics_kappa'     => { 'RESULT-statistics' => 'kappa' },
		 'Statistics_lambda'    => { 'RESULT-statistics' => 'lambda' },
		 'Statistics_entropy'   => { 'RESULT-statistics' => 'entropy'},
		 );
    eval {  require Time::HiRes };	
    if( $@ ) { $DEBUG = 0; }
}


use base qw(Bio::SearchIO);

=head2 new

 Title   : new
 Usage   : my $searchio = new Bio::SearchIO(-format => 'blastxml',
					    -file   => 'filename',
					    -tempfile => 1);
 Function: Initializes the object - this is chained through new in SearchIO
 Returns : Bio::SearchIO::blastxml object
 Args    : One additional argument from the format and file/fh parameters.
           -tempfile => boolean.  Defaults to false.  Write out XML data
                                  to a temporary filehandle to send to 
                                  PerlSAX parser.
=cut

=head2 _initialize

 Title   : _initialize
 Usage   : private
 Function: Initializes the object - this is chained through new in SearchIO

=cut

sub _initialize{
    my ($self,@args) = @_;   
    $self->SUPER::_initialize(@args);
    my ($usetempfile) = $self->_rearrange([qw(TEMPFILE)],@args);
    defined $usetempfile && $self->use_tempfile($usetempfile);
    $self->{'_xmlparser'} = XML::SAX::ParserFactory->parser(Handler => $self);
    my $local_parser = ref($self->{'_xmlparser'});
    if ($local_parser eq 'XML::SAX::Expat') {
        $self->warn('XML::SAX::Expat not currently supported; '.
                    'must have local copies of NCBI DTD docs!');
    }    
    $DEBUG = 1 if( ! defined $DEBUG && $self->verbose > 0);
}

=head2 next_result

 Title   : next_result
 Usage   : my $hit = $searchio->next_result;
 Function: Returns the next Result from a search
 Returns : Bio::Search::Result::ResultI object
 Args    : none

=cut

sub next_result {
    my ($self) = @_;
    local $/ = "\n";
    local $_;
 
    my $data = '';
    my $firstline = 1;
    my ($tfh);
    if( $self->use_tempfile ) {
	$tfh = IO::File->new_tmpfile or $self->throw("Unable to open temp file: $!");	
	$tfh->autoflush(1);
    }
   
    my ($sawxmlheader,$okaytoprocess,$sawdoctype);
    while( defined( $_ = $self->_readline) ) {
	if( /^RPS-BLAST/i ) {
	    $self->{'_type'} = 'RPS-BLAST';
	    next;
	}
	if( /^<\?xml version/ ) {
	    if( ! $firstline ) {
		$self->_pushback($_);
		last;
	    }
	    $sawxmlheader = 1;
	} 
	# for the non xml version prefixed in each section
	if( /DOCTYPE/ ) { #|| /<BlastOutput>/
	    if(  $sawdoctype ) {
		if( ! $sawxmlheader ) { 
		    $self->_pushback("<?xml version=\"1.0\"?>\n");
		}
		$self->_pushback($_);
		last;
	    }
	    $sawdoctype = 1;
	    unless( $sawxmlheader ) {
		$self->debug( "matched here\n");
		$self->_pushback("<?xml version=\"1.0\"?>\n");
		$self->_pushback($_);
		next;
	    }
	}
	$okaytoprocess = 1;
	if( defined $tfh ) {
	    print $tfh $_;
	} else {
	    $data .= $_;
	}
	$firstline = 0;
    }
    return unless( $okaytoprocess);
    
    my %parser_args;
    if( defined $tfh ) {
	seek($tfh,0,0);
	%parser_args = ('Source' => { 'ByteStream' => $tfh });
    } else {
	%parser_args = ('Source' => { 'String' => $data });
    }
    my $result;
    my $starttime;
    #if(  $DEBUG ) {  $starttime = [ Time::HiRes::gettimeofday() ]; }

    eval { 
	$result = $self->{'_xmlparser'}->parse(%parser_args);
        $self->{'_result_count'}++;
    };
    if( $@ ) {
	$self->warn("error in parsing a report:\n $@");
	$result = undef;
    }    
    #if( $DEBUG ) {
	#$self->debug( sprintf("parsing took %f seconds\n", Time::HiRes::tv_interval($starttime)));
    #}
    # parsing magic here - but we call event handlers rather than 
    # instantiating things 
    return $result;
}

=head2 SAX methods

=cut

=head2 start_document

 Title   : start_document
 Usage   : $parser->start_document;
 Function: SAX method to indicate starting to parse a new document
 Returns : none
 Args    : none


=cut

sub start_document{
    my ($self) = @_;
    $self->{'_lasttype'} = '';
    $self->{'_values'} = {};
    $self->{'_result'}= undef;
}

=head2 end_document

 Title   : end_document
 Usage   : $parser->end_document;
 Function: SAX method to indicate finishing parsing a new document
 Returns : Bio::Search::Result::ResultI object
 Args    : none

=cut

sub end_document{
   my ($self,@args) = @_;
   return $self->{'_result'};
}

=head2 start_element

 Title   : start_element
 Usage   : $parser->start_element($data)
 Function: SAX method to indicate starting a new element
 Returns : none
 Args    : hash ref for data

=cut

sub start_element{
    my ($self,$data) = @_;
    # we currently don't care about attributes
    my $nm = $data->{'Name'};    

    if( my $type = $MODEMAP{$nm} ) {
	if( $self->_eventHandler->will_handle($type) ) {
	    my $func = sprintf("start_%s",lc $type);
	    $self->_eventHandler->$func($data->{'Attributes'});
	}						     
    }

    if($nm eq 'BlastOutput') {
	$self->{'_values'} = {};
	$self->{'_result'}= undef;
    }
}

=head2 end_element

 Title   : end_element
 Usage   : $parser->end_element($data)
 Function: Signals finishing an element
 Returns : Bio::Search object dpending on what type of element
 Args    : hash ref for data

=cut

sub end_element{
    my ($self,$data) = @_;

    my $nm = $data->{'Name'};
    my $rc;
    if($nm eq 'BlastOutput_program' &&
       $self->{'_last_data'} =~ /(t?blast[npx])/i ) {
	$self->{'_type'} = uc $1; 
    }

    if( my $type = $MODEMAP{$nm} ) {
	if( $self->_eventHandler->will_handle($type) ) {
	    my $func = sprintf("end_%s",lc $type);
	    $rc = $self->_eventHandler->$func($self->{'_type'},
					      $self->{'_values'});
	}
    } elsif( $MAPPING{$nm} ) { 
	if ( ref($MAPPING{$nm}) =~ /hash/i ) {
	    my $key = (keys %{$MAPPING{$nm}})[0];
	    $self->{'_values'}->{$key}->{$MAPPING{$nm}->{$key}} = $self->{'_last_data'};
	} else {
	    $self->{'_values'}->{$MAPPING{$nm}} = $self->{'_last_data'};
	}
    } elsif( $nm eq 'Iteration' || $nm eq 'Hit_hsps' || $nm eq 'Parameters' ||
	     $nm eq 'BlastOutput_param' || $nm eq 'Iteration_hits' || 
	     $nm eq 'Statistics' || $nm eq 'BlastOutput_iterations' ){
        # ignores these elements for now; no iteration parsing
    } else { 	
	
	$self->debug("ignoring unrecognized element type $nm\n");
    }
    $self->{'_last_data'} = ''; # remove read data if we are at 
				# end of an element
    $self->{'_result'} = $rc if( $nm eq 'BlastOutput' );
    return $rc;
}

=head2 characters

 Title   : characters
 Usage   : $parser->characters($data)
 Function: Signals new characters to be processed
 Returns : characters read
 Args    : hash ref with the key 'Data'


=cut

sub characters{
   my ($self,$data) = @_;   
   return unless ( defined $data->{'Data'} && $data->{'Data'} !~ /^\s+$/ );
   $self->{'_last_data'} = &decode_entities($data->{'Data'}); 
}

=head2 use_tempfile

 Title   : use_tempfile
 Usage   : $obj->use_tempfile($newval)
 Function: Get/Set boolean flag on whether or not use a tempfile
 Example : 
 Returns : value of use_tempfile
 Args    : newvalue (optional)


=cut

sub use_tempfile{
   my ($self,$value) = @_;
   if( defined $value) {
      $self->{'_use_tempfile'} = $value;
    }
    return $self->{'_use_tempfile'};
}

sub result_count {
    my $self = shift;
    return $self->{'_result_count'};
}

1;
