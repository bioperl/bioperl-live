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

    # do not use this object directly.  See Bio::SearchIO documentation

=head1 DESCRIPTION

This object implements a NCBI Blast XML parser. 

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

package Bio::SearchIO::blastxml;
use vars qw(@ISA $DTD %MAPPING %MODEMAP);
use strict;

$DTD = 'ftp://ftp.ncbi.nlm.nih.gov/blast/documents/NCBI_BlastOutput.dtd';
# Object preamble - inherits from Bio::Root::Root

use Bio::Root::Root;
use Bio::SearchIO::EventGeneratorI;
use Bio::SearchIO;
use XML::Parser::PerlSAX;
use XML::Handler::Subs;

BEGIN { 
    # mapping of NCBI Blast terms to Bioperl hash keys
    %MODEMAP = ('BlastOutput' => 'result',
		'Hit'         => 'hit',
		'Hsp'         => 'hsp'
		);

    %MAPPING = ( 
		 'Hsp_bit-score' => 'bits',
		 'Hsp_score'     => 'score',
		 'Hsp_evalue'    => 'evalue',
		 'Hsp_query-from'=> 'querystart',
		 'Hsp_query-to'  => 'queryend',
		 'Hsp_hit-from'  => 'hitstart',
		 'Hsp_hit-to'    => 'hitend',
		 'Hsp_positive'  => 'conserved',
		 'Hsp_identity'  => 'identical',
		 'Hsp_gaps'      => 'gaps',
		 'Hsp_hitgaps'   => 'hitgaps',
		 'Hsp_querygaps' => 'querygaps',
		 'Hsp_qseq'      => 'queryseq',
		 'Hsp_hseq'      => 'hitseq',
		 'Hsp_midline'   => 'homolseq',
		 'Hsp_align-len' => 'hsplen',
		 'Hsp_query-frame'=> 'queryframe',
		 'Hsp_hit-frame'  => 'hitframe',

		 'Hit_id'        => 'hitname',
		 'Hit_len'       => 'hitlen',
		 'Hit_accession' => 'hitacc',
		 'Hit_def'       => 'hitdesc',
		 
		 'BlastOutput_program'  => 'programname',
		 'BlastOutput_version'  => 'programver',
		 'BlastOutput_query-def'=> 'queryname',
		 'BlastOutput_query-len'=> 'querylen',
		 'BlastOutput_db'       => 'dbname',
		 
		 'Iteration_iter-num'   => 'iternum',
		 'Parameters_matrix'    => { 'param' => 'matrix'},
		 'Parameters_expect'    => { 'param' => 'expect'},
		 'Parameters_include'   => { 'param' => 'include'},
		 'Parameters_sc-match'  => { 'param' => 'match'},
		 'Parameters_sc-mismatch' => { 'param' => 'mismatch'},
		 'Parameters_gap-open'  => { 'param' => 'gapopen'},
		 'Parameters_gap-extend'=> { 'param' => 'gapext'},
		 'Parameters_filter'    => {'param' => 'filter'},
		 'Statistics_db-num'    => 'dbsize',
		 'Statistics_db-len'    => 'dblets',
		 'Statistics_hsp-len'   => { 'stat' => 'hsplength'},
		 'Statistics_eff-space' => { 'stat' => 'effectivespace'},
		 'Statistics_kappa'     => { 'stat' => 'kappa' },
		 'Statistics_lambda'    => { 'stat' => 'lambda' },
		 'Statistics_entropy'   => { 'stat' => 'entropy'},
		 );
}

@ISA = qw(Bio::SearchIO );

=head2 _initialize

 Title   : _initialize
 Usage   : private
 Function: Initializes the object - this is chained through new in SearchIO

=cut

sub _initialize{
   my ($self,@args) = @_;   
   $self->SUPER::_initialize(@args);
   $self->{'_xmlparser'} = new XML::Parser::PerlSAX();
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
 
    my $data = '';
    my $firstline = 1;
    while( defined( $_ = $self->_readline) ) {
	if( /^RPS-BLAST/i ) {
	    $self->{'_type'} = 'RPSBLAST';
	    next;
	}
	if( /^<\?xml version/ && ! $firstline) { 
	    $self->_pushback($_);
	    last;
	}
	s/&apos;/\`/g;	
	$data .= $_;
	$firstline = 0;
    }

    return undef unless( $data);
    my %parser_args = ('Source' => { 'String' => $data },
		       'Handler' => $self);
    
    my $result = $self->{'_xmlparser'}->parse(%parser_args);

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
    $self->{'_mode'} = '';
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
	$self->_mode($type);
	if( $self->_eventHandler->will_handle($type) ) {
	    my $func = sprintf("start_%s",lc $type);
	    $self->_eventHandler->$func($data->{'Attributes'});
	}						     
    }

    if($nm eq 'BlastOutput') {
	$self->{'_values'} = {};
	$self->{'_result'}= undef;
	$self->{'_mode'} = '';
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
    } else { 
	
	$self->debug("ignoring unrecognized element type $nm");
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
   
   $self->{'_last_data'} = $data->{'Data'}; 
}


=head2 _mode

 Title   : _mode
 Usage   : $obj->_mode($newval)
 Function: private methods do not use outside of object
 Example : 
 Returns : value of _mode
 Args    : newvalue (optional)


=cut

sub _mode{
   my ($self,$value) = @_;
   if( defined $value) {
      $self->{'_mode'} = $value;
    }
    return $self->{'_mode'};

}
1;
