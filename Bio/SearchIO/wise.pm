#
# BioPerl module for Bio::SearchIO::wise
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

Bio::SearchIO::wise - Parsing of wise output as alignments

=head1 SYNOPSIS

  use Bio::SearchIO;
  my $parser = Bio::SearchIO->new(-file    => 'file.genewise', 
                                 -format  => 'wise',
                                 -wisetype=> 'genewise');

  while( my $result = $parser->next_result ) {}

=head1 DESCRIPTION

This object parsers Wise output using Bio::Tools::Genewise or
Bio::Tools::Genomewise as a helper.

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


package Bio::SearchIO::wise;
use vars qw(%MAPPING %MODEMAP $DEFAULT_WRITER_CLASS);
use strict;

# Object preamble - inherits from Bio::Root::Root

use base qw(Bio::SearchIO);

%MODEMAP = ('WiseOutput' => 'result',
	    'Hit'             => 'hit',
	    'Hsp'             => 'hsp'
    );
%MAPPING =
    (
     'Hsp_query-from'=>  'HSP-query_start',
     'Hsp_query-to'  =>  'HSP-query_end',
     'Hsp_hit-from'  =>  'HSP-hit_start',
     'Hsp_hit-to'    =>  'HSP-hit_end',
     'Hsp_qseq'      =>  'HSP-query_seq',
     'Hsp_hseq'      =>  'HSP-hit_seq',
     'Hsp_midline'   =>  'HSP-homology_seq',
     'Hsp_score'     =>  'HSP-score',
     'Hsp_qlength'   =>  'HSP-query_length',
     'Hsp_hlength'   =>  'HSP-hit_length',
     'Hsp_align-len' =>  'HSP-hsp_length',
     'Hsp_positive'   => 'HSP-conserved',
     'Hsp_identity'   => 'HSP-identical',
     #'Hsp_gaps'       => 'HSP-hsp_gaps',
     #'Hsp_hitgaps'    => 'HSP-hit_gaps',
     #'Hsp_querygaps'  => 'HSP-query_gaps',
     
     'Hit_id'        => 'HIT-name',
#    'Hit_desc'      => 'HIT-description',
#    'Hit_len'       => 'HIT-length',
     'Hit_score'     => 'HIT-score',

     'WiseOutput_program'   => 'RESULT-algorithm_name',
     'WiseOutput_query-def' => 'RESULT-query_name',
     'WiseOutput_query-desc'=> 'RESULT-query_description',
     'WiseOutput_query-len' => 'RESULT-query_length',
    );

$DEFAULT_WRITER_CLASS = 'Bio::SearchIO::Writer::HitTableWriter';


use Bio::Tools::Genewise;
use Bio::Tools::Genomewise;

=head2 new

 Title   : new
 Usage   : my $obj = Bio::SearchIO::wise->new();
 Function: Builds a new Bio::SearchIO::wise object 
 Returns : an instance of Bio::SearchIO::wise
 Args    : -wise => a Bio::Tools::Genewise or Bio::Tools::Genomewise object


=cut

sub _initialize {
    my ($self,@args) = @_;
    my ( $wisetype, $file,$fh ) =
	$self->_rearrange([qw(WISETYPE FILE FH)], @args);
    my @newargs;
    while( @args ) {
	my $a = shift @args;
	if( $a =~ /FILE|FH/i ) {
	    shift @args;
	    next;
	}
	push @newargs, $a, shift @args;
    }
    $self->SUPER::_initialize(@newargs);

    # Optimization: caching the EventHandler 
    # since it's use a lot during the parse.
    $self->{'_handler_cache'} = $self->_eventHandler;

    $self->wisetype($wisetype);
    my @ioargs;
    if( $fh ) { 
	push @ioargs, ('-fh' => $fh);
    } elsif( $file ) {
	push @ioargs, ('-file' => $file);
    }

    if( $wisetype =~ /genewise/i ) {
	$self->wise(Bio::Tools::Genewise->new(@ioargs));
    } elsif( $wisetype =~ /genomewise/i ) {
	$self->wise(Bio::Tools::Genomewise->new(@ioargs));
    } else { 
	$self->throw("Must supply a -wisetype to ".ref($self)." which is one of 'genomewise' 'genewise'\n");
    }
    return $self;
}


=head2 next_result

 Title   : next_result
 Usage   : my $hit = $searchio->next_result;
 Function: Returns the next Result from a search
 Returns : Bio::Search::Result::ResultI object
 Args    : none

=cut

sub next_result{
   my ($self) = @_;
   local $/ = "\n";
   local $_;

   return unless $self->wise;
   my $prediction = $self->wise->next_prediction;
   return unless $prediction;
   $self->{'_reporttype'} = uc $self->wisetype;
   $self->start_element({'Name' => 'WiseOutput'});
   $self->element({'Name' => 'WiseOutput_program',
		   'Data' => $self->wisetype});
   $self->element({'Name' => 'WiseOutput_query-def',
		   'Data' => $self->wise->_prot_id});
   my @transcripts = $prediction->transcripts;

   foreach my $transcript ( @transcripts ) {
       my @exons =  $transcript->exons;
       my $protid;
       $self->start_element({'Name' => 'Hit'});
       
       if( $exons[0]->has_tag('supporting_feature') ) {
	   my ($supporting_feature) = $exons[0]->get_tag_values('supporting_feature');
	   $protid = $supporting_feature->feature2->seq_id;
	   $self->element({'Name' => 'Hit_id',
			   'Data' => $self->wise->_target_id});       
       } 
       $self->element({'Name' => 'Hit_score',
		       'Data' => $self->wise->_score});
       foreach my $exon ( @exons ) {
	   $self->start_element({'Name' => 'Hsp'});
	   if( $exon->strand < 0 ) { 
	       $self->element({'Name' => 'Hsp_query-from',
			       'Data' => $exon->end});
	       $self->element({'Name' => 'Hsp_query-to',
			       'Data' => $exon->start});
	   } else { 
	       $self->element({'Name' => 'Hsp_query-from',
			       'Data' => $exon->start});
	       $self->element({'Name' => 'Hsp_query-to',
			       'Data' => $exon->end});
	   }
	   $self->element({'Name' => 'Hsp_score',
			   'Data' => $self->wise->_score});
	   if( $exon->has_tag('supporting_feature') ) {
	       my ($sf) = $exon->get_tag_values('supporting_feature');
	       my $protein = $sf->feature2;
	       if( $protein->strand < 0 ) {
		   $self->element({'Name' => 'Hsp_hit-from',
				   'Data' => $protein->end});
		   $self->element({'Name' => 'Hsp_hit-to',
				   'Data' => $protein->start});
	       } else { 
		   $self->element({'Name' => 'Hsp_hit-from',
				   'Data' => $protein->start});
		   $self->element({'Name' => 'Hsp_hit-to',
				   'Data' => $protein->end});
	       }
	   }
	   $self->element({'Name' => 'Hsp_identity',
			   'Data' => 0});
	   $self->element({'Name' => 'Hsp_positive',
			   'Data' => 0});
	   $self->end_element({'Name' => 'Hsp'});
       }
       $self->end_element({'Name' => 'Hit'});
   }
   $self->end_element({'Name' => 'WiseOutput'});
   return $self->end_document();
}

=head2 start_element

 Title   : start_element
 Usage   : $eventgenerator->start_element
 Function: Handles a start element event
 Returns : none
 Args    : hashref with at least 2 keys 'Data' and 'Name'


=cut

sub start_element{
   my ($self,$data) = @_;
   # we currently don't care about attributes
   my $nm = $data->{'Name'};
   my $type = $MODEMAP{$nm};

   if( $type ) {
       if( $self->_eventHandler->will_handle($type) ) {
	   my $func = sprintf("start_%s",lc $type);
	   $self->_eventHandler->$func($data->{'Attributes'});
       }
       unshift @{$self->{'_elements'}}, $type;

       if($type eq 'result') {
	   $self->{'_values'} = {};
	   $self->{'_result'}= undef;
       }
   }

}

=head2 end_element

 Title   : start_element
 Usage   : $eventgenerator->end_element
 Function: Handles an end element event
 Returns : none
 Args    : hashref with at least 2 keys 'Data' and 'Name'


=cut

sub end_element {
    my ($self,$data) = @_;
    my $nm = $data->{'Name'};
    my $type = $MODEMAP{$nm};
    my $rc;

    if( $type = $MODEMAP{$nm} ) {
	if( $self->_eventHandler->will_handle($type) ) {
	    my $func = sprintf("end_%s",lc $type);
	    $rc = $self->_eventHandler->$func($self->{'_reporttype'},
					      $self->{'_values'});
	}
	shift @{$self->{'_elements'}};

    } elsif( $MAPPING{$nm} ) {

	if ( ref($MAPPING{$nm}) =~ /hash/i ) {
	    my $key = (keys %{$MAPPING{$nm}})[0];
	    $self->{'_values'}->{$key}->{$MAPPING{$nm}->{$key}} = $self->{'_last_data'};
	} else {
	    $self->{'_values'}->{$MAPPING{$nm}} = $self->{'_last_data'};
	}
    } else {
	$self->debug( "unknown nm $nm, ignoring\n");
    }
    $self->{'_last_data'} = ''; # remove read data if we are at
				# end of an element
    $self->{'_result'} = $rc if( defined $type && $type eq 'result' );
    return $rc;
}

=head2 element

 Title   : element
 Usage   : $eventhandler->element({'Name' => $name, 'Data' => $str});
 Function: Convience method that calls start_element, characters, end_element
 Returns : none
 Args    : Hash ref with the keys 'Name' and 'Data'


=cut

sub element{
   my ($self,$data) = @_;
   $self->start_element($data);
   $self->characters($data);
   $self->end_element($data);
}

=head2 characters

 Title   : characters
 Usage   : $eventgenerator->characters($str)
 Function: Send a character events
 Returns : none
 Args    : string


=cut

sub characters{
   my ($self,$data) = @_;

   return unless ( defined $data->{'Data'} && $data->{'Data'} !~ /^\s+$/ );

   $self->{'_last_data'} = $data->{'Data'};
}

=head2 within_element

 Title   : within_element
 Usage   : if( $eventgenerator->within_element($element) ) {}
 Function: Test if we are within a particular element
           This is different than 'in' because within can be tested
           for a whole block.
 Returns : boolean
 Args    : string element name


=cut

sub within_element{
   my ($self,$name) = @_;
   return 0 if ( ! defined $name &&
		 ! defined  $self->{'_elements'} ||
		 scalar @{$self->{'_elements'}} == 0) ;
   foreach (  @{$self->{'_elements'}} ) {
       if( $_ eq $name  ) {
	   return 1;
       }
   }
   return 0;
}


=head2 in_element

 Title   : in_element
 Usage   : if( $eventgenerator->in_element($element) ) {}
 Function: Test if we are in a particular element
           This is different than 'in' because within can be tested
           for a whole block.
 Returns : boolean
 Args    : string element name


=cut

sub in_element{
   my ($self,$name) = @_;
   return 0 if ! defined $self->{'_elements'}->[0];
   return ( $self->{'_elements'}->[0] eq $name)
}

=head2 start_document

 Title   : start_document
 Usage   : $eventgenerator->start_document
 Function: Handle a start document event
 Returns : none
 Args    : none


=cut

sub start_document{
    my ($self) = @_;
    $self->{'_lasttype'} = '';
    $self->{'_values'} = {};
    $self->{'_result'}= undef;
    $self->{'_elements'} = [];
    $self->{'_reporttype'} = 'exonerate';
}


=head2 end_document

 Title   : end_document
 Usage   : $eventgenerator->end_document
 Function: Handles an end document event
 Returns : Bio::Search::Result::ResultI object
 Args    : none


=cut

sub end_document{
   my ($self,@args) = @_;
   return $self->{'_result'};
}


sub write_result {
   my ($self, $blast, @args) = @_;

   if( not defined($self->writer) ) {
       $self->warn("Writer not defined. Using a $DEFAULT_WRITER_CLASS");
       $self->writer( $DEFAULT_WRITER_CLASS->new() );
   }
   $self->SUPER::write_result( $blast, @args );
}

sub result_count {
    my $self = shift;
    return $self->{'_result_count'};
}

sub report_count { shift->result_count }


=head2 wise

 Title   : wise
 Usage   : $obj->wise($newval)
 Function: Get/Set the Wise object parser
 Returns : value of wise (a scalar)
 Args    : on set, new value (a scalar or undef, optional)


=cut

sub wise{
    my $self = shift;
    return $self->{'wise'} = shift if @_;
    return $self->{'wise'};
}

=head2 wisetype

 Title   : wisetype
 Usage   : $obj->wisetype($newval)
 Function: Wise program type
 Returns : value of wisetype (a scalar)
 Args    : on set, new value (a scalar or undef, optional)


=cut

sub wisetype{
    my $self = shift;

    return $self->{'wisetype'} = shift if @_;
    return $self->{'wisetype'};
}

1;
