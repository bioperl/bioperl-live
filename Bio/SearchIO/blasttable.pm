# $Id$
#
# BioPerl module for Bio::SearchIO::blasttable
#
# Cared for by Jason Stajich <jason-at-bioperl-dot-org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::SearchIO::blasttable - Driver module for SearchIO for parsing NCBI -m 8/9 format

=head1 SYNOPSIS

  # do not use this module directly
  use Bio::SearchIO;
  my $parser = new Bio::SearchIO(-file   => $file,
                                 -format => 'blasttable');

  while( my $result = $parser->next_result ) {
  }

=head1 DESCRIPTION

This module will support parsing NCBI -m 8 or -m 9 tabular output.

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


package Bio::SearchIO::blasttable;
use vars qw(@ISA %MAPPING %MODEMAP $DEFAULT_WRITER_CLASS $DefaultProgramName);
use strict;
use Bio::SearchIO;
use Bio::Search::Result::ResultFactory;
use Bio::Search::Hit::HitFactory;
use Bio::Search::HSP::HSPFactory;

$DefaultProgramName = 'BLASTN';
$DEFAULT_WRITER_CLASS = 'Bio::Search::Writer::HitTableWriter';

# mapping of terms to Bioperl hash keys
%MODEMAP = (
	    'Result'             => 'result',
	    'Hit'                => 'hit',
	    'Hsp'                => 'hsp'
	    );

%MAPPING = ( 
	     'Hsp_bit-score'  => 'HSP-bits',
	     'Hsp_score'      => 'HSP-score',
	     'Hsp_evalue'     => 'HSP-evalue',
	     'Hsp_query-from' => 'HSP-query_start',
	     'Hsp_query-to'   => 'HSP-query_end',
	     'Hsp_hit-from'   => 'HSP-hit_start',
	     'Hsp_hit-to'     => 'HSP-hit_end',
	     'Hsp_positive'   => 'HSP-conserved',
	     'Hsp_identity'   => 'HSP-identical',
	     'Hsp_mismatches' => 'HSP-mismatches',
	     'Hsp_qgapblocks' => 'HSP-query_gapblocks',
	     'Hsp_hgapblocks' => 'HSP-hit_gapblocks',
	     'Hsp_gaps'       => 'HSP-hsp_gaps',
	     'Hsp_hitgaps'    => 'HSP-hit_gaps',
	     'Hsp_querygaps'  => 'HSP-query_gaps',
	     'Hsp_align-len'  => 'HSP-hsp_length',
	     'Hsp_query-frame'=> 'HSP-query_frame',
	     'Hsp_hit-frame'  => 'HSP-hit_frame',

	     'Hit_id'        => 'HIT-name',
	     'Hit_len'       => 'HIT-length',
	     'Hit_accession' => 'HIT-accession',
	     'Hit_def'       => 'HIT-description',
	     'Hit_signif'    => 'HIT-significance',
	     'Hit_score'     => 'HIT-score',
	     'Hit_bits'      => 'HIT-bits',

	     'Result_program'  => 'RESULT-algorithm_name',
	     'Result_version'  => 'RESULT-algorithm_version',
	     'Result_query-def'=> 'RESULT-query_name',
	     'Result_query-len'=> 'RESULT-query_length',
	     'Result_query-acc'=> 'RESULT-query_accession',
	     'Result_querydesc'=> 'RESULT-query_description',
	     'Result_db'       => 'RESULT-database_name',
	     'Result_db-len'   => 'RESULT-database_entries',
	     'Result_db-let'   => 'RESULT-database_letters',
	     );

@ISA = qw(Bio::SearchIO );

=head2 new

 Title   : new
 Usage   : my $obj = new Bio::SearchIO::blasttable();
 Function: Builds a new Bio::SearchIO::blasttable object 
 Returns : an instance of Bio::SearchIO::blasttable
 Args    :


=cut

sub _initialize {
    my ($self,@args) = @_;
    $self->SUPER::_initialize(@args);

    my ($pname) = $self->_rearrange([qw(PROGRAM_NAME)],
				    @args);
    $self->program_name($pname || $DefaultProgramName);
    $self->_eventHandler->register_factory('result', Bio::Search::Result::ResultFactory->new(-type => 'Bio::Search::Result::GenericResult'));
    $self->_eventHandler->register_factory('hit', Bio::Search::Hit::HitFactory->new(-type => 'Bio::Search::Hit::GenericHit'));
    $self->_eventHandler->register_factory('hsp', Bio::Search::HSP::HSPFactory->new(-type => 'Bio::Search::HSP::GenericHSP'));
}


=head2 next_result

 Title   : next_result
 Usage   : my $result = $parser->next_result
 Function: Parse the next result from the data stream
 Returns : L<Bio::Search::Result::ResultI>
 Args    : none


=cut

sub next_result{
   my ($self) = @_;
   my ($lastquery,$lasthit) = undef;
   local ($_);
   while( defined ($_ = $self->_readline) ) {
       next if /^\#/ || /^\s+$/;
       my ($qname,$hname, $percent_id, $hsp_len, $mismatches,$gapsm,
            $qstart,$qend,$hstart,$hend,$evalue,$bits) = split;
       
       # Remember Jim's code is 0 based
       if( defined $lastquery && 
	   $lastquery ne $qname ) {
	   $self->end_element({'Name' => 'Hit'});
	   $self->end_element({'Name' => 'Result'});
	   $self->_pushback($_);
	   return $self->end_document;
       } elsif( ! defined $lastquery ) {
	   $self->{'_result_count'}++;
	   $self->start_element({'Name' => 'Result'});
	   $self->element({'Name' => 'Result_program',
			   'Data' => $self->program_name});
	   $self->element({'Name' => 'Result_query-def',
			   'Data' => $qname});
	   $self->start_element({'Name' => 'Hit'});
	   $self->element({'Name' => 'Hit_id',
			   'Data' => $hname});
	   # we'll store the 1st hsp bits as the hit bits
	   $self->element({'Name' => 'Hit_bits',			   
			   'Data' => $bits});	   
           # we'll store the 1st hsp value as the hit evalue
	   $self->element({'Name' => 'Hit_signif',			   
			   'Data' => $evalue});
	   
       } elsif( $lasthit ne $hname ) {
	   if( $self->in_element('hit') ) {	       
	       $self->end_element({'Name' => 'Hit'});
	   }
	   $self->start_element({'Name' => 'Hit'});
	   $self->element({'Name' => 'Hit_id',
			   'Data' => $hname});
	   # we'll store the 1st hsp bits as the hit bits
	   $self->element({'Name' => 'Hit_bits',			   
			   'Data' => $bits});	   
           # we'll store the 1st hsp value as the hit evalue
	   $self->element({'Name' => 'Hit_signif',			   
			   'Data' => $evalue});
       }
       my $identical = $hsp_len - $mismatches - $gapsm;
       $self->start_element({'Name' => 'Hsp'});
       $self->element({'Name' => 'Hsp_evalue',			   
		       'Data' => $evalue});       
       $self->element({'Name' => 'Hsp_bit-score',
		       'Data' => $bits});
       $self->element({'Name' => 'Hsp_identity',
		       'Data' => $identical});
       $self->element({'Name' => 'Hsp_positive',
		       'Data' => $identical});
       $self->element({'Name' => 'Hsp_gaps',
		       'Data' => $gapsm});
       $self->element({'Name' => 'Hsp_query-from',
		       'Data' => $qstart});
       $self->element({'Name' => 'Hsp_query-to',
		       'Data' => $qend});

       $self->element({'Name' => 'Hsp_hit-from',
		       'Data' => $hstart });
       $self->element({'Name' => 'Hsp_hit-to',
		       'Data' => $hend });
       $self->element({'Name' => 'Hsp_align-len',
		       'Data' => $hsp_len});
       $self->end_element({'Name' => 'Hsp'});
       $lastquery = $qname;
       $lasthit   = $hname;
   }
   # fencepost
   if( defined $lasthit && defined $lastquery ) {
       if( $self->in_element('hit') ) {
	   $self->end_element({'Name' => 'Hit'});
       }
       $self->end_element({'Name' => 'Result'});
       return $self->end_document;
   }
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
   if( my $type = $MODEMAP{$nm} ) {
	$self->_mode($type);
	if( $self->_will_handle($type) ) {
	    my $func = sprintf("start_%s",lc $type);
	    $self->_eventHandler->$func($data->{'Attributes'});
	}						 
	unshift @{$self->{'_elements'}}, $type;
    }
    if($nm eq 'Result') {
	$self->{'_values'} = {};
	$self->{'_result'}= undef;
	$self->{'_mode'} = '';
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
    my $rc;
    # Hsp are sort of weird, in that they end when another
    # object begins so have to detect this in end_element for now
        
    if( my $type = $MODEMAP{$nm} ) {
	if( $self->_will_handle($type) ) {
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
	$self->warn( "unknown nm $nm ignoring\n");
    }
    $self->{'_last_data'} = ''; # remove read data if we are at 
				# end of an element
    $self->{'_result'} = $rc if( $nm eq 'Result' );
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

   return unless ( defined $data->{'Data'} );
   if( $data->{'Data'} =~ /^\s+$/ ) {
       return unless $data->{'Name'} =~ /Hsp\_(midline|qseq|hseq)/;
   }

   if( $self->in_element('hsp') && 
       $data->{'Name'} =~ /Hsp\_(qseq|hseq|midline)/ ) {
       
       $self->{'_last_hspdata'}->{$data->{'Name'}} .= $data->{'Data'};
   }  
   
   $self->{'_last_data'} = $data->{'Data'}; 
}

=head2 _mode

 Title   : _mode
 Usage   : $obj->_mode($newval)
 Function: 
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
 Function: Handles a start document event
 Returns : none
 Args    : none


=cut

sub start_document{
    my ($self) = @_;
    $self->{'_lasttype'} = '';
    $self->{'_values'} = {};
    $self->{'_result'}= undef;
    $self->{'_mode'} = '';
    $self->{'_elements'} = [];
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

=head2 result_count

 Title   : result_count
 Usage   : my $count = $searchio->result_count
 Function: Returns the number of results we have processed
 Returns : integer
 Args    : none


=cut

sub result_count {
    my $self = shift;
    return $self->{'_result_count'};
}

sub report_count { shift->result_count }


=head2 program_name

 Title   : program_name
 Usage   : $obj->program_name($newval)
 Function: Get/Set the program name
 Returns : value of program_name (a scalar)
 Args    : on set, new value (a scalar or undef, optional)


=cut

sub program_name{
    my $self = shift;

    $self->{'program_name'} = shift if @_;
    return $self->{'program_name'} || $DefaultProgramName;
}


=head2 _will_handle

 Title   : _will_handle
 Usage   : Private method. For internal use only.
              if( $self->_will_handle($type) ) { ... }
 Function: Provides an optimized way to check whether or not an element of a 
           given type is to be handled.
 Returns : Reference to EventHandler object if the element type is to be handled.
           undef if the element type is not to be handled.
 Args    : string containing type of element.

Optimizations:

=over 2

=item 1

Using the cached pointer to the EventHandler to minimize repeated
lookups.

=item 2

Caching the will_handle status for each type that is encountered so
that it only need be checked by calling
handler-E<gt>will_handle($type) once.

=back

This does not lead to a major savings by itself (only 5-10%).  In
combination with other optimizations, or for large parse jobs, the
savings good be significant.

To test against the unoptimized version, remove the parentheses from
around the third term in the ternary " ? : " operator and add two
calls to $self-E<gt>_eventHandler().

=cut

sub _will_handle {
    my ($self,$type) = @_;
    my $handler = $self->{'_handler'};
    my $will_handle = defined($self->{'_will_handle_cache'}->{$type})
                             ? $self->{'_will_handle_cache'}->{$type}
                             : ($self->{'_will_handle_cache'}->{$type} =
                               $handler->will_handle($type));

    return $will_handle ? $handler : undef;
}

1;
