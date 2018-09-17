#
# BioPerl module for Bio::SearchIO::megablast
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

Bio::SearchIO::megablast - a driver module for Bio::SearchIO to parse
megablast reports (format 0)

=head1 SYNOPSIS

# do not use this module directly

  use Bio::SearchIO;
  # for default format output from megablast
  my $in = Bio::SearchIO->new(-file   => 'file.mbl',
                             -format => 'megablast',
                             -report_format => 0);

  while( my $r = $in->next_result ) {
    while( my $hit = $r->next_hit ) {
      while( my $hsp = $hit->next_hsp ) {
      }
    }
  }

=head1 DESCRIPTION

Beware!

Because of the way megablast report format 0 is coded, realize that score
means # gap characters + # mismatches for a HSP.

The docs from NCBI regarding FORMAT 0
#   0: Produce one-line output for each alignment, in the form
#
#   'subject-id'=='[+-]query-id' (s_off q_off s_end q_end) score
#
#   Here subject(query)-id is a gi number, an accession or some other type of
#   identifier found in the FASTA definition line of the respective sequence.
#
#   + or - corresponds to same or different strand alignment.
#
#   Score for non-affine gapping parameters means the total number of
#   differences (mismatches + gap characters). For affine case it is the
#   actual (raw) score of the alignment.

FORMAT 1 parsing has not been implemented
FORMAT 2 parsing should work with the SearchIO 'blast' parser

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

Email jason-at-bioperl.org

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::SearchIO::megablast;
use strict;
use vars qw(%MAPPING %MODEMAP $DEFAULT_BLAST_WRITER_CLASS);

use base qw(Bio::SearchIO);

BEGIN {
    # mapping of MegaBlast terms to Bioperl hash keys
    %MODEMAP = ('MegaBlastOutput' => 'result',
		'Hit'         => 'hit',
		'Hsp'         => 'hsp'
		);

    # This should really be done more intelligently, like with
    # XSLT

    %MAPPING =
	(
	  'Hsp_query-from' => 'HSP-query_start',
	  'Hsp_query-to'   => 'HSP-query_end',
	  'Hsp_hit-from'   => 'HSP-hit_start',
	  'Hsp_hit-to'     => 'HSP-hit_end',
	  'Hit_score'      => 'HIT-score',
	  'Hsp_score'      => 'HSP-score',
	
	  'Hsp_identity'   => 'HSP-identical',
	  'Hsp_positive'   => 'HSP-conserved',

	  'Hit_id'         => 'HIT-name',
	
	  'MegaBlastOutput_program'  => 'RESULT-algorithm_name',
	  'MegaBlastOutput_query-def'=> 'RESULT-query_name',
	  );


    $DEFAULT_BLAST_WRITER_CLASS = 'Bio::SearchIO::Writer::HitTableWriter';
}


=head2 new

 Title   : new
 Usage   : my $obj = Bio::SearchIO::blast->new();
 Function: Builds a new Bio::SearchIO::blast object
 Returns : Bio::SearchIO::blast
 Args    : -fh/-file => filehandle/filename to BLAST file
           -format   => 'blast'

=cut

sub _initialize {
    my ($self,@args) = @_;
    $self->SUPER::_initialize(@args);
    my ($fmt) = $self->_rearrange([qw(REPORT_FORMAT)], @args);

    $self->throw("Must provide a value for -report_format when initializing a megablast parser") unless defined $fmt ;
    $self->report_format($fmt);
    return 1;
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

   my $fmt = $self->report_format;
   my ($lastquery,$lasthit);
   while( defined($_ = $self->_readline) ) {
       if( $fmt == 0 ) {
	   if( /^\'(\S+)\'\=\=\'(\+|\-)(\S+)\'\s+
	       \((\d+)\s+(\d+)\s+(\d+)\s+(\d+)\)\s+
	       (\d+)/ox )
	   {
	       my ($hit,$strand,$query,
		   $h_start,$q_start,$h_end,$q_end,
		   $score) = ($1,$2,$3,$4,$5,$6,$7,$8);
	       if( ! defined $lastquery ) {
		   $self->start_element({'Name' => 'MegaBlastOutput'});
		   $self->element({'Name' => 'MegaBlastOutput_program',
				   'Data' => 'MEGABLAST'});
		   $self->element({'Name' => 'MegaBlastOutput_query-def',
				   'Data' => $query});
	       } elsif( $lastquery ne $query ) {
		   $self->_pushback($_);
		   $self->end_element({'Name' => 'Hit'}) if( defined $lasthit);
		   $self->end_element({ 'Name' => 'MegaBlastOutput'});
		   $lasthit = undef;
		   $lastquery = undef;
		   return $self->end_document();
	       }

	       if( ! defined $lasthit || $lasthit ne $hit  ) {
		   $self->end_element({'Name' => 'Hit'}) if( defined $lasthit);
		   $self->start_element({'Name' => 'Hit'});
		   $self->element({'Name' => 'Hit_id',
				   'Data' => $hit});
	       }
	       $self->start_element({'Name' => 'Hsp'});
	       $self->element({'Name' => 'Hsp_score',
			       'Data' => $score});

	       # flip flop start/end if strand is < 0
	       # since strandedness is inferred from the query
	       # because of the way it is coded all queries will
	       # be on the forward strand and hits will be either
	       # +/-

	       # also the NCBI docs state:
#   0: Produce one-line output for each alignment, in the form
#
#   'subject-id'=='[+-]query-id' (s_off q_off s_end q_end) score
#
#   Here subject(query)-id is a gi number, an accession or some other type of
#   identifier found in the FASTA definition line of the respective sequence.
#
#   + or - corresponds to same or different strand alignment.
#
#   Score for non-affine gapping parameters means the total number of
#   differences (mismatches + gap characters). For affine case it is the
#   actual (raw) score of the alignment.
	
	       # and yet when rev strand hits are made I see
	       # (MBL 2.2.4)
	       # 'Contig634'=='-503384' (1 7941 321 7620) 19
	       # so the query is on the rev strand and the
	       # subject is on the fwd strand
	       # so I am flip-flopping everything when I see a '-'
	       if( $strand eq '-' ) {
		   ($h_start,$h_end) = ( $h_end,$h_start);
		   ($q_start,$q_end) = ( $q_end,$q_start);
	       }
	       $self->element({'Name' => 'Hsp_hit-from',
			       'Data' => $h_start});
	       $self->element({'Name' => 'Hsp_hit-to',
			       'Data' => $h_end});
	       $self->element({'Name' => 'Hsp_query-from',
			       'Data' => $q_start});
	       $self->element({'Name' => 'Hsp_query-to',
			       'Data' => $q_end});
	       # might not be quite right -- need to know length of the HSP
	       my $numid = (abs($q_end - $q_start) - $score);

	       $self->element({'Name' => 'Hsp_identity',
			       'Data' => $numid});
	       $self->element({'Name' => 'Hsp_positive',
			       'Data' => $numid});

	       $self->end_element({'Name' => 'Hsp'});
	       $lasthit   = $hit;
	       $lastquery = $query;
	   } else {
	       $self->debug("Unknown line in fmt0 parsing: $_");
	   }
       }
   }
   if( defined $lastquery && $fmt == 0 ) {
       $self->end_element({'Name' => 'Hit'}) if( defined $lasthit);
       $self->end_element({ 'Name' => 'MegaBlastOutput'});
       return $self->end_document();
   }
   return;
}

=head2 report_format

 Title   : report_format
 Usage   : $obj->report_format($newval)
 Function: Get/Set the report_format value
 Returns : value of report_format (a scalar)
 Args    : on set, new value (a scalar or undef, optional)


=cut

sub report_format{
    my $self = shift;
    return $self->{'_report_format'} = shift if @_;
    return $self->{'_report_format'};
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
    # we currently do not care about attributes
    my $nm = $data->{'Name'};
   if( my $type = $MODEMAP{$nm} ) {
	$self->_mode($type);
	if( $self->_eventHandler->will_handle($type) ) {
	    my $func = sprintf("start_%s",lc $type);
	    $self->_eventHandler->$func($data->{'Attributes'});
	}
	unshift @{$self->{'_elements'}}, $type;
    }

    if($nm eq 'MegaBlastOutput') {
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

    if( my $type = $MODEMAP{$nm} ) {
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
	$self->warn( "unknown nm $nm ignoring\n");
    }
    $self->{'_last_data'} = ''; # remove read data if we are at
				# end of an element
    $self->{'_result'} = $rc if( $nm eq 'MegaBlastOutput' );
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
   return unless defined $data->{'Data'};
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

1;
