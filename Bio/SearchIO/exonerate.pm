# $Id$
#
# BioPerl module for Bio::SearchIO::exonerate
#
# Cared for by Jason Stajich <jason@bioperl.org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::SearchIO::exonerate - parser for Exonerate 

=head1 SYNOPSIS

# do not use this module directly, it is a driver for SearchIO

use Bio::SearchIO;
my $searchio = new Bio::SearchIO(-file => 'file.exonerate',
                                 -format => 'exonerate');

while( my $r = $searchio->next_result ) {
  print $r->query_name, "\n";
}

=head1 DESCRIPTION

This is a driver for the SearchIO system for parsing Exonerate (Guy
Slater) output.  You can get Exonerate at
http://cvsweb.sanger.ac.uk/cgi-bin/cvsweb.cgi/exonerate/?cvsroot=Ensembl
[until Guy puts up a Web reference,publication for it.]).

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


package Bio::SearchIO::exonerate;
use strict;
use vars qw(@ISA @STATES %MAPPING %MODEMAP $DEFAULT_WRITER_CLASS);
use Bio::SearchIO;

@ISA = qw(Bio::SearchIO );

use POSIX;

BEGIN { 
    %MODEMAP = ('ExonerateOutput' => 'result',
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
	  'Hsp_identical' =>  'HSP-identical',
	  
	  'Hit_id'        => 'HIT-name',
	  'Hit_desc'      => 'HIT-description',
	  
	  'ExonerateOutput_program'  => 'RESULT-algorithm_name',
	  'ExonerateOutput_query-def'=> 'RESULT-query_name',
	  'ExonerateOutput_query-desc'=> 'RESULT-query_description',
	  
	  );
    $DEFAULT_WRITER_CLASS = 'Bio::Search::Writer::HitTableWriter';
}

=head2 new

 Title   : new
 Usage   : my $obj = new Bio::SearchIO::exonerate();
 Function: Builds a new Bio::SearchIO::exonerate object 
 Returns : an instance of Bio::SearchIO::exonerate
 Args    :


=cut

=head2 next_result

 Title   : next_result
 Usage   : my $hit = $searchio->next_result;
 Function: Returns the next Result from a search
 Returns : Bio::Search::Result::ResultI object
 Args    : none

=cut

sub next_result{
   my ($self) = @_;
   $self->{'_last_data'} = '';
   my ($reporttype,$seenquery,$reportline);
   $self->start_document();
   my @hit_signifs;
   my $seentop;
   while( defined($_ = $self->_readline) ) {  
       if( /^Query:\s+(\S+)(\s+(.+))?/ ) {
	   $seentop = 1;
	   if( $self->within_element('hsp') ) {
	       $self->end_element({'Name' => 'Hsp'});
	       $self->end_element({'Name' => 'Hit'});
	       $self->end_element({'Name' => 'ExonerateOutput'});
	       $self->_pushback($_);
	       return $self->end_document();
	   }
	   my ($nm,$desc) = ($1,$2);
	   chomp($desc) if defined $desc;
	   $self->{'_result_count'}++;
	   $self->start_element({'Name' => 'ExonerateOutput'});
	   $self->element({'Name' => 'ExonerateOutput_query-def',
			   'Data' => $nm });
	   $self->element({'Name' => 'ExonerateOutput_query-desc',
			   'Data' => $desc });
	   $self->element({'Name' => 'ExonerateOutput_program', 
			    'Data' => 'Exonerate' });
	   	   
       } elsif ( /^Target:\s+(\S+)(\s+(.+))?/ ) {
	    my ($nm,$desc) = ($1,$2);
	   chomp($desc) if defined $desc;
	   $self->start_element({'Name' => 'Hit'});
	   $self->element({'Name' => 'Hit_id',
			   'Data' => $nm});
	   $self->element({'Name' => 'Hit_desc',
			   'Data' => $desc});
	   $self->start_element({'Name' => 'Hsp'});
       } elsif( /^(\s+(\d+)\s+\:\s+)(.+)\s+\:\s+(\d+)\s*$/ ) {
	   my ($qstart,$qmatch,$qend) = ($2,$3,$4);
#	   Hsp_qseq Hsp_midline Hsp_hseq
	   $self->element({'Name'  => 'Hsp_qseq',
			   'Data'  => $qmatch});
	   
	   unless( defined( $_ = $self->_readline) ) {
	       $self->warn("Unexpected ending to Exonerate file");
	   }
	   my $midline = substr($_,length($1),length($qmatch));
	   $self->element({'Name'  => 'Hsp_midline',
			   'Data'  => $midline});
	   unless( defined( $_ = $self->_readline) &&
		   /^\s+(\d+)\s+\:\s+(.+)\s+\:\s+(\d+)\s*$/ ) {
	       $self->warn("Unexpected ending to Exonerate file");
	   }
	   my ($hstart,$hmatch,$hend) = ( $1,$2,$3 );
	   $self->element({'Name'  => 'Hsp_hseq',
			   'Data'  => $hmatch})
       } elsif(  s/^cigar:\s+(\S+)\s+          # query sequence id
		 (\d+)\s+(\d+)\s+([\-\+])\s+   # query start-end-strand
		 (\S+)\s+                      # target sequence id
		 (\d+)\s+(\d+)\s+([\-\+])\s+   # target start-end-strand
		 (\d+)\s+                      # score
		 //ox ) {	   
	   my ($qs,$qe,$qstrand) = ($2,$3,$4);
	   my ($hs,$he,$hstrand) = ($6,$7,$8);
	   
	   if( $qstrand eq '-' ) {
	       ($qs,$qe) = ($qe,$qs); # flip-flop if we're on opp strand
	   }
	   if( $hstrand eq '-' ) {
	       ($hs,$he) = ($he,$hs); # flip-flop if we're on opp strand
	   }
	   
	   $self->element({'Name' => 'Hsp_query-from',
			   'Data' => $qs+1});
	   $self->element({'Name' => 'Hsp_query-to',
			   'Data' => $qe+1});
	   
	   $self->element({'Name' => 'Hsp_hit-from',
			   'Data' => $hs+1});
	   $self->element({'Name' => 'Hsp_hit-to',
			   'Data' => $he+1});

	   $self->element({'Name' => 'Hsp_score',
			   'Data' => $9});
	   my (@rest) = split;
	   my ($state,$len) = ( shift @rest, shift @rest);
	   $self->element({'Name' => 'Hsp_align-len',
			   'Data' => $len});	   
	   $self->element({'Name' => 'Hsp_qlength',
			   'Data' => abs($3-$2)});	   
	   $self->element({'Name' => 'Hsp_hlength',
			   'Data' => abs($7-$2)});
	   my ($i,$lastindex,$count) = (0,-1,0);
	   while( ($i = index ($self->{'_last_hspdata'}->{'Hsp_midline'},
			       '|',$lastindex+1)) >= $lastindex ) {
	       $count++;
	       $lastindex = $i;
	   }
	   $self->element({'Name' => 'Hsp_identical',
			   'Data' => $count});	   
	   $self->end_element({'Name' => 'Hsp'});	   
	   $self->end_element({'Name' => 'Hit'});
	   $self->end_element({'Name' => 'ExonerateOutput'});
	   return $self->end_document();	  
       } else { 
	   
       }       
   }
   return $self->end_document() if( $seentop );
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
    # Hsp are sort of weird, in that they end when another
    # object begins so have to detect this in end_element for now
    if( $nm eq 'Hsp' ) {
	foreach ( qw(Hsp_qseq Hsp_midline Hsp_hseq) ) {
	    $self->element({'Name' => $_,
			    'Data' => $self->{'_last_hspdata'}->{$_}});
	}
	$self->{'_last_hspdata'} = {};
    }
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

   if( $self->in_element('hsp') && 
       $data->{'Name'} =~ /Hsp\_(qseq|hseq|midline)/ ) {       
       $self->{'_last_hspdata'}->{$data->{'Name'}} .= $data->{'Data'};
   }  
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

1;
