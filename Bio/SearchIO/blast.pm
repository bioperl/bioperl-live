# $Id$
#
# BioPerl module for Bio::SearchIO::blast
#
# Cared for by Jason Stajich <jason@bioperl.org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::SearchIO::blast - Event generator for event based parsing of blast reports 

=head1 SYNOPSIS

   # Do not use this object directly - it is used as part of the
   # Bio::SearchIO system.

    use Bio::SearchIO;
    my $searchio = new Bio::SearchIO(-format => 'blast',
				     -file   => 't/data/ecolitst.bls');
    while( my $result = $searchio->next_result ) {
	while( my $hit = $result->next_hit ) {
	    while( my $hsp = $hit->next_hsp ) {
		# ...
	    }
	}
    }

=head1 DESCRIPTION

This object encapsulated the necessary methods for generating events
suitable for building Bio::Search objects from a BLAST report file.

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



package Bio::SearchIO::blast;
use strict;
use vars qw(@ISA %MAPPING %MODEMAP);
use Bio::SearchIO;

@ISA = qw(Bio::SearchIO );

BEGIN { 
    # mapping of NCBI Blast terms to Bioperl hash keys
    %MODEMAP = ('BlastOutput' => 'result',
		'Hit'         => 'hit',
		'Hsp'         => 'hsp'
		);

    # This should really be done more intelligently, like with
    # XSLT

    %MAPPING = ( 
		 'Hsp_bit-score'  => 'bits',
		 'Hsp_score'      => 'score',
		 'Hsp_evalue'     => 'evalue',
		 'Hsp_pvalue'     => 'pvalue',
		 'Hsp_query-from' => 'querystart',
		 'Hsp_query-to'   => 'queryend',
		 'Hsp_hit-from'   => 'hitstart',
		 'Hsp_hit-to'     => 'hitend',
		 'Hsp_positive'   => 'conserved',
		 'Hsp_identity'   => 'identical',
		 'Hsp_gaps'      => 'gaps',
		 'Hsp_hitgaps'   => 'hitgaps',
		 'Hsp_querygaps' => 'querygaps',
		 'Hsp_qseq'       => 'queryseq',
		 'Hsp_hseq'       => 'hitseq',
		 'Hsp_midline'    => 'homolseq',
		 'Hsp_align-len'  => 'hsplen',
		 'Hsp_query-frame'=> 'queryframe',
		 'Hsp_hit-frame'  => 'hitframe',

		 'Hit_id'        => 'hitname',
		 'Hit_len'       => 'hitlen',
		 'Hit_accession' => 'hitacc',
		 'Hit_def'       => 'hitdesc',
		 'Hit_signif'    => 'hitsignif',
		 'Hit_score'     => 'hitscore',
		 
		 'BlastOutput_program'  => 'programname',
		 'BlastOutput_version'  => 'programver',
		 'BlastOutput_query-def'=> 'queryname',
		 'BlastOutput_query-len'=> 'querylen',
		 'BlastOutput_query-acc'=> 'queryacc',
		 'BlastOutput_db'       => 'dbname',
		 'BlastOutput_db-len'   => 'dbsize',
		 'BlastOutput_db-let'   => 'dblets',
		 'BlastOutput_query-acc'=> 'queryacc',

		 'Iteration_iter-num'   => 'iternum',
		 'Parameters_matrix'    => { 'param' => 'matrix'},
		 'Parameters_expect'    => { 'param' => 'expect'},
		 'Parameters_include'   => { 'param' => 'include'},
		 'Parameters_sc-match'  => { 'param' => 'match'},
		 'Parameters_sc-mismatch' => { 'param' => 'mismatch'},
		 'Parameters_gap-open'  =>   { 'param' => 'gapopen'},
		 'Parameters_gap-extend'=>   { 'param' => 'gapext'},
		 'Parameters_filter'    =>  {'param' => 'filter'},
		 'Statistics_db-len'    => {'stat' => 'dbentries'},
		 'Statistics_db-let'    => { 'stat' => 'dbletters'},
		 'Statistics_hsp-len'   => { 'stat' => 'hsplength'},
		 'Statistics_query-len'   => { 'stat' => 'querylength'},
		 'Statistics_eff-space' => { 'stat' => 'effectivespace'},
		 'Statistics_eff-spaceused' => { 'stat' => 'effectivespaceused'},
		 'Statistics_eff-dblen'    => { 'stat' => 'effectivedblength'},
		 'Statistics_kappa'     => { 'stat' => 'kappa' },
		 'Statistics_lambda'    => { 'stat' => 'lambda' },
		 'Statistics_entropy'   => { 'stat' => 'entropy'},
		 'Statistics_framewindow'=> { 'stat' => 'frameshiftwindow'},
		 'Statistics_decay'=> { 'stat' => 'decayconst'},

		 'Statistics_T'=> { 'stat' => 'T'},
		 'Statistics_A'=> { 'stat' => 'A'},
		 'Statistics_X1'=> { 'stat' => 'X1'},
		 'Statistics_X2'=> { 'stat' => 'X2'},
		 'Statistics_S1'=> { 'stat' => 'S1'},
		 'Statistics_S2'=> { 'stat' => 'S2'},
		 
		 );
}


=head2 new

 Title   : new
 Usage   : my $obj = new Bio::SearchIO::blast();
 Function: Builds a new Bio::SearchIO::blast object 
 Returns : Bio::SearchIO::blast
 Args    :

=cut

sub new {
  my($class,@args) = @_;

  my $self = $class->SUPER::new(@args);
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
   
   my $data = '';
   my $seentop = 0;
   my $reporttype;
   $self->start_document();
   my @hit_signifs;
   while( defined ($_ = $self->_readline )) {
       next if( /^\s+$/); # skip empty lines
       if( /^([T]?BLAST[NPX])\s*(\S+)/i ) {
	   if( $seentop ) {
	       $self->_pushback($_);
	       $self->end_element({ 'Name' => 'BlastOutput'});
	       return $self->end_document();
	   }
	   $self->start_element({ 'Name' => 'BlastOutput' } );
	   $seentop = 1;
	   $reporttype = $1;
	   $self->element({ 'Name' => 'BlastOutput_program',
			    'Data' => $reporttype});
	   
	   $self->element({ 'Name' => 'BlastOutput_version',
			    'Data' => $2});
       } elsif ( /^Query=\s*(.+)$/ ) {
	   my $q = $1;
	   my $size = 0;      
	   $_ = $self->_readline;
	   while( defined ($_) && $_ !~ /^\s+$/ ) {	
	       chomp;
	       if( /\(([\d,]+)\s+letters\)/ ) {		   
		   $size = $1;
		   $size =~ s/,//g;
		   last;
	       } else { 
		   $q .= " $_";
		   $q =~ s/ +/ /g;
		   $q =~ s/^ | $//g;
	       }
	       $_ = $self->_readline;
	   }
	   chomp($q);
	   $self->element({ 'Name' => 'BlastOutput_query-def',
			    'Data' => $q});
	   $self->element({ 'Name' => 'BlastOutput_query-len', 
			    'Data' => $size});
	   my ($firstpart) = split(/\s+/,$q);
	   my @pieces = split(/\|/,$firstpart);
	   my $acc = pop @pieces;
	   $self->element({ 'Name' =>  'BlastOutput_query-acc',
			    'Data'  => $acc});	   
       } elsif( /Sequences producing significant alignments:/ ) {
	   # skip the next whitespace line
	   $_ = $self->_readline();
	   while( defined ($_ = $self->_readline() ) && 
		  ! /^\s+$/ ) {	       
	       my @line = split;
	       push @hit_signifs, [ pop @line, pop @line];
	   }
       } elsif( /Sequences producing High-scoring Segment Pairs:/ ) {
	   # skip the next line
	   $_ = $self->_readline();
	   
	    while( defined ($_ = $self->_readline() ) && 
		  ! /^\s+$/ ) {	
	       my @line = split;
	       pop @line; # throw away first number which is for 'N'col
	       push @hit_signifs, [ pop @line, pop @line];
	   }
       } elsif ( /^Database:\s*(.+)$/ ) {
	   my $db = $1;

	   while( defined($_ = $self->_readline) ) {
	       if( /^\s+([\d\,]+)\s+sequences\;\s+([\d,]+)\s+total\s+letters/){
		   my ($s,$l) = ($1,$2);
		   $s =~ s/,//g;
		   $l =~ s/,//g;
		   $self->element({'Name' => 'BlastOutput_db-len',
				   'Data' => $s});	       
		   $self->element({'Name' => 'BlastOutput_db-let',
				   'Data' => $l});
		   last;
	       } else {
		   chomp;
		   $db .= $_;
	       }
	   }	       
	   $self->element({'Name' => 'BlastOutput_db',
			   'Data' => $db});
       } elsif( /^>(\S+)\s*(.*)?/ ) {
	   $self->in_element('hsp') && $self->end_element({ 'Name' => 'Hsp'});
	   $self->in_element('hit') && $self->end_element({ 'Name' => 'Hit'});
	   
	   $self->start_element({ 'Name' => 'Hit'});
	   my $id = $1;	  
	   $self->element({ 'Name' => 'Hit_id',
			    'Data' => $id});
	   my @pieces = split(/\|/,$id);
	   my $acc = pop @pieces;
	   $self->element({ 'Name' =>  'Hit_accession',
			    'Data'  => $acc});	   
	   $self->element({ 'Name' => 'Hit_def',
			    'Data' => $2});
	   my $v = shift @hit_signifs;
	   if( defined $v ) {
	       $self->element({'Name' => 'Hit_signif',
			       'Data' => $v->[0]});
	       $self->element({'Name' => 'Hit_score',
			       'Data' => $v->[1]});
	   }

      } elsif( /\s+(Plus|Minus) Strand HSPs:/i ) {
	   next;
       } elsif(  /Length\s*=\s*([\d,]+)/ ) {
	   $self->element({ 'Name' => 'Hit_len',
			    'Data' => $1 });
       } elsif( ($self->in_element('hit') || $self->in_element('hsp')) && # wublast
		/Score\s*=\s*(\d+)\s*\(([\d\.]+)\s*bits\),\s*Expect\s*=\s*([^,\s]+),\s*P\s*=\s*([^,\s]+)/ ) {
	   $self->in_element('hsp') && $self->end_element({'Name' => 'Hsp'});
	   $self->start_element({'Name' => 'Hsp'});
       	   $self->element( { 'Name' => 'Hsp_score',
			     'Data' => $1});
	   $self->element( { 'Name' => 'Hsp_bit-score',
			     'Data' => $2});
	   $self->element( { 'Name' => 'Hsp_evalue',
			     'Data' => $3});
	   $self->element( {'Name'  => 'Hsp_pvalue',
			    'Data'  =>$4});
       } elsif( ($self->in_element('hit') || $self->in_element('hsp')) && # ncbi blast
		/Score\s*=\s*(\S+)\s*bits\s*\((\d+)\),\s*Expect(\(\d+\))?\s*=\s*(\S+)/) {
	   $self->in_element('hsp') && $self->end_element({ 'Name' => 'Hsp'});
	   
	   $self->start_element({'Name' => 'Hsp'});
	   $self->element( { 'Name' => 'Hsp_score',
			     'Data' => $2});
	   $self->element( { 'Name' => 'Hsp_bit-score',
			     'Data' => $1});
	   $self->element( { 'Name' => 'Hsp_evalue',
			     'Data' => $4});
       } elsif( $self->in_element('hsp') &&
		/Identities\s*=\s*(\d+)\s*\/\s*(\d+)\s*[\d\%\(\)]+\s*(,\s*Positives\s*=\s*(\d+)\/(\d+)\s*[\d\%\(\)]+\s*)?(\,\s*Gaps\s*=\s*(\d+)\/(\d+))?/i ) {
	   
	   $self->element( { 'Name' => 'Hsp_identity',
			     'Data' => $1});
	   $self->element( {'Name' => 'Hsp_align-len',
			    'Data' => $2});
	   if( defined $3 ) {
	       $self->element( { 'Name' => 'Hsp_positive',
				 'Data' => $4});
	   }
	   if( defined $6 ) { 
	       $self->element( { 'Name' => 'Hsp_gaps',
				 'Data' => $7});	   
	   } else { 
	       $self->element( { 'Name' => 'Hsp_gaps',
				 'Data' => 0});
	   }
	   $self->{'_Query'} = {'begin' => 0, 'end' => 0};
	   $self->{'_Sbjct'} = { 'begin' => 0, 'end' => 0};
       } elsif( $self->in_element('hsp') &&
		/Strand\s*=\s*(Plus|Minus)\s*\/\s*(Plus|Minus)/i ) {
	   next;
       } elsif( $self->in_element('hsp') &&
		/Frame\s*=\s*([\+\-][1-3])\s*(\/\s*([\+\-][1-3]))?/ ){

	   my ($queryframe,$hitframe);
	   if( $reporttype eq 'TBLASTX' ) {
	       ($queryframe,$hitframe) = ($1,$2);
	       $hitframe =~ s/\/\s*//g;
	   } elsif( $reporttype eq 'TBLASTN' ) {
	       ($hitframe,$queryframe) = ($1,0);
	   } else { 
	       ($queryframe,$hitframe) = ($1,0);
	   }
	   $self->element({'Name' => 'Hsp_query-frame',
			   'Data' => $queryframe});
	   	   
	   $self->element({'Name' => 'Hsp_hit-frame',
			   'Data' => $hitframe});
       } elsif(  /^Parameters:/ || /\s+Database:/ ) {
	   $self->in_element('hsp') && $self->end_element({'Name' => 'Hsp'});
	   $self->in_element('hit') && $self->end_element({'Name' => 'Hit'});
	   my $blast = ( /Parameters/ ) ? 'wublast' : 'ncbi'; 
	   my $last = '';
	   while( defined ($_ = $self->_readline ) ) {
	       if( /^([T]?BLAST[NPX])\s*([\d\.]+)/i ) {
		   $self->_pushback($_);
		   # let's handle this in the loop
		   last;
	       }
	       # here is where difference between wublast and ncbiblast
	       # is better handled by different logic
	       if( /Number of Sequences:\s+([\d\,]+)/i ||
			/of sequences in database:\s+([\d,]+)/i) {
		   my $c = $1;
		   $c =~ s/\,//g;
		   $self->element({'Name' => 'Statistics_db-len',
				   'Data' => $c});
	       } elsif ( /letters in database:\s+([\d,]+)/i) {	   
		   my $s = $1;
		   $s =~ s/,//g;
		   $self->element({'Name' => 'Statistics_db-let',
				   'Data' => $s});
	       } elsif( $blast eq 'wublast' ) {
		   if( /E=(\S+)/ ) {
		       $self->element({'Name' => 'Parameters_expect',
				       'Data' => $1});
		   } elsif( $last =~ /Frame\s+MatID\s+Matrix name/i ) {
		       s/^\s+//;
                       #throw away first two slots
		       my @vals = split;
		       splice(@vals, 0,2);
		       my ($matrix,$lambda,$kappa,$entropy) = @vals;
		       $self->element({'Name' => 'Parameters_matrix',
				       'Data' => $matrix});
		       $self->element({'Name' => 'Statistics_lambda',
				       'Data' => $lambda});
		       $self->element({'Name' => 'Statistics_kappa',
				       'Data' => $kappa});
		       $self->element({'Name' => 'Statistics_entropy',
				       'Data' => $entropy});
		   } 
	       } elsif ( $blast eq 'ncbi' ) {
		   if( /^Matrix:\s+(\S+)/i ) {
		       $self->element({'Name' => 'Parameters_matrix',
				       'Data' => $1});		       
		   } elsif( /Lambda/ ) {
		       $_ = $self->_readline;
		       s/^\s+//;
		       my ($lambda, $kappa, $entropy) = split;
		       $self->element({'Name' => 'Statistics_lambda',
				       'Data' => $lambda});
		       $self->element({'Name' => 'Statistics_kappa',
				       'Data' => $kappa});
		       $self->element({'Name' => 'Statistics_entropy',
				       'Data' => $entropy});
		   } elsif( /effective\s+search\s+space\s+used:\s+(\d+)/ ) {
		       $self->element({'Name' => 'Statistics_eff-spaceused',
				       'Data' => $1});		       
		   } elsif( /effective\s+search\s+space:\s+(\d+)/ ) {
		       $self->element({'Name' => 'Statistics_eff-space',
				       'Data' => $1});
		   } elsif( /Gap\s+Penalties:\s+Existence:\s+(\d+)\,\s+Extension:\s+(\d+)/) {
		       $self->element({'Name' => 'Parameters_gap-open',
				       'Data' => $1});
		       $self->element({'Name' => 'Parameters_gap-extend',
				       'Data' => $2});
		   } elsif( /effective\s+HSP\s+length:\s+(\d+)/ ) {
		        $self->element({'Name' => 'Statistics_hsp-len',
					'Data' => $1});
		   } elsif( /effective\s+length\s+of\s+query:\s+([\d\,]+)/ ) {
		       my $c = $1;
		       $c =~ s/\,//g;
		        $self->element({'Name' => 'Statistics_query-len',
					'Data' => $c});
		   } elsif( /effective\s+length\s+of\s+database:\s+([\d\,]+)/){
		       my $c = $1;
		       $c =~ s/\,//g;
		       $self->element({'Name' => 'Statistics_eff-dblen',
				       'Data' => $c});
		   } elsif( /^(T|A|X1|X2|S1|S2):\s+(\d+)/ ) {
		       $self->element({'Name' => "Statistics_$1",
				       'Data' => $2})
		       } elsif( /frameshift\s+window\,\s+decay\s+const:\s+(\d+)\,\s+([\.\d]+)/ ) {
			   $self->element({'Name'=> 'Statistics_framewindow',
					   'Data' => $1});
			   $self->element({'Name'=> 'Statistics_decay',
					   'Data' => $2});
		       }
	       }
	       $last = $_;
	   }
       } elsif( $self->in_element('hsp') ) {
           # let's read 3 lines at a time;
	   my %data = ( 'Query' => '',
			'Mid' => '',
			'Hit' => '' );
	   my $len;
	   for( my $i = 0; 
		defined($_) && $i < 3; 
		$i++ ){
	       chomp;	       
	       if( /^((Query|Sbjct):\s+(\d+)\s+)(\S+)\s+(\d+)/ ) {
		   $data{$2} = $4;
		   $len = length($1);
		   $self->{"\_$2"}->{'begin'} = $3 unless $self->{"_$2"}->{'begin'};
		   $self->{"\_$2"}->{'end'} = $5;
	       } else { 
		   $self->throw("no data for midline $_") unless defined $_ && defined $len;
		   $data{'Mid'} = substr($_,$len);
	       }
	       $_ = $self->_readline();	       
	   }
	   $self->characters({'Name' => 'Hsp_qseq',
			      'Data' => $data{'Query'} });
	   $self->characters({'Name' => 'Hsp_hseq',
			      'Data' => $data{'Sbjct'}});
	   $self->characters({'Name' => 'Hsp_midline',
			      'Data' => $data{'Mid'} });
       } else { 
	  # print "unrecognized line $_";
       }
   } 
   $self->end_element({'Name' => 'BlastOutput'}) unless ! $seentop;
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
    if( my $type = $MODEMAP{$nm} ) {
	if( $self->_eventHandler->will_handle($type) ) {
	    my $func = sprintf("start_%s",lc $type);
	    $self->_eventHandler->$func($data->{'Attributes'});
	}						 
	unshift @{$self->{'_elements'}}, $type;
    }
    if($nm eq 'BlastOutput') {
	$self->{'_values'} = {};
	$self->{'_result'}= undef;
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
    if($nm eq 'BlastOutput_program' &&
       $self->{'_last_data'} =~ /(t?blast[npx])/i ) {
	$self->{'_reporttype'} = uc $1; 	    
    }   

    # Hsp are sort of weird, in that they end when another
    # object begins so have to detect this in end_element for now
    if( $nm eq 'Hsp' ) {
	foreach ( qw(Hsp_qseq Hsp_midline Hsp_hseq) ) {
	    $self->element({'Name' => $_,
			    'Data' => $self->{'_last_hspdata'}->{$_}});
	}
	$self->{'_last_hspdata'} = {};
	$self->element({'Name' => 'Hsp_query-from',
			'Data' => $self->{'_Query'}->{'begin'}});
	$self->element({'Name' => 'Hsp_query-to',
			'Data' => $self->{'_Query'}->{'end'}});
	
	$self->element({'Name' => 'Hsp_hit-from',
			'Data' => $self->{'_Sbjct'}->{'begin'}});
	$self->element({'Name' => 'Hsp_hit-to',
			'Data' => $self->{'_Sbjct'}->{'end'}});
    }
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
	print "unknown nm $nm, ignoring\n";
    }
    $self->{'_last_data'} = ''; # remove read data if we are at 
				# end of an element
    $self->{'_result'} = $rc if( $nm eq 'BlastOutput' );
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

1;
