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
Read the L<Bio::SearchIO> for more information about how to use this.

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

    %MAPPING = 
	( 
	  'Hsp_bit-score'  => 'HSP-bits',
	  'Hsp_score'      => 'HSP-score',
	  'Hsp_evalue'     => 'HSP-evalue',
	  'Hsp_pvalue'     => 'HSP-pvalue',
	  'Hsp_query-from' => 'HSP-query_start',
	  'Hsp_query-to'   => 'HSP-query_end',
	  'Hsp_hit-from'   => 'HSP-hit_start',
	  'Hsp_hit-to'     => 'HSP-hit_end',
	  'Hsp_positive'   => 'HSP-conserved',
	  'Hsp_identity'   => 'HSP-identical',
	  'Hsp_gaps'       => 'HSP-hsp_gaps',
	  'Hsp_hitgaps'    => 'HSP-hit_gaps',
	  'Hsp_querygaps'  => 'HSP-query_gaps',
	  'Hsp_qseq'       => 'HSP-query_seq',
	  'Hsp_hseq'       => 'HSP-hit_seq',
	  'Hsp_midline'    => 'HSP-homology_seq',
	  'Hsp_align-len'  => 'HSP-hsp_length',
	  'Hsp_query-frame'=> 'HSP-query_frame',
	  'Hsp_hit-frame'  => 'HSP-hit_frame',

	  'Hit_id'        => 'HIT-name',
	  'Hit_len'       => 'HIT-length',
	  'Hit_accession' => 'HIT-accession',
	  'Hit_def'       => 'HIT-description',
	  'Hit_signif'    => 'HIT-significance',
	  'Hit_score'     => 'HIT-score',
	  'Iteration_iter-num'   => 'HIT-iteration',

	  'BlastOutput_program'  => 'RESULT-algorithm_name',
	  'BlastOutput_version'  => 'RESULT-algorithm_version',
	  'BlastOutput_query-def'=> 'RESULT-query_name',
	  'BlastOutput_query-len'=> 'RESULT-query_length',
	  'BlastOutput_query-acc'=> 'RESULT-query_accession',
	  'BlastOutput_querydesc'=> 'RESULT-query_description',
	  'BlastOutput_db'       => 'RESULT-database_name',
	  'BlastOutput_db-len'   => 'RESULT-database_entries',
	  'BlastOutput_db-let'   => 'RESULT-database_letters',

	  'Parameters_matrix'    => { 'RESULT-parameters' => 'matrix'},
	  'Parameters_expect'    => { 'RESULT-parameters' => 'expect'},
	  'Parameters_include'   => { 'RESULT-parameters' => 'include'},
	  'Parameters_sc-match'  => { 'RESULT-parameters' => 'match'},
	  'Parameters_sc-mismatch' => { 'RESULT-parameters' => 'mismatch'},
	  'Parameters_gap-open'  =>   { 'RESULT-parameters' => 'gapopen'},
	  'Parameters_gap-extend'=>   { 'RESULT-parameters' => 'gapext'},
	  'Parameters_filter'    =>  {'RESULT-parameters' => 'filter'},
	  'Parameters_allowgaps' =>   { 'RESULT-parameters' => 'allowgaps'},

	  'Statistics_db-len'    => {'RESULT-statistics' => 'dbentries'},
	  'Statistics_db-let'    => { 'RESULT-statistics' => 'dbletters'},
	  'Statistics_hsp-len'   => { 'RESULT-statistics' => 'hsplength'},
	  'Statistics_query-len'   => { 'RESULT-statistics' => 'querylength'},
	  'Statistics_eff-space' => { 'RESULT-statistics' => 'effectivespace'},
	  'Statistics_eff-spaceused' => { 'RESULT-statistics' => 'effectivespaceused'},
	  'Statistics_eff-dblen' => { 'RESULT-statistics' => 'effectivedblength'},
	  'Statistics_kappa'     => { 'RESULT-statistics' => 'kappa' },
	  'Statistics_lambda'    => { 'RESULT-statistics' => 'lambda' },
	  'Statistics_entropy'   => { 'RESULT-statistics' => 'entropy'},
	  'Statistics_framewindow'=> { 'RESULT-statistics' => 'frameshiftwindow'},
	  'Statistics_decay'=> { 'RESULT-statistics' => 'decayconst'},

	  'Statistics_T'=> { 'RESULT-statistics' => 'T'},
	  'Statistics_A'=> { 'RESULT-statistics' => 'A'},
	  'Statistics_X1'=> { 'RESULT-statistics' => 'X1'},
	  'Statistics_X2'=> { 'RESULT-statistics' => 'X2'},
	  'Statistics_S1'=> { 'RESULT-statistics' => 'S1'},
	  'Statistics_S2'=> { 'RESULT-statistics' => 'S2'},

	  # WU-BLAST stats
	  'Statistics_DFA_states'=> { 'RESULT-statistics' => 'num_dfa_states'},
	  'Statistics_DFA_size'=> { 'RESULT-statistics' => 'dfa_size'},

	  'Statistics_search_cputime' => { 'RESULT-statistics' => 'search_cputime'},
	  'Statistics_total_cputime' => { 'RESULT-statistics' => 'total_cputime'},
	  'Statistics_search_actualtime' => { 'RESULT-statistics' => 'search_actualtime'},
	  'Statistics_total_actualtime' => { 'RESULT-statistics' => 'total_actualtime'},

	  'Statistics_noprocessors' => { 'RESULT-statistics' => 'no_of_processors'},
	  'Statistics_neighbortime' => { 'RESULT-statistics' => 'neighborhood_generate_time'},
	  'Statistics_starttime' => { 'RESULT-statistics' => 'start_time'},
	  'Statistics_endtime' => { 'RESULT-statistics' => 'end_time'},
	  );
}


=head2 new

 Title   : new
 Usage   : my $obj = new Bio::SearchIO::blast();
 Function: Builds a new Bio::SearchIO::blast object 
 Returns : Bio::SearchIO::blast
 Args    : -fh/-file => filehandle/filename to BLAST file
           -format   => 'blast'

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
   
   my $data = '';
   my $seentop = 0;
   my $reporttype;
   $self->start_document();
   my @hit_signifs;
   while( defined ($_ = $self->_readline )) {       
       next if( /^\s+$/); # skip empty lines
       next if( /CPU time:/);
       next if( /^>\s*$/);
       
       if( /^([T]?BLAST[NPX])\s*(.+)$/i ) {
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
	   my ($nm,$desc) = split(/\s+/,$q,2);
	   $self->element({ 'Name' => 'BlastOutput_query-def',
			    'Data' => $nm});
	   $self->element({ 'Name' => 'BlastOutput_query-len', 
			    'Data' => $size});
	   defined $desc && $desc =~ s/\s+$//;
	   $self->element({ 'Name' => 'BlastOutput_querydesc', 
			    'Data' => $desc});
	   
	   if( my @pieces = split(/\|/,$nm) ) {
	       my $acc = pop @pieces;
	       $acc = pop @pieces if( ! defined $acc || $acc =~ /^\s+$/);
	       $self->element({ 'Name' =>  'BlastOutput_query-acc',
				'Data'  => $acc});
	   }
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
	   chomp;
	   $self->in_element('hsp') && $self->end_element({ 'Name' => 'Hsp'});
	   $self->in_element('hit') && $self->end_element({ 'Name' => 'Hit'});
	   
	   $self->start_element({ 'Name' => 'Hit'});
	   my $id = $1;	  
	   my $restofline = $2;
	   $self->element({ 'Name' => 'Hit_id',
			    'Data' => $id});
	   my @pieces = split(/\|/,$id);
	   my $acc = pop @pieces;
	   $self->element({ 'Name' =>  'Hit_accession',
			    'Data'  => $acc});	   

	   my $v = shift @hit_signifs;
	   if( defined $v ) {
	       $self->element({'Name' => 'Hit_signif',
			       'Data' => $v->[0]});
	       $self->element({'Name' => 'Hit_score',
			       'Data' => $v->[1]});
	   }
	   while(defined($_ = $self->_readline()) ) {
	       next if( /^\s+$/ );
	       chomp;
	       if(  /Length\s*=\s*([\d,]+)/ ) {
		   my $l = $1;
		   $l =~ s/\,//g;
		   $self->element({ 'Name' => 'Hit_len',
				    'Data' => $l });
		   last;
	       }  elsif ( /Score/ ) {
		   $self->_pushback($_);
		   last;
	       } else { 
		   $restofline .= $_;
	       }
	   }
	   $restofline =~ s/\s+/ /g;
	   $self->element({ 'Name' => 'Hit_def',
			    'Data' => $restofline});       
      } elsif( /\s+(Plus|Minus) Strand HSPs:/i ) {
	   next;
       } elsif( ($self->in_element('hit') || 
		 $self->in_element('hsp')) && # wublast
	       /Score\s*=\s*(\S+)\s*\(([\d\.]+)\s*bits\),\s*Expect\s*=\s*([^,\s]+),\s*(Sum)?\s*P(\(\d+\))?\s*=\s*([^,\s]+)/ 
		  ) {
	   $self->in_element('hsp') && $self->end_element({'Name' => 'Hsp'});
	   $self->start_element({'Name' => 'Hsp'});
       	   $self->element( { 'Name' => 'Hsp_score',
			     'Data' => $1});
	   $self->element( { 'Name' => 'Hsp_bit-score',
			     'Data' => $2});
	   $self->element( { 'Name' => 'Hsp_evalue',			     
			     'Data' => $3});
	   $self->element( {'Name'  => 'Hsp_pvalue',
			    'Data'  =>$6});       
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
	   } 
	   $self->{'_Query'} = { 'begin' => 0, 'end' => 0};
	   $self->{'_Sbjct'} = { 'begin' => 0, 'end' => 0};

	   if( /(Frame\s*=\s*.+)$/ ) {
	       # handle wu-blast Frame listing on same line
	       $self->_pushback($1);
	   }	   
       } elsif( $self->in_element('hsp') &&
		/Strand\s*=\s*(Plus|Minus)\s*\/\s*(Plus|Minus)/i ) {
	   # consume this event
	   next;
       } elsif( $self->in_element('hsp') &&
		/Frame\s*=\s*([\+\-][1-3])\s*(\/\s*([\+\-][1-3]))?/ ){
	   my ($queryframe,$hitframe);
	   if( $reporttype eq 'TBLASTX' ) {
	       ($queryframe,$hitframe) = ($1,$2);
	       $hitframe =~ s/\/\s*//g;
	   } elsif( $reporttype eq 'TBLASTN' ) {
	       ($hitframe,$queryframe) = ($1,0);	       
	   } elsif( $reporttype eq 'BLASTX' ) {	       
	       ($queryframe,$hitframe) = ($1,0);
	   } 
	   $self->element({'Name' => 'Hsp_query-frame',
			   'Data' => $queryframe});
	   	   
	   $self->element({'Name' => 'Hsp_hit-frame',
			   'Data' => $hitframe});
       } elsif(  /^Parameters:/ || /^\s+Database:\s+?/ || 
		 ( $self->in_element('hsp') && (/WARNING/ || /NOTE/)) ) {
	   $self->in_element('hsp') && $self->end_element({'Name' => 'Hsp'});
	   $self->in_element('hit') && $self->end_element({'Name' => 'Hit'});
	   my $blast = ( /Parameters\:/ ) ? 'wublast' : 'ncbi'; 
	   my $last = '';
	   # default is that gaps are allowed
	   $self->element({'Name' => 'Parameters_allowgaps',
			   'Data' => 'yes'});
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
		   } elsif( /nogaps/ ) {
		       $self->element({'Name' => 'Parameters_allowgaps',
				       'Data' => 'no'});
		   } elsif( $last =~ /(Frame|Strand)\s+MatID\s+Matrix name/i ) {
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
		   } elsif( /(\S+\s+\S+)\s+DFA:\s+(\S+)\s+\((.+)\)/ ) {
		       if( $1 eq 'states in') { 
			   $self->element({'Name' => 'Statistics_DFA_states',
					   'Data' => "$2 $3"});
		       } elsif( $1 eq 'size of') {
			   $self->element({'Name' => 'Statistics_DFA_size',
					   'Data' => "$2 $3"});
		       }
		   } elsif( /^\s+Time to generate neighborhood:\s+(\S+\s+\S+\s+\S+)/ ) { 
		       $self->element({'Name' => 'Statistics_neighbortime',
				       'Data' => $1});
		   } elsif( /processors\s+used:\s+(\d+)/ ) {
		          $self->element({'Name' => 'Statistics_noprocessors',
					   'Data' => $1});
		   } elsif( /^\s+(\S+)\s+cpu\s+time:\s+(\S+\s+\S+\s+\S+)\s+Elapsed:\s+(\S+)/ ) {
		       my $cputype = lc($1);
		       $self->element({'Name' => "Statistics_$cputype\_cputime",
				       'Data' => $2});
		       $self->element({'Name' => "Statistics_$cputype\_actualtime",
				       'Data' => $3});
		   } elsif( /^\s+Start:/ ) {
		       my ($junk,$start,$stime,$end,$etime) = split(/\s+(Start|End)\:\s+/,$_);
		       chomp($stime);
		       $self->element({'Name' => 'Statistics_starttime',
				       'Data' => $stime});
		       chomp($etime);
		       $self->element({'Name' => 'Statistics_endtime',
				       'Data' => $etime});
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
	       if( /^((Query|Sbjct):\s+(\d+)\s*)(\S+)\s+(\d+)/ ) {
		   $data{$2} = $4;
		   $len = length($1);
		   $self->{"\_$2"}->{'begin'} = $3 unless $self->{"_$2"}->{'begin'};
		   $self->{"\_$2"}->{'end'} = $5;
	       } else { 
		   $self->throw("no data for midline $_") 
		       unless (defined $_ && defined $len);
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
	   $self->debug( "unrecognized line $_");
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
