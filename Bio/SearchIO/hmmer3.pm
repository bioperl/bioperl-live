# $Id: bioperl.lisp 15559 2009-02-23 12:11:20Z maj $
#
# BioPerl module for Bio::SearchIO::hmmer3
#
# Please direct questions and support issues to <bioperl-l@bioperl.org>
#
# Cared for by Thomas Sharpton <thomas.sharpton@gmail.com>
#
# Copyright Thomas Sharpton
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::SearchIO::hmmer3 - DESCRIPTION of Object

=head1 SYNOPSIS

Give standard usage here

=head1 DESCRIPTION

Describe the object here

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

  bioperl-l@bioperl.org                  - General discussion
  http://bioperl.org/wiki/Mailing_lists  - About the mailing lists

=head2 Support

Please direct usage questions or support issues to the mailing list:

L<bioperl-l@bioperl.org>

rather than to the module maintainer directly. Many experienced and
reponsive experts will be able look at the problem and quickly
address it. Please include a thorough description of the problem
with code and data examples if at all possible.

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution. Bug reports can be submitted via
the web:

  https://redmine.open-bio.org/projects/bioperl/

=head1 AUTHOR - Thomas Sharpton

Email thomas.sharpton@gmail.com

Describe contact details here

=head1 CONTRIBUTORS

Additional contributors names and emails here

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...

package Bio::SearchIO::hmmer3;

use strict;
use Data::Dumper;
use Bio::Factory::ObjectFactory;
use vars qw(%MAPPING %MODEMAP);
use base qw(Bio::SearchIO::hmmer);

BEGIN {

    # mapping of HMMER items to Bioperl hash keys
    %MODEMAP = (
        'HMMER_Output' => 'result',
        'Hit'          => 'hit',
        'Hsp'          => 'hsp'
    );

    %MAPPING = (
        'Hsp_bit-score'   => 'HSP-bits',
        'Hsp_score'       => 'HSP-score',
        'Hsp_evalue'      => 'HSP-evalue',
        'Hsp_query-from'  => 'HSP-query_start',
        'Hsp_query-to'    => 'HSP-query_end',
        'Hsp_hit-from'    => 'HSP-hit_start',
        'Hsp_hit-to'      => 'HSP-hit_end',
        'Hsp_positive'    => 'HSP-conserved',
        'Hsp_identity'    => 'HSP-identical',
        'Hsp_gaps'        => 'HSP-hsp_gaps',
        'Hsp_hitgaps'     => 'HSP-hit_gaps',
        'Hsp_querygaps'   => 'HSP-query_gaps',
        'Hsp_qseq'        => 'HSP-query_seq',
        'Hsp_hseq'        => 'HSP-hit_seq',
        'Hsp_midline'     => 'HSP-homology_seq',
        'Hsp_align-len'   => 'HSP-hsp_length',
        'Hsp_query-frame' => 'HSP-query_frame',
        'Hsp_hit-frame'   => 'HSP-hit_frame',

        'Hit_id'        => 'HIT-name',
        'Hit_len'       => 'HIT-length',
        'Hit_accession' => 'HIT-accession',
        'Hit_desc'      => 'HIT-description',
        'Hit_signif'    => 'HIT-significance',
        'Hit_score'     => 'HIT-score',

        'HMMER_program'   => 'RESULT-algorithm_name',
        'HMMER_version'   => 'RESULT-algorithm_version',
        'HMMER_query-def' => 'RESULT-query_name',
        'HMMER_query-len' => 'RESULT-query_length',
        'HMMER_query-acc' => 'RESULT-query_accession',
        'HMMER_querydesc' => 'RESULT-query_description',
        'HMMER_hmm'       => 'RESULT-hmm_name',
        'HMMER_seqfile'   => 'RESULT-sequence_file',
        'HMMER_db'        => 'RESULT-database_name',
    );
}

=head2 new

 Title   : new
 Usage   : my $obj = new Bio::SearchIO::Hmmer3->new();
 Function: Builds a new Bio::SearchIO::Hmmer3 object
 Returns : an instance of Bio::SearchIO::Hmmer3
 Args    : -fh/-file => HMMER filename
           -format   => 'hmmer3'

=cut

sub _initialize {
  my( $self,@args ) = @_;
  $self->SUPER::_initialize(@args);
  $self->{'_hmmidline'} = 'HMMER 3.0b placeholder';
  $self->{'_alnreport'} = 1; #does report include alignments
}

=head2 next_result

 Title   : next_result
 Usage   : my $hit = $searchio->next_result;
 Function: Returns the next Result from a search
 Returns : Bio::Search::Result::ResultI object
 Args    : none

=cut

sub next_result{
   my ($self)  = @_;
   my $seentop = 0; #Placeholder for when we deal with multi-query reports
   my $reporttype;
   my ( $last, @hit_list, @hsp_list, %hspinfo, %hitinfo, %domaincounter );
   local $/ = "\n";
   local $_;

   my $verbose = $self->verbose;  # cache for speed? JES's idea in hmmer.pm
   $self->start_document();
   local ($_);
   #This is here to ensure that next_result doesn't produce infinite loop
   if(!defined( $_ = $self->_readline) ) {
       return undef;
   }
   else{
       $self->_pushback($_);
   }
   #Regex goes here for HMMER3
   #Start with hmmsearch processing
   while ( defined( $_ = $self->_readline ) ) {
       my $lineorig = $_;
       chomp;
       #Grab the program name.
       if ( $_ =~ m/^\#\s(\S+)\s\:\:\s/ ){
	   my $prog = $1;
	   #TO DO LATER: customize the above regex to adapt to other
	   #program types!!! (hmmscan, etc)
	   $self->start_element( { 'Name' => 'HMMER_Output' } );
	   $self->{'_result_count'}++; #Might need to move to another block
	   $self->element(
	       {
		   'Name' => 'HMMER_program',
		   'Data' => uc($prog)
	       }
	   );
       }
       #Get the HMMER package version and release date
       elsif ( $_ =~ m/^\#\sHMMER\s+(\S+)\s+\((.+)\)/ ) {
	   my $version     = $1;
	   my $versiondate = $2;
	   $self->{'_hmmidline'} = $_;
	   $self->element(
	       {
		   'Name' => 'HMMER_version',
		   'Data' => $version
	       }
           );
       }
       #Get the query info
       elsif( $_ =~ /^\#\squery \w+ file\:\s+(\S+)/ ){
	   if( $self->{'_reporttype'} eq 'HMMSEARCH') {
	       $self->{'_hmmfileline'} = $lineorig;
	       $self->element(
		   {
		       'Name' => 'HMMER_hmm',
		       'Data' => $1
		   }
	       );
	   }
	   elsif( $self->{'_reporttype'} eq 'HMMSCAN' ) {
	       $self->{'_hmmseqline'} = $lineorig;
	       $self->element(
		   {
		       'Name' => 'HMMER_seqfile',
		       'Data' => $1
		   }
	      );
	   }
       }
       #If this is a report without alignments
       elsif( $_ =~ m/^\#\sshow\salignments\sin\soutput/ ){
	   $self->{'_alnreport'} = 0;
       }
       #Get the database info
       elsif( $_ =~ m/^\#\starget\s\S+\sdatabase\:\s+(\S+)/ ){
#	   $self->{'_hmmseqline'} = $lineorig;
#	   $self->element(
#	       {
#		   'Name' => 'HMMER_seqfile',
#		   'Data' => $1
#	       }
#	   );
	   if( $self->{'_reporttype'} eq 'HMMSEARCH') {
	       $self->{'_hmmseqline'} = $lineorig;
	       $self->element(
		   {
		       'Name' => 'HMMER_seqfile',
		       'Data' => $1
		   }
	       );
	   }
	   elsif( $self->{'_reporttype'} eq 'HMMSCAN' ) {
	       $self->{'_hmmfileline'} = $lineorig;
	       $self->element(
		   {
		       'Name' => 'HMMER_hmm',
		       'Data' => $1
		   }
	      );
	   }
       }
       #Get query data
       elsif( $_ =~ s/^Query:\s+// ) {
           #TODO Code to deal with multi query report
	   unless( s/\s+\[[L|M]\=(\d+)\]$// ){
	       warn "Error parsing length for query, offending line $_\n";
	       exit(0);
	   }
	   my $querylen = $1;
	   $self->element(
	       {
		   'Name' => 'HMMER_query-len',
		   'Data' => $querylen
	       }
	   );
	   $self->element(
	       {
		   'Name' => 'HMMER_query-def',
		   'Data' => $_
	       }
	   );
       }
       #Get Accession data
       elsif( $_ =~ s/^Accession:\s+// ){
	   s/\s+$//;
	   $self->element(
	       {
		   'Name' => 'HMMER_query-acc',
		   'Data' => $_
	       }
	   );
       }
       #Get description data
       elsif( $_ =~ s/^Description:\s+// ){
	   s/\s+$//;
	   $self->element(
	       {
		   'Name' => 'HMMER_querydesc',
		   'Data' => $_
	       }
	   );
       }
       #PROCESS HMMSEARCH AND HMMSCAN RESULTS SPECIFIC FORMATTING HERE
       elsif( defined $self->{'_reporttype'} &&
	      (
	       $self->{'_reporttype'} eq 'HMMSEARCH' ||
	       $self->{'_reporttype'} eq 'HMMSCAN'
	      )
	   ){
	   #Complete sequence table data above inclusion threshold
	   if( $_ =~ m/Scores for complete sequence/){
	       while (defined( $_ = $self->_readline ) ) {
		   if ($_ =~ m/inclusion threshold/ || m/Domain( and alignment)? annotation for each/ ||
                       m/\[No hits detected/ || m!^//! ){
		       $self->_pushback($_);
		       last;
		   }
		   #grab table data
		   next if ( m/\-\-\-/ || m/^\s+E\-value\s+score/ || m/^$/);
		   my (
		       $eval_full, $score_full, $bias_full,
		       $eval_best, $score_best, $bias_best,
		       $exp, $n, $hitid, $desc, @hitline
		       );
                   @hitline = split(" ", $_);
                   $eval_full  = shift @hitline;
                   $score_full = shift @hitline;
                   $bias_full  = shift @hitline;
                   $eval_best  = shift @hitline;
                   $score_best = shift @hitline;
                   $bias_best  = shift @hitline;
                   $exp        = shift @hitline;
                   $n          = shift @hitline;
                   $hitid      = shift @hitline;
                   $desc       = join " ", @hitline;

		   if( !defined( $desc ) ){
		       $desc = "";
		   }
		   push @hit_list, [ $hitid, $desc, $eval_full, $score_full ];
		   $hitinfo{$hitid}= $#hit_list;
	       }
	   }
	   #Complete sequence table data below inclusion threshold
	   #not currently fully implemented
	   elsif( $_ =~ m/inclusion threshold/ ){
	       while( defined( $_ = $self->_readline ) ) {
		   if( $_ =~ m/Domain( and alignment)? annotation for each/ ||
                       m/Internal pipeline statistics summary/ ){
		       $self->_pushback($_);
		       last;
		   }
		   next if( $_ =~ m/^$/ );
		   my (
		       $eval_full, $score_full, $bias_full,
		       $eval_best, $score_best, $bias_best,
		       $exp, $n, $hitid, $desc, @hitline
		       );
                   @hitline = split(" ", $_);
                   $eval_full  = shift @hitline;
                   $score_full = shift @hitline;
                   $bias_full  = shift @hitline;
                   $eval_best  = shift @hitline;
                   $score_best = shift @hitline;
                   $bias_best  = shift @hitline;
                   $exp        = shift @hitline;
                   $n          = shift @hitline;
                   $hitid      = shift @hitline;
                   $desc       = join " ", @hitline;

		   $hitinfo{$hitid} = "below_inclusion";
	       }
	   }
	   #Domain annotation for each sequence table data
	   elsif( $_ =~ m/Domain( and alignment)? annotation for each/){
	       @hsp_list = (); #here for multi-query reports
	       my $name;

	       while( defined( $_ = $self->_readline ) ) {
 		   if ($_ =~ m/Internal pipeline statistics/ || m/\[No targets detected/ ){
		       $self->_pushback($_);
		       last;
		   }
		   if( $_ =~ m/^\>\>\s(.*?)\s+/ ) {
 		       $name = $1;
		       #skip hits below inclusion threshold
		       next if( $hitinfo{$name} eq "below_inclusion");
                       $domaincounter{$name} = 0;

		       while( defined( $_ = $self->_readline ) ) {
			   #grab table data for sequence
			   if ($_ =~ m/Internal pipeline statistics/ ||
			       $_ =~ m/^\>\>/                        ){
			       $self->_pushback($_);
			       last;
			   }
			   if ( $_ =~ m/Alignments for each domain/ ) {
			       $self->_pushback($_);
			       last;
			   }
			   if ( $_ =~ m/^\s+\#\s+score/ ||
			        $_ =~ m/^\s\-\-\-\s+/   ||
#			        $_ =~ m/^\>\>/          ||
			        $_ =~ m/^$/             ){
			       next;
			   }

#			   grab hsp data from table, push into @hsp;
			   if(
			       my ($domain_num, $score, $bias, $ceval,
				   $ieval, $hmmstart, $hmmstop,
				   $qalistart, $qalistop, $envstart,
				   $envstop, $envbound, $acc) =
			       m!^\s+(\d+)\s\!*\?*\s+       #domain num
                                   (\S+)\s+(\S+)\s+             #score, bias
                                   (\S+)\s+(\S+)\s+             #c-eval, i-eval
                                   (\d+)\s+(\d+).+?             #hmm start, stop
                                   (\d+)\s+(\d+).+?             #query start, stop
                                   (\d+)\s+(\d+).+?             #env start, stop
                                   (\S+)                        #acc
                                   \s*$!ox
			       ){
			       #keeping simple for now. let's customize later
			       my @vals = ($hmmstart, $hmmstop, $qalistart, $qalistop, $score, $ceval, '', '', '');
			       my $info = $hit_list[ $hitinfo{$name} ];
			       if( !defined $info ){
				   $self->warn(
				       "Incomplete sequence information; can't find $name, hitinfo says $hitinfo{$name}\n"
				       );
				   next;
			       }
			       $domaincounter{$name}++;
                               my $hsp_key = $name . "_" . $domaincounter{$name};
			       push @hsp_list, [ $name, @vals ];
                               $hspinfo{$hsp_key} = $#hsp_list;
			   }
			   else{
			       print "missed this line: $_\n";
			   }
		       }
		   }
		   elsif ($_ =~ m/Alignments for each domain/ ) {
		       my $domain_count = 0;
                       #line counter
                       my $count = 0;
                       # There's an optional block, so we sometimes need to
                       # count to 3, and sometimes to 4.
                       my $max_count = 3;
		       my $lastdomain;
                       my $hsp;
                       my ($hline, $midline, $qline);

		       while( defined( $_ = $self->_readline ) ) {
			   if( $_ =~ m/^\>\>/ ||
			       $_ =~ m/Internal pipeline statistics/){
			       $self->_pushback($_);
			       last;
			   }
			   elsif( $hitinfo{$name} eq "below_inclusion" ||
			           $_ =~ m/^$/ ) {
			       next;
			   }
			   elsif( $_ =~ /\s\s\=\=\sdomain\s(\d+)\s+/){
			       my $domainnum = $1;
			       $count = 0;
                               my $key = $name . "_" . $domainnum;
                               $hsp = $hsp_list[ $hspinfo{$key} ];
                               $hline = $$hsp[-3];
                               $midline = $$hsp[-2];
                               $qline = $$hsp[-1];
			       $lastdomain = $name;
			   }
                           # model data track, some reports don't have
                           elsif( $_ =~ m/\s+\S+\sCS$/ ){
			       my $modeltrack = $_;
                               $max_count++;
			       $count++;
			       next;
			   }
			   elsif( $count == $max_count - 3 ){
			       #hit sequence
			       my @data = split(" ", $_);
			       my $seq = $data[-2];
                               $hline .= $seq;
			       $count++;
			       next;
			   }
			   elsif( $count == $max_count - 2 ){
			       #conservation track
			       #storage isn't quite right - need to remove
			       #leading/lagging whitespace while preserving
			       #gap data (latter isn't done, former is)
			       $_ =~ s/^\s+//;
			       $_ =~ s/\s+$//;
                               $midline .= $_;
			       $count++;
			       next;
			   }
			   elsif( $count == $max_count - 1 ){
			       #query track
			       my @data = split(" ", $_);
			       my $seq = $data[-2];
                               $qline .= $seq;
			       $count++;
			       next;
			   }
			   elsif( $count == $max_count ){
			       #pval track
			       my $pvals = $_;
			       $count = 0;
                               $max_count = 3;
                               $$hsp[-3] = $hline;
                               $$hsp[-2] = $midline;
                               $$hsp[-1] = $qline;
			       next;
			   }
			   else{
                               print "missed $_\n";
			   }
		       }
		   }
	       }
	   }
	   elsif( m/Internal pipeline statistics/ || m!^//! ){
#	       if within hit, hsp close;
	       if ( $self->within_element('hit') ) {
		   if ( $self->within_element('hsp') ) {
		       $self->end_element( { 'Name' => 'Hsp' } );
		   }
		   $self->end_element( { 'Name' => 'Hit' } );
	       }
	       #grab summary statistics of run
	       while( defined( $_ = $self->_readline ) ) {
                   last if ( $_ =~ m/^\/\/$/ );
	       }

	       #Jason does a lot of processing of hits/hsps here;
	       while( my $hit = shift @hit_list ) {
                   my $hit_name = shift @$hit;
                   my $hit_desc = shift @$hit;
                   my $hit_signif = shift @$hit;
                   my $hit_score = shift @$hit;
                   my $num_domains = $domaincounter{$hit_name} || 0;

                   $self->start_element( { 'Name' => 'Hit' } );
                   $self->element(
                       {
                           'Name' => 'Hit_id',
                           'Data' => $hit_name
                       }
                   );
                   $self->element(
                       {
                           'Name' => 'Hit_desc',
                           'Data' => $hit_desc
                       }
                   );
                   $self->element(
                       {
                           'Name' => 'Hit_signif',
                           'Data' => $hit_signif
                       }
                   );
                   $self->element(
                       {
                           'Name' => 'Hit_score',
                           'Data' => $hit_score
                       }
                   );
                   for my $i (1..$num_domains) {
                       my $key = $hit_name . "_" . $i;
                       my $hsp = $hsp_list[ $hspinfo{$key} ];
                       if(defined $hsp) {
                           my $hsp_name = shift @$hsp;
                           $self->start_element( { 'Name' => 'Hsp' } );
                           $self->element( {
                                   'Name' => 'Hsp_identity',
                                   'Data' => 0
                               } );
                           $self->element( {
                                   'Name' => 'Hsp_positive',
                                   'Data' => 0
                               } );
                           $self->element( {
                                   'Name' => 'Hsp_hit-from',
                                   'Data' => shift @$hsp
                               } );
                           $self->element( {
                                   'Name' => 'Hsp_hit-to',
                                   'Data' => shift @$hsp
                               } );
                           $self->element( {
                                   'Name' => 'Hsp_query-from',
                                   'Data' => shift @$hsp
                               } );
                           $self->element( {
                                   'Name' => 'Hsp_query-to',
                                   'Data' => shift @$hsp
                               } );
                           $self->element( {
                                   'Name' => 'Hsp_score',
                                   'Data' => shift @$hsp
                               } );
                           $self->element( {
                                   'Name' => 'Hsp_evalue',
                                   'Data' => shift @$hsp
                               } );
                           $self->element( {
                                   'Name' => 'Hsp_hseq',
                                   'Data' => shift @$hsp
                               } );
                           $self->element( {
                                   'Name' => 'Hsp_midline',
                                   'Data' => shift @$hsp
                               } );
                           $self->element( {
                                   'Name' => 'Hsp_qseq',
                                   'Data' => shift @$hsp
                               } );
                           $self->end_element( { 'Name' => 'Hsp' } );
                       }
                   }
                   $self->end_element( { 'Name' => 'Hit' } );
	       }
	       @hit_list = ();
	       %hitinfo = ();
	       last;
	   }
       }
       else{
	   print "missed: $_\n";
	   $self->debug($_);
       }
       $last = $_;
   }
   $self->end_element( { 'Name' => 'HMMER_Output' } );
   my $result = $self->end_document();
   return $result;
}



=head2 start_element

 Title   : start_element
 Usage   : $eventgenerator->start_element
 Function: Handles a start event
 Returns : none
 Args    : hashref with at least 2 keys 'Data' and 'Name'

=cut

sub start_element{

    my ( $self, $data ) = @_;

    # we currently don't care about attributes
    my $nm   = $data->{'Name'};
    my $type = $MODEMAP{$nm};
    if ($type) {
        if ( $self->_eventHandler->will_handle($type) ) {
            my $func = sprintf( "start_%s", lc $type );
            $self->_eventHandler->$func( $data->{'Attributes'} );
        }
        unshift @{ $self->{'_elements'} }, $type;
    }
    if ( defined $type
        && $type eq 'result' )
    {
        $self->{'_values'} = {};
        $self->{'_result'} = undef;
    }
}

=head2 end_element

 Title   : end_element
 Usage   : $eventgeneartor->end_element
 Function: Handles and end element event
 Returns : none
 Args    : hashref with at least 2 keys 'Data' and 'Name'

=cut

sub end_element{

    my ( $self, $data ) = @_;
    my $nm   = $data->{'Name'};
    my $type = $MODEMAP{$nm};
    my $rc;

    if ( $nm eq 'HMMER_program' ) {
        if ( $self->{'_last_data'} =~ /(HMM\S+)/i ) {
            $self->{'_reporttype'} = uc $1;
        }
    }

    # Hsp are sort of weird, in that they end when another
    # object begins so have to detect this in end_element for now
    if ( $nm eq 'Hsp' ) {
        foreach (qw(Hsp_qseq Hsp_midline Hsp_hseq)) {
            my $data = $self->{'_last_hspdata'}->{$_};
            if ($data && $_ eq 'Hsp_hseq') {
                # replace hmm '.' gap symbol by '-'
                $data =~ s/\./-/g;
            }
            $self->element(
                {
                    'Name' => $_,
                    'Data' => $data
                }
            );
        }
        $self->{'_last_hspdata'} = {};
    }
    if ($type) {
        if ( $self->_eventHandler->will_handle($type) ) {
            my $func = sprintf( "end_%s", lc $type );
            $rc = $self->_eventHandler->$func( $self->{'_reporttype'},
                $self->{'_values'} );
	}
        my $lastelem = shift @{ $self->{'_elements'} };
    }
    elsif ( $MAPPING{$nm} ) {
       if ( ref( $MAPPING{$nm} ) =~ /hash/i ) {
            my $key = ( keys %{ $MAPPING{$nm} } )[0];
            $self->{'_values'}->{$key}->{ $MAPPING{$nm}->{$key} } =
              $self->{'_last_data'};
        }
        else {
            $self->{'_values'}->{ $MAPPING{$nm} } = $self->{'_last_data'};
#	    print "lastdata is " . $self->{'_last_data'} . "\n";
        }
    }
    else {
        $self->debug("unknown nm $nm, ignoring\n");
    }
    $self->{'_last_data'} = '';    # remove read data if we are at
                                   # end of an element
    $self->{'_result'} = $rc if ( defined $type && $type eq 'result' );
    return $rc;
}

=head2 element

 Title   : element
 Usage   : $eventhandler->element({'Name' => $name, 'Data' => $str});
 Function: Convienence method that calls start_element, characters, end_element
 Returns : none
 Args    : Hash ref with the keys 'Name' and 'Data'

=cut

sub element{
    my ( $self, $data ) = @_;
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
    my ( $self, $data ) = @_;

    if (   $self->in_element('hsp')
        && $data->{'Name'} =~ /Hsp\_(qseq|hseq|midline)/o
        && defined $data->{'Data'} )
    {
        $self->{'_last_hspdata'}->{ $data->{'Name'} } .= $data->{'Data'};
    }
    return unless ( defined $data->{'Data'} && $data->{'Data'} !~ /^\s+$/o );

    $self->{'_last_data'} = $data->{'Data'};
}

=head2 within_element

 Title   : within_element
 Usage   : if( $eventgenerator->within_element( $element ) ) {}
 Function: Test if we are within a particular element
           This is different than 'in' because within can be tested for
           a whole block
 Returns : boolean
 Args    : string element name

=cut

sub within_element{
    my ( $self, $name ) = @_;
    return 0
      if ( !defined $name
        || !defined $self->{'_elements'}
        || scalar @{ $self->{'_elements'} } == 0 );
    foreach ( @{ $self->{'_elements'} } ) {
        return 1 if ( $_ eq $name );
    }
    return 0;
}

=head2 in_element

 Title   : in_element
 Usage   : if( $eventgenerator->in_element( $element ) ) {}
 Function: Test if we are in a particular element
           This is different than 'within' because 'in' only
           tests its immediate parent
 Returns : boolean
 Args    : string element name

=cut

sub in_element{
    my ( $self, $name ) = @_;
    return 0 if !defined $self->{'_elements'}->[0];
    return ( $self->{'_elements'}->[0] eq $name );
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
   $self->{'_values'}   = {};
   $self->{'_result'}   = undef;
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
   my ($self) = @_;
   return $self->{'_result'};
}

=head2 result_count

 Title   : result_count
 Usage   : my $count = $searchio->result_count
 Function: Returns the number of results processed
 Returns : interger
 Args    : none

=cut

sub result_count{
   my $self = shift;
   return $self->{'_result_count'};
}


1;
