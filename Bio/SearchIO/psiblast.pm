#------------------------------------------------------------------
# $Id$
#
# BioPerl module Bio::SearchIO::psiblast
#
# Cared for by Steve Chervitz <sac@bioperl.org>
#
# You may distribute this module under the same terms as perl itself
#------------------------------------------------------------------

=head1 NAME

Bio::SearchIO::psiblast - Parser for traditional BLAST and PSI-BLAST reports

=head1 SYNOPSIS

    use Bio::SearchIO;

    my $in = Bio::SearchIO->new( -format  => 'psiblast',
                                 -file    => 'report.blastp' );

    while ( my $blast = $in->next_result() ) {
	foreach my $hit ( $blast->hits ) {
            print "Hit: $hit\n";
        }
   }

    # BLAST hit filtering function. All hits of each BLAST report must satisfy 
    # this criteria to be retained. If a hit fails this test, it is ignored.
    # If all hits of a report fail, the report will be considered hitless.
    # But we can distinguish this from the case where there were no
    # hits in the report by testing the function $blast->no_hits_found().

    my $filt_func = sub{ my $hit=shift; 
    			 $hit->frac_identical('query') >= 0.5 
    			     && $hit->frac_aligned_query >= 0.50
    			 };

    # Not supplying a -file or -fh parameter means read from STDIN

    my $in2 = Bio::SearchIO->new( -format  => 'psiblast',
                                  -hit_filter => $filt_func
                                 );


=head1 DESCRIPTION

This module parses BLAST and PSI-BLAST reports and acts as a factory for
objects that encapsulate BLAST results:
L<Bio::Search::Result::BlastResult>, L<Bio::Search::Hit::BlastHit>,
L<Bio::Search::HSP::BlastHSP>.

This module does not parse XML-formatted BLAST reports.
See L<Bio::SearchIO::blastxml|Bio::SearchIO::blastxml> if you need to do that.

To use this module, the only module you need to C<use> is
Bio::SearchIO.pm. SearchIO knows how to load this module when you
supply a C<-format =E<gt> 'psiblast'> parameters to its C<new>() 
function. For more information about the SearchIO system, see
documentation in Bio::SearchIO.pm.

=head2 PSI-BLAST Support

In addition to BLAST1 and BLAST2 reports, this module can also handle
PSI-BLAST reports. When accessing the set of Hits in a result, hits
from different iterations are lumped together but can be distinguished by
interrogating L<Bio::Search::Hit::BlastHit::iteration> and 
L<Bio::Search::Hit::BlastHit::found_again>.

If you want to collect hits only from a certain iteration during parsing,
supply a function using the C<-HIT_FILTER> parameter.

=head1 EXAMPLES

To get a feel for how to use this, have look at scripts in the
B<examples/searchio> and B<examples/searchio/writer> directory of the Bioperl 
distribution as well as the test script B<t/SearchIO.t>.

=head1 SEE ALSO

For more documentation about working with Blast result objects that are
produced by this parser, see L<Bio::Search::Result::BlastResult>, 
L<Bio::Search::Hit::BlastHit>, L<Bio::Search::HSP::BlastHSP>.

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
  http://bugzilla.bioperl.org/

=head1 AUTHOR 

Steve Chervitz E<lt>sac@bioperl.orgE<gt>

See L<the FEEDBACK section | FEEDBACK> for where to send bug reports and comments.

=head1 ACKNOWLEDGEMENTS

I would like to acknowledge my colleagues at Affymetrix for useful
feedback.

=head1 COPYRIGHT

Copyright (c) 2001 Steve Chervitz. All Rights Reserved.

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.

=head1 DISCLAIMER

This software is provided "as is" without warranty of any kind.

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::SearchIO::psiblast;

use strict;
use vars qw( @ISA 
	     $MAX_HSP_OVERLAP 
	     $DEFAULT_MATRIX  
	     $DEFAULT_SIGNIF  
	     $DEFAULT_SCORE
	     $DEFAULT_BLAST_WRITER_CLASS 
	     $DEFAULT_HIT_FACTORY_CLASS 
	     $DEFAULT_RESULT_FACTORY_CLASS  
	     );

use Bio::SearchIO;
use Bio::Search::Result::BlastResult;
use Bio::Factory::BlastHitFactory;
use Bio::Factory::BlastResultFactory;
use Bio::Tools::StateMachine::IOStateMachine qw($INITIAL_STATE $FINAL_STATE);

@ISA = qw( Bio::SearchIO
           Bio::Tools::StateMachine::IOStateMachine  );

$MAX_HSP_OVERLAP  = 2;  # Used when tiling multiple HSPs.
$DEFAULT_MATRIX   = 'BLOSUM62';
$DEFAULT_SIGNIF   = 999;# Value used as significance cutoff if none supplied.
$DEFAULT_SCORE    = 0;  # Value used as score cutoff if none supplied.
$DEFAULT_BLAST_WRITER_CLASS    = 'Bio::Search::Writer::HitTableWriter';
$DEFAULT_HIT_FACTORY_CLASS     = 'Bio::Factory::BlastHitFactory';
$DEFAULT_RESULT_FACTORY_CLASS  = 'Bio::Factory::BlastResultFactory';

my %state = (
              'header'    => 'Header',
              'descr'     => 'Descriptions',
              'align'     => 'Alignment',
              'footer'    => 'Footer',
              'iteration' => 'Iteration',      # psiblast
              'nohits'    => 'No Hits'
             );

# These are the expected transitions assuming a "normal" report (Blast2 or PSI-Blast).
my @state_transitions = (  [ $INITIAL_STATE, $state{'header'}],
                           [ $state{'header'}, $state{'descr'} ],
                           [ $state{'header'}, $state{'iteration'} ],
                           [ $state{'iteration'},  $state{'descr'} ],
                           [ $state{'iteration'}, $state{'nohits'} ],
                           [ $state{'descr'}, $state{'align'} ],
                           [ $state{'align'}, $state{'align'} ],
                           [ $state{'align'}, $state{'footer'} ],
                           [ $state{'align'}, $state{'iteration'} ],   # psiblast
                           [ $state{'nohits'}, $state{'iteration'} ],  # psiblast
                           [ $state{'nohits'}, $state{'footer'} ],
                           [ $state{'footer'}, $state{'header'} ],
                           [ $state{'footer'}, $FINAL_STATE]
                        );

my $current_iteration;  # psiblast

=head2 new

 Usage     : Bio::SearchIO::psiblast->new( %named_parameters )
 Purpose   : Parse traditional BLAST or PSI-BLAST data a file or input stream.
           : Handles Blast1, Blast2, NCBI and WU Blast reports.
           : Populates Bio::Search:: objects with data extracted from the report.
           : (The exact type of Bio::Search objects depends on the type of
           : Bio::Factory::ResultFactory and Bio::Factory::HitFactory you hook up
           : to the SearchIO object.)
 Returns   : Bio::SearchIO::psiblast object.
 Argument  : Named parameters:  (PARAMETER TAGS CAN BE UPPER OR LOWER CASE).
           : These are in addition to those specified by Bio::SearchIO::new() (see).
	   : -SIGNIF     => number (float or scientific notation number to be used
	   :                         as a P- or Expect value cutoff; default = 999.)
	   : -SCORE     => number (integer or scientific notation number to be used
	   :                         as a score value cutoff; default = 0.)
	   : -HIT_FILTER  => func_ref (reference to a function to be used for
           :                          filtering out hits based on arbitrary criteria.
           :                          This function should take a
           :                          Bio::Search::Hit::BlastHit.pm object as its first
           :                          argument and return a boolean value,
	   :                          true if the hit should be filtered out).
           :                          Sample filter function:
           :                          -HIT_FILTER => sub { $hit = shift;
	   :				                  $hit->gaps == 0; },
           :                 Historical note: This parameter was formerly
                             called -FILT_FUNC in the older
                             Bio::Tools::Blast::parse method. Using
                             -FILT_FUNC will still work for backward
                             compatibility.
           : -CHECK_ALL_HITS => boolean (check all hits for significance against
           :                             significance criteria.  Default = false.
	   :			         If false, stops processing hits after the first
           :                             non-significant hit or the first hit that fails
           :                             the hit_filter call. This speeds parsing,
           :                             taking advantage of the fact that the hits
           :                             are processed in the order they are ranked.)
           : -MIN_LEN     => integer (to be used as a minimum query sequence length
           :                          sequences below this length will not be processed).
           :                          default = no minimum length).
	   : -STATS       => boolean (collect key statistical parameters for the report: 
           :                          matrix, filters, etc. default = false). 
           :                          This requires extra parsing
           :                          so if you aren't interested in this info, don't
           :                          set this parameter. Note that the unparsed 
           :                          parameter section of a Blast report is always
           :                          accessible via $blast->raw_statistics().
	   : -BEST        => boolean (only process the best hit of each report;
           :                          default = false).
           : -OVERLAP     => integer (the amount of overlap to permit between
           :                          adjacent HSPs when tiling HSPs. A reasonable value is 2.
           :                          Default = $Bio::SearchIO::psiblast::MAX_HSP_OVERLAP)
           : -HOLD_RAW_DATA => boolean (store the raw alignment sections for each hit.
           :                            used with the -SHALLOW_PARSE option).
           : -SHALLOW_PARSE => boolean (only minimal parsing; does not parse HSPs.
           :                            Hit data is limited to what can be obtained
           :                            from the description line.
           :                            Replaces the older NO_ALIGNS option.)
	   :		
           :
 Comments  : Do NOT remove the HTML from an HTML-formatted Blast report by using the
           : "Save As" option of a web browser to save it as text. This renders the
           : report unparsable.
 Throws    : An exception will be thrown if a BLAST report contains a FATAL: error.

=cut

sub new {
  my ($class, %args) = @_; 

  # TODO: Resolve this issue:
  # Temporary hack to allow factory-based and non-factory based 
  # SearchIO objects co-exist.
  # steve --- Sat Dec 22 04:41:20 2001
  $args{-USE_FACTORIES} = 1;

  my $self = $class->Bio::SearchIO::new(%args);

  $self->_init_state_machine( %args, -transition_table => \@state_transitions);

  $self->_init_parse_params( %args );

  $self->pause_between_reports( 1 );

  $self->{'_result_count'} = 0;

  return $self;
}


sub default_result_factory_class {
  my $self = shift;
  return $DEFAULT_RESULT_FACTORY_CLASS;
}

sub default_hit_factory_class {
  my $self = shift;
  return $DEFAULT_HIT_FACTORY_CLASS;
}

sub check_for_new_state {
    my ($self) = @_;

    # Ignore blank lines
    my $chunk = $self->SUPER::check_for_new_state(1);

    my $newstate = undef;

    # End the run if there's no more input.
    if( ! $chunk ) {
	return $self->final_state;
    }
    $self->clear_state_change_cache;
    my $curr_state = $self->current_state;

    if( $chunk =~ /^(<.*>)?T?BLAST[NPX]/ ) {
	$newstate = $state{header};
	$self->state_change_cache( $chunk );
    }

    elsif ($chunk =~ /^Sequences producing/ ) {
	$newstate = $state{descr};
	$self->state_change_cache( $chunk );
    }

    elsif ($chunk =~ /No hits found/i ) {
	$newstate = $state{nohits};
	$self->state_change_cache( $chunk );
    }

    elsif ($chunk =~ /^\s*Searching/ ) {
	$newstate = $state{iteration};
    }

    elsif ($chunk =~ /^>(.*)/ ) {
	$newstate = $state{align};
	$self->state_change_cache( "$1\n" );
    }

    elsif ($chunk =~ /^(CPU time|Parameters):/ ) {
	$newstate = $state{footer};
	$self->state_change_cache( $chunk );
    }

    # Necessary to distinguish "  Database:" lines that start a footer section
    # from those that are internal to a footer section.
    elsif ($chunk =~ /^\s+Database:/ && $curr_state ne $state{'footer'}) {
	$newstate = $state{footer};
	$self->state_change_cache( $chunk );
    }

    if( $curr_state ne $INITIAL_STATE and not $newstate  ) {
#        print "Appending input cache with ($curr_state): $chunk\n";
	$self->append_input_cache( $chunk );
    }

    return $newstate;
}


sub change_state {
    my ($self, $state) = @_;

    my $from_state = $self->current_state;
    my $verbose = $self->verbose;
    $verbose and print STDERR ">>>>> Changing state from $from_state to $state\n";

    if ( $self->validate_transition( $from_state, $state ) ) {

        # Now we know the current state is complete 
        # and all data from it is now in the input cache.
        my @data = $self->get_input_cache();

#	 if($from_state eq $state{iteration} ) {
#	   do{
#	     print STDERR "State change cache: ", $self->state_change_cache, "\n";
#	     print STDERR "Input cache ($from_state):\n@data\n\n";
#	 };
#	 }

        # Now we need to process the input cache data.
        # Remember, at this point, the current state is the "from" state
        # of the state transition. The "to" state doesn't get set until
        # the finalize_state_change() call at the end of this method.

        if($from_state eq $state{header} ) {
            $self->_process_header( @data );
        }
        elsif($from_state eq $state{descr} ) {
            $self->_process_descriptions( @data );
        }
        elsif($from_state eq $state{iteration} ) {
            $self->_process_iteration( @data, $self->state_change_cache() );
        }
        elsif($from_state eq $state{align} ) {
            $self->_process_alignment( @data );
        }
        elsif($from_state eq $state{footer} ) {
	  my $ok_to_pause = not $state eq $self->final_state;
	  $self->_process_footer( $ok_to_pause, @data );
        }

        $self->finalize_state_change( $state, 1 );
    }
}


sub _add_error {
    my ($self, $err) = @_;
    if( $err ) {
      push @{$self->{'_blast_errs'}}, $err;
    }
}

sub _clear_errors {
    my $self = shift;
    $self->{'_blast_errs'} = undef;
}

#---------
sub errors {
#---------
    my $self = shift;
    my @errs = ();
    @errs = @{$self->{'_blast_errs'}} if ref $self->{'_blast_errs'};
    return @errs;
}


#----------------------
sub _init_parse_params {
#----------------------
#Initializes parameters used during parsing of Blast reports.

    my ( $self, @param ) = @_;

    # -FILT_FUNC has been replaced by -HIT_FILTER.
    # Leaving -FILT_FUNC in place for backward compatibility
    my($signif, $filt_func, $hit_filter, $min_len, $check_all, $gapped, $score, 
       $overlap, $stats, $best, $shallow_parse, $hold_raw) =
	$self->_rearrange([qw(SIGNIF FILT_FUNC HIT_FILTER MIN_LEN CHECK_ALL_HITS GAPPED SCORE
			      OVERLAP STATS BEST SHALLOW_PARSE HOLD_RAW_DATA)], @param);

    $self->{'_hit_filter'} = $hit_filter || $filt_func || 0;
    $self->{'_check_all'} = $check_all || 0;
    $self->{'_shallow_parse'} = $shallow_parse || 0;
    $self->{'_hold_raw_data'} = $hold_raw || 0;

    $self->_set_signif($signif, $min_len, $self->{'_hit_filter'}, $score);
    $self->best_hit_only($best) if $best;
    $self->{'_blast_count'} = 0;

    $self->{'_collect_stats'} = defined($stats) ? $stats : 0;

    # TODO: Automatically determine whether gapping was used.
    # e.g., based on version number. Otherwise, have to read params.
    $self->{'_gapped'} = $gapped || 1;

    # Clear any errors from previous parse.
    $self->_clear_errors;
    undef $self->{'_hit_count'};
    undef $self->{'_num_hits_significant'};
}

#=head2 _set_signif
#
# Usage     : n/a; called automatically by _init_parse_params()
# Purpose   : Sets significance criteria for the BLAST object.
# Argument  : Obligatory three arguments:
#           :   $signif = float or sci-notation number or undef
#           :   $min_len = integer or undef
#           :   $hit_filter = function reference or undef
#           :
#           :   If $signif is undefined, a default value is set
#           :   (see $DEFAULT_SIGNIF; min_length = not set).
# Throws    : Exception if significance value is defined but appears
#           :   out of range or invalid.
#           : Exception if $hit_filter if defined and is not a func ref.
# Comments  : The significance of a BLAST report can be based on
#           : the P (or Expect) value and/or the length of the query sequence.
#           : P (or Expect) values GREATER than '_max_significance' are not significant.
#           : Query sequence lengths LESS than '_min_length' are not significant.
#           :
#           : Hits can also be screened using arbitrary significance criteria
#           : as discussed in the parse() method.
#           :
#           : If no $signif is defined, the '_max_significance' level is set to
#           : $DEFAULT_SIGNIF (999).
#
#See Also   : L<signif>(), L<min_length>(), L<_init_parse_params>(), L<parse>()
#
#=cut

#-----------------
sub _set_signif {
#-----------------
    my( $self, $sig, $len, $func, $score ) = @_;

    if(defined $sig) {
	$self->{'_confirm_significance'} = 1;
	if( $sig =~ /[^\d.e-]/ or $sig <= 0) {
	    $self->throw(-class => 'Bio::Root::BadParameter',
                         -text => "Invalid significance value: $sig\n".
			 "Must be a number greater than zero.");
	}
	$self->{'_max_significance'} = $sig;
    } else {
	$self->{'_max_significance'}   = $DEFAULT_SIGNIF;
    }

    if(defined $score) {
	$self->{'_confirm_significance'} = 1;
	if( $score =~ /[^\de+]/ or $score <= 0) {
	    $self->throw(-class => 'Bio::Root::BadParameter',
                         -text => "Invalid score value: $score\n".
			 "Must be an integer greater than zero.");
	}
	$self->{'_min_score'} = $score;
    } else {
	$self->{'_min_score'}  = $DEFAULT_SCORE;
    }

    if(defined $len) {
	if($len =~ /\D/ or $len <= 0) {
	    $self->warn("Invalid minimum length value: $len",
			"Value must be an integer > 0. Value not set.");
	} else {
	    $self->{'_min_length'} = $len;
	}
    }

    if(defined $func) {
        $self->{'_check_all'} = 1;
	$self->{'_confirm_significance'} = 1;
	if($func and not ref $func eq 'CODE') {
	    $self->throw("Not a function reference: $func",
			  "The -hit_filter parameter must be function reference.");
	  }
    }
}

=head2 signif

Synonym for L<max_significance()|max_significance>

=cut

#-----------
sub signif { shift->max_significance }


=head2 max_significance

 Usage     : $obj->max_significance();
 Purpose   : Gets the P or Expect value used as significance screening cutoff.
             This is the value of the -signif parameter supplied to new().
             Hits with P or E-value above this are skipped.
 Returns   : Scientific notation number with this format: 1.0e-05.
 Argument  : n/a
 Comments  : Screening of significant hits uses the data provided on the
           : description line. For NCBI BLAST1 and WU-BLAST, this data 
           : is P-value. for NCBI BLAST2 it is an Expect value.

=cut

#-----------
sub max_significance {
#-----------
    my $self = shift;
    my $sig = $self->{'_max_significance'};
    sprintf "%.1e", $sig;
}

=head2 min_score

 Usage     : $obj->min_score();
 Purpose   : Gets the Blast score used as screening cutoff.
             This is the value of the -score parameter supplied to new().
             Hits with scores below this are skipped.
 Returns   : Integer or scientific notation number.
 Argument  : n/a
 Comments  : Screening of significant hits uses the data provided on the
           : description line. 

=cut

#-----------
sub min_score {
#-----------
    my $self = shift;
    return $self->{'_min_score'};
}

=head2 min_length

 Usage     : $obj->min_length();
 Purpose   : Gets the query sequence length used as screening criteria.
             This is the value of the -min_len parameter supplied to new().
             Hits with sequence length below this are skipped.
 Returns   : Integer
 Argument  : n/a

See Also   : L<signif()|signif>

=cut

#--------------
sub min_length {
#--------------
    my $self = shift;
    $self->{'_min_length'};
}

=head2 highest_signif

 Usage     : $value = $obj->highest_signif();
 Purpose   : Gets the largest significance (P- or E-value) observed in
             the report.
           : For NCBI BLAST1 and WU-BLAST, this is a P-value. 
           : For NCBI BLAST2 it is an Expect value.
 Returns   : Float or sci notation number
 Argument  : n/a

=cut

sub highest_signif { shift->{'_highestSignif'} }

=head2 lowest_signif

 Usage     : $value = $obj->lowest_signif();
 Purpose   : Gets the largest significance (P- or E-value) observed in
             the report.
           : For NCBI BLAST1 and WU-BLAST, this is a P-value. 
           : For NCBI BLAST2 it is an Expect value.
 Returns   : Float or sci notation number
 Argument  : n/a

=cut

sub lowest_signif { shift->{'_lowestSignif'} }

=head2 highest_score

 Usage     : $value = $obj->highest_score();
 Purpose   : Gets the largest BLAST score observed in the report.
 Returns   : Integer or sci notation number
 Argument  : n/a

=cut

sub highest_score { shift->{'_highestScore'} }

=head2 lowest_score

 Usage     : $value = $obj->lowest_score();
 Purpose   : Gets the largest BLAST score observed in the report.
 Returns   : Integer or sci notation number
 Argument  : n/a

=cut

sub lowest_score { shift->{'_lowestScore'} }


# Throws : Exception if BLAST report contains a FATAL: error.
sub _process_header {
    my ($self, @data) = @_;

#    print STDERR "Processing Header...\n";

    $current_iteration = 0;
    $self->{'_result_count'}++;
    # Finish off the current Blast object, if any
    my $blast = $self->{'_current_blast'} = $self->result_factory->create_result();

    my ($set_method, $set_query, $set_db, $set_length);
    my ($query_start, $query_desc);
    
    foreach my $line (@data) {
        if( $line =~ /^(<.*>)?(T?BLAST[NPX])\s+(.*)$/ ) {
            $blast->analysis_method( $2 );
            $blast->analysis_method_version( $3 );
            $set_method = 1;
        }
        elsif ($line =~ /^Query= *(.+)$/ ) {
            $query_start = 1;
            my $info = $1;
            $info =~ s/TITLE //;
            # Split the query line into two parts.
            # Using \s instead of ' '
            $info =~ /(\S+?)\s+(.*)/;
            # set name of Blast object and return.
            $blast->query_name($1 || 'UNKNOWN');
            $query_desc = $2 || '';
            $set_query = 1;
        }
        elsif ($line =~ /^Database:\s+(.+)$/ ) {
            require Bio::Search::GenericDatabase;
            my $blastdb = Bio::Search::GenericDatabase->new( -name => $1 );
            $blast->analysis_subject( $blastdb );
            $set_db = 1;
        }
        elsif( $line =~ m/^\s+\(([\d|,]+) letters\)/ ) {
            my $length = $1;
            $length =~ s/,//g;
            $self->_set_query_length( $length );
            $set_length = 1;
            $blast->query_description( $query_desc );
            $query_start = 0;
        }
        elsif( $line =~ /WARNING: (.+?)/ ) {
            $self->warn( $1 );
        }
        elsif( $line =~ /FATAL: (.+?)/ ) {
            $self->throw("FATAL BLAST Report Error: $1");
        }
        # This needs to be the last elsif block.
        elsif( $query_start ) {
            # Handling multi-line query descriptions.
            chomp( $line );
            $query_desc .= " $line";
        }
    }
    if (!$set_method) {
        $self->throw("Can't determine type of BLAST program.");
    }
    if (!$set_query) {
        $self->throw("Can't determine name of query sequence.");
    }
    if(!$set_db) {
        $self->throw("Can't determine name of database.");
    }
    if(!$set_length) {
        $self->throw("Can't determine sequence length from BLAST report.");
    }

}

sub _process_descriptions {
    my ($self, @data) = @_;
#    print STDERR "Processing Descriptions...\n";

    # Step through each line parsing out the P/Expect value
    # All we really need to do is check the first one, if it doesn't
    # meet the significance requirement, we can skip the report.
    # BUT: we want to collect data for all hits anyway to get min/max signif.

    my $max_signif = $self->max_significance;
    my $min_score  = $self->min_score;
    my $layout_set = $self->{'_layout'} || 0;
    my ($layout, $sig, $hitid, $score, $is_p_value);

    if( $data[0] =~ /^\s*Sequences producing.*Score\s+P/i ) {
        $is_p_value = 1;
    } else {
        $is_p_value = 0;
    }

    my $hit_found_again;

    desc_loop:
  foreach my $line (@data) {
      last desc_loop if $line =~ / NONE |End of List|Results from round/;
      next desc_loop if $line =~ /^\.\./;

      if($line =~ /^Sequences used in model/ ) {
	#Sequences used in model and found again:
	$hit_found_again = 1;
	next;
      }
      elsif($line =~ /^Sequences not found previously/ ) {
	#Sequences not found previously or not previously below threshold:
	$hit_found_again  = 0;
	next;
      }

      ## Checking the significance value (P- or Expect value) of the hit
      ## in the description line.

      next desc_loop unless $line =~ /\d/;

      # TODO: These regexps need testing on a variety of reports.
      if ( $line =~ /^(\S+)\s+.*\s+([\de.+-]+)\s{1,5}[\de.-]+\s*$/) {
          $hitid = $1;
          $score = $2;
          $layout = 2;
      } elsif( $line =~ /^(\S+)\s+.*\s+([\de.+-]+)\s{1,5}[\de.-]+\s{1,}\d+\s*$/) {
          $hitid = $1;
          $score = $2;
          $layout = 1;
      } else {
	$self->warn("Can't parse description line\n $line");
	next desc_loop;
      }
      not $layout_set and ($self->_layout($layout), $layout_set = 1);

      $sig = &_parse_signif( $line, $layout, $self->{'_gapped'} );

#      print STDERR "  Parsed signif for $hitid: $sig (layout=$layout)\n";

      $self->{'_hit_hash'}->{$hitid}->{'signif'} = $sig;
      $self->{'_hit_hash'}->{$hitid}->{'score'} = $score;
      $self->{'_hit_hash'}->{$hitid}->{'found_again'} = $hit_found_again;
      $self->{'_hit_hash'}->{'is_pval'} = $is_p_value;

      last desc_loop if (not $self->{'_check_all'} and 
                         ($sig > $max_signif or $score < $min_score));

      $self->_process_significance($sig, $score);
    }

#  printf "\n%d SIGNIFICANT HITS.\nDONE PARSING DESCRIPTIONS.\n", $self->{'_num_hits_significant'};
}


#=head2 _set_query_length
#
# Usage     : n/a; called automatically during Blast report parsing.
# Purpose   : Sets the length of the query sequence (extracted from report).
# Returns   : integer (length of the query sequence)
# Throws    : Exception if cannot determine the query sequence length from
#           :           the BLAST report.
#           : Exception if the length is below the min_length cutoff (if any).
# Comments  : The logic here is a bit different from the other _set_XXX()
#           : methods since the significance of the BLAST report is assessed
#           : if MIN_LENGTH is set.
#
#=cut

#---------------
sub _set_query_length {
#---------------
    my ($self, $length) = @_;

    my($sig_len);
    if(defined($self->{'_min_length'})) {
      local $^W = 0;
      if($length < $self->{'_min_len'}) {
	$self->throw("Query sequence too short (Query= ${\$self->{'_current_blast'}->query_name}, length= $length)",
		     "Minimum  length is $self->{'_min_len'}");
      }
    }

    $self->{'_current_blast'}->query_length($length); 
}


# Records the highest and lowest significance (P- or E-value) and
# score encountered in a given report. 
sub _set_hi_low_signif_and_score {
    my($self, $sig, $score) = @_;

    my $hiSig = $self->{'_highestSignif'} || 0;
    my $lowSig = $self->{'_lowestSignif'} || $DEFAULT_SIGNIF;
    my $hiScore = $self->{'_highestScore'} || 0;
    my $lowScore = $self->{'_lowestScore'} || $DEFAULT_SIGNIF;

    $self->{'_highestSignif'} = ($sig > $hiSig)
   	                        ? $sig : $hiSig;

    $self->{'_lowestSignif'} = ($sig < $lowSig)
                                 ? $sig : $lowSig;

    $self->{'_highestScore'} = ($score > $hiScore)
   	                        ? $score : $hiScore;

    $self->{'_lowestScore'} = ($score < $lowScore)
                                 ? $score : $lowScore;
}


sub _process_significance {
    my($self, $sig, $score) = @_;

    $self->_set_hi_low_signif_and_score($sig, $score);

    # Significance value assessment.
    if($sig <= $self->{'_max_significance'} and $score >= $self->{'_min_score'}) {
        $self->{'_num_hits_significant'}++;
    }
    $self->{'_num_hits'}++;

    $self->{'_is_significant'} = 1 if $self->{'_num_hits_significant'};
}

#=head2 _layout
#
# Usage     : n/a; internal method.
# Purpose   : Set/Get indicator for the layout of the report.
# Returns   : Integer (1 | 2)
#           : Defaults to 2 if not set.
# Comments  : Blast1 and WashU-Blast2 have a layout = 1.
#           : This is intended for internal use by this and closely
#           : allied modules like BlastHit.pm and BlastHSP.pm.
#
#=cut

#------------
sub _layout {
#------------
    my $self = shift;
    if(@_) {
	$self->{'_layout'} = shift;
    }
    $self->{'_layout'} || 2;
}

#=head2 _parse_signif
#
# Usage     : $signif = _parse_signif(string, layout, gapped);
#           : This is a class function.
# Purpose   : Extracts the P- or Expect value from a single line of a BLAST description section.
# Example   : _parse_signif("PDB_UNIQUEP:3HSC_  heat-shock cognate ...   799  4.0e-206  2", 1);
#           : _parse_signif("gi|758803  (U23828) peritrophin-95 precurs   38  0.19", 2);
# Argument  : string = line from BLAST description section
#           : layout = integer (1 or 2)
#           : gapped = boolean (true if gapped Blast).
# Returns   : Float (0.001 or 1e-03)
# Status    : Static
#
#=cut

#------------------
sub _parse_signif {
#------------------
    my ($line, $layout, $gapped) = @_;

    local $_ = $line;
    my @linedat = split();

    # When processing both Blast1 and Blast2 reports
    # in the same run, offset needs to be configured each time.
    # NOTE: We likely won't be supporting mixed report streams. Not needed.

    my $offset  = 0;
    $offset  = 1 if $layout == 1 or not $gapped;

    my $signif = $linedat[ $#linedat - $offset ];

    # fail-safe check
    if(not $signif =~ /[.-]/) {
	$offset = ($offset == 0 ? 1 : 0);
	$signif = $linedat[ $#linedat - $offset ];
    }

    $signif = "1$signif" if $signif =~ /^e/i;
    return $signif;
}


=head2 best_hit_only

 Usage     : print "only getting best hit.\n" if $obj->best_hit_only();
 Purpose   : Set/Get the indicator for whether or not to processing only 
           : the best BlastHit.
 Returns   : Boolean (1 | 0)
 Argument  : n/a

=cut

#----------
sub best_hit_only {
#----------
    my $self = shift;
    if(@_) { $self->{'_best'} = shift; }
    $self->{'_best'};
}

sub _process_alignment {
    my ($self, @data) = @_;
#    print STDERR "Processing Alignment...\n";

    # If all of the significant hits have been parsed,
    # return if we're not checking all or if we don't need to get
    # the Blast stats (parameters at footer of report).
    if(defined $self->{'_hit_count'} and
      defined $self->{'_num_hits_significant'}) {
      return if $self->{'_hit_count'} >= $self->{'_num_hits_significant'} and
	not ($self->{'_check_all'} or $self->{'_collect_stats'});
    }

    # Return if we're only interested in the best hit.
    # This has to occur after checking for the parameters section
    # in the footer (since we may still be interested in them).
    return if $self->best_hit_only and ( defined $self->{'_hit_count'} and $self->{'_hit_count'} >=1);

    push @data, 'end';

#    print STDERR "\nALIGNMENT DATA:\n@data\n";

    my $max_signif  = $self->max_significance;
    my $min_score   = $self->min_score;

    my ($hitid, $score, $signif, $is_pval, $found_again);
    if( $data[0] =~ /^(\S+)\s+/ ) {
        $hitid = $1;
        return unless defined $self->{'_hit_hash'}->{$hitid};
        $score  = $self->{'_hit_hash'}->{$hitid}->{'score'};
        $signif = $self->{'_hit_hash'}->{$hitid}->{'signif'};
        $found_again = $self->{'_hit_hash'}->{$hitid}->{'found_again'};
        $is_pval =  $self->{'_hit_hash'}->{'is_pval'};
#        print STDERR "  Got hitid: $hitid ($signif, $score, P?=$is_pval)\n";
    }
    
    # Now construct the BlastHit objects from the alignment section

#	debug(1);

    $self->{'_hit_count'}++;

    # If not confirming significance, _process_descriptions will not have been run,
    # so we need to count the total number of hits here.
    if( not $self->{'_confirm_significance'}) {
      $self->{'_num_hits'}++;
    }

    my %hit_params = ( -RESULT     => $self->{'_current_blast'},
		       -RAW_DATA   =>\@data,
		       -SIGNIF     => $signif,
		       -IS_PVAL    => $is_pval,
		       -SCORE      => $score,
		       -RANK       => $self->{'_hit_count'},
		       -RANK_BY    => 'order',
		       -OVERLAP    => $self->{'_overlap'} || $MAX_HSP_OVERLAP,
		       -FOUND_AGAIN => $found_again,
		       -SHALLOW_PARSE  => $self->{'_shallow_parse'},
		       -HOLD_RAW_DATA  => $self->{'_hold_raw_data'},
		     );

    my $hit; 
    $hit = $self->hit_factory->create_hit( %hit_params );

    #printf STDERR "NEW HIT: %s, SIGNIFICANCE = %g\n", $hit->name, $hit->expect;  <STDIN>;
    # The BLAST report may have not had a description section.
    if(not $self->{'_has_descriptions'}) {
      $self->_process_significance($hit->signif, $score);
    }
    
    # Collect overall signif data if we don't already have it,
    # (as occurs if no -signif or -score parameter are supplied).
    my $hit_signif = $hit->signif;
    
    if (not $self->{'_confirm_significance'} ) {
	$self->_set_hi_low_signif_and_score($hit_signif, $score);
    }

    # Test significance using custom function (if supplied)
    my $add_hit = 0;

    my $hit_filter  = $self->{'_hit_filter'} || 0;

    if($hit_filter) {
        if(&$hit_filter($hit)) {
            $add_hit = 1;
        }
    } elsif($hit_signif <= $max_signif and $score >= $min_score) {
        $add_hit = 1;
    }

    $add_hit && $self->{'_current_blast'}->add_hit( $hit );
}


sub _process_footer {
    my ($self, $ok_to_pause, @data) = @_;
#    print STDERR "Processing Footer...\n";

    $self->{'_current_blast'}->raw_statistics( [@data] );

    if($self->{'_collect_stats'}) {
        foreach my $line (@data) {
           if( $line =~ /^\s*Matrix:\s*(\S+)/i ) {
                $self->{'_current_blast'}->matrix( $1 );
            }
            elsif( $line =~ /^\s*Number of Sequences:\s*(\d+)/i ) {
                $self->{'_current_blast'}->analysis_subject->entries( $1 );
            }
            elsif( $line =~ /^\s*length of database:\s*(\d+)/i ) {
                $self->{'_current_blast'}->analysis_subject->letters( $1 );
            }
            elsif( $line =~ /^\s*Posted date:\s*(.+)$/i ) {
                $self->{'_current_blast'}->analysis_subject->date( $1 );
            }
        }
    }

    if( $self->errors ) {
        my $num_err = scalar($self->errors);
        $self->warn( "$num_err Blast parsing errors occurred.");
	foreach( $self->errors ) { print STDERR "$_\n"; };
    }

    if( $self->{'_pause_between_reports'} and $ok_to_pause ) {
        $self->pause;
    }

}

sub _process_nohits {
    my $self = shift;
#    print STDERR "Processing No Hits (iteration = $current_iteration)\n";
    $self->{'_current_blast'}->set_no_hits_found( $current_iteration );
}


sub _process_iteration {
    my ($self, @data) = @_;
#    print STDERR "Processing Iteration\n";
#    print STDERR "   Incrementing current iteration (was=$current_iteration)\n";
    $current_iteration++;
    $self->{'_current_blast'}->iterations( $current_iteration );

    foreach( @data ) {
        if( /Results from round \d+/i ) {
	  $self->{'_current_blast'}->psiblast( 1 );
	}
        elsif( /No hits found/i ) {
            $self->_process_nohits();
            last;
        }
        elsif( /^\s*Sequences/i ) {
            $self->_process_descriptions( @data );
            last;
        }
    }
}

sub pause_between_reports {
    my ($self, $setting) = @_;
    if( defined $setting ) {
        $self->{'_pause_between_reports'} = $setting;
    }
    $self->{'_pause_between_reports'};
}

sub result_count {
    my $self = shift;
    return $self->{'_result_count'};
}

# For backward compatibility:
sub report_count { shift->result_count }

sub next_result {
   my ($self) = @_;
  # print STDERR "next_result() called\n";
   if( not $self->running ) {
       $self->run;
   }
   else {
       $self->resume;
   }
   my $blast = $self->{'_current_blast'};
   $self->{'_current_blast'} = undef;
   return $blast;
}

=head2 write_result

 Title   : write_result
 Usage   : $stream->write_result($result_result, @other_args)
 Function: Writes data from the $result_result object into the stream.
         : Delegates to the to_string() method of the associated 
         : WriterI object.
 Returns : 1 for success and 0 for error
 Args    : Bio::Search:Result::ResultI object,
         : plus any other arguments for the Writer
 Throws  : Bio::Root::Exception if a Writer has not been set.

See L<Bio::Root::Exception>

=cut

sub write_result {
   my ($self, $blast, @args) = @_;

   if( not defined($self->writer) ) {
       $self->warn("Writer not defined. Using a $DEFAULT_BLAST_WRITER_CLASS");
       $self->writer( $DEFAULT_BLAST_WRITER_CLASS->new() );
   }
   $self->SUPER::write_result( $blast, @args );
}



1;
__END__

