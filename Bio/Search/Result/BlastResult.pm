#-----------------------------------------------------------------
# $Id$
#
# BioPerl module Bio::Search::Result::BlastResult
#
# Cared for by Steve Chervitz <sac@bioperl.org>
#
# You may distribute this module under the same terms as perl itself
#-----------------------------------------------------------------

# POD documentation - main docs before the code

=head1 NAME

Bio::Search::Result::BlastResult - A top-level BLAST Report object

=head1 SYNOPSIS

The construction of BlastResult objects is performed by
by the B<Bio::SearchIO::psiblast> parser.
Therefore, you do not need to
use B<Bio::Search::Result::BlastResult>) directly. If you need to construct
BlastHits directly, see the new() function for details.

For B<Bio::SearchIO> BLAST parsing usage examples, see the
B<examples/search-blast> directory of the Bioperl distribution.

=head1 DESCRIPTION

This module supports BLAST versions 1.x and 2.x, gapped and ungapped,
and PSI-BLAST.

=head1 DEPENDENCIES

Bio::Search::Result::BlastResult.pm is a concrete class that inherits from B<Bio::Root::Root> and B<Bio::Search::Result::ResultI>. It  relies on two other modules:

=over 4

=item B<Bio::Search::Hit::BlastHit> 

Encapsulates a single a single BLAST hit.

=item B<Bio::Search::GenericDatabase>

Provides an interface to a blast database metadata.

=back

=head1 FEEDBACK

=head2 Mailing Lists 

User feedback is an integral part of the evolution of this and other
Bioperl modules.  Send your comments and suggestions preferably to one
of the Bioperl mailing lists.  Your participation is much appreciated.

    bioperl-l@bioperl.org              - General discussion
    http://bio.perl.org/MailList.html  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution. Bug reports can be submitted via email
or the web:

    bioperl-bugs@bio.perl.org                   
    http://bugzilla.bioperl.org/           

=head1 AUTHOR 

Steve Chervitz E<lt>sac@bioperl.orgE<gt>

See L<the FEEDBACK section | FEEDBACK> for where to send bug reports and comments.

=head1 ACKNOWLEDGEMENTS

This software was originally developed in the Department of Genetics
at Stanford University. I would also like to acknowledge my
colleagues at Affymetrix for useful feedback.

=head1 COPYRIGHT

Copyright (c) 2001 Steve Chervitz. All Rights Reserved.

=cut

=head1 DISCLAIMER

This software is provided "as is" without warranty of any kind.

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...

package Bio::Search::Result::BlastResult;

use strict;

use Bio::Search::Result::ResultI;
use Bio::Root::Root;

use overload 
    '""' => \&to_string;

use vars qw(@ISA $Revision );

$Revision = '$Id$';  #'
@ISA = qw( Bio::Root::Root Bio::Search::Result::ResultI);

#----------------
sub new {
#----------------
    my ($class, @args) = @_; 
    my $self = $class->SUPER::new(@args);
    return $self;
}

#sub DESTROY {
#    my $self = shift;
#    print STDERR "->DESTROYING $self\n";
#}


#=================================================
# Begin Bio::Search::Result::ResultI implementation
#=================================================

=head2 next_hit

See L<Bio::Search::Result::ResultI::next_hit()|Bio::Search::Result::ResultI> for documentation

=cut

#----------------
sub next_hit {
#----------------
    my ($self) = @_;

    unless($self->{'_hit_queue_started'}) {
        $self->{'_hit_queue'} = [$self->hits()];
        $self->{'_hit_queue_started'} = 1;
    }

    pop @{$self->{'_hit_queue'}};
}

=head2 query_name

See L<Bio::Search::Result::ResultI::query_name()|Bio::Search::Result::ResultI> for documentation

=cut

#----------------
sub query_name {
#----------------
    my $self = shift;
    if (@_) { 
        my $name = shift;
        $name =~ s/^\s+|(\s+|,)$//g;
        $self->{'_query_name'} = $name;
    }
    return $self->{'_query_name'};
}

=head2 query_length

See L<Bio::Search::Result::ResultI::query_length()|Bio::Search::Result::ResultI> for documentation

=cut

#----------------
sub query_length {
#----------------
    my $self = shift;
    if(@_) { $self->{'_query_length'} = shift; }
    return $self->{'_query_length'};
}

=head2 query_description

See L<Bio::Search::Result::ResultI::query_description()|Bio::Search::Result::ResultI> for documentation

=cut

#----------------
sub query_description {
#----------------
    my $self = shift;
    if(@_) { 
        my $desc = shift;
        defined $desc && $desc =~ s/(^\s+|\s+$)//g;
        # Remove duplicated ID at beginning of description string
        defined $desc && $desc =~ s/^$self->{'_query_name'}//o;
        $self->{'_query_query_desc'} = $desc || '';
    }
    return $self->{'_query_query_desc'};
}


=head2 analysis_method

See L<Bio::AnalysisResultI::analysis_method()|Bio::AnalysisResultI> for documentation

This implementation ensures that the name matches /blast/i.

=cut

#----------------
sub analysis_method { 
#----------------
    my ($self, $method) = @_;  
    if($method ) {
      if( $method =~ /blast/i) {
	$self->{'_analysis_prog'} = $method;
      } else {
	$self->throw("method $method not supported in " . ref($self));
      }
    }
    return $self->{'_analysis_prog'}; 
}

=head2 analysis_method_version

See L<Bio::AnalysisResultI::analysis_method_version()|Bio::AnalysisResultI> for documentation

=cut

#----------------
sub analysis_method_version {
#----------------
    my ($self, $version) = @_; 
    if($version) {
	$self->{'_analysis_progVersion'} = $version;
    }
    return $self->{'_analysis_progVersion'}; 
}


=head2 analysis_query

See L<Bio::AnalysisResultI::analysis_query()|Bio::AnalysisResultI> for documentation

=cut

#----------------
sub analysis_query {
#----------------

    my ($self) = @_;
    if(not defined $self->{'_analysis_query'}) {
        require Bio::PrimarySeq;
        my $moltype = $self->analysis_method =~ /blastp|tblastn/i ? 'protein' : 'dna';
	$self->{'_analysis_query'} =  Bio::PrimarySeq->new( -display_id => $self->query_name,
                                                            -desc => $self->query_description,
                                                            -moltype => $moltype
                                                          );
        $self->{'_analysis_query'}->length( $self->query_length );
    }
    return $self->{'_analysis_query'};
}

=head2 analysis_subject

 Usage     : $blastdb = $result->analyis_subject();
 Purpose   : Get a Bio::Search::DatabaseI object containing
             information about the database used in the BLAST analysis.
 Returns   : Bio::Search::DatabaseI object.
 Argument  : n/a

=cut

#---------------
sub analysis_subject { 
#---------------
    my ($self, $blastdb) = @_; 
    if($blastdb) {
        if( ref $blastdb and $blastdb->isa('Bio::Search::DatabaseI')) {
            $self->{'_analysis_sbjct'} = $blastdb;
        }
        else {
            $self->throw(-class =>'Bio::Root::BadParameter',
                         -text => "Can't set BlastDB: not a Bio::Search::DatabaseI $blastdb"
                         );
        }
    }
    return $self->{'_analysis_sbjct'};
}

=head2 next_feature

 Title   : next_feature
 Usage   : while( my $feat = $blast_result->next_feature ) { # do something }
 Function: Returns the next feature available in the analysis result, or
           undef if there are no more features.
 Example :
 Returns : A Bio::SeqFeatureI compliant object, in this case, 
           each Bio::Search::HSP::BlastHSP object within each BlastHit.
 Args    : None

=cut

#---------------
sub next_feature{
#---------------
   my ($self) = @_;
   my ($hit, $hsp);
   $hit = $self->{'_current_hit'};
   unless( defined $hit ) {
       $hit = $self->{'_current_hit'} = $self->next_hit;
       return undef unless defined $hit;
   }
   $hsp = $hit->next_hsp;
   unless( defined $hsp ) {
       $self->{'_current_hit'} = undef;
       return $self->next_feature;
   }
   return $hsp || undef;
}


sub algorithm { shift->analysis_method( @_ ); }
sub algorithm_version { shift->analysis_method_version( @_ ); }

# No-ops for now...
sub available_parameters{ return ''; }
sub get_parameter{ return ''; }
sub available_statistics{ return ''; }
sub get_statistic{ return ''; }

#=================================================
# End Bio::Search::Result::ResultI implementation
#=================================================


=head2 to_string

 Title   : to_string
 Usage   : print $blast->to_string;
 Function: Returns a string representation for the Blast result. 
           Primarily intended for debugging purposes.
 Example : see usage
 Returns : A string of the form:
           [BlastResult] <analysis_method> query=<name> <description> db=<database
           e.g.:
           [BlastResult] BLASTP query=YEL060C vacuolar protease B, db=PDBUNIQ 
 Args    : None

=cut

#---------------
sub to_string {
#---------------
    my $self = shift;
    my $str = "[BlastResult] " . $self->analysis_method . " query=" . $self->query_name . " " . $self->query_description .", db=" . $self->database_name;
    return $str;
}

#---------------
sub database_name {
#---------------
    my $self = shift;
    my $dbname = '';
    if( ref $self->analysis_subject) {
      $dbname = $self->analysis_subject->name;
    } 
    return $dbname;
}

#---------------
sub database_entries {
#---------------
    my $self = shift;
    my $dbentries = '';
    if( ref $self->analysis_subject) {
      $dbentries = $self->analysis_subject->entries;
    } 
    return $dbentries;
}

#---------------
sub database_letters {
#---------------
    my $self = shift;
    my $dbletters = '';
    if( ref $self->analysis_subject) {
      $dbletters = $self->analysis_subject->letters;
    } 
    return $dbletters;
}

#---------------
sub hits {
#---------------
    my $self = shift;
    my @hits = ();
    if( ref $self->{'_hits'}) {
        @hits = @{$self->{'_hits'}};
    }
    return @hits;
}

=head2 add_hit

 Usage     : $blast->add_hit( $hit );
 Purpose   : Adds a hit object to the collection of hits in this BLAST result.
 Returns   : n/a
 Argument  : A Bio::Search::Hit::HitI object
 Comments  : For PSI-BLAST, hits from all iterations are lumped together.
             For any given hit, you can determine the iteration in which it was
             found by checking $hit->iteration().

=cut

#---------------
sub add_hit {
#---------------
    my ($self, $hit) = @_;
    my $add_it = 1;

    unless( ref $hit and $hit->isa('Bio::Search::Hit::HitI')) {
        $add_it = 0;
        $self->throw(-class =>'Bio::Root::BadParameter',
                     -text => "Can't add hit: not a Bio::Search::Hit::HitI: $hit"
                    );
    }

    # Avoid adding duplicate hits if we're doing multiple iterations (PSI-BLAST)
#    if( $self->iterations > 1 ) {
#        my $hit_name = $hit->name;
#        if( grep $hit_name eq $_, @{$self->{'_hit_names'}}) {
#            $add_it = 0;
#        }
#    }

    if( $add_it ) {
        push @{$self->{'_hits'}}, $hit;
        push @{$self->{'_hit_names'}}, $hit->name;
    }
}


=head2 is_signif

 Usage     : $blast->is_signif();
 Purpose   : Determine if the BLAST report contains significant hits.
 Returns   : Boolean
 Argument  : n/a
 Comments  : BLAST reports without significant hits but with defined
           : significance criteria will throw exceptions during construction.
           : This obviates the need to check significant() for
           : such objects.

=cut

#------------
sub is_signif { my $self = shift; return $self->{'_is_significant'}; }
#------------


=head2 matrix

 Usage     : $blast_object->matrix();
 Purpose   : Get the name of the scoring matrix used.
           : This is extracted from the report.
 Argument  : n/a
 Returns   : string or undef if not defined
 Comments  : TODO: Deprecate this and implement get_parameter('matrix').

=cut

#------------
sub matrix { 
#------------
    my $self = shift; 
    if(@_) {
         $self->{'_matrix'} = shift;
    }
    $self->{'_matrix'};
}


=head2 raw_statistics

 Usage     : @stats = $blast_result->raw_statistics();
 Purpose   : Get the raw, unparsed statistical parameter section of the Blast report.
             This is the section at the end after the last HSP alignment.
 Argument  : n/a
 Returns   : Array of strings

=cut

#------------
sub raw_statistics { 
#------------
    my $self = shift; 
    if(@_) {
         my $params = shift;
         if( ref $params eq 'ARRAY') {
              $self->{'_raw_statistics'} = $params;
         }
        else {
            $self->throw(-class =>'Bio::Root::BadParameter',
                         -text => "Can't set statistical params: not an ARRAY ref: $params"
                         );
        }
    }
    if(not defined $self->{'_raw_statistics'}) {
              $self->{'_raw_statistics'} = [];
    }

    @{$self->{'_raw_statistics'}};
}



=head2 no_hits_found

 Usage     : $nohits = $blast->no_hits_found( [iteration_number] ); 
 Purpose   : Get boolean indicator indicating whether or not any hits
             were present in the report.

             This is NOT the same as deteriming the number of hits via
             the hits() method, which will return zero hits if there were no
             hits in the report or if all hits were filtered out during the parse.

             Thus, this method can be used to distinguish these possibilities
             for hitless reports generated when filtering.

 Returns   : Boolean
 Argument  : (optional) integer indicating the iteration number (PSI-BLAST)
             If iteration number is not specified and this is a PSI-BLAST result,
             then this method will return true only if all iterations had
             no hits found.

=cut

#-----------
sub no_hits_found {
#-----------
    my ($self, $round) = @_;

    my $result = 0;   # final return value of this method.
    # Watch the double negative! 
    # result = 0 means "yes hits were found"
    # result = 1 means "no hits were found" (for the indicated iteration or all iterations)

    # If a iteration was not specified and there were multiple iterations,
    # this method should return true only if all iterations had no hits found.
    if( not defined $round ) {
        if( $self->{'_iterations'} > 1) {
            $result = 1;
            foreach my $i( 1..$self->{'_iterations'} ) {
                if( not defined $self->{"_iteration_$i"}->{'_no_hits_found'} ) {
                    $result = 0;
                    last;
                }
            }
        }
        else {
            $result = $self->{"_iteration_1"}->{'_no_hits_found'};
        }
    }
    else {
        $result = $self->{"_iteration_$round"}->{'_no_hits_found'};
    }

    return $result;
}


=head2 set_no_hits_found

 Usage     : $blast->set_no_hits_found( [iteration_number] ); 
 Purpose   : Set boolean indicator indicating whether or not any hits
             were present in the report.
 Returns   : n/a
 Argument  : (optional) integer indicating the iteration number (PSI-BLAST)

=cut

#-----------
sub set_no_hits_found {
#-----------
    my ($self, $round) = @_;
    $round ||= 1;
    $self->{"_iteration_$round"}->{'_no_hits_found'} = 1;
}


=head2 iterations

 Usage     : $num_iterations = $blast->iterations;  (get)
             $blast->iterations($num_iterations);   (set)
 Purpose   : Set/get the number of iterations in the Blast Report (PSI-BLAST).
 Returns   : Total number of iterations in the report
 Argument  : integer  (when setting)

=cut

#----------------
sub iterations {
#----------------
    my ($self, $num ) = @_;
    if( defined $num ) {
        $self->{'_iterations'} = $num;
    }
    return $self->{'_iterations'};
}


=head2 psiblast

 Usage     : if( $blast->psiblast ) { ... }
 Purpose   : Set/get a boolean indicator whether or not the report 
             is a PSI-BLAST report.
 Returns   : 1 if PSI-BLAST, undef if not.
 Argument  : 1 (when setting)

=cut

#----------------
sub psiblast {
#----------------
    my ($self, $val ) = @_;
    if( $val ) {
        $self->{'_psiblast'} = 1;
    }
    return $self->{'_psiblast'};
}


1;
__END__
