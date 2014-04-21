#-----------------------------------------------------------------
#
# BioPerl module Bio::SearchIO::Writer::GbrowseGFF.pm
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Mark Wilkinson <markw@illuminae.com>
#
# You may distribute this module under the same terms as perl itself
#-----------------------------------------------------------------

=head1 NAME

Bio::SearchIO::Writer::GbrowseGFF - Interface for outputting parsed search results in Gbrowse GFF format

=head1 SYNOPSIS

  use Bio::SearchIO;
  my $in = Bio::SearchIO->new(-file   => 'result.blast',      
                             -format => 'blast');
  my $out = Bio::SearchIO->new(-output_format  => 'GbrowseGFF',
                              -output_cigar   => 1,
                              -output_signif  => 1,
                              -file           => ">result.gff");
  while( my $r = $in->next_result ) {
    $out->write_result($r);
  }

=head1 DESCRIPTION

This writer produces Gbrowse flavour GFF from a Search::Result object.

=head1 AUTHOR  Mark Wilkinson

Email markw-at-illuminae-dot-com

=head1 CONTRIBUTORS

Susan Miller sjmiller at email-DOT-arizon-DOT-edu
Jason Stajich jason at bioperl-dot-org

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

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


package Bio::SearchIO::Writer::GbrowseGFF;
use vars qw(%Defaults);
use strict;

$Defaults{'Prefix'}   = 'EST';
$Defaults{'HSPTag'}   = 'HSP';
$Defaults{'MatchTag'} = 'match';

use base qw(Bio::Root::Root Bio::SearchIO::SearchWriterI);


=head2 new

 Title   : new
 Usage   : my $obj = Bio::SearchIO::Writer::GbrowseGFF->new(@args);
 Function: Builds a new Bio::SearchIO::Writer::GbrowseGFF object 
 Returns : an instance of Bio::SearchIO::Writer::GbrowseGFF
 Args    :  -e_value => 10   : set e_value parsing cutoff (default undef)
            (note the -e_value flag is deprecated.)

=cut

sub new {
    my($class,@args) = @_;

    my $self = $class->SUPER::new(@args);
    ($self->{'_evalue'},
     $self->{'_cigar'},
     $self->{'_prefix'},
     $self->{'_signif'} ) = $self->_rearrange([qw(E_VALUE OUTPUT_CIGAR PREFIX
						 OUTPUT_SIGNIF)], @args);
    $self->{'_evalue'} && warn( "Use of the -e_value argument is deprecated.\nIn future, use \$writer->filter(\"type\", \&code) instead.\n\tparsing will proceed correctly with this e_value\n");
    $self->{Gbrowse_HSPID} = 0;
    $self->{Gbrowse_HITID} = 0;
    $self->{'_prefix'} ||= $Defaults{'Prefix'};
    return $self;
}

sub _incrementHSP {
    my ($self) = @_;
    return ++$self->{Gbrowse_HSPID};
}

sub _incrementHIT {
    my ($self) = @_;
    return ++$self->{Gbrowse_HITID}
}
# according to the GFF3 spec:
#"match".  In addition to the generic "match"
#type, there are the subclasses "cDNA_match," "EST_match,"
#"translated_nucleotide_match," "nucleotide_to_protein_match," and
#"nucleotide_motif."

=head2 to_string

 Purpose   : Produce the Gbrowse format GFF lines for a Result
 Usage     : print $writer->to_string( $result_obj, @args);
 Argument  : $result_obj = A Bio::Search::Result::ResultI object
             -version => 1|2|2.5|3  ; the GFF format you want to output (default 3)
             -match_tag => match|cDNA_match|EST_match|translated_nucleotide_match
                           nucleotide_to_protein_match|nucleotide_motif
                           This is the SO term to be placed in GFF column 3.
             -prefix => String to prefix the group by, default is EST 
                        (see %Defaults class variable) A default can also
                        be set on object init
 Returns   : String containing data for each search Result or any of its
           : sub-objects (Hits and HSPs).
 Throws    : n/a

=cut

             #-reference => 'hit'|'query' ; whether the hit sequence name or the
             #                              query sequence name is used as the
             #                              'reference' sequence (GFF column 1)

sub to_string {
    my ($self, $result, @args) = @_;
    my ($format, $reference, 
	$match_tag,$hsp_tag,
	$prefix) = $self->_rearrange([qw
				      (VERSION 
				       REFERENCE 
				       MATCH_TAG HSP_TAG
				       PREFIX)], @args);
    $self->warn($reference) if $reference; 
    $reference ||='hit'; # default is that the hit sequence (db sequence) becomes the reference sequence.  I think this is fairly typical...
    $match_tag ||= $Defaults{'MatchTag'}; # default is the generic 'match' tag.
    $hsp_tag   ||= $Defaults{'HSPTag'}; # default is the generic 'hsp' tag.
    $prefix    ||= $self->{'_prefix'};
    $self->throw("$reference must be one of 'query', or 'hit'\n") unless $reference;
    
    #*************  THIS IS WHERE I STOPPED  ****************   
    # *****************************************************
    #*************************************************
    
    $format ||='3';
    my $gffio = Bio::Tools::GFF->new(-gff_version => $format); # try to set it
    
    # just in case that behaviour changes (at the moment, an invalid format throws an exception, but it might return undef in the future
    return "" unless defined $gffio;  # be kind and don't return undef in case the person is putting teh output directly into a printstatement without testing it
    # now $gffio is either false, or a valid GFF formatter

    my ($GFF,$cigar,$score);
    my ($resultfilter,$hitfilter,$hspfilter) = (
        $self->filter('RESULT'),
        $self->filter('HIT'),
        $self->filter('HSP'));
	$result->can('rewind') &&  $result->rewind(); # ensure we're at the beginning
    next if (defined $resultfilter && ! (&{$resultfilter}($result)) );

    while( my $hit = $result->next_hit ) {
        
        if (defined $self->{_evalue}){
            next unless ($hit->significance < $self->{_evalue});
        }
        next if( defined $hitfilter && ! &{$hitfilter}($hit) ); # test against filter code

        my $refseq  = $reference eq 'hit' ? $hit->name : $result->query_name;
        my $seqname = $reference eq 'hit' ? $result->query_name : $hit->name;  # hopefully this will be a simple identifier without a full description line!!
	if ($self->{_signif}) {
	    $score = $hit->significance;
	} else {
	    $score = $hit->raw_score;
	}
	$self->throw("No reference sequence name found in hit; required for GFF (this may not be your fault if your report type does not include reference sequence names)\n") unless $refseq;
        my $source = $hit->algorithm;
        $self->throw("No algorithm name found in hit; required for GFF (this may not be your fault if your report type does not include algorithm names)\n") unless $refseq;
        $self->throw("This module only works on BLASTN reports at this time.  Sorry.\n") unless $source eq "BLASTN";
        
        my @plus_hsps;
        my @minus_hsps;
        
        # pre-process the HSP's because we later need to know
        # the extents of the plus and munus strand
        # on both the subject and query strands individually
        my ($qpmin, $qpmax, $qmmin, $qmmax, $spmin, $spmax, $smmin, $smmax); # variables for the plus/minus strand min start and max end to know the full extents of the hit
        while( my $hsp = $hit->next_hsp ) {
            if ( defined $self->{_evalue} ) {  
                # for backward compatibility only
                next unless ($hsp->significance < $self->{_evalue});
	    }
            next if( defined $hspfilter && ! &{$hspfilter}($hsp) ); # test against HSP filter
            if ($hsp->hit->strand >= 0 ){
                push @plus_hsps, $hsp;
                if (defined $qpmin){  # set or reset the minimum and maximum extent of the plus-strand hit
                    $qpmin = $hsp->query->start if $hsp->query->start < $qpmin;
                    $qpmax = $hsp->query->end if $hsp->query->end > $qpmax;
                    $spmin = $hsp->hit->start if $hsp->hit->start < $spmin;
                    $spmax = $hsp->hit->end if $hsp->hit->end > $spmax;
                } else {
                    $qpmin = $hsp->query->start;
                    $qpmax = $hsp->query->end;
                    $spmin = $hsp->hit->start;
                    $spmax = $hsp->hit->end;
                }
            } 
            if ($hsp->hit->strand < 0 ){
                push @minus_hsps, $hsp;
                if (defined $qmmin){ # set or reset the minimum and maximum extent of the minus-strand hit
                    $qmmin = $hsp->query->start if $hsp->query->start < $qmmin;
                    $qmmax = $hsp->query->end if $hsp->query->end > $qmmax;
                    $smmin = $hsp->hit->start if $hsp->hit->start < $smmin;
                    $smmax = $hsp->hit->end if $hsp->hit->end > $smmax;
                } else {
                    $qmmin = $hsp->query->start;
                    $qmmax = $hsp->query->end;
                    $smmin = $hsp->hit->start;
                    $smmax = $hsp->hit->end;
                }
            }
            #else next if there is no strand, but that makes no sense..??
        }
        next unless (scalar(@plus_hsps) + scalar(@minus_hsps));  # next if no hsps (??)
        my $ID = $self->_incrementHIT();
        # okay, write out the index line for the entire hit before 
	# processing HSP's
        # unfortunately (or not??), HitI objects do not implement 
	# SeqFeatureI, so we can't just call ->gff_string
        # as a result, this module is quite brittle to changes 
	# in the GFF format since we are hard-coding the GFF output here :-(
	
        if (scalar(@plus_hsps)){
	    my %tags = ( 'ID' => "match_sequence$ID");

	    if ($format==2.5) {
		$tags{'Target'} = "$prefix:$seqname";
		$tags{'tstart'} = $qmmin;
		$tags{'tend'}   = $qmmax;
	    } else {
		$tags{'Target'} = "$prefix:$seqname $qpmin $qpmax";
	    }
	    if ( $self->{'_cigar'} ) {
		$tags{'Gap'} = $cigar;
	    }
            my $feat = Bio::SeqFeature::Generic->new(
                -seq_id      => $refseq,
                -source_tag  => $source,
                -primary_tag => $match_tag,
                -start       => $spmin,
                -end         => $spmax,
                -score       => $score,
                -strand      => '+',
                -frame       => '.',
                -tag         => \%tags 
            );


            my $formatter = Bio::Tools::GFF->new(-gff_version => $format);
            $GFF .= $feat->gff_string($formatter)."\n";
        }
        if (scalar(@minus_hsps)){
	    my %tags  = ( 'ID' => "match_sequence$ID");

            if ($format==2.5) {
                $tags{'Target'} = "$prefix:$seqname";
                $tags{'tstart'} = $qpmax;
                $tags{'tend'}   = $qpmin;
            }
            else {
                $tags{'Target'} = "$prefix:$seqname $qpmax $qpmin";
            }
            my $feat = Bio::SeqFeature::Generic->new(
                -seq_id      => $refseq,
                -source_tag  => $source,
                -primary_tag => $match_tag,
                -start       => $smmin,
                -end         => $smmax,
                -score       => $score,
                -strand      => '-',
                -frame       => '.',
                -tag         => \%tags 
	    );

            my $formatter = Bio::Tools::GFF->new(-gff_version => $format);
            $GFF .= $feat->gff_string($formatter)."\n";
        }
        
        # process + strand hsps
        foreach my $hsp (@plus_hsps){
            my $hspID  = $self->_incrementHSP();
            my $qstart = $hsp->query->start;
            my $qend   = $hsp->query->end;
            my $sstart = $hsp->hit->start;
            my $send   = $hsp->hit->end;
            my $score  = $hsp->score;
	    
	    my %tags  = ( 'ID'     => "match_hsp$hspID",
		          'Parent' => "match_sequence$ID" );
	    
            if ($format==2.5) {
                $tags{'Target'} = "$prefix:$seqname";
                $tags{'tstart'} = $qstart;
                $tags{'tend'}   = $qend;
            }
            else {
                $tags{'Target'} = "$prefix:$seqname $qstart $qend";
            }
	    if ( $self->{'_cigar'} ) {
		$tags{'Gap'} = $hsp->cigar_string;
	    }

            my $feat = Bio::SeqFeature::Generic->new(
                -seq_id      => $refseq,
                -source_tag  => $source,
                -primary_tag => $hsp_tag,
                -start       => $sstart,
                -end         => $send,
                -score       => $score,
                -strand      => '+',
                -frame       => '.',
                -tag         => \%tags 
            );

            my $formatter = Bio::Tools::GFF->new(-gff_version => $format);
            $GFF .= $feat->gff_string($formatter)."\n";
        }

        foreach my $hsp (@minus_hsps) {
            my $hspID  = $self->_incrementHSP();
            my $qstart = $hsp->query->start;
            my $qend   = $hsp->query->end;
            my $sstart = $hsp->hit->start;
            my $send   = $hsp->hit->end;
            my $score  = $hsp->score;

            my %tags  = ( 'ID'     => "match_hsp$hspID",
                          'Parent' => "match_sequence$ID" );

            if ($format==2.5) {
                $tags{'Target'} = "$prefix:$seqname";
                $tags{'tstart'} = $qend;
                $tags{'tend'}   = $qstart;
            }
            else {
                $tags{'Target'} = "$prefix:$seqname $qend $qstart";
            }
	    if ( $self->{'_cigar'} ) {
		$tags{'Gap'} = $hsp->cigar_string;
	    }

            my $feat = Bio::SeqFeature::Generic->new(
                -seq_id      => $refseq,
                -source_tag  => $source,
                -primary_tag => $hsp_tag,
                -start       => $sstart,
                -end         => $send,
                -score       => $score,
                -strand      => '-',
                -frame       => '.',
		-tag         => \%tags 
            );

            my $formatter = Bio::Tools::GFF->new(-gff_version => $format);
            $GFF .= $feat->gff_string($formatter) ."\n";
        }
    }
    return $GFF;
}

sub significance_filter {
    my ($self,$method,$code) = @_;    
    return unless $method;
    $method = uc($method);
    if( $method ne 'HSP' &&
	$method ne 'HIT' &&
	$method ne 'RESULT' ) {
	$self->warn("Unknown method $method");
	return;
    }
    if( $code )  {
	$self->throw("Must provide a valid code reference") unless ref($code) =~ /CODE/;
	$self->{$method} = $code;
    }
    return $self->{$method};
}

=head2 start_report

 Title   : start_report
 Usage   : $self->start_report()
 Function: has no function, returns nothing
 Returns : empty string
 Args    : none

=cut

sub start_report { return '' }

=head2 end_report

 Title   : end_report
 Usage   : $self->end_report()
 Function: has no function, returns nothing
 Returns : empty string
 Args    : none


=cut

sub end_report {  return '' }

=head2 filter

 Title   : filter
 Usage   : $writer->filter('hsp', \&hsp_filter);
 Function: Filter out either at HSP,Hit,or Result level
 Returns : none
 Args    : string => data type,
           CODE reference
 Note    : GbrowseGFF.pm makes no changes to the default filter code


=cut

1;


