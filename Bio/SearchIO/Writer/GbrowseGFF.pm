#-----------------------------------------------------------------
# $Id$
#
# BioPerl module Bio::SearchIO::Writer::GbrowseGFF.pm
#
# Cared for by Mark Wilkinson <markw@illuminae.com>
#
# You may distribute this module under the same terms as perl itself
#-----------------------------------------------------------------

=head1 NAME

Bio::SearchIO::Writer::GbrowseGFF - Interface for outputting parsed search results in Gbrowse GFF format

=head1 SYNOPSIS

  use Bio::SearchIO;
  my $in = new Bio::SearchIO(-file   => 'result.blast',
                             -format => 'blast');
  my $out = new Bio::SearchIO(-output_format  => 'GbrowseGFF',
                              -file           => ">result.gff");
  while( my $r = $in->next_result ) {
    $out->write_result($r);
  }

=head1 DESCRIPTION

This writer produces Gbrowse flavour GFF from a Search::Result object.

=head1 AUTHOR  Mark Wilkinson

Email markw-at-illuminae-dot-com

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

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


package Bio::SearchIO::Writer::GbrowseGFF;
use Bio::SearchIO::SearchWriterI;
use Bio::Root::RootI;
use vars qw(@ISA);
use strict;

@ISA = qw(Bio::Root::Root Bio::SearchIO::SearchWriterI);


=head2 new

 Title   : new
 Usage   : my $obj = new Bio::SearchIO::Writer::GbrowseGFF(@args);
 Function: Builds a new Bio::SearchIO::Writer::GbrowseGFF object 
 Returns : an instance of Bio::SearchIO::Writer::GbrowseGFF
 Args    :  -e_value => 10   : set e_value parsing cutoff (default undef)
            (note the -e_value flag is deprecated.)

=cut

sub new {
    my($class,@args) = @_;

    my $self = $class->SUPER::new(@args);
    my ($evalue) = $self->_rearrange(["E_VALUE"], @args);
    $self->{_evalue} = $evalue;
    $evalue && print STDERR 'use of the -e_value argument is deprecated.  In future, use $writer->filter("type", \&code)  instead.\n\tparsing will proceed correctly with this e_value\n';
    $self->{Gbrowse_HSPID} = 0;
    $self->{Gbrowse_HITID} = 0;

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
 Returns   : String containing data for each search Result or any of its
           : sub-objects (Hits and HSPs).
 Throws    : n/a

=cut

             #-reference => 'hit'|'query' ; whether the hit sequence name or the
             #                              query sequence name is used as the
             #                              'reference' sequence (GFF column 1)

sub to_string {
    my ($self, $result, @args) = @_;
    my ($format, $reference, $match_tag) = $self->_rearrange(["VERSION", "REFERENCE", "MATCH_TAG"], @args);
    $reference ||='hit'; # default is that the hit sequence (db sequence) becomes the reference sequence.  I think this is fairly typical...
    $match_tag ||='match'; # default is the generic 'match' tag.
    $self->throw("$reference must be one of 'query', or 'hit'\n") unless $reference;
    
    #*************  THIS IS WHERE I STOPPED  ****************   
    # *****************************************************
    #*************************************************
    
    $format ||='3';
    my $gffio = Bio::Tools::GFF->new(-gff_version => $format); # try to set it
    
    # just in case that behaviour changes (at the moment, an invalid format throws an exception, but it might return undef in the future
    return "" unless defined $gffio;  # be kind and don't return undef in case the person is putting teh output directly into a printstatement without testing it
    # now $gffio is either false, or a valid GFF formatter

    my $GFF;
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

        my $refseq = $hit->name;
        my $seqname = $result->query_name;  # hopefully this will be a simple identifier without a full description line!!
        my $score = $hit->raw_score;
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
            if (defined $self->{_evalue}){  # for backward compatibility only
                next unless ($hsp->significance < $self->{_evalue});
            }
            next if( defined $hspfilter && ! &{$hspfilter}($hsp) ); # test against HSP filter
            if ($hsp->strand('subject') eq "1"){
                push @plus_hsps, $hsp;
                if (defined $qpmin){  # set or reset the minimum and maximum extent of the plus-strand hit
                    $qpmin = $hsp->start('query') if $hsp->start('query') < $qpmin;
                    $qpmax = $hsp->end('query') if $hsp->end('query') > $qpmax;
                    $spmin = $hsp->start('subject') if $hsp->start('subject') < $spmin;
                    $spmax = $hsp->end('subject') if $hsp->end('subject') > $spmax;                    
                } else {
                    $qpmin = $hsp->start('query');
                    $qpmax = $hsp->end('query');
                    $spmin = $hsp->start('subject');
                    $spmax = $hsp->end('subject');
                }
            } 
            if ($hsp->strand('subject') eq "-1"){
                push @minus_hsps, $hsp;
                if (defined $qmmin){ # set or reset the minimum and maximum extent of the minus-strand hit
                    $qmmin = $hsp->start('query') if $hsp->start('query') < $qmmin;
                    $qmmax = $hsp->end('query') if $hsp->end('query') > $qmmax;
                    $smmin = $hsp->start('subject') if $hsp->start('subject') < $smmin;
                    $smmax = $hsp->end('subject') if $hsp->end('subject') > $smmax;                    
                } else {
                    $qmmin = $hsp->start('query');
                    $qmmax = $hsp->end('query');
                    $smmin = $hsp->start('subject');
                    $smmax = $hsp->end('subject');
                }
            }
            #else next if there is no strand, but that makes no sense..??
        }
        next unless (scalar(@plus_hsps) + scalar(@minus_hsps));  # next if no hsps (??)
        my $ID = $self->_incrementHIT();
        # okay, write out the index line for the entire hit before processing HSP's
        # unfortunately (or not??), HitI objects do not implement SeqFeatureI, so we can't just call ->gff_string
        # as a result, this module is quite brittle to changes in the GFF format since we are hard-coding the GFF output here :-(
        if (scalar(@plus_hsps)){
            my $feat = Bio::SeqFeature::Generic->new(
                -seq_id      => $refseq,
                -source_tag  => $source,
                -primary_tag => $match_tag,
                -start       => $spmin,
                -end         => $spmax,
                -score      => $score,
                -strand     => '+',
                -frame      => '.',
                -tag         => {
                    'ID'    => "match_sequence$ID",
                    'Target' => (($format==2.5)?"EST:$seqname":"EST:$seqname+$qpmin+$qpmax"),
                    'tstart' => $qpmin,
                    'tend' => $qpmax,
                }
                                         );
            my $formatter = Bio::Tools::GFF->new(-gff_version => $format);
            $GFF .= $feat->gff_string($formatter)."\n";
            #    
            #    $GFF .= _GFF25($refseq,$source,$match_tag,$spmin,$spmax,$score,'+','.',$seqname,"match_sequence$ID",$qpmin,$qpmax);
            #} elsif ($format == 2){
            #    $GFF .= _GFF2($refseq,$source,$match_tag,$spmin,$spmax,$score,'+','.',$seqname,"match_sequence$ID",$qpmin,$qpmax);
            #} elsif ($format == 1){
            #    $GFF .= _GFF1($refseq,$source,$match_tag,$spmin,$spmax,$score,'+','.',$seqname,"match_sequence$ID",$qpmin,$qpmax);
            #} else { #$format == 3
            #    $GFF .= _GFF3($refseq,$source,$match_tag,$spmin,$spmax,$score,'+','.',$seqname,"match_sequence$ID",$qpmin,$qpmax);
            #}                
        }
        if (scalar(@minus_hsps)){
            my $feat = Bio::SeqFeature::Generic->new(
                -seq_id      => $refseq,
                -source_tag  => $source,
                -primary_tag => $match_tag,
                -start       => $smmin,
                -end         => $smmax,
                -score      => $score,
                -strand     => '-',
                -frame      => '.',
                -tag         => {
                    'ID'    => "match_sequence$ID",
                    'Target' => (($format==2.5)?"EST:$seqname":"EST:$seqname+$qmmax+$qmmin"),
                    'tstart' => $qmmax,
                    'tend' => $qmmin,
                }
                                         );
            my $formatter = Bio::Tools::GFF->new(-gff_version => $format);
            $GFF .= $feat->gff_string($formatter)."\n";
#
#            if ($format == 2.5){
#                $GFF .= _GFF25($refseq,$source,$match_tag,$smmin,$smmax,$score,'-','.',$seqname,"match_sequence$ID",$qmmax,$qmmin);
#            } elsif ($format ==2){
#                $GFF .= _GFF2($refseq,$source,$match_tag,$smmin,$smmax,$score,'-','.',$seqname,"match_sequence$ID",$qmmax,$qmmin);
#            } elsif ($format == 1){
#                $GFF .= _GFF1($refseq,$source,$match_tag,$smmin,$smmax,$score,'-','.',$seqname,"match_sequence$ID",$qmmax,$qmmin);
#            } else { #$format == 3
#                $GFF .= _GFF3($refseq,$source,$match_tag,$smmin,$smmax,$score,'-','.',$seqname,"match_sequence$ID",$qmmax,$qmmin);
#            }
        }
        # process + strand hsps
        foreach my $hsp(@plus_hsps){
            my $hspID = $self->_incrementHSP();
            my $qstart = $hsp->start('query');
            my $qend = $hsp->end('query');
            my $sstart = $hsp->start('subject');
            my $send = $hsp->end('subject');
            my $score = $hsp->score;
            my $feat = Bio::SeqFeature::Generic->new(
                -seq_id      => $refseq,
                -source_tag  => $source,
                -primary_tag => $match_tag,
                -start       => $sstart,
                -end         => $send,
                -score      => $score,
                -strand     => '+',
                -frame      => '.',
                -tag         => {
                    'ID'    => "match_hsp$hspID",
                    'Target' => (($format==2.5)?"EST:$seqname":"EST:$seqname+$qstart+$qend"),
                    'tstart' => $qstart,
                    'tend' => $qend,
                    'Parent' => "match_sequence$ID",
                }
                                         );
            my $formatter = Bio::Tools::GFF->new(-gff_version => $format);
            $GFF .= $feat->gff_string($formatter)."\n";
            #if ($format == 2.5){
            #    $GFF .= _GFF25($refseq,$source,$match_tag,$sstart,$send,$score,'+','.',$seqname,"match_hsp$hspID",$qstart,$qend,"match_sequence$ID");
            #} elsif ($format ==2){
            #    $GFF .= _GFF2($refseq,$source,$match_tag,$sstart,$send,$score,'+','.',$seqname,"match_hsp$hspID",$qstart,$qend,"match_sequence$ID");
            #} elsif ($format == 1){
            #    $GFF .= _GFF1($refseq,$source,$match_tag,$sstart,$send,$score,'+','.',$seqname,"match_hsp$hspID",$qstart,$qend,"match_sequence$ID");
            #} else { #$format == 3
            #    $GFF .= _GFF3($refseq,$source,$match_tag,$sstart,$send,$score,'+','.',$seqname,"match_hsp$hspID",$qstart,$qend,"match_sequence$ID");
            #}
        }
        foreach my $hsp(@minus_hsps){
            my $hspID = $self->_incrementHSP();
            my $qstart = $hsp->start('query');
            my $qend = $hsp->end('query');
            my $sstart = $hsp->start('subject');
            my $send = $hsp->end('subject');
            my $score = $hsp->score;
            my $feat = Bio::SeqFeature::Generic->new(
                -seq_id      => $refseq,
                -source_tag  => $source,
                -primary_tag => $match_tag,
                -start       => $sstart,
                -end         => $send,
                -score      => $score,
                -strand     => '-',
                -frame      => '.',
                -tag         => {
                    'ID'    => "match_hsp$hspID",
                    'Target' => (($format==2.5)?"EST:$seqname":"EST:$seqname+$qend+$qstart"),
                    'tstart' => $qend,
                    'tend' => $qstart,
                    'Parent' => "match_sequence$ID",
                }
                                         );
            my $formatter = Bio::Tools::GFF->new(-gff_version => $format);
            $GFF .= $feat->gff_string($formatter) ."\n";
            #if ($format == 2.5){
            #    $GFF .= _GFF25($refseq,$source,$match_tag,$sstart,$send,$score,'-','.',$seqname,"match_hsp$hspID",$qend,$qstart,"match_sequence$ID");
            #} elsif ($format ==2){
            #    $GFF .= _GFF2($refseq,$source,$match_tag,$sstart,$send,$score,'-','.',$seqname,"match_hsp$hspID",$qend,$qstart,"match_sequence$ID");
            #} elsif ($format == 1){
            #    $GFF .= _GFF1($refseq,$source,$match_tag,$sstart,$send,$score,'-','.',$seqname,"match_hsp$hspID",$qend,$qstart,"match_sequence$ID");
            #} else { #$format == 3
            #    $GFF .= _GFF3($refseq,$source,$match_tag,$sstart,$send,$score,'-','.',$seqname,"match_hsp$hspID",$qend,$qstart,"match_sequence$ID");
            #}
        }
    }
    return $GFF;
}

#sub _GFF1 {
#    my ($refseq, $source, $match_tag, $start, $stop, $score, $strand, $frame, $EST, $ID, $hitmin, $hitmax, $parentID) = @_;
#    return "$refseq\t$source\t$match_tag\t$start\t$stop\t$score\tstrand\t$frame";
#    
#}
#
#sub _GFF2 {
#    my ($refseq, $source, $match_tag, $start, $stop, $score, $strand, $frame, $EST, $ID, $hitmin, $hitmax, $parentID) = @_;
#    return"$refseq\t$source\t$match_tag\t$start\t$stop\t$score\t$strand\t$frame\tTarget EST:$EST ; tstart $hitmin ; tend $hitmax\n";
#    
#}
#
#sub _GFF25 {
#    my ($refseq, $source, $match_tag, $start, $stop, $score, $strand, $frame, $EST, $ID, $hitmin, $hitmax, $parentID) = @_;
#    return "$refseq\t$source\t$match_tag\t$start\t$stop\t$score\t$strand\t$frame\tTarget EST:$EST ; tstart $hitmin ; tend $hitmax\n";
#
#}
#
#sub _GFF3 {
#    my ($refseq, $source, $match_tag, $start, $stop, $score, $strand, $frame, $EST, $ID, $hitmin, $hitmax, $parentID) = @_;
#    return"$refseq\t$source\t$match_tag\t$start\t$stop\t$score\t$strand\t$frame\tID=$ID ; Target=EST:$EST+$hitmin+$hitmax".($parentID?" ; Parent=$parentID":"")."\n";
#    
#}
#
sub significance_filter {
    my ($self,$method,$code) = @_;    
    return undef unless $method;
    $method = uc($method);
    if( $method ne 'HSP' &&
	$method ne 'HIT' &&
	$method ne 'RESULT' ) {
	$self->warn("Unknown method $method");
	return undef;
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


