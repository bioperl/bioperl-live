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

@ISA = qw(Bio::Root::Root Bio::SearchIO::SearchWriterI);


=head2 new

 Title   : new
 Usage   : my $obj = new Bio::SearchIO::Writer::GbrowseGFF(@args);
 Function: Builds a new Bio::SearchIO::Writer::GbrowseGFF object 
 Returns : an instance of Bio::SearchIO::Writer::GbrowseGFF
 Args    :  -e_value => 10   : set e_value parsing cutoff (default 10)


=cut

sub new {
  my($class,@args) = @_;

  my $self = $class->SUPER::new(@args);
  my ($evalue) = $self->_rearrange(["E_VALUE"], @args);
  $self->{_evalue} = $evalue > 0?$evalue:10;
  return $self;
}


=head2 to_string

 Purpose   : Produce the Gbrowse format GFF lines for a Result
 Usage     : print $writer->to_string( $result_obj, @args);
 Argument  : $result_obj = A Bio::Search::Result::ResultI object
           : @args = none at present...
 Returns   : String containing data for each search Result or any of its
           : sub-objects (Hits and HSPs).
 Throws    : n/a

=cut

sub to_string {
    my ($self, $result, @args) = @_;
    my $GFF;
    while( my $hit = $result->next_hit ) {
        my $significance = $hit->significance;
        next unless ($significance < $self->{_evalue} &&  ($self->{_evalue} > 0)); 
        my $refseq = $hit->name;
        my $seqname = $hit->locus?$hit->locus:$hit->description;  # hopefully this will be a simple identifier without a full description line!!
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
            if ($hsp->strand('subject') eq "1"){
                push @plus_hsps, $hsp;
                if (defined $qpmin){  # set or reset the minimum and maximum extent of the plus-strand hit
                    $qpmin = $hsp->start if $hsp->start('query') < $qpmin;
                    $qpmax = $hsp->end if $hsp->end('query') > $qpmax;
                    $spmin = $hsp->start if $hsp->start('subject') < $spmin;
                    $spmax = $hsp->end if $hsp->end('subject') > $spmax;                    
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
                    $qmmin = $hsp->start if $hsp->start('query') < $qmmin;
                    $qmmax = $hsp->end if $hsp->end('query') > $qmmax;
                    $smmin = $hsp->start if $hsp->start('subject') < $smmin;
                    $smmax = $hsp->end if $hsp->end('subject') > $smmax;                    
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
        # okay, write out the index line for the entire hit before processing HSP's
        if (scalar(@plus_hsps)){
            $GFF .= "$refseq\t$source\tmatch\t$spmin\t$spmax\t$score\t+\t.\tTarget EST:$seqname $qpmin $qpmax\n";
        }
        if (scalar(@minus_hsps)){
            $GFF .= "$refseq\t$source\tmatch\t$smmin\t$smmax\t$score\t+\t.\tTarget EST:$seqname $qmmax $qmmin\n";  # note reversal of max and min in column 9, as per the spec
        }
        # process + strand hsps
        my $strand = "+";
        foreach my $hsp(@plus_hsps){
            my $qstart = $hsp->start('query');
            my $qend = $hsp->end('query');
            my $sstart = $hsp->start('subject');
            my $send = $hsp->end('subject');
            my $score = $hsp->score;
            $GFF .= "$refseq\t$source\tHSP\t$sstart\t$send\t$score\t+\t.\tTarget EST:$seqname $qstart $qend\n";
        }
        foreach my $hsp(@minus_hsps){
            my $qstart = $hsp->start('query');
            my $qend = $hsp->end('query');
            my $sstart = $hsp->start('subject');
            my $send = $hsp->end('subject');
            my $score = $hsp->score;
            $GFF .= "$refseq\t$source\tHSP\t$sstart\t$send\t$score\t-\t.\tTarget EST:$seqname $qend $qstart\n";  # note reversal of qstart/qend as per spec
        }
    }
    return $GFF;
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


