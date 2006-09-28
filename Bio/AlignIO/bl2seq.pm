# $Id$
#
# BioPerl module for Bio::AlignIO::bl2seq

#   based on the Bio::SeqIO modules
#       by Ewan Birney <birney@ebi.ac.uk>
#       and Lincoln Stein  <lstein@cshl.org>
#
#   the Bio::Tools::BPlite modules by
#   Ian Korf (ikorf@sapiens.wustl.edu, http://sapiens.wustl.edu/~ikorf),
#   Lorenz Pollak (lorenz@ist.org, bioperl port)
#
#       and the SimpleAlign.pm module of Ewan Birney
#
# Copyright Peter Schattner
#
# You may distribute this module under the same terms as perl itself
# _history
# September 5, 2000
# POD documentation - main docs before the code

=head1 NAME

Bio::AlignIO::bl2seq - bl2seq sequence input/output stream

=head1 SYNOPSIS

Do not use this module directly.  Use it via the L<Bio::AlignIO> class, as in:

    use Bio::AlignIO;

    $in  = Bio::AlignIO->new(-file   => "inputfilename" ,
                             -format => "bl2seq",
                             -report_type => "blastn");
    $aln = $in->next_aln();


=head1 DESCRIPTION

This object can create L<Bio::SimpleAlign> sequence alignment objects (of
two sequences) from C<bl2seq> BLAST reports.

A nice feature of this module is that - in combination with
L<Bio::Tools::Run::StandAloneBlast.pm> or a remote BLAST - it can be used to
align two sequences and make a L<Bio::SimpleAlign> object from them which
can then be manipulated using any L<Bio::SimpleAlign> methods, eg:

   # Get two sequences
   $str = Bio::SeqIO->new(-file=>'t/amino.fa' , '-format' => 'Fasta', );
   my $seq3 = $str->next_seq();
   my $seq4 = $str->next_seq();

   # Run bl2seq on them
   $factory = Bio::Tools::StandAloneBlast->new('program' => 'blastp',
                                               'outfile' => 'bl2seq.out');
   my $bl2seq_report = $factory->bl2seq($seq3, $seq4);
   # Note that report is a Bio::SearchIO object

   # Use AlignIO.pm to create a SimpleAlign object from the bl2seq report
   $str = Bio::AlignIO->new(-file=> 'bl2seq.out','-format' => 'bl2seq');
   $aln = $str->next_aln();

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists.  Your participation is much appreciated.

  bioperl-l@bioperl.org                  - General discussion
  http://bioperl.org/wiki/Mailing_lists  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via the
web:

  http://bugzilla.open-bio.org/

=head1 AUTHOR - Peter Schattner

Email: schattner@alum.mit.edu

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::AlignIO::bl2seq;
use strict;

use Bio::SearchIO;

use base qw(Bio::AlignIO);

=head2 new

 Title   : new
 Usage   : my $alignio = Bio::SimpleAlign->new(-format => 'bl2seq',
                                               -file   => 'filename',
                                               -report_type => 'blastx');
 Function: Get a L<Bio::SimpleAlign>
 Returns : L<Bio::SimpleAlign> object
 Args    : -report_type => report type (blastn,blastx,tblastx,tblastn,blastp)


=cut


sub new {
    my ($class) = shift;
    my $self = $class->SUPER::new(@_);
    my ($rt) = $self->_rearrange([qw(REPORT_TYPE)],@_);
    defined $rt && $self->report_type($rt);
    return $self;
}

=head2 next_aln

 Title   : next_aln
 Usage   : $aln = $stream->next_aln()
 Function: returns the next alignment in the stream.
 Returns : L<Bio::Align::AlignI> object on success,
           undef on error or end of file
 Args    : none

=cut

sub next_aln {
    my $self = shift;
    my $aln =  Bio::SimpleAlign->new(-source => 'bl2seq');
    $self->{'bl2seqobj'} = $self->{'bl2seqobj'} ||
	Bio::SearchIO->new(-fh => $self->_fh,
			   -format => 'blast');
    my $bl2seqobj = $self->{'bl2seqobj'};
    my $result = $self->{'_result'} || $bl2seqobj->next_result;
    $self->{'result'} = undef, return unless defined $result;

    my $hit = $self->{'_hit'} || $result->next_hit;
    $self->{'_hit'} = undef, return unless defined $hit;

    my $hsp  = $hit->next_hsp;
    return unless defined $hsp;
    return $hsp->get_aln;

# much easier above, eh?
#     my ($start,$end,$name,$seqname,$seq,$seqchar,$strand);
#     $seqchar = $hsp->query_string;
#     $start   = $hsp->query->start;
#     $end     = $hsp->query->end;
#      # Query name typically not present in bl2seq report
#     $seqname = $hsp->query->seq_id || 'Query-sequence';
#     $strand  = $hsp->query->strand;

#     #    unless ($seqchar && $start && $end  && $seqname) {return 0} ;
#     unless ($seqchar && $start && $end ) {return 0} ;

#     $seq = new Bio::LocatableSeq('-seq'   =>$seqchar,
# 				 '-id'    =>$seqname,
# 				 '-start' =>$start,
# 				 '-end'   =>$end,
# 				 '-strand'=>$strand,
# 				 );

#     $aln->add_seq($seq);

#     $seqchar  = $hsp->hit_string;
#     $start    = $hsp->hit->start;
#     $end      = $hsp->hit->end;
#     $seqname  = $hsp->hit->seq_id;
#     $strand   = $hsp->hit->strand;

#     unless ($seqchar && $start && $end  && $seqname) {return 0} ;

#     $seq = new Bio::LocatableSeq('-seq'   =>$seqchar,
# 				 '-id'    =>$seqname,
# 				 '-start' =>$start,
# 				 '-end'   =>$end,
# 				 '-strand'=>$strand,
# 				 );
#     $aln->add_seq($seq);
#     return $aln;
}


=head2 write_aln (NOT IMPLEMENTED)

 Title   : write_aln
 Usage   : $stream->write_aln(@aln)
 Function: writes the $aln object into the stream in bl2seq format
 Returns : 1 for success and 0 for error
 Args    : L<Bio::Align::AlignI> object


=cut

sub write_aln {
    my ($self,@aln) = @_;
    $self->throw_not_implemented();
}

=head2 report_type

 Title   : report_type
 Usage   : $obj->report_type($newval)
 Function: Sets the report type (blastn, blastp...)
 Returns : value of report_type (a scalar)
 Args    : on set, new value (a scalar or undef, optional)


=cut

sub report_type{
    my $self = shift;
    return $self->{'report_type'} = shift if @_;
    return $self->{'report_type'};
}

1;
