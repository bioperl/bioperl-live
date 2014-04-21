#
#
# BioPerl module for Bio::DB::EMBL
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Heikki Lehvaslaiho <heikki-at-bioperl-dot-org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::DB::RefSeq - Database object interface for RefSeq retrieval

=head1 SYNOPSIS

  use Bio::DB::RefSeq;

  $db = Bio::DB::RefSeq->new();

  # most of the time RefSeq_ID eq RefSeq acc
  $seq = $db->get_Seq_by_id('NM_006732'); # RefSeq ID
  print "accession is ", $seq->accession_number, "\n";

  # or changeing to accession number and Fasta format ...
  $db->request_format('fasta');
  $seq = $db->get_Seq_by_acc('NM_006732'); # RefSeq ACC
  print "seq is ", $seq->seq, "\n";

  # especially when using versions, you better be prepared
  # in not getting what what want
  eval {
      $seq = $db->get_Seq_by_version('NM_006732.1'); # RefSeq VERSION
  };
  print "accesion is ", $seq->accession_number, "\n" unless $@;

  # or ... best when downloading very large files, prevents
  # keeping all of the file in memory

  # also don't want features, just sequence so let's save bandwith
  # and request Fasta sequence
  $db = Bio::DB::RefSeq->new(-retrievaltype => 'tempfile' ,
 			       -format => 'fasta');
  my $seqio = $db->get_Stream_by_id(['NM_006732', 'NM_005252'] );
  while( my $seq  =  $seqio->next_seq ) {
 	print "seqid is ", $seq->id, "\n";
  }

=head1 DESCRIPTION

Allows the dynamic retrieval of sequence objects L<Bio::Seq> from the
RefSeq database using the dbfetch script at EBI:

http://www.ebi.ac.uk/Tools/dbfetch/dbfetch

In order to make changes transparent we have host type (currently only
ebi) and location (defaults to ebi) separated out.  This allows later
additions of more servers in different geographical locations.

The functionality of this module is inherited from L<Bio::DB::DBFetch>
which implements L<Bio::DB::WebDBSeqI>.

This module retrieves entries from EBI although it
retrives database entries produced at NCBI. When read into bioperl
objects, the parser for GenBank format it used. RefSeq is a
NONSTANDARD GenBank file so be ready for surprises.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists.  Your participation is much appreciated.

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
the bugs and their resolution.  Bug reports can be submitted via the
web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Heikki Lehvaslaiho

Email Heikki Lehvaslaiho E<lt>heikki-at-bioperl-dot-orgE<gt>

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::DB::RefSeq;
use strict;
use vars qw($MODVERSION %HOSTS %FORMATMAP $DEFAULTFORMAT);

$MODVERSION = '0.1';

use base qw(Bio::DB::DBFetch);

BEGIN {
    # you can add your own here theoretically.
    %HOSTS = (
	       'dbfetch' => {
		   baseurl => 'http://%s/Tools/dbfetch/dbfetch?db=refseq&style=raw',
		   hosts   => {
		       'ebi'  => 'www.ebi.ac.uk'
		       }
	       }
	      );
    %FORMATMAP = ( 'embl'    => 'embl',
		   'genbank' => 'genbank',
		   'fasta' => 'fasta'
		   );
    $DEFAULTFORMAT = 'genbank';
}

sub new {
    my ($class, @args ) = @_;
    my $self = $class->SUPER::new(@args);

    $self->{ '_hosts' } = {};
    $self->{ '_formatmap' } = {};

    $self->hosts(\%HOSTS);
    $self->formatmap(\%FORMATMAP);
    $self->{'_default_format'} = $DEFAULTFORMAT;

    return $self;
}


sub get_seq_stream {
   my ($self,%qualifiers) = @_;
   if( exists $qualifiers{'-uids'} ) {
       if( ref($qualifiers{'-uids'}) =~ /ARRAY/i ) {
	   foreach my $u ( @{$qualifiers{'-uids'}} ) {
	       $u =~ s/^(\S+)\|//;
	   }
       } else { 
	   $qualifiers{'-uids'} =~ s/^(\S+)\|//;
       }
   }
   $self->SUPER::get_seq_stream(%qualifiers);
}

1;
