#
# $Id$
#
# BioPerl module for Bio::DB::EMBL
#
# Cared for by Heikki Lehvaslaiho <Heikki@ebi.ac.uk>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::DB::RefSeq - Database object interface for RefSeq retrieval

=head1 SYNOPSIS

  use Bio::DB::RefSeq;

  $db = new Bio::DB::RefSeq;

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
  }
  print "accesion is ", $seq->accession_number, "\n" unless $@;

  # or ... best when downloading very large files, prevents
  # keeping all of the file in memory

  # also don't want features, just sequence so let's save bandwith
  # and request Fasta sequence
  $db = new Bio::DB::RefSeq(-retrievaltype => 'tempfile' ,
 			       -format => 'fasta');
   my $seqio = $embl->get_Stream_by_batch(['NM_006732', 'NM_005252'] );
  while( my $seq  =  $seqio->next_seq ) {
 	print "seqid is ", $clone->id, "\n";
  }

=head1 DESCRIPTION


Allows the dynamic retrieval of sequence objects L<Bio::Seq> from the
RefSeq database using the dbfetch script at EBI:
L<http:E<sol>E<sol>www.ebi.ac.ukE<sol>cgi-binE<sol>dbfetch>.

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

  bioperl-l@bioperl.org                         - General discussion
  http://bio.perl.org/MailList.html             - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
 the bugs and their resolution.
 Bug reports can be submitted via email or the web:

  bioperl-bugs@bio.perl.org
  http://bio.perl.org/bioperl-bugs/

=head1 AUTHOR - Heikki Lehvaslaiho

Email Heikki Lehvaslaiho E<lt>Heikki@ebi.ac.ukE<gt>

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::DB::RefSeq;
use strict;
use vars qw(@ISA $MODVERSION %HOSTS  %FORMATMAP $DEFAULTFORMAT);

$MODVERSION = '0.1';
use Bio::DB::DBFetch;

@ISA = qw(Bio::DB::DBFetch);

BEGIN {
    # you can add your own here theoretically.
    %HOSTS = (
	       'dbfetch' => {
		   baseurl => 'http://%s/cgi-bin/dbfetch?db=refseq&style=raw',
		   hosts   => {
		       'ebi'  => 'www.ebi.ac.uk'
		       }
	       }
	      );
    %FORMATMAP = ( 'genbank' => 'genbank',
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

1;
