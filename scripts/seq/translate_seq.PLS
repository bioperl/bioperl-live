#!perl -w
use strict;
# $Id: translate_seq.PLS,v 1.1 2003-05-12 14:58:16 birney Exp $


=head1 NAME

translate_seq - translates a sequence

=head1 SYNOPSIS

translate_seq < cdna_cds.fa > protein.fa

=head1 DESCRIPTION 

The script will translate one fasta file (on stdin) to protein on stdout

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

 http://bugzilla.bioperl.org/

=head1 AUTHOR

  Ewan Birney <birney@ebi.ac.uk>

=cut

use Bio::SeqIO;
use Getopt::Long;

my ($format);


GetOptions(
	   'format:s'  => \$format,
	);

my $oformat = 'fasta';

# this implicity uses the <> file stream
my $seqin = Bio::SeqIO->new( -format => $format); 
my $seqout = Bio::SeqIO->new( -format => $oformat, -file => ">-" );


while( (my $seq = $seqin->next_seq()) ) {
	my $pseq = $seq->translate();
	$seqout->write_seq($pseq);
}
