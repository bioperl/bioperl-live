#!perl -w
use strict;
# $Id: split_seq.PLS,v 1.1 2003-05-12 14:58:16 birney Exp $


=head1 NAME

split_seq - splits a sequence into equal sized chunks

=head1 SYNOPSIS

split_seq -chunk 10000 seq.in

=head1 DESCRIPTION 

The script will split sequences into chunks

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

my ($format,$chunk,$do_file);

$chunk = 10000;
$do_file =0;

GetOptions(
	   'chunk:i' => \$chunk,
	   'format:s'  => \$format,
	   'file' => \$do_file
	);

my $oformat = 'fasta';

# this implicity uses the <> file stream
my $seqin = Bio::SeqIO->new( -format => $format); 
my $seqout = Bio::SeqIO->new( -format => $oformat, -file => ">-" );


while( (my $seq = $seqin->next_seq()) ) {
   for(my $i = 1;$i < $seq->length;$i+=$chunk) {
       my $id = $seq->id;
       $id .= ".$i";
       if( $do_file ) {
	   $seqout = Bio::SeqIO->new( -format => $oformat, -file => ">$id.seq" );
       }
       my $seqtrunc = $seq->trunc($i,$i+$chunk-1);
       $seqtrunc->id($id);
       $seqout->write_seq($seqtrunc);
   }
}
