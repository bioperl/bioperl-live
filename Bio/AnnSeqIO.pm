
#
# BioPerl module for Bio::SeqIO
#
# Cared for by Ewan Birney <birney@sanger.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::AnnSeqIO - Handler for AnnSeqIO Formats

=head1 SYNOPSIS

    use Bio::AnnSeqIO;

    $in  = Bio::AnnSeqIO->new(-file => "inputfilename" , -format => 'Fasta');
    $out = Bio::AnnSeqIO->new(-file => ">outputfilename" , -format => 'EMBL');

    while ($seq = $in->next_annseq() ) {
	$out->write_annseq($seq);
    }

=head1 DESCRIPTION

Bio::AnnSeqIO is a handler module for the formats in the AnnSeqIO set
(eg, Bio::AnnSeqIO::EMBL). It is the officially sanctioned way of
getting at the format objects, which most people should use.


The idea is that you request a stream object for a particular format.
All the stream objects have a notion of an internal file that is read
from or written to (the same object handles both input and output).
A physical example of a stream object is the Bio::AnnSeqIO::EMBL object.

Each stream object has functions

   $stream->next_annseq();

and

   $stream->write_annseq($seq);

also 

   $stream->type() # returns 'INPUT' or 'OUTPUT'


=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this
and other Bioperl modules. Send your comments and suggestions preferably
 to one of the Bioperl mailing lists.
Your participation is much appreciated.

  vsns-bcd-perl@lists.uni-bielefeld.de          - General discussion
  vsns-bcd-perl-guts@lists.uni-bielefeld.de     - Technically-oriented discussion
  http://bio.perl.org/MailList.html             - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
 the bugs and their resolution.
 Bug reports can be submitted via email or the web:

  bioperl-bugs@bio.perl.org
  http://bio.perl.org/bioperl-bugs/

=head1 AUTHOR - Ewan Birney

Email birney@sanger.ac.uk

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::AnnSeqIO;
use Exporter;
use vars qw($AUTOLOAD @ISA);
use strict;

@ISA = qw(Exporter);

# _initialize is where the heavy stuff will happen when new is called

=head2 new

 Title   : new
 Usage   : $stream = Bio::AnnSeqIO->new(-file => $filename, -format => 'Format')
 Function: Returns a new seqstream
 Returns : A Bio::AnnSeqIO::* initialsed with the appropiate format
 Args    : -file => $filename 
           -format => format
           -fh => filehandle to attach to
           


=cut

my $entry = 0;

sub new {
   my ($class,%param) = @_;
   my ($format);
   my ($handler,$stream);
   # we can ignore class!

   $format = $param{'-format'} || $param{'-FORMAT'};

   if( &_load_format_module($format) == 0 ) {
       return undef;
   }

   $stream = "Bio::AnnSeqIO::$format"->new(%param);
   return $stream;
}


=head2 _load_format_module

 Title   : _load_format_module
 Usage   : *INTERNAL AnnSeqIO stuff*
 Function: Loads up (like use) a module at run time on demand
 Example :
 Returns : 
 Args    :


=cut

sub _load_format_module{
   my ($format) = @_;
   my ($module,$load,$m);

   $module = "_<Bio/AnnSeqIO/$format.pm";
   $load = "Bio/AnnSeqIO/$format.pm";

#   foreach $m ( keys %main:: ) {
#       print "Got [$m] vs [$module]\n";
#       if( $m eq $module ) {
#	   print STDERR "Found it!\n";
#       }
#   }

   if( $main::{$module} ) {
     #  print STDERR "$module Module loaded already..\n";
       return 1;
   } else {
       eval {
	   require $load;
       };
       if( $@ ) {
	   print STDERR "$load: $format cannot be found\nException $@\nFor more information about the AnnSeqIO system please see the AnnSeqIO docs.\nThis includes ways of checking for formats at compile time, not run time\n";
	   return 0;
       }
       return 1;
   }

}


1;







