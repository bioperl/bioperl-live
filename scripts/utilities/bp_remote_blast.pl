#!/usr/bin/perl
#
# BioPerl module for bp_remote_blast.pl
#
# Revived by Evan Weaver for bioperl-1.5.1
# 3/14/2006
#
# Copyright Jason Stajich, Evan Weaver
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs after the code

use strict;
use warnings;
use vars qw($USAGE);

use Bio::Tools::Run::RemoteBlast;
use Bio::SeqIO;
use Getopt::Long;

$USAGE = "remote_blast.pl [-h] [-p prog] [-d db] [-e expect] [-mod
Blast] [-f seqformat] -z=\"entrez query\" -v 1 -t output_format -i
seqfile\n";

my ($prog, $db, $expect,$method) = ( 'blastp', 'nr', '10', 'Blast');

my ($sequencefile,$sequenceformat,$help, $entrez, $outputformat,
$verbose) = (undef, 'fasta',undef, undef, undef, 1);

&GetOptions('prog|p=s'               => \$prog,
				'db|d=s'                 => \$db,
				'expect|e=s'             => \$expect,
				'blsmod|module|method=s' => \$method,
				'input|i=s'              => \$sequencefile,  	
				'format|f=s'             => \$sequenceformat,
				'help|h'                 => \$help,
				'entrez|z=s'             => \$entrez,
				'output_format|t=s'      => \$outputformat,
				'verbose|v=s'            => \$verbose
			  );

if( $help ) {
    exec('perldoc', $0);
    die;
}

if( !defined $prog ) {
    die($USAGE . "\n\tMust specify a valid program name ([t]blast[pxn])\n");
}
if( !defined $db ) {
    die($USAGE . "\n\tMust specify a db to search\n");
}
if( !defined $sequencefile ) {
    die($USAGE . "\n\tMust specify an input file\n");
}

my $blastfactory = new Bio::Tools::Run::RemoteBlast ('-prog' => $prog,
																	  '-data'      => $db,
																	  '-expect'    => $expect,
																	  'readmethod' => $method,
																	 );

if ($entrez) {
  if ($verbose) {
    print "Entrez query (submission side): $entrez\n";
  }
  #$Bio::Tools::Run::RemoteBlast::RETRIEVALHEADER{ FORMAT_ENTREZ_QUERY} = $entrez;
  $Bio::Tools::Run::RemoteBlast::HEADER{ ENTREZ_QUERY } = $entrez;
}
if ($outputformat) {
  print "Don't use output format type; it doesn't work.\n";
  $Bio::Tools::Run::RemoteBlast::RETRIEVALHEADER{ FORMAT_TYPE } = $outputformat;
}

# submit_blast can only currenly handle fasta format files so I'll
# preprocess outside of the module but I'd rather be sure here

my $input;
if( $sequenceformat !~ /fasta/ ) {
	my @seqs;
	my $seqio = new Bio::SeqIO('-format' => $sequenceformat,
										'-file'   => $sequencefile );
	while( my $seq = $seqio->next_seq() ) {
		push @seqs, $seq;
	}
	$input = \@seqs;
} else {
	$input = $sequencefile;
}

my $r = $blastfactory->submit_blast($input);
#my $r = $factory->submit_blast(?amino.fa?);

print STDERR "waiting...\n" if( $verbose > 0 );
while ( my @rids = $blastfactory->each_rid ) {
  foreach my $rid ( @rids ) {
    my $rc = $blastfactory->retrieve_blast($rid);
    if( !ref($rc) ) {
      if( $rc < 0 ) {
        $blastfactory->remove_rid($rid);
      }
      print STDERR "   checking $rid\n" if ( $verbose > 0 );
      sleep 5;
    } else {
      my $result = $rc->next_result();
      #save the output
      my $filename = $result->query_name()."\.out";
      $blastfactory->save_output($filename);
      $blastfactory->remove_rid($rid);
      print "\nQuery Name: ", $result->query_name(), "\n";
      while ( my $hit = $result->next_hit ) {
        next unless ( $verbose > 0);
        print "\thit name is ", $hit->name, "\n";
        while( my $hsp = $hit->next_hsp ) {
          print "\t\tscore is ", $hsp->score, "\n";
        }
      }
    }
  }
  print STDERR scalar(@rids) . " left\n";
}
	

__END__

#
# BioPerl module for bp_remote_blast.pl
#
# Cared for by Jason Stajich <jason-at-bioperl-dot-org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

bp_remote_blast.pl - script for submitting jobs to a remote blast server
(ncbi blast queue at this time)

=head1 SYNOPSIS

% bp_remote_blast.pl -p blastp -d ecoli -e 1e-5 -i myseqs.fa

=head1 DESCRIPTION

This module will run a remote blast on a set of sequences by
submitting them to the NCBI blast queue and printing the output of the
request.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

  bioperl-l@bioperl.org                  - General discussion
  http://bioperl.org/wiki/Mailing_lists  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via the
web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Jason Stajich

Email jason-at-bioperl-dot-org

=cut
