# BioPerl module for remote_blast.pl
#
# Cared for by Jason Stajich <jason@chg.mc.duke.edu>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs after the code

#!/usr/local/bin/perl -w

use strict;
use vars qw($USAGE);

use Bio::Tools::Run::RemoteBlast;
use Bio::SeqIO;
use Getopt::Long;

$USAGE = "remote_blast.pl [-h] [-p prog] [-d db] [-e expect] [-mod BPlite|Blast] [-f seqformat] -i seqfile\n";

my ($prog, $db, $expect,$method) = ( 'blastp', 'nr', '10', 'BPlite');

my ($sequencefile,$sequenceformat,$help) = (undef, 'fasta',undef);

&GetOptions('prog|p=s'               => \$prog,
	    'db|d=s'                 => \$db,
	    'expect|e=s'             => \$expect,
	    'blsmod|module|method=s' => \$method,
	    'input|i=s'              => \$sequencefile,  	    
	    'format|f=s'             => \$sequenceformat,
	    'help|h'                 => \$help, 
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

my $blastfactory = new Bio::Tools::Run::RemoteBlast ('-prog'      => $prog,
						     '-data'      => $db,
						     '-expect'    => $expect,
						     'readmethod' => $method,
						     );

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
$blastfactory->submit_blast($input);
print "retrieving blasts..";
while( my @rids = $blastfactory->each_rid ) {
    foreach my $rid ( @rids ) {
	my $rc = $blastfactory->retrieve_blast($rid);
	if( !ref($rc) ) {
	    if( $rc < 0 ) { 
		# there was an error and this RID is not coming back
		$blastfactory->remove_rid($rid);
	    }
	    print ".";
	    sleep 5;
	} else { 
	    print "\n";
	    $blastfactory->remove_rid($rid);
	    while( my $sbjct = $rc->nextSbjct ) {		
		print "sbjct name is ", $sbjct->name, "\n";
		while( my $hsp = $sbjct->nextHSP ) {
		    print "score is ", $hsp->score, "\n";
		} 
	    }
	}
    }
}

__END__

#
# BioPerl module for remote_blast.pl
#
# Cared for by Jason Stajich <jason@chg.mc.duke.edu>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

remote_blast.pl - script for submitting jobs to a remote blast server
(ncbi blast queue at this time)

=head1 SYNOPSIS

% remote_blast.pl -p blastp -d ecoli -e 1e-5 -i myseqs.fa

=head1 DESCRIPTION

This module will run a remote blast on a set of sequences by
submitting them to the NCBI blast queue and printing the output of the
request.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

  bioperl-l@bioperl.org               - General discussion
  http://bioperl.org/MailList.shtml   - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via email
or the web:

  bioperl-bugs@bioperl.org
  http://bioperl.org/bioperl-bugs/

=head1 AUTHOR - Jason Stajich

Email jason@chg.mc.duke.edu

=cut
