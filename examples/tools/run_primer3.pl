#!/usr/bin/perl

=head1 NAME

run_primer3.pl - run primer3 and parse its output

=head1 SYNOPSIS

  ./run_primer3.pl -i test.fa

  #or

  ./run_primer3.pl --input=test.fa

=head1 DESCRIPTION

Example of how to run primer3 and parse its output, essentially taken from an
email written by Paul Wiersma to bioperl-l.

=head1 FEEDBACK

User feedback is an integral part of the evolution of this and other
Bioperl scripts. Send your comments and suggestions to the Bioperl 
mailing list.  Your participation is much appreciated.

  bioperl-l@bioperl.org                  - General discussion
  http://bioperl.org/wiki/Mailing_lists  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution. Bug reports can be submitted via the
web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR

Brian Osborne, bosborne at alum.mit.edu

=cut

use strict;
use Getopt::Long;
use Bio::Tools::Run::Primer3;
use Bio::SeqIO;

my $in_file;

GetOptions("i|input:s" => \$in_file );

usage() unless $in_file;

my $seqio = Bio::SeqIO->new(-file => $in_file);

while (my $seq = $seqio->next_seq) {
	my $primer3 = Bio::Tools::Run::Primer3->new(-seq => $seq);
	$primer3->program_name('primer3_core') unless $primer3->executable;

	$primer3->add_targets('PRIMER_MIN_TM' => 56, 'PRIMER_MAX_TM' => 90);

	my $results = $primer3->run;

	unless ($results->number_of_results) {
		print "No results for ",$seq->display_id;
		next;
	}

	my @out_keys_part = qw(START
								    LENGTH
								    TM
									 GC_PERCENT
									 SELF_ANY
									 SELF_END
									 SEQUENCE );

	print "\n", $seq->display_id, "\n";

	for (my $i = 0 ; $i < $results->number_of_results ; $i++){
		my $result = $results->primer_results($i);
 
		print "\n", $i + 1;	
		for my $key qw(PRIMER_LEFT PRIMER_RIGHT){	
			my ($start, $length) = split /,/, $result->{$key};
			$result->{$key . "_START"} = $start;
			$result->{$key . "_LENGTH"} = $length;
			foreach my $partkey (@out_keys_part) {
				print "\t", $result->{$key . "_" . $partkey};
			} 
			print "\n";
		}
		print "\tPRODUCT SIZE: ", $result->{'PRIMER_PRODUCT_SIZE'}, ", PAIR ANY COMPL: ",
		  $result->{'PRIMER_PAIR_COMPL_ANY'};
		print ", PAIR 3\' COMPL: ", $result->{'PRIMER_PAIR_COMPL_END'}, "\n";
	}
}

sub usage {
   exec('perldoc',$0);
   exit(0);
}

__END__
