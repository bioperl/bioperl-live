# $Id$
# BioPerl module for Bio::ClusterIO::dbsnp
#
# Copyright Allen Day <allenday@ucla.edu>, Stan Nelson <snelson@ucla.edu>
# Human Genetics, UCLA Medical School, University of California, Los Angeles

# POD documentation - main docs before the code

=head1 NAME

Bio::ClusterIO::dbsnp - dbSNP input stream

=head1 SYNOPSIS

Do not use this module directly.  Use it via the Bio::ClusterIO class.

=head1 DESCRIPTION

Parse dbSNP files, one refSNP at a time.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists.  Your participation is much appreciated.

  bioperl-l@bioperl.org            - General discussion
  http://bioperl.org/MailList.shtml - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
 the bugs and their resolution.
 Bug reports can be submitted via email or the web:

  bioperl-bugs@bio.perl.org
  http://bio.perl.org/bioperl-bugs/

=head1 AUTHOR

Allen Day <allenday@ucla.edu>

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...
package Bio::ClusterIO::dbsnp;
use vars qw(@ISA);
use strict;

use Bio::ClusterIO;

#$::RD_TRACE = 1;
$::RD_WARN = 1;
$::RD_HINT = 1;

@ISA = qw(Bio::ClusterIO);

=head2 next_cluster

 Title   : next_cluster
 Usage   : $unigene = $stream->next_cluster()
 Function: returns the next refSNP in the stream
 Returns : Bio::Variation::SNP object
 Args    : NONE

=cut

sub next_cluster {
	my $self = shift;
	my($b,$e,$xml_rs);

	while($b = $self->_readline){next unless $b =~ m!<NBE-rs>!}
	$xml_rs = $b;
	while($e = $self->_readline){$xml_rs .= $e; last if $e =~ m!</NBE-rs>!}

	my @xml_ss = $xml_rs =~ s!(<NSE-ss>.+?</NSE-ss>)!!gsi;

	my $refsnp = Bio::Cluster::dbSNP->newSNP($xml_rs);
	foreach my $xml_ss (@xml_ss){
		my $subsnp = Bio::Cluster::dbSNP->newSNP($xml_ss);
		$refsnp->add_subsnp($subsnp);
	}

	return $refsnp;
}

1;

