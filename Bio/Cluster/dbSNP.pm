# $Id$
#
# BioPerl module for Bio::Cluster::dbSNP.pm
#
# Copyright Allen Day <allenday@ucla.edu>, Stan Nelson <snelson@ucla.edu>
# Human Genetics, UCLA Medical School, University of California, Los Angeles

=head1 NAME

Bio::Cluster::dbSNP - create Bio::Variation::SNP objects from dbSNP entry

=head1 SYNOPSIS

  #

=head1 DESCRIPTION

  #

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists. Your participation is much appreciated.

  bioperl-l@bioperl.org              - General discussion
  http://bio.perl.org/MailList.html  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via email
or the web:

  bioperl-bugs@bioperl.org
  http://bio.perl.org/bioperl-bugs/

=head1 AUTHOR

Allen Day E<lt>allenday@ucla.eduE<gt>

=head1 APPENDIX


The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a "_".

=cut

# Let the code begin...

package Bio::Cluster::dbSNP;
use vars qw(@ISA $VERSION);
use strict;


use Bio::Root::Root;
use Bio::Variation::SNP;
use XML::Parser;

$VERSION = '1.0';
@ISA = qw(Bio::Root::Root);

my($snp,$up,$dn,%in,$meth,$char);

=head2 newSNP

 Argument: a dbSNP NSE-rs or NSE-ss entry
 Title   : newSNP
 Usage   : used by ClusterIO
 Returns : a new Bio::Variation::SNP object

=cut

sub newSNP {
	my($class,$xml) = @_;

	$snp = Bio::Variation::SNP->new();
	
	my $parser = XML::Parser->new	(
					Handlers => 	{
							Start => \&_start,
							End   => \&_end,
							Char  => \&_char,
							}
					);
	$parser->parse($xml);
}

sub _start {
	my($p,$tag,%attr) = @_;
	$up = $tag;
}

sub _end {
	my($p,$tag,%attr) = @_;
	if($up eq $tag){
		print "<$up>$char</$up>\n";
		undef $char;
	}
}

sub _char {
	my($p,$str) = @_;
	$char .= $str;
}

1;
