# $Id$
#
# BioPerl module for Bio::DB::Align
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Jun Yin <jun dot yin at ucd dot ie>
#
# Copyright Jun Yin
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code
# 
# Interface from Bio::DB::Align packages

=head1 NAME

Bio::DB::Align - interface fro webagent packages retrieves alignment sequences 
from online databases, e.g. Pfam

=head1 SYNOPSIS

  # ...To be added!

=head1 DESCRIPTION

	# ...To be added!

=head1 TODO

=over 3

=item * Finish documentation

HOWTOs (both standard and Cookbook).

=item * Cookbook tests

Set up dev-only tests for Cookbook examples to make sure they are consistently
updated.

=back

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the 
evolution of this and other Bioperl modules. Send
your comments and suggestions preferably to one
of the Bioperl mailing lists. Your participation
is much appreciated.

  bioperl-l@lists.open-bio.org               - General discussion
  http://www.bioperl.org/wiki/Mailing_lists  - About the mailing lists

=head2 Support 

Please direct usage questions or support issues to the mailing list:

I<bioperl-l@bioperl.org>

rather than to the module maintainer directly. Many experienced and 
reponsive experts will be able look at the problem and quickly 
address it. Please include a thorough description of the problem 
with code and data examples if at all possible.

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to
help us keep track the bugs and their resolution.
Bug reports can be submitted via the web.

  http://bugzilla.open-bio.org/

=head1 AUTHOR 

Email jun dot yin at ucd dot ie

=head1 APPENDIX

The rest of the documentation details each of the
object methods. Internal methods are usually
preceded with a _

=cut

# Let the code begin...

package Bio::DB::Align;
use strict;
use warnings;
use LWP::UserAgent;
use HTTP::Request;
use Bio::AlignIO;
use Bio::Root::IO;
use vars qw(%FORMATS %ALNTYPE $HOSTBASE);

use base qw(Bio::Root::Root Bio::Root::IO Bio::DB::GenericWebAgent);

BEGIN {
	$HOSTBASE = 'http://pfam.sanger.ac.uk';
	%FORMATS=qw(fasta 1 stockholm 1 selex 1 msf 1); #supported formats in pfam
	%ALNTYPE=qw(seed 1 full 1 ncbi 1 metagenomics 1); #supported alignment types
}

sub new {
	my($class,@args) = @_;
	my $self = $class->SUPER::new(@args);
   my $ua = new LWP::UserAgent(env_proxy => 1);
   #$ua->agent(ref($self) ."/$MODVERSION");
   $self->ua($ua);  
   $self->{'_authentication'} = [];	
	return $self;	
}


=head2 get_Aln_by_id

	Title   : get_Aln_by_id
	Usage   : $aln = $db->get_Aln_by_id('Piwi')
	Function: This method uses Pfam id conversion service id2acc to convert id 
	          to accession. Then, it gets a Bio::SimpleAlign object 
	          using get_Aln_by_acc
	Returns : a Bio::SimpleAlign object
	Args    : -id  the id as a string
	Note    : 
=cut

sub get_Aln_by_id {
    my ($self,$aln) = @_;
    $self->throw("Sorry, you cannot read from a generic Bio::DB::Align object.");
}


=head2 get_Aln_by_acc

  Title   : get_Aln_by_acc
  Usage   : $seq = $db->get_Aln_by_acc($acc);
  Function: Gets a Bio::SimpleAlign object by accession numbers
  Returns : a Bio::SimpleAlign object
  Args    : -accession  the accession number as a string
  Note    : 
=cut

sub get_Aln_by_acc {
    my ($self,$aln) = @_;
    $self->throw("Sorry, you cannot read from a generic Bio::DB::Align object.");
}

1;


                  