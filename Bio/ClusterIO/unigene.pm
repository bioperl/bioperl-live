# $Id$
# BioPerl module for Bio::ClusterIO::unigene
#
# Cared for by Andrew Macgregor <andrew@anatomy.otago.ac.nz>
#
# Copyright Andrew Macgregor, Jo-Ann Stanton, David Green
# Molecular Embryology Group, Anatomy & Structural Biology, University of Otago
# http://anatomy.otago.ac.nz/meg
#
# You may distribute this module under the same terms as perl itself
#
# _history
# April 17, 2002 - Initial implementation by Andrew Macgregor

# POD documentation - main docs before the code

=head1 NAME

Bio::ClusterIO::unigene - UniGene input stream

=head1 SYNOPSIS

Do not use this module directly.  Use it via the Bio::ClusterIO class.

=head1 DESCRIPTION

This object reads from Unigene *.data files downloaded from ftp://ncbi.nlm.nih.gov/repository/UniGene/.
It doesn't download and decompress the file, you have to do that yourself.

This module requires that you have the Parse::RecDescent module.

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

=head1 AUTHORS - Andrew Macgregor

Email: andrew@anatomy.otago.ac.nz


=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...
package Bio::ClusterIO::unigene;
use vars qw(@ISA);
use strict;

use Bio::ClusterIO;
use Parse::RecDescent;

#$::RD_TRACE = 1;
$::RD_WARN = 1;
$::RD_HINT = 1;

@ISA = qw(Bio::ClusterIO);


# this is the guts of it, the grammar to parse the unigene records
my $grammar = <<'EOGRAMMAR';

{
my $unigene_id;
my @sequence;
my @sts;
my @txmap;
my @protsim;
}

record: <rulevar: local $UGobj = $arg[0]>

record			:	id
					title
					gene(?)
					cytoband(?)
					mgi(?)
					locuslink(?)
					express(?)
					gnm_terminus(?)
					chromosome(?)
					sts(s?)
					txmap(s?)
					protsim(s?)
					scount
					sequence(s)
					delimiter
														{
															
															$UGobj->unigene_id($item{id});
															$UGobj->title($item{title});
															if (defined $item{gene}->[0]) { $UGobj->gene($item{gene}->[0]) };
															if (defined $item{cytoband}->[0]) { $UGobj->cytoband($item{cytoband}->[0]) };
															if (defined $item{locuslink}->[0]) { $UGobj->locuslink($item{locuslink}->[0]) };
															if (defined $item{gnm_terminus}->[0]) { $UGobj->gnm_terminus($item{gnm_terminus}->[0]) };
															if (defined $item{chromosome}->[0]) { $UGobj->chromosome($item{chromosome}->[0]) };
															$UGobj->scount($item{scount});
															$UGobj->sts(\@sts);
															$UGobj->txmap(\@txmap);
															$UGobj->protsim(\@protsim);
															$UGobj->sequence(\@sequence);
														}
					| <error>
					
		


id				:	'ID' organism '.' unigene_no		{
															$unigene_id = "$item{organism}.$item{unigene_no}";
															$return = "$item{organism}.$item{unigene_no}";
														}

title			:	'TITLE' /.+/

gene			:	'GENE' /.+/

cytoband		:	'CYTOBAND' /.+/

mgi				:	'MGI' /.+/ 

express			:	'EXPRESS' /.+/						{
															$item[2] =~ s/^;//;
															my @express = split /;/ , $item[2];
															$UGobj->express(\@express);
														}

gnm_terminus	:	'GNM_TERMINUS' /.+/ 

locuslink		:	'LOCUSLINK' /[0-9]+/

chromosome		:	'CHROMOSOME' /[0-9XY|Un]+/

sts				:	'STS' /.+/							{ push @sts, $item[2]; }

txmap			:	'TXMAP'	/.+/						{ push @txmap, $item[2]; }

protsim			:	'PROTSIM' /.+/						{ push @protsim, $item[2]; }

scount			:	'SCOUNT' /[0-9]+/

sequence		:	'SEQUENCE'
					acc(?)
					nid(?)
					pid(?)
					clone(?)
					end(?)
					lid(?)
					mgc(?)					{ 
												my $seq = {};
												$seq->{unigene_id} = $unigene_id;
												$seq->{acc} = $item{acc}->[0] if defined $item{acc}->[0];
												$seq->{nid} = $item{nid}->[0] if defined $item{nid}->[0];
												$seq->{pid} = $item{pid}->[0] if defined $item{pid}->[0];
												$seq->{clone} = $item{clone}->[0] if defined $item{clone}->[0];
												$seq->{end} = $item{end}->[0] if defined $item{end}->[0];
												$seq->{lid} = $item{lid}->[0] if defined $item{lid}->[0];
												$seq->{mgc} = $item{mgc}->[0] if defined $item{mgc}->[0];
												push @sequence, $seq;
											}			


organism		:	/At|Bt|Dr|Hs|Hv|Mm|Os|Rn|Ta|Xl|Zm/

unigene_no		:	/[0-9]+/

acc				:	'ACC=' /\w+/ seq_delimiter(s?)		{ $return =  $item[2] }

nid				:	'NID=' /\w+/ seq_delimiter(s?)		{ $return =  $item[2] }

pid				:	'PID=' /\w+/ seq_delimiter(s?)		{ $return =  $item[2] }

clone			:	'CLONE=' /[^;\n]+/ seq_delimiter(s?)	{ $return =  $item[2] }

end				:	'END=' /5'|3'/ seq_delimiter(s?)		{ $return =  $item[2] }

lid				:	'LID=' /\w+/ seq_delimiter(s?)		{ $return =  $item[2] }

mgc				:	'MGC=' /\w+/ seq_delimiter(s?)		{ $return =  $item[2] }



seq_delimiter 	:	';'

delimiter		:	/\/\/\Z/
EOGRAMMAR


my $parser = new Parse::RecDescent ($grammar);

=head2 next_unigene

 Title   : next_unigene
 Usage   : $unigene = $stream->next_unigene()
 Function: returns the next unigene in the stream
 Returns : Bio::Cluster::UniGene object
 Args    : NONE

=cut

sub next_unigene {
	my( $self) = @_;
	local $/ = "//";
	return unless my $entry = $self->_readline;
	
	my $UGobj = Bio::Cluster::UniGene->new();
	$parser->record($entry,1,$UGobj);
	
	return $UGobj;
}

1;

