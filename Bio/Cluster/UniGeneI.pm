# $Id$
#
# BioPerl module for Bio::Cluster::UniGeneI.pm
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
# April 31, 2002 - Initial implementation by Andrew Macgregor
# POD documentation - main docs before the code

=head1 NAME

Bio::Cluster::UniGeneI - abstract interface of UniGene object

=head1 SYNOPSIS

  #

=head1 DESCRIPTION

  #

=head1 FEEDBACK

  #

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

=head1 AUTHOR - Andrew Macgregor

Email andrew@anatomy.otago.ac.nz


=head1 APPENDIX


The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a "_".

=cut

# Let the code begin...


package Bio::Cluster::UniGeneI;
use vars qw(@ISA $VERSION);
use strict;


use Bio::Root::Root;
use Bio::Seq;
use Bio::PrimarySeq;

$VERSION = '1.0';
@ISA = qw(Bio::Root::Root);


=head2 new

 Title   : new
 Usage   : used by ClusterIO
 Returns : a new Bio::Cluster::Unigene object

=cut

sub new {
	my ($self) = @_;
	$self->throw_not_implemented;
}


=head2 unigene_id

 Title   : unigene_id
 Usage   : unigene_id();
 Function: Returns the unigene_id associated with the object.
 Example : $id = $unigene->unigene_id or $unigene->unigene_id($id)
 Returns : A string
 Args    : None or an id


=cut

sub unigene_id {
	my ($self) = @_;
	$self->throw_not_implemented;
}



=head2 title

 Title   : title
 Usage   : title();
 Function: Returns the title associated with the object.
 Example : $title = $unigene->title or $unigene->title($title)
 Returns : A string
 Args    : None or a title


=cut

sub title {
	my ($self) = @_;
	$self->throw_not_implemented;
}


=head2 gene

 Title   : gene
 Usage   : gene();
 Function: Returns the gene associated with the object.
 Example : $gene = $unigene->gene or $unigene->gene($gene)
 Returns : A string
 Args    : None or a gene


=cut

sub gene {
	my ($self) = @_;
	$self->throw_not_implemented;
}


=head2 cytoband

 Title   : cytoband
 Usage   : cytoband();
 Function: Returns the cytoband associated with the object.
 Example : $cytoband = $unigene->cytoband or $unigene->cytoband($cytoband)
 Returns : A string
 Args    : None or a cytoband


=cut

sub cytoband {
	my ($self) = @_;
	$self->throw_not_implemented;
}


=head2 locuslink

 Title   : locuslink
 Usage   : locuslink();
 Function: Returns the locuslink associated with the object.
 Example : $locuslink = $unigene->locuslink or $unigene->locuslink($locuslink)
 Returns : A string
 Args    : None or a locuslink

=cut

sub locuslink {
	my ($self) = @_;
	$self->throw_not_implemented;
}


=head2 gnm_terminus

 Title   : gnm_terminus
 Usage   : gnm_terminus();
 Function: Returns the gnm_terminus associated with the object.
 Example : $gnm_terminus = $unigene->gnm_terminus or $unigene->gnm_terminus($gnm_terminus)
 Returns : A string
 Args    : None or a gnm_terminus

=cut

sub gnm_terminus {
	my ($self) = @_;
	$self->throw_not_implemented;
}


=head2 scount

 Title   : scount
 Usage   : scount();
 Function: Returns the scount associated with the object.
 Example : $scount = $unigene->scount or $unigene->scount($scount)
 Returns : A string
 Args    : None or a scount

=cut

sub scount {
	my ($self) = @_;
	$self->throw_not_implemented;
}



=head2 express

 Title   : express
 Usage   : express();
 Function: Returns or stores a reference to an array containing tissue expression data
 		   This should really only be used by ClusterIO, not directly
 Returns : An array reference
 Args    : None or an array reference

=cut

sub express {
	my ($self) = @_;
	$self->throw_not_implemented;
}


=head2 next_express

 Title   : next_express
 Usage   : next_express();
 Function: Returns the next tissue from an array referred to using $obj->{'express'}
 Example : 	while ( my $express = $in->next_express() ) {
				print "$express\n";
			}
 Returns : String
 Args    : None

=cut

sub next_express {
	my ($self) = @_;
	$self->throw_not_implemented;
}


=head2 chromosome

 Title   : chromosome
 Usage   : chromosome();
 Function: Returns or stores a reference to an array containing chromosome lines
 		   This should really only be used by ClusterIO, not directly
 Returns : An array reference
 Args    : None or an array reference

=cut

sub chromosome {
	my ($self) = @_;
	$self->throw_not_implemented;
}


=head2 next_chromosome

 Title   : next_chromosome
 Usage   : next_chromosome();
 Function: Returns the next chromosome line from an array referred to using $obj->{'chromosome'}
 Example : 	while ( my $chromosome = $in->next_chromosome() ) {
				print "$chromosome\n";
			}
 Returns : String
 Args    : None

=cut

sub next_chromosome {
	my ($self) = @_;
	$self->throw_not_implemented;
}


=head2 sts

 Title   : sts
 Usage   : sts();
 Function: Returns or stores a reference to an array containing sts lines
 		   This should really only be used by ClusterIO, not directly
 Returns : An array reference
 Args    : None or an array reference

=cut

sub sts {
	my ($self) = @_;
	$self->throw_not_implemented;
}


=head2 next_sts

 Title   : next_sts
 Usage   : next_sts();
 Function: Returns the next sts line from an array referred to using $obj->{'sts'}
 Example : 	while ( my $sts = $in->next_sts() ) {
				print "$sts\n";
			}
 Returns : String
 Args    : None

=cut

sub next_sts {
	my ($self) = @_;
	$self->throw_not_implemented;
}


=head2 txmap

 Title   : txmap
 Usage   : txmap();
 Function: Returns or stores a reference to an array containing txmap lines
 		   This should really only be used by ClusterIO, not directly
 Returns : An array reference
 Args    : None or an array reference

=cut

sub txmap {
	my ($self) = @_;
	$self->throw_not_implemented;
}


=head2 next_txmap

 Title   : next_txmap
 Usage   : next_txmap();
 Function: Returns the next txmap line from an array referred to using $obj->{'txmap'}
 Example : 	while ( my $tsmap = $in->next_txmap() ) {
				print "$txmap\n";
			}
 Returns : String
 Args    : None

=cut

sub next_txmap {
	my ($obj) = @_;
	shift @{$obj->{'txmap'}};
}


=head2 protsim

 Title   : protsim
 Usage   : protsim();
 Function: Returns or stores a reference to an array containing protsim lines
 		   This should really only be used by ClusterIO, not directly
 Returns : An array reference
 Args    : None or an array reference

=cut

sub protsim {
	my ($self) = @_;
	$self->throw_not_implemented;
}


=head2 next_protsim

 Title   : next_protsim
 Usage   : next_protsim();
 Function: Returns the next protsim line from an array referred to using $obj->{'protsim'}
 Example : 	while ( my $protsim = $in->next_protsim() ) {
				print "$protsim\n";
			}
 Returns : String
 Args    : None

=cut

sub next_protsim {
	my ($self) = @_;
	$self->throw_not_implemented;
}



=head2 sequence

 Title   : sequence
 Usage   : sequence();
 Function: Returns or stores a reference to an array containing sequence data
 		   This should really only be used by ClusterIO, not directly
 Returns : An array reference
 Args    : None or an array reference

=cut

sub sequence {
	my ($self) = @_;
	$self->throw_not_implemented;
}


=head2 next_seq

 Title   : next_seq
 Usage   : next_seq();
 Function: Returns the next seq as a Seq object, at present with just the accession_number
 Example : 		while ( my $sequence = $in->next_seq() ) {
					print $sequence->accession_number() . "\n";
				}	

 Returns : String
 Args    : None

=cut

sub next_seq {
	my ($self) = @_;
	$self->throw_not_implemented;
}


1;
