#
# BioPerl module for Bio::Cluster::UniGeneI.pm
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Andrew Macgregor <andrew at cbbc.murdoch.edu.au>
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

This is the general interface for a UniGene cluster representation in
Bioperl. You cannot use this module directly, use an implementation
instead.

You can create UniGene cluster objects yourself by instantiating
L<Bio::Cluster::UniGene>. If you read UniGene clusters from a
ClusterIO parser, you will get objects implementing this interface,
most likely instances of said UniGene class.

L<Bio::Cluster::UniGeneI> inherits from L<Bio::ClusterI>, so you can
use it wherever a cluster object is expected.

=head1 FEEDBACK

  #

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists. Your participation is much appreciated.

  bioperl-l@bioperl.org                  - General discussion
  http://bioperl.org/wiki/Mailing_lists  - About the mailing lists

=head2 Support 

Please direct usage questions or support issues to the mailing list:

I<bioperl-l@bioperl.org>

rather than to the module maintainer directly. Many experienced and 
reponsive experts will be able look at the problem and quickly 
address it. Please include a thorough description of the problem 
with code and data examples if at all possible.

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via the
web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Andrew Macgregor

Email andrew at cbbc.murdoch.edu.au


=head1 APPENDIX


The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a "_".

=cut

# Let the code begin...


package Bio::Cluster::UniGeneI;
use strict;


use base qw(Bio::ClusterI);


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


=head2 mgi

 Title   : mgi
 Usage   : mgi();
 Function: Returns the mgi associated with the object.
 Example : $mgi = $unigene->mgi or $unigene->mgi($mgi)
 Returns : A string
 Args    : None or a mgi


=cut

sub mgi {
    my ($self) = @_;
    $self->throw_not_implemented;
}


=head2 locuslink

 Title   : locuslink
 Usage   : locuslink();
 Function: Returns or stores a reference to an array containing locuslink data.
           This should really only be used by ClusterIO, not directly
 Returns : An array reference
 Args    : None or an array reference

=cut

sub locuslink {
    my ($self) = @_;
    $self->throw_not_implemented;
}


=head2 homol

 Title   : homol
 Usage   : homol();
 Function: Returns the homol entry associated with the object.
 Example : $homol = $unigene->homol or $unigene->homol($homol)
 Returns : A string
 Args    : None or a homol entry

=cut

sub homol {
    my ($self) = @_;
    $self->throw_not_implemented;
}


=head2 restr_expr

 Title   : restr_expr
 Usage   : restr_expr();
 Function: Returns the restr_expr entry associated with the object.
 Example : $restr_expr = $unigene->restr_expr or $unigene->restr_expr($restr_expr)
 Returns : A string
 Args    : None or a restr_expr entry

=cut

sub restr_expr {
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
 Function: Returns or stores a reference to an array containing tissue expression data.
           This should really only be used by ClusterIO, not directly
 Returns : An array reference
 Args    : None or an array reference

=cut

sub express {
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


=head2 txmap

 Title   : txmap
 Usage   : txmap();
 Function: Returns or stores a reference to an array containing txmap lines
 Returns : An array reference
 Args    : None or an array reference

=cut

sub txmap {
    my ($self) = @_;
    $self->throw_not_implemented;
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

=head2 species

 Title   : species
 Usage   : $obj->species($newval)
 Function: Get the species object for this Unigene cluster.
 Example : 
 Returns : value of species (a L<Bio::Species> object)
 Args    : 


=cut

sub species{
    shift->throw_not_implemented();
}

=head1 Methods inherited from L<Bio::ClusterI>

=cut

=head2 display_id

 Title   : display_id
 Usage   : 
 Function: Get/set the display name or identifier for the cluster
 Returns : a string
 Args    : optional, on set the display ID ( a string)

=cut

=head2 description

 Title   : description
 Usage   : Bio::ClusterI->description("POLYUBIQUITIN")
 Function: get/set for the consensus description of the cluster
 Returns : the description string 
 Args    : Optional the description string 

=cut

=head2 size

 Title   : size
 Usage   : Bio::ClusterI->size();
 Function: get/set for the size of the family, 
           calculated from the number of members
 Returns : the size of the family 
 Args    : 

=cut

=head2 cluster_score

 Title   : cluster_score
 Usage   : $cluster ->cluster_score(100);
 Function: get/set for cluster_score which
           represent the score in which the clustering
           algorithm assigns to this cluster.
 Returns : a number

=cut

=head2 get_members

 Title   : get_members
 Usage   : Bio::ClusterI->get_members(($seq1, $seq2));
 Function: retrieve the members of the family by some criteria, for
           example :
           $cluster->get_members(-species => 'homo sapiens'); 

           Will return all members if no criteria are provided.

 Returns : the array of members
 Args    : 

=cut

1;
