# $Id$
#
# bioperl module for Bio::LiveSeq::Repeat_Unit
#
# Cared for by Joseph Insana <insana@ebi.ac.uk> <jinsana@gmx.net>
#
# Copyright Joseph Insana
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=head1 NAME

  Bio::LiveSeq::Repeat_Unit - Repeat_Unit class for LiveSeq

=head1 SYNOPSIS


=head1 DESCRIPTION

Class for REPEAT_UNIT objects. They consist of a beginlabel, an endlabel (both
referring to a LiveSeq DNA object) and a strand.
The strand could be 1 (forward strand, default), -1 (reverse strand).

=head1 AUTHOR - Joseph A.L. Insana

Email:  Insana@ebi.ac.uk, jinsana@gmx.net

Address: 

     EMBL Outstation, European Bioinformatics Institute
     Wellcome Trust Genome Campus, Hinxton
     Cambs. CB10 1SD, United Kingdom 

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::LiveSeq::Repeat_Unit;
$VERSION=1.0;

# Version history:
# Tue Apr  4 18:11:31 BST 2000 v 1.0 created

use strict;
use vars qw($VERSION @ISA);
use Bio::LiveSeq::Repeat_Region 1.0; # uses Repeat_Region, inherits from it
@ISA=qw(Bio::LiveSeq::Repeat_Region);

=head1 new

  Title   : new
  Usage   : $intron1=Bio::LiveSeq::Repeat_Unit->new(-seq => $objref,
					      -start => $startlabel,
					      -end => $endlabel, -strand => 1);

  Function: generates a new Bio::LiveSeq::Repeat_Unit
  Returns : reference to a new object of class Repeat_Unit
  Errorcode -1
  Args    : two labels and an integer

=cut

1;
