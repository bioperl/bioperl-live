#-----------------------------------------------------------------------------
# PACKAGE : Bio::Root::Xref.pm
# AUTHOR  : Steve A. Chervitz (sac@genome.stanford.edu)
# CREATED : 8 May 1997
# REVISION: $Id$
# STATUS  : Pre-Alpha 
#
# WARNING: This is considered an experimental module.
#
# Copyright (c) 1997-8 Steve A. Chervitz. All Rights Reserved.
#           This module is free software; you can redistribute it and/or 
#           modify it under the same terms as Perl itself.
#-----------------------------------------------------------------------------

package Bio::Root::Xref;

use Bio::Root::Global;
use Bio::Root::Object ();
use Bio::Root::Vector ();

@Bio::Root::Xref::ISA = qw( Bio::Root::Vector Bio::Root::Object );

use vars qw($ID $VERSION);
$ID = 'Bio::Root::Xref';
$VERSION = 0.01;

## POD Documentation:

=head1 NAME

Bio::Root::Xref - A generic cross-reference object.

B<WARNING: This module is still in the experimental phase and has not been tested.>

=head1 SYNOPSIS

=head2 Object Creation
  
 use Bio::Root::Object;

 $myObj->xref($object_ref);
 
=head2 Object Manipulation
 
 Accessors
 ---------------------------------------------------------------------
 obj()         - Get the cross-referenced object.
 desc()        - Description of the nature of the cross-reference.
 set_desc()    - Set description.
 type()        - Symmetric or assymetric.
 
 Methods
 ---------------------------------------------------------------------
 clear()       - remove all cross-references within the Xref object (not implemented).

=head1 DESCRIPTION

An instance of B<Bio::Root::Xref.pm> manages sets of objects not
necessarily related by inheritance or composition, but by an arbitrary
criterion defined by the client. Currently, Bio::Root::Xref inherits
from both B<Bio::Root::Object.pm> and B<Bio::Root::Vector.pm>. An Xref
object is an example of a heterogeneous Vector object since different
objects in the vector need not all derive from the same base class.

The two objects involved in the cross-reference typically involve a
symmetrical relationship in which each will have a Xref object relating it
to the other object. This relationship is not necessarily transitive,
however: if A is an xref of B and B is an xref of C, A is not
necessarily an xref of C. Assymetric Xrefs are also possible.

The establishment of cross-references is managed by B<Bio::Root::Object.pm>.
See the xref() method in that module.

B<The API for this module is not complete since the module is under development. Caveat emptor.>


=head1 SEE ALSO

  Bio::Root::Object.pm       - Core object
  Bio::Root::Global.pm       - Manages global variables/constants

  http://bio.perl.org/Projects/modules.html  - Online module documentation
  http://bio.perl.org/                       - Bioperl Project Homepage 
 
=head1 FEEDBACK

=head2 Mailing Lists 

User feedback is an integral part of the evolution of this and other Bioperl modules.
Send your comments and suggestions preferably to one of the Bioperl mailing lists.
Your participation is much appreciated.

    vsns-bcd-perl@lists.uni-bielefeld.de          - General discussion
    vsns-bcd-perl-guts@lists.uni-bielefeld.de     - Technically-oriented discussion
    http://bio.perl.org/MailList.html             - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track the bugs and 
their resolution. Bug reports can be submitted via email or the web:

    bioperl-bugs@bio.perl.org                   
    http://bio.perl.org/bioperl-bugs/           

=head1 AUTHOR

Steve A. Chervitz, sac@genome.stanford.edu

See the L<FEEDBACK> section for where to send bug reports and comments.

=head1 VERSION

Bio::Root::Xref.pm, 0.01 pre-alpha

=head1 COPYRIGHT

Copyright (c) 1997-8 Steve A. Chervitz. All Rights Reserved.
This module is free software; you can redistribute it and/or 
modify it under the same terms as Perl itself.


=head1 TODO

Update documentation to work with pod2html from Perl 5.004.

=cut

#
##
###
#### END of main POD documentation.
###
##
#



#####################################################################################
##                                 CONSTRUCTOR                                     ##
#####################################################################################

sub _initialize {
    my( $self, %param ) = @_;

    $self->SUPER::_initialize(%param);
    
    $self->{'_obj'} = ($param{-OBJ} || undef);

    ## By default, all Xrefs are symmetric.
    ## Create symmetric cross-reference in obj.
    if(!$param{-ASYM}) {
	$self->{'_obj'}->xref(-OBJ=>$param{-PARENT});
	$self->{'_type'} = 'sym';
    } else {
	$self->{'_type'} = 'asym';
    }	
}


#####################################################################################
##                                  ACCESSORS                                      ##
#####################################################################################

sub obj {my ($self) = shift; return $self->{'_obj'}; }
sub desc {my ($self) = shift; return $self->{'_desc'}; }
sub type {my ($self) = shift; return $self->{'_type'}; }
    
sub set_desc {my ($self,$desc) = @_; 
	     $self->{'_desc'} = $desc;
	 }

sub clear {
## Not implemented. Need to do this carefully.
## Not sure if this method is needed.    
    my ($self) = @_;
}

1;
__END__

#####################################################################################
#                                  END OF CLASS                                     #
#####################################################################################

=head1 DATA MEMBERS 

 _obj   : The object being cross-referenced to the parent.
 _type  : Symmetric or asymmetric
 _desc  : Description associated with the cross-reference
 
 INHERITED DATA MEMBERS (from Bio::Root::Object)

 _parent : The object receiving the cross-reference.
 _name   : Descriptive nature of the cross-reference.

=cut


