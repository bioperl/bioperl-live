
#
# BioPerl module for Bio::DBLinkContainerI
#
# Cared for by Heikki Lehvaslaiho <heikki@ebi.ac.uk>
#
# Copyright Heikki Lehvaslaiho
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::DBLinkContainerI - Abstract interface for any object wanting to use  
                        database cross references
=head1 SYNOPSIS

    # get an objects containing database cross reference

        foreach $obj ( @objs ) {
                if( $obj->isa('Bio::DBLinkContainerI') ) {
                        foreach $dblink ( $obj->each_DBLink() ) {
                                # do stuff
                        }
                }
        }
             

=head1 DESCRIPTION

This interface defines the functions one can expect for any object wanting
to use database cross-references. This class doesn't actually provide
any implemention, it just provides the definitions of what methods one
can call.

The database cross-references are implemented as L<Bio::Annotation::DBLink>
objects.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists.  Your participation is much appreciated.

  vsns-bcd-perl@lists.uni-bielefeld.de          - General discussion
  vsns-bcd-perl-guts@lists.uni-bielefeld.de     - Technically-oriented discussion
  http://bio.perl.org/MailList.html             - About the mailing lists

=head2 Reporting Bugs

report bugs to the Bioperl bug tracking system to help us keep track
 the bugs and their resolution.  Bug reports can be submitted via
 email or the web:

  bioperl-bugs@bio.perl.org
  http://bio.perl.org/bioperl-bugs/

=head1 AUTHOR - Heikki Lehvaslaiho

Email:  heikki@ebi.ac.uk
Address: 

     EMBL Outstation, European Bioinformatics Institute
     Wellcome Trust Genome Campus, Hinxton
     Cambs. CB10 1SD, United Kingdom 

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::DBLinkContainerI;
use vars qw(@ISA);
use strict;

use Carp;


=head2 each_DBLink

 Title   : each_DBLink
 Usage   : foreach $ref ( $self->each_DBlink() )
 Function: gets an array of DBlink of objects
 Example :
 Returns : an array of Bio::Annotation::DBLink objects
 Args    : none


=cut

sub each_DBLink{
   my ($self) = @_;

   $self->_abstractDeath();

}


sub _abstractDeath {
  my $self = shift;
  my $package = ref $self;
  my $caller = (caller)[1];
  
  confess "Abstract method '$caller' defined in interface 
          Bio::DBLinkContainerI not implemented by package $package";
}


1;

