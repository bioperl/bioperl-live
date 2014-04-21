#
# BioPerl module for Bio::Ontology::RelationshipFactory
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Hilmar Lapp <hlapp at gmx.net>
#
# Copyright Hilmar Lapp
#
# You may distribute this module under the same terms as perl itself

#
# (c) Hilmar Lapp, hlapp at gmx.net, 2002.
# (c) GNF, Genomics Institute of the Novartis Research Foundation, 2002.
#
# You may distribute this module under the same terms as perl itself.
# Refer to the Perl Artistic License (see the license accompanying this
# software package, or see http://www.perl.com/language/misc/Artistic.html)
# for the terms under which you may use, modify, and redistribute this module.
# 
# THIS PACKAGE IS PROVIDED "AS IS" AND WITHOUT ANY EXPRESS OR IMPLIED
# WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED WARRANTIES OF
# MERCHANTIBILITY AND FITNESS FOR A PARTICULAR PURPOSE.
#

# POD documentation - main docs before the code

=head1 NAME

Bio::Ontology::RelationshipFactory - Instantiates a new
Bio::Ontology::RelationshipI (or derived class) through a factory

=head1 SYNOPSIS

    use Bio::Ontology::RelationshipFactory;

    # the default type is Bio::Ontology::Relationship
    my $factory = Bio::Ontology::RelationshipFactory->new(
                                 -type => 'Bio::Ontology::GOterm');
    my $clu = $factory->create_object(-name => 'peroxisome',
                                      -ontology => 'Gene Ontology',
                                      -identifier => 'GO:0005777');


=head1 DESCRIPTION

This object will build L<Bio::Ontology::RelationshipI> objects generically.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

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
of the bugs and their resolution. Bug reports can be submitted via
the web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Hilmar Lapp

Email hlapp at gmx.net


=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Ontology::RelationshipFactory;
use strict;

use Bio::Root::Root;

use base qw(Bio::Factory::ObjectFactory);

=head2 new

 Title   : new
 Usage   : my $obj = Bio::Ontology::RelationshipFactory->new();
 Function: Builds a new Bio::Ontology::RelationshipFactory object 
 Returns : Bio::Ontology::RelationshipFactory
 Args    : -type => string, name of a Bio::Ontology::RelationshipI
                    derived class.
                    The default is Bio::Ontology::Relationship.

See L<Bio::Ontology::Relationship>, L<Bio::Ontology::RelationshipI>.

=cut

sub new {
    my($class,@args) = @_;

    my $self = $class->SUPER::new(@args);
  
    # make sure this matches our requirements
    $self->interface("Bio::Ontology::RelationshipI");
    $self->type($self->type() || "Bio::Ontology::Relationship");

    return $self;
}

1;
