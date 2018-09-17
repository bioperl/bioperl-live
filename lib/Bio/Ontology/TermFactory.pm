#
# BioPerl module for Bio::Ontology::TermFactory
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

Bio::Ontology::TermFactory - Instantiates a new 
Bio::Ontology::TermI (or derived class) through a factory

=head1 SYNOPSIS

    use Bio::Ontology::TermFactory;

    # the default type is Bio::Ontology::Term
    my $factory = Bio::Ontology::TermFactory->new(
                        -type => 'Bio::Ontology::GOterm');
    my $term = $factory->create_object(-name => 'peroxisome',
                                       -ontology => 'Gene Ontology',
                                       -identifier => 'GO:0005777');


=head1 DESCRIPTION

This object will build L<Bio::Ontology::TermI> objects generically.

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


package Bio::Ontology::TermFactory;
use strict;

use Bio::Root::Root;

use base qw(Bio::Factory::ObjectFactory);

=head2 new

 Title   : new
 Usage   : my $obj = Bio::Ontology::TermFactory->new();
 Function: Builds a new Bio::Ontology::TermFactory object 
 Returns : Bio::Ontology::TermFactory
 Args    : -type => string, name of a Bio::Ontology::TermI derived class.
                    The default is Bio::Ontology::Term.

See L<Bio::Ontology::TermI>, L<Bio::Ontology::Term>.

=cut

sub new {
    my($class,@args) = @_;

    my $self = $class->SUPER::new(@args);
  
    # make sure this matches our requirements
    $self->interface("Bio::Ontology::TermI");
    $self->type($self->type() || "Bio::Ontology::Term");

    return $self;
}


=head2 create_object

 Title   : create_object
 Usage   : my $term = $factory->create_object(<named parameters>);
 Function: Instantiates new Bio::Ontology::TermI (or one of its child classes)

           This object allows us to genericize the instantiation of
           Term objects.

 Returns : Bio::Ontology::TermI compliant object
           The return type is configurable using new(-type =>"...").
 Args    : initialization parameters specific to the type of term
           object we want.  Typically 
           -name        => $name
           -identifier  => identifier for the term
           -ontology    => ontology for the term

See L<Bio::Ontology::TermI>.

=cut

1;
