#
# BioPerl module for Bio::OntologyIO::goflat
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Christian M. Zmasek <czmasek-at-burnham.org> or <cmzmasek@yahoo.com>
#
# (c) Christian M. Zmasek, czmasek-at-burnham.org, 2002.
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
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::OntologyIO::goflat - a parser for the Gene Ontology flat-file format

=head1 SYNOPSIS

  use Bio::OntologyIO;

  # do not use directly -- use via Bio::OntologyIO
  my $parser = Bio::OntologyIO->new
	( -format       => "go",
     -defs_file    => "/home/czmasek/GO/GO.defs",
	  -files        => ["/home/czmasek/GO/component.ontology",
	                    "/home/czmasek/GO/function.ontology",
	                    "/home/czmasek/GO/process.ontology"] );

  my $go_ontology = $parser->next_ontology();

  my $IS_A    = Bio::Ontology::RelationshipType->get_instance( "IS_A" );
  my $PART_OF = Bio::Ontology::RelationshipType->get_instance( "PART_OF" );

=head1 DESCRIPTION

Needs Graph.pm from CPAN.

This is essentially a very thin derivation of the dagflat parser.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to the 
Bioperl mailing lists  Your participation is much appreciated.

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

=head1 AUTHOR

Christian M. Zmasek

Email: czmasek-at-burnham.org  or  cmzmasek@yahoo.com

WWW:   http://monochrome-effect.net/

Address: 

  Genomics Institute of the Novartis Research Foundation
  10675 John Jay Hopkins Drive
  San Diego, CA 92121

=head2 CONTRIBUTOR

 Hilmar Lapp, hlapp at gmx.net

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package  Bio::OntologyIO::goflat;

use strict;

use Bio::Ontology::TermFactory;

use constant TRUE         => 1;
use constant FALSE        => 0;


use base qw(Bio::OntologyIO::dagflat);


=head2 new

 Title   : new
 Usage   : $parser = Bio::OntologyIO->new(
                             -format => "go",
                             -defs_file => "/path/to/GO.defs",
                             -files => ["/path/to/component.ontology",
                                        "/path/to/function.ontology",
                                        "/path/to/process.ontology"] );
 Function: Creates a new goflat parser.
 Returns : A new goflat parser object, implementing Bio::OntologyIO.
 Args    : -defs_file  => the name of the file holding the term
                          definitions
           -files      => a single ontology flat file holding the
                          term relationships, or an array ref holding
                          the file names (for GO, there will usually be
                          3 files: component.ontology, function.ontology,
                          process.ontology)
           -file       => if there is only a single flat file, it may
                          also be specified via the -file parameter
           -ontology_name => the name of the ontology; if not specified the
                          parser will auto-discover it by using the term
                          that starts with a $, and converting underscores
                          to spaces
           -engine     => the Bio::Ontology::OntologyEngineI object
                          to be reused (will be created otherwise); note
                          that every Bio::Ontology::OntologyI will
                          qualify as well since that one inherits from the
                          former.

See L<Bio::OntologyIO>.

=cut

# in reality, we let OntologyIO::new do the instantiation, and override
# _initialize for all initialization work
sub _initialize {
    my ($self, @args) = @_;
    
    $self->SUPER::_initialize( @args );

    # default term object factory
    $self->term_factory(Bio::Ontology::TermFactory->new(
					  -type => "Bio::Ontology::GOterm"))
	unless $self->term_factory();

} # _initialize

  
1;
