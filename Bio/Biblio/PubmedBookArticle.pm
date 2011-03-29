#
# BioPerl module for Bio::Biblio::PubmedBookArticle
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Martin Senger <senger@ebi.ac.uk>
# For copyright and disclaimer see below.

# POD documentation - main docs before the code

=head1 NAME

Bio::Biblio::PubmedBookArticle - Representation of a PUBMED book article

=head1 SYNOPSIS

    $obj = Bio::Biblio::PubmedBookArticle->new
                  (-title => 'Still getting started'.
                   -book => Bio::Biblio::MedlineBook->new());
    # note that there is no specialised class PubmedBook

  #--- OR ---

    $obj = Bio::Biblio::PubmedBookArticle->new();
    $obj->title ('Still getting started');

=head1 DESCRIPTION

A storage object for a PUBMED book article.
See its place in the class hierarchy in
http://www.ebi.ac.uk/~senger/openbqs/images/bibobjects_perl.gif

=head2 Attributes

There are no specific attributes in this class
(however, you can set and get all attributes defined in the parent classes).

=head1 SEE ALSO

=over 4

=item *

OpenBQS home page: http://www.ebi.ac.uk/~senger/openbqs/

=item *

Comments to the Perl client: http://www.ebi.ac.uk/~senger/openbqs/Client_perl.html

=back

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
of the bugs and their resolution. Bug reports can be submitted via the
web:

  https://redmine.open-bio.org/projects/bioperl/

=head1 AUTHOR

Martin Senger (senger@ebi.ac.uk)

=head1 COPYRIGHT

Copyright (c) 2002 European Bioinformatics Institute. All Rights Reserved.

This module is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.

=head1 DISCLAIMER

This software is provided "as is" without warranty of any kind.

=cut


# Let the code begin...


package Bio::Biblio::PubmedBookArticle;
use strict;
use vars qw(@ISA);

use base qw(Bio::Biblio::PubmedArticle Bio::Biblio::MedlineBookArticle);

#
# a closure with a list of allowed attribute names (these names
# correspond with the allowed 'get' and 'set' methods); each name also
# keep what type the attribute should be (use 'undef' if it is a
# simple scalar)
#
{
    my %_allowed =
	(
	 );

    # return 1 if $attr is allowed to be set/get in this class
    sub _accessible {
	my ($self, $attr) = @_;
	return 1 if exists $_allowed{$attr};
        foreach my $parent (@ISA) {
	    return 1 if $parent->_accessible ($attr);
	}
    }

    # return an expected type of given $attr
    sub _attr_type {
	my ($self, $attr) = @_;
	if (exists $_allowed{$attr}) {
	    return $_allowed{$attr};
	} else {
	    foreach my $parent (@ISA) {
		if ($parent->_accessible ($attr)) {
		    return $parent->_attr_type ($attr);
		}
	    }
	}
	return 'unknown';
    }
}


1;
__END__
