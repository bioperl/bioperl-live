#
# BioPerl module for Bio::Biblio::MedlineArticle
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Martin Senger <senger@ebi.ac.uk>
# For copyright and disclaimer see below.

# POD documentation - main docs before the code

=head1 NAME

Bio::Biblio::MedlineArticle - Representation of a MEDLINE article

=head1 SYNOPSIS

  $obj = Bio::Biblio::MedlineArticle->new(-mesh_headings =>
                                            #array ref of hashes
                                         );

  # how are Mesh terms stored:
  use Data::Dumper;
  print Data::Dumper->Dump ( [$obj->mesh_headings], ['MeshHeadings']);

  #It produces (something like) this:
  #'MeshHeadings' => [
  #       { 'descriptorName' => 'Adult' },
  #       { 'descriptorName' => 'Cardiovascular Diseases',
  #         'subHeadings'    => [ { 'subHeading' => 'etiology' },
  #                               { 'majorTopic' => 'Y',
  #                                 'subHeading' => 'mortality' } ] },
  #       { 'descriptorName' => 'Child Development',
  #         'subHeadings'    => [ { 'majorTopic' => 'Y',
  #                                 'subHeading' => 'physiology' } ] },
  #       { 'descriptorName' => 'Human' },
  #      ]

=head1 DESCRIPTION

A storage object for a MEDLINE article.
See its place in the class hierarchy in
http://www.ebi.ac.uk/~senger/openbqs/images/bibobjects_perl.gif

=head2 Attributes

The following attributes are specific to this class
(however, you can also set and get all attributes defined in the parent classes):

  affiliation
  chemicals                      type: array ref of hashes
  citation_owner
  comment_ins                    type: array ref of hashes
  comment_ons                    type: array ref of hashes
  date_of_electronic_publication
  erratum_fors                   type: array ref of hashes
  erratum_in                     type: array ref of hashes
  gene_symbols
  general_notes                  type: array ref of hashes
  grant_list_complete
  grants                         type: array ref of hashes
  medline_date
  medline_id
  medline_page
  mesh_headings                  type: array ref of hashes
  number_of_references
  original_report_ins            type: array ref of hashes
  other_abstracts                type: array ref of hashes
  other_ids                      type: array ref of hashes
  other_languages
  pmid
  republished_froms              type: array ref of hashes
  republished_ins                type: array ref of hashes
  retraction_ins                 type: array ref of hashes
  retraction_ofs                 type: array ref of hashes
  season
  status
  summary_for_patients_ins       type: array ref of hashes
  update_ins                     type: array ref of hashes
  update_ofs                     type: array ref of hashes
  vernacular_title

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

=head1 AUTHORS

Heikki Lehvaslaiho (heikki-at-bioperl-dot-org),
Martin Senger (senger@ebi.ac.uk)

=head1 COPYRIGHT

Copyright (c) 2002 European Bioinformatics Institute. All Rights Reserved.

This module is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.

=head1 DISCLAIMER

This software is provided "as is" without warranty of any kind.

=cut


# Let the code begin...


package Bio::Biblio::MedlineArticle;
use strict;


use base qw(Bio::Biblio::Article);

#
# a closure with a list of allowed attribute names (these names
# correspond with the allowed 'get' and 'set' methods); each name also
# keep what type the attribute should be (use 'undef' if it is a
# simple scalar)
#
{
    my %_allowed =
	(
	 _affiliation => undef,
	 _chemicals => 'ARRAY',
	 _citation_owner => undef,
	 _comment_ins => 'ARRAY',
	 _comment_ons => 'ARRAY',
	 _date_of_electronic_publication => undef,
	 _erratum_fors => 'ARRAY',
	 _erratum_ins => 'ARRAY',
	 _gene_symbols => undef,
	 _general_notes => 'ARRAY',
	 _grant_list_complete => undef,
	 _grants => 'ARRAY',
	 _medline_date => undef,
	 _medline_id => undef,
	 _medline_page => undef,
	 _mesh_headings => 'ARRAY',
	 _number_of_references => undef,
	 _original_report_ins => 'ARRAY',
	 _other_abstracts => 'ARRAY',
	 _other_ids => 'ARRAY',
	 _other_languages => undef,
	 _pmid => undef,
	 _republished_froms => 'ARRAY',
	 _republished_ins => 'ARRAY',
	 _retraction_ins => 'ARRAY',
	 _retraction_ofs => 'ARRAY',
	 _season => undef,
	 _status => undef,
	 _summary_for_patients_ins => 'ARRAY',
	 _update_ins => 'ARRAY',
	 _update_ofs => 'ARRAY',
	 _vernacular_title => undef,
	 );

    # return 1 if $attr is allowed to be set/get in this class
    sub _accessible {
	my ($self, $attr) = @_;
	exists $_allowed{$attr} or $self->SUPER::_accessible ($attr);
    }

    # return an expected type of given $attr
    sub _attr_type {
	my ($self, $attr) = @_;
	if (exists $_allowed{$attr}) {
	    return $_allowed{$attr};
	} else {
	    return $self->SUPER::_attr_type ($attr);
	}
    }
}


1;
__END__
