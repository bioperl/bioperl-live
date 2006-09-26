# $Id$
#------------------------------------------------------------------
#
# BioPerl module Bio::Restriction::EnzymeI
#
# Cared for by Heikki Lehvaslaiho, heikki-at-bioperl-dot-org
#
# You may distribute this module under the same terms as perl itself
#------------------------------------------------------------------

## POD Documentation:

=head1 NAME

Bio::Restriction::EnzymeI - Interface class for restriction endonuclease

=head1 SYNOPSIS

  # do not run this class directly

=head1 DESCRIPTION

This module defines methods for a single restriction endonuclease.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists. Your participation is much appreciated.

  bioperl-l@bioperl.org                  - General discussion
  http://bioperl.org/wiki/Mailing_lists  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution. Bug reports can be submitted via the
web:

  http://bugzilla.open-bio.org/

=head1 AUTHOR

Heikki Lehvaslaiho, heikki-at-bioperl-dot-org

=head1 CONTRIBUTORS

Rob Edwards, redwards@utmem.edu

=head1 SEE ALSO

L<Bio::Restriction::Enzyme>

=head1 APPENDIX

Methods beginning with a leading underscore are considered private and
are intended for internal use by this module. They are not considered
part of the public interface and are described here for documentation
purposes only.

=cut

package Bio::Restriction::EnzymeI;
use strict;



use base qw(Bio::Root::RootI);

sub name {  shift->throw_not_implemented; }
sub site {  shift->throw_not_implemented; }
sub cuts_after {  shift->throw_not_implemented; }
sub cut {  shift->throw_not_implemented; }
sub complementary_cut {  shift->throw_not_implemented; }
sub type {  shift->throw_not_implemented; }
sub seq {  shift->throw_not_implemented; }
sub string {  shift->throw_not_implemented; }
sub revcom {  shift->throw_not_implemented; }
sub recognition_length {  shift->throw_not_implemented; }
sub non_ambiguous_length {  shift->throw_not_implemented; }
sub cutter {  shift->throw_not_implemented; }
sub palindromic {  shift->throw_not_implemented; }
sub overhang {  shift->throw_not_implemented; }
sub overhang_seq {  shift->throw_not_implemented; }
sub is_ambiguous {  shift->throw_not_implemented; }
sub is_prototype {  shift->throw_not_implemented; }
sub prototype_name {  shift->throw_not_implemented; }
sub isoschizomers {  shift->throw_not_implemented; }
sub purge_isoschizomers {  shift->throw_not_implemented; }
sub methylation_sites {  shift->throw_not_implemented; }
sub purge_methylation_sites {  shift->throw_not_implemented; }
sub microbe {  shift->throw_not_implemented; }
sub source {  shift->throw_not_implemented; }
sub vendors {  shift->throw_not_implemented; }
sub purge_vendors {  shift->throw_not_implemented; }
sub vendor {  shift->throw_not_implemented; }
sub references {  shift->throw_not_implemented; }
sub purge_references {  shift->throw_not_implemented; }
sub xxxxxxx {  shift->throw_not_implemented; }

1;

