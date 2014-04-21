#------------------------------------------------------------------
#
# BioPerl module Bio::Restriction::Enzyme::MultiSite
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Heikki Lehvaslaiho, heikki-at-bioperl-dot-org
#
# You may distribute this module under the same terms as perl itself
#------------------------------------------------------------------

## POD Documentation:

=head1 NAME

Bio::Restriction::Enzyme::MultiSite - A single restriction endonuclease

=head1 SYNOPSIS

  # set up a single restriction enzyme. This contains lots of
  # information about the enzyme that is generally parsed from a
  # rebase file and can then be read back

  use Bio::Restriction::Enzyme;


=head1 DESCRIPTION

This module is used for restriction enzymes that recogonize more than
one site. There are some enzymes that recognize sites that cannot be
represented by the ambiguous genetic code. For example, M.PhiBssHII
recognizes the sites: ACGCGT,CCGCGG,RGCGCY,RCCGGY, and GCGCGC

Each site gets its own object that Bio::Restriction::Enzyme will
refer to. Each also correlates with the other sites using the 
method L<others|others> which stores references to other objects 
with alternative sites.

In this schema each object within an EnzymeCollection can be checked
for matching a sequence.

=head1 FEEDBACK

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
the bugs and their resolution. Bug reports can be submitted via the
web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR

Heikki Lehvaslaiho, heikki-at-bioperl-dot-org

=head1 CONTRIBUTORS

Rob Edwards, redwards@utmem.edu

=head1 COPYRIGHT

Copyright (c) 2003 Rob Edwards.

Some of this work is Copyright (c) 1997-2002 Steve A. Chervitz. All
Rights Reserved.  This module is free software; you can redistribute
it and/or modify it under the same terms as Perl itself.

=head1 SEE ALSO

L<Bio::Restriction::Enzyme>, L<Bio::Restriction::Analysis>, 
L<Bio::Restriction::EnzymeCollection>

=head1 APPENDIX

Methods beginning with a leading underscore are considered private and
are intended for internal use by this module. They are not considered
part of the public interface and are described here for documentation
purposes only.

=cut

package Bio::Restriction::Enzyme::MultiSite;
use strict;

use Data::Dumper;

use vars qw ();
use base qw(Bio::Restriction::Enzyme);

=head2 new

 Title     : new
 Function
 Function  : Initializes the enzyme object
 Returns   : The Restriction::Enzyme::MultiSite object
 Argument  : 

=cut

sub new {
    my($class, @args) = @_;
    my $self = $class->SUPER::new(@args);

    my ($others) =
            $self->_rearrange([qw(
                                  OTHERS
                                 )], @args);

    $others && $self->others($others);
    return $self;
}

=head2 others

 Title     : others
 Usage     : $re->others(@others);
 Function  : Gets/Sets the a list of other sites that this enzyme recoginizes
 Arguments : An array containing the other Bio::Restriction::Enzyme::MultiSite
             objects.
 Returns   : An array containing the other Bio::Restriction::Enzyme::MultiSite
             objects.

=cut

sub others {
    my $self = shift;
    push @{$self->{_others}}, @_ if @_;
    return unless $self->{_others};
    return @{$self->{'_others'}};
}


=head2 purge_others

 Title     : purge_others
 Usage     : $re->purge_references();
 Function  : Purges the set of references for this enzyme
 Arguments : 
 Returns   : 

=cut

sub purge_others {
    my ($self) = shift;
    $self->{_others} = [];

}


1;

