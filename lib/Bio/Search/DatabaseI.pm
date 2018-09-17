#-----------------------------------------------------------------
#
# BioPerl module Bio::Search::DatabaseI
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Steve Chervitz <sac@bioperl.org>
#
# You may distribute this module under the same terms as perl itself
#-----------------------------------------------------------------

# POD documentation - main docs before the code

=head1 NAME

Bio::Search::DatabaseI - Interface for a database used in a sequence search

=head1 SYNOPSIS

Bio::Search::DatabaseI objects should not be instantiated since this
module defines a pure interface.

Given an object that implements the Bio::Search::DatabaseI  interface,
you can do the following things with it:

    $name = $db->name();

    $date = $db->date();

    $num_letters = $db->letters();

    $num_entries = $db->entries();

=head1 DESCRIPTION

This module defines methods for an object that provides metadata
information about a database used for sequence searching.

=head1 FEEDBACK

=head2 Mailing Lists 

User feedback is an integral part of the evolution of this and other
Bioperl modules.  Send your comments and suggestions preferably to one
of the Bioperl mailing lists.  Your participation is much appreciated.

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

Steve Chervitz E<lt>sac@bioperl.orgE<gt>

See L<the FEEDBACK section | FEEDBACK> for where to send bug reports and comments.

=head1 COPYRIGHT

Copyright (c) 2001 Steve Chervitz. All Rights Reserved.

=head1 DISCLAIMER

This software is provided "as is" without warranty of any kind.

=cut

=head1 APPENDIX

The rest of the documentation details each of the object methods.

=cut

# Let the code begin...

package Bio::Search::DatabaseI;

use strict;

use base qw(Bio::Root::RootI);


=head2 name

 Usage     : $name = $db->name();
 Purpose   : Get the name of the database searched.
 Returns   : String
 Argument  : n/a

=cut

sub name {
    my $self = shift;
    $self->throw_not_implemented;
}

=head2 date

 Usage     : $date = $db->date();
 Purpose   : Get the creation date of the queried database.
 Returns   : String
 Argument  : n/a

=cut

sub date {
    my $self = shift;
    $self->throw_not_implemented;
}


=head2 letters

 Usage     : $num_letters = $db->letters();
 Purpose   : Get the number of letters in the queried database.
 Returns   : Integer
 Argument  : n/a

=cut

sub letters {
    my $self = shift;
    $self->throw_not_implemented;
}


=head2 entries

 Usage     : $num_entries = $db->entries();
 Purpose   : Get the number of entries in the queried database.
 Returns   : Integer
 Argument  : n/a

=cut

sub entries {
    my $self = shift;
    $self->throw_not_implemented;
}

1;
