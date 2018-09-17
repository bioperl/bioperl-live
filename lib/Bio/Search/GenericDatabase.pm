#-----------------------------------------------------------------
#
# BioPerl module Bio::Search::GenericDatabase
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Steve Chervitz <sac@bioperl.org>
#
# You may distribute this module under the same terms as perl itself
#-----------------------------------------------------------------

# POD documentation - main docs before the code

=head1 NAME

Bio::Search::GenericDatabase - Generic implementation of Bio::Search::DatabaseI

=head1 SYNOPSIS

    use Bio::Search::GenericDatabase;

    $db = Bio::Search::GenericDatabase->new( -name => 'my Blast db',
					     -date => '2001-03-13',
					     -length => 2352511,
					     -entries => 250000 );

    $name = $db->name();
    $date = $db->date();
    $num_letters = $db->letters();
    $num_entries = $db->entries();

=head1 DESCRIPTION

This module provides a basic implementation of L<Bio::Search::DatabaseI>.
See documentation in that module for more information.

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

package Bio::Search::GenericDatabase;

use strict;

use base qw(Bio::Root::Root Bio::Search::DatabaseI);

sub new {
    my ($class, @args) = @_; 
    my $self = $class->SUPER::new(@args);
    my ($name, $date, $length, $ents) = 
        $self->_rearrange( [qw(NAME DATE LENGTH ENTRIES)], @args);

    $name    && $self->name($name);
    $date    && $self->date($date);
    $length  && $self->letters($length);
    $ents    && $self->entries($ents);

    return $self;
}

=head2 name

See L<Bio::Search::DatabaseI::name>() for documentation

This implementation is a combined set/get accessor.

=cut

#---------------
sub name {
#---------------
    my $self = shift;
    if(@_) { 
        my $name = shift;
        $name =~ s/(^\s+|\s+$)//g;
        $self->{'_db'} = $name;
    }
    $self->{'_db'};
}

=head2 date

See L<Bio::Search::DatabaseI::date>() for documentation

This implementation is a combined set/get accessor.

=cut

#-----------------------
sub date {
#-----------------------
    my $self = shift;
    if(@_) { $self->{'_dbDate'} = shift; }
    $self->{'_dbDate'};
}


=head2 letters

See L<Bio::Search::DatabaseI::letters>() for documentation

This implementation is a combined set/get accessor.

=cut

#----------------------
sub letters {
#----------------------
    my $self = shift;
    if(@_) { $self->{'_dbLetters'} = shift; }
    $self->{'_dbLetters'};
}


=head2 entries

See L<Bio::Search::DatabaseI::entries>() for documentation

This implementation is a combined set/get accessor.

=cut

#------------------
sub entries {
#------------------
    my $self = shift;
    if(@_) { $self->{'_dbEntries'} = shift; }
    $self->{'_dbEntries'};
}

1;
