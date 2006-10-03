# $Id$
#
# BioPerl module for Bio::DB::EUtilities::DocSum
#
# Cared for by Chris Fields
#
# Copyright Chris Fields
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code
# 
# Part of the EUtilities BioPerl package

=head1 NAME

Bio::DB::EUtilities::DocSum 

=head1 SYNOPSIS

*** Give standard usage here

=head1 DESCRIPTION

*** Describe the object here

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the
evolution of this and other Bioperl modules. Send
your comments and suggestions preferably to one
of the Bioperl mailing lists. Your participation
is much appreciated.

  bioperl-l@lists.open-bio.org               - General discussion
  http://www.bioperl.org/wiki/Mailing_lists  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to
help us keep track the bugs and their resolution.
Bug reports can be submitted via the web.

  http://bugzilla.open-bio.org/

=head1 AUTHOR 

Email cjfields at uiuc dot edu

=head1 APPENDIX

The rest of the documentation details each of the
object methods. Internal methods are usually
preceded with a _

=cut

# Let the code begin...

package Bio::DB::EUtilities::DocSum;
use strict;
use warnings;
use Data::Dumper;

#use Data::Dumper;

use base qw(Bio::Root::Root);

sub new {
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(@args);
    return $self;
}

# this should make a deep copy of the XML data for each docsum

sub _add_data {
    my ($self, $docsum) = @_;
    if (!$docsum || ref($docsum) !~ /array/i) {
        $self->throw("Bad ESummary DocSum");
    }
    #$self->debug("Docsum: ".Dumper($docsum));
    for my $item (@ {$docsum} ) {
        my ($name, $type, $content) = ($item->{Name},
                                       $item->{Type},
                                       $item->{content});
        # layered docsum data
        if ($type eq 'List' && exists $item->{Item}) {
            my $ds = $self->new(-verbose => $self->verbose);
            $ds->_add_data($item->{Item});
            $content = $ds;
        }
        push @{$self->{'_docdata'}}, [$name, $type, $content];
    }
}

sub esummary_id {
    my $self = shift;
    return $self->{'_esum_id'} = shift if @_;
    return $self->{'_esum_id'};
}

1;

__END__