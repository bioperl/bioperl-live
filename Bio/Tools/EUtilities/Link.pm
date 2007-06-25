# $Id$
#
# BioPerl module for Bio::Tools::EUtilities::Link
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

Bio::Tools::EUtilities::Link - general API for accessing data retrieved from elink queries

=head1 SYNOPSIS

...

=head1 DESCRIPTION

...

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

package Bio::Tools::EUtilities::Link;
use strict;
use warnings;
use Data::Dumper;

use base qw(Bio::Tools::EUtilities Bio::Tools::EUtilities::EUtilDataI);

use Bio::Tools::EUtilities::Link::LinkSet;

# private EUtilDataI method

{
    my %SUBCLASS = (
                    'LinkSetDb' => 'dblink',
                    'LinkSetDbHistory' => 'history',
                    'IdUrlList' => 'urllink',
                    'IdCheckList' => 'idcheck',
                    'NoLinks' => 'nolinks',
                    );

sub _add_data {
    my ($self, $data) = @_;
    # divide up per linkset
    if (!exists $data->{LinkSet}) {
        $self->warn("No linksets returned");
        return;
    }
    for my $ls (@{ $data->{LinkSet} }) {
        my $subclass;
        # caching for efficiency; no need to recheck
        if (!exists $self->{'_subclass_type'}) {
            ($subclass) = grep { exists $ls->{$_} } qw(LinkSetDb LinkSetDbHistory IdUrlList IdCheckList);
            $subclass ||= 'NoLinks';
            $self->{'_subclass_type'} = $subclass; 
        } else {
            $subclass = $self->{'_subclass_type'};
        }
        # split these up by ID, since using correspondence() clobbers them...
        if ($subclass eq 'IdUrlList' || $subclass eq 'IdCheckList') {
            my $list = $subclass eq 'IdUrlList' ? 'IdUrlSet' :
                $subclass eq 'IdCheckList' && exists $ls->{$subclass}->{IdLinkSet} ? 'IdLinkSet' :
                'Id';
            $ls->{$subclass} = $ls->{$subclass}->{$list};
        }
        # divide up linkset per link
        for my $ls_sub (@{ $ls->{$subclass} }) {
            for my $key (qw(WebEnv DbFrom IdList)) {
                $ls_sub->{$key} = $ls->{$key} if exists $ls->{$key};
            }
            my $obj = Bio::Tools::EUtilities::Link::LinkSet->new(-eutil => 'elink',
                                                    -datatype => $SUBCLASS{$subclass},
                                                    -verbose => $self->verbose);
            $obj->_add_data($ls_sub);
            if ($self->is_lazy) {
                return $obj
            } else {
                push @{$self->{'_linksets'}}, $obj;
            }
        }
    }
}

}

1;

__END__