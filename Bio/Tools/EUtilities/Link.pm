#
# BioPerl module for Bio::Tools::EUtilities::Link
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
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

  ...TODO

=head1 DESCRIPTION

Bio::Tools::EUtilities::Link is a loadable plugin for Bio::Tools::EUtilities
that specifically handles NCBI elink-related data.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the
evolution of this and other Bioperl modules. Send
your comments and suggestions preferably to one
of the Bioperl mailing lists. Your participation
is much appreciated.

  bioperl-l@lists.open-bio.org               - General discussion
  http://www.bioperl.org/wiki/Mailing_lists  - About the mailing lists

=head2 Support 

Please direct usage questions or support issues to the mailing list:

I<bioperl-l@bioperl.org>

rather than to the module maintainer directly. Many experienced and 
reponsive experts will be able look at the problem and quickly 
address it. Please include a thorough description of the problem 
with code and data examples if at all possible.

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to
help us keep track the bugs and their resolution.
Bug reports can be submitted via the web.

  https://redmine.open-bio.org/projects/bioperl/

=head1 AUTHOR 

Email cjfields at bioperl dot org

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
        # attempt to catch linkset errors
        if (exists $ls->{ERROR}) {
            my ($error, $dbfrom) = ($ls->{ERROR},$ls->{DbFrom});
            $self->warn("NCBI LinkSet error: $dbfrom: $error\n");
            # try to save the rest of the data, if any
            next;
        }
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
            push @{$self->{'_linksets'}}, $obj;
            # push only potential history-carrying objects into history queue
            if ($subclass eq 'LinkSetDbHistory') {
                push @{$self->{'_histories'}}, $obj;
            }
        }
    }
}

}

=head2 to_string

 Title    : to_string
 Usage    : $foo->to_string()
 Function : converts current object to string
 Returns  : none
 Args     : (optional) simple data for text formatting
 Note     : Used generally for debugging and for various print methods

=cut

sub to_string {
    my $self = shift;
    my $string = $self->SUPER::to_string;
    while (my $ls = $self->next_LinkSet) {
        $string .= $ls->to_string;
    }
    return $string;
}

1;
