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

=head2 get_ids

 Title    : get_ids
 Usage    : my @ids = $eutil->get_ids
 Function : returns list of retrieved IDs
 Returns  : array or array ref of IDs
 Args     : [optional] a single string (database name) or a callback (code ref)
            which is passed the LinkSet.
            
            In most cases no arg is needed. However when multiple databases are
            queried each search returns a LinkSet with it's own database, IDs,
            linkname, etc., with some databases reporting back more than one
            group of IDs (eg more than one linkset).

=cut

sub get_ids {
    my ($self, $request) = @_;
    $self->parse_data unless $self->data_parsed;
    my @ids;
    if ($request) {
        if (ref $request eq 'CODE') {
            push @ids, map {$_->get_ids }
                grep { $request->($_) } $self->get_LinkSets;
        } else {
            push @ids, map {$_->get_ids }
                grep {$_->get_dbto eq $request} $self->get_LinkSets;
        }
    } else {
        $self->warn('Multiple database present, IDs will be globbed together')
            if $self->get_dbs > 1;
        push @ids, map {$_->get_ids } $self->get_LinkSets;
    }
    return wantarray ? @ids : \@ids;
}

=head2 get_dbs

 Title    : get_dbs
 Usage    : my @dbs = $eutil->get_dbs
 Function : returns list of databases linked to in linksets
 Returns  : array of databases
 Args     : none

=cut

sub get_dbs {
    my $self = shift;
    $self->parse_data unless $self->data_parsed;
    unless (exists $self->{'_db'}) {
        my %temp;
        # make sure unique db is returned
        # do the linksets have a db? (URLs, db checks do not)
        
        push @{$self->{'_db'}}, map {$_->get_dbto}
            grep { $_->get_dbto ? !$temp{$_->get_dbto}++: 0 } $self->get_LinkSets;
    }
    return @{$self->{'_db'}};
}

=head2 next_LinkSet

 Title    : next_LinkSet
 Usage    : 
 Function : 
 Returns  : 
 Args     : 

=cut

sub next_LinkSet {
    my $self = shift;
    $self->parse_data unless $self->data_parsed;
    unless ($self->{"_linksets_it"}) {
        my @ls = $self->get_LinkSets;
        $self->{"_linksets_it"} = sub {return shift @ls}
    }
    $self->{'_linksets_it'}->();
}

=head2 get_LinkSets

 Title    : get_LinkSets
 Usage    : 
 Function : 
 Returns  : 
 Args     : 

=cut

# add support for retrieval of data if lazy parsing is enacted

sub get_LinkSets {
    my $self = shift;
    $self->parse_data unless $self->data_parsed;
    return ref $self->{'_linksets'} ? @{ $self->{'_linksets'} } : return;
}

=head2 rewind

 Title    : rewind
 Usage    : 
 Function : 
 Returns  : 
 Args     : 

=cut

sub rewind {
    my $self = shift;
    delete $self->{'_linksets_it'}
}

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
    if (exists $data->{LinkSet}) {
        for my $ls (@{ $data->{LinkSet} }) {
            my $subclass;
            # caching for efficiency
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
            }
        }
    }
    $self->{'_parsed'} = 1;
}

}

1;

__END__