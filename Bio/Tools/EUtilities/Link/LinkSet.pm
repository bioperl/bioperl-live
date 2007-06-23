# $Id$
#
# BioPerl module for Bio::Tools::EUtilities::Link::LinkSet
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

Bio::Tools::EUtilities::Link::LinkSet

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

package Bio::Tools::EUtilities::Link::LinkSet;
use base qw(Bio::Root::Root Bio::Tools::EUtilities::HistoryI);
use Bio::Tools::EUtilities::Link::UrlLink;
use Bio::Tools::EUtilities::Info::LinkInfo;
use Data::Dumper;

sub new {
    my ($class,@args) = @_;
    my $self = $class->SUPER::new(@args);
    my ($type) = $self->_rearrange([qw(DATATYPE)],@args);
    $type ||= 'linkset';
    $self->eutil('elink');
    $self->datatype($type);
    return $self;
}

=head2 get_ids

 Title    : get_ids
 Usage    : my @ids = $linkset->get_ids
 Function : returns list of retrieved IDs
 Returns  : array or array ref of IDs
 Args     : none

=cut

sub get_ids {
    my $self = shift;
    unless ($self->{'_sorted_id'}) {
        @{$self->{'_sorted_id'}} =
            sort {
                $self->{'_id'}->{$a}->[0] <=>
                $self->{'_id'}->{$b}->[0]
            } keys %{$self->{'_id'}};
    }
    return wantarray ? @{$self->{'_sorted_id'}} : $self->{'_sorted_id'};
}

=head2 get_dbto

 Title    : get_dbto
 Usage    : my $string = $linkset->get_dbto;
 Function : retrieve database referred to for this linkset
 Returns  : string
 Args     : none

=cut

sub get_dbto { return shift->{'_dbto'} }

=head2 get_dbfrom

 Title    : get_dbfrom
 Usage    : my $string = $linkset->get_dbfrom;
 Function : retrieve originating database for this linkset
 Returns  : string
 Args     : none

=cut

sub get_dbfrom { return shift->{'_dbfrom'} }

=head2 get_linkname

 Title    : get_linkname
 Usage    : my $string = $linkset->get_linkname;
 Function : retrieve eutil-compatible link name
 Returns  : string
 Args     : none

=cut

sub get_linkname { return shift->{'_linkname'} }

=head2 get_submitted_ids

 Title    : get_submitted_ids
 Usage    : my $string = $linkset->get_submitted_ids;
 Function : retrieve original ID list
 Returns  : string
 Args     : none

=cut

sub get_submitted_ids {
    my $self = shift;
    my $datatype = $self->datatype;
    if ($datatype eq 'idcheck' || $datatype eq 'urllink') {
        return wantarray ? $self->get_ids : [$self->get_ids];
    } elsif ($self->{'_submitted_ids'}) {
        return wantarray ? @{$self->{'_submitted_ids'}} : $self->{'_submitted_ids'};
    } else {
        return wantarray ? () : undef;
    }
}

=head2 has_scores

 Title    : has_scores
 Usage    : if (my $linkset->has_scores) {...}
 Function : returns TRUE if score data is present 
 Returns  : Boolean 
 Args     : none

=cut

sub has_scores {
    my $self = shift;
    return exists $self->{'_has_scores'} ? $self->{'_has_scores'} : 0;
}

=head2 get_scores

 Title    : get_scores
 Usage    : %scores = $linkset->get_scores;
 Function : returns flattened list or hash ref containing ID => score pairs
 Returns  : hash or hash ref (based on list or scalar context)
 Args     : none

=cut

sub get_scores {
    my $self = shift;
    # do we want to cache this or generate only when needed?  Likely won't be
    # called more than once...
    return unless $self->has_scores;
    my %scores = map {$_ => $self->{'_id'}->{$_}->[1]} keys %{$self->{'_id'}};
    return wantarray ? %scores : \%scores; 
}

=head2 get_score_by_id

 Title    : get_score_by_id
 Usage    : $score = $linkset->get_score_by_id($id);
 Function : returns the score for a particular primary ID
 Returns  : integer
 Args     : [REQUIRED] Primary ID for the score lookup

=cut

sub get_score_by_id {
    my ($self, $id) = @_;
    ($id && exists $self->{'_id'}->{$id}) ? return $self->{'_id'}->{$id}->[1] :
        return;
}

=head2 has_linkout

 Title    : has_linkout
 Usage    : if ($linkset->has_linkout) {...}
 Function : returns TRUE if the single ID present in this linkset has a linkout
 Returns  : boolean
 Args     : none (uses the ID in get_ids(), which for these cases is always 
            only one)

=cut

sub has_linkout {
    my $self = shift;
    if (exists $self->{'_haslinkout'}) {
        return $self->{'_haslinkout'} eq 'Y' ? 1 : 0;
    } else {
        $self->warn('No data present; did you use cmd lcheck?');
        return;
    }
}

=head2 has_neighbor

 Title    : has_neighbor
 Usage    : if ($linkset->has_neighbor) {...}
 Function : returns TRUE if the single ID present in this linkset has a neighbor
            in the same database
 Returns  : boolean
 Args     : none (uses the ID in get_ids(), which for these cases is always 
            only one)

=cut

sub has_neighbor {
    my $self = shift;
    if (exists $self->{'_hasneighbor'}) {
        return $self->{'_hasneighbor'} eq 'Y' ? 1 : 0;
    } else {
        $self->warn('No data present; did you use cmd ncheck?');
        return;
    }
}

=head2 next_UrlLink

 Title    : next_UrlLink
 Usage    : while (my $url = $linkset->next_UrlLink) {...}
 Function : iterate through UrlLink objects
 Returns  : Bio::Tools::EUtilities::Link::UrlLink
 Args     : 

=cut

sub next_UrlLink {
    my $self = shift;
    unless ($self->{"_urllinks_it"}) {
        my @ul = $self->get_UrlLinks;
        $self->{"_urllinks_it"} = sub {return shift @ul}
    }
    $self->{'_urllinks_it'}->();
}

=head2 get_UrlLinks

 Title    : get_UrlLinks
 Usage    : my @urls = $linkset->get_UrlLinks
 Function : returns all UrlLink objects
 Returns  : list of Bio::Tools::EUtilities::Link::UrlLink
 Args     : 

=cut

sub get_UrlLinks {
    my $self = shift;
    return ref $self->{'_urllinks'} ? @{ $self->{'_urllinks'} } : return;
}

=head2 next_LinkInfo

 Title    : next_LinkInfo
 Usage    : while (my $info = $linkset->next_LinkInfo) {...}
 Function : iterate through LinkInfo objects
 Returns  : Bio::Tools::EUtilities::Link::LinkInfo
 Args     : 

=cut

sub next_LinkInfo {
    my $self = shift;
    unless ($self->{"_linkinfo_it"}) {
        my @li = $self->get_LinkInfo;
        $self->{"_linkinfo_it"} = sub {return shift @li}
    }
    $self->{'_linkinfo_it'}->();
}

=head2 get_LinkInfo

 Title    : get_LinkInfo
 Usage    : my @urls = $linkset->get_LinkInfo
 Function : returns all LinkInfo objects
 Returns  : list of Bio::Tools::EUtilities::Link::LinkInfo
 Args     : 

=cut

sub get_LinkInfo {
    return ref $self->{'_linkinfo'} ? @{ $self->{'_linkinfo'} } : return;
}

# private methods and handlers

{
    my %DATA_HANDLER = (
        'IdList' => \&_add_submitted_ids,
        'Id'     => \&_add_retrieved_ids,
        'LinkInfo' => \&_add_linkinfo,
        'Link'   => \&_add_retrieved_ids,
        'ObjUrl' => \&_add_objurls,
        );

sub _add_data {
    my ($self, $data) = @_;
    for my $key (qw(IdList Link Id ObjUrl LinkInfo)) {
        next if !exists $data->{$key};
        my $handler = $DATA_HANDLER{$key};
        $self->$handler($data);
        delete $data->{$key};
    }
    # map the rest
    if ($self->datatype eq 'idcheck' && exists $data->{content}) {
        %{$self->{'_id'} } = ($data->{content} => [1]);
        delete $data->{content}
    }
    map {$self->{'_'.lc $_} = $data->{$_}} keys %$data;
}

}

sub _add_submitted_ids {
    my ($self, $data) = @_;
    @{$self->{'_submitted_ids'}} = @{$data->{IdList}->{Id}} ;
}

sub _add_retrieved_ids {
    my ($self, $data) = @_;
    # map all IDs to deal with possible scores
    # ID => {'count' = POSITION, 'score' => SCORE}
    if (exists $data->{Link}) {
        my $ct = 0;
        for my $link (@{$data->{Link}}) {
            if (exists $link->{Score}) {
                $self->{'_has_scores'}++;
                $self->{'_id'}->{$link->{Id}->[0]} = [ $ct++,$link->{Score}];
            } else {
                $self->{'_id'}->{$link->{Id}->[0]} = [ $ct++ ];
            }
        }
    }
    elsif (exists $data->{Id}) { # urls
        %{$self->{'_id'} } = ($data->{Id}->[0] => [1]);
    }
}

sub _add_objurls {
    my ($self, $data) = @_;
    for my $urldata (@{$data->{ObjUrl}}) {
        my $obj = Bio::Tools::EUtilities::Link::UrlLink->new(-eutil => 'elink',
                                                             -datatype => 'urldata',
                                                             -verbose => $self->verbose
                                                             );
        $obj->_add_data($urldata);
        push @{$self->{'_urllinks'}}, $obj;
    }
}

sub _add_linkinfo {
    my ($self, $data) = @_;
    for my $linkinfo (@{$data->{LinkInfo}}) {
        my $obj = Bio::Tools::EUtilities::Info::LinkInfo->new(-eutil => 'elink',
                                                             -datatype => 'linkinfo',
                                                             -verbose => $self->verbose
                                                             );
        $obj->_add_data($linkinfo);
        push @{$self->{'_linkinfo'}}, $obj;
    }
}

1;

__END__ 