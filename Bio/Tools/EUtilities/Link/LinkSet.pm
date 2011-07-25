#
# BioPerl module for Bio::Tools::EUtilities::Link::LinkSet
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

Bio::Tools::EUtilities::Link::LinkSet - class for EUtils LinkSets

=head1 SYNOPSIS

  # ...

=head1 DESCRIPTION

  # ...

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

package Bio::Tools::EUtilities::Link::LinkSet;
use strict;
use warnings;

use base qw(Bio::Root::Root Bio::Tools::EUtilities::HistoryI);
use Bio::Tools::EUtilities::Link::UrlLink;
use Bio::Tools::EUtilities::Info::LinkInfo;

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
 Returns  : array of IDs
 Args     : none
 Notes    : Cmd                   Description
            acheck                same as get_submitted_ids
            lcheck                same as get_submitted_ids
            ncheck                same as get_submitted_ids
            prlinks               same as get_submitted_ids
            llinks                same as get_submitted_ids
            llinkslib             same as get_submitted_ids
            neighbor              linked IDs for database in get_database
            neighbor_history      linked IDs for database in get_database

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
    return @{$self->{'_sorted_id'}};
}

=head2 get_database

 Title    : get_database
 Usage    : my $db = $info->get_database;
 Function : returns single database name (eutil-compatible).  This is the
            queried database. For elinks (which have 'db' and 'dbfrom')
            this is equivalent to db/dbto (use get_dbfrom() to for the latter).
            Note that this only returns the first db; in some cases this may
            not be what you want (when multiple dbs are queried, for instance)
 Returns  : string
 Args     : none
 Notes    : with all elink cmd arguments

=cut

sub get_database {
    return ($_[0]->get_databases)[0];
}

=head2 get_db (alias for get_database)

=cut

sub get_db {
    return shift->get_database;
}

=head2 get_dbto (alias for get_database)

=cut

sub get_dbto {
    return shift->get_database;
}

=head2 get_databases

 Title    : get_databases
 Usage    : my $string = $linkset->get_databases;
 Function : retrieve databases referred to for this linkset
            these may be present as a single database or embedded in 
 Returns  : array of strings
 Args     : none

=cut

sub get_databases {
    my $self = shift;
    my %tmp;
    my @dbs = sort map {$_->get_database} 
        grep {!$tmp{$_->get_database}++} ($self->get_LinkInfo);
    unshift @dbs, $self->{'_dbto'} if $self->{'_dbto'} && !$tmp{$self->{'_dbto'}}++;
    return @dbs;
}

=head2 get_dbs (alias for get_databases)

=cut

sub get_dbs {
    return shift->get_databases;
}

=head2 get_dbfrom

 Title    : get_dbfrom
 Usage    : my $string = $linkset->get_dbfrom;
 Function : retrieve originating database for this linkset
 Returns  : string
 Args     : none

=cut

sub get_dbfrom { return shift->{'_dbfrom'} }

=head2 get_link_names

 Title    : get_link_names
 Usage    : my $string = $linkset->get_link_names;
 Function : retrieve eutil-compatible link names
 Returns  : array of strings
 Args     : none
 Notes    : Each LinkSet can hold multiple LinkInfo objects (each containing
            a link name). Also, some LinkSets define a single link name. This
            returns an array with all unique linknames globbed both sources, if
            present and defined

=cut

sub get_link_names {
    my ($self) = shift;
    my %tmps;
    my @lns;
    if ($self->{'_linkname'}) {
        push @lns, $self->{'_linkname'};
        $tmps{$self->{'_linkname'}}++;
    }
    push @lns, map {$_->get_link_name} $self->get_LinkInfo;
    return @lns;
}

=head2 get_link_name

 Title    : get_link_name
 Usage    : my $string = $linkset->get_link_name;
 Function : retrieve eutil-compatible link name
 Returns  : single link name
 Args     : none

=cut

sub get_link_name {
    return ($_[0]->get_link_names)[0];
}

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
        return $self->get_ids;
    } elsif ($self->{'_submitted_ids'}) {
        return @{$self->{'_submitted_ids'}};
    } else {
        return ();
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
    return exists $self->{'_has_scores'} ? 1 : 0;
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
    return %scores; 
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
 Args     : none
 Notes    : this checks cmd=lcheck (boolean for a linkout) and also backchecks
            cmd=acheck for databases with name 'LinkOut'

=cut

sub has_linkout {
    my $self = shift;
    if (exists $self->{'_haslinkout'}) {
        return $self->{'_haslinkout'} eq 'Y' ? 1 : 0;
    } else  {
        return (grep {$_ eq 'LinkOut'} $self->get_databases) ? 1 : 0;
    } 
}

=head2 has_neighbor

 Title    : has_neighbor
 Usage    : if ($linkset->has_neighbor) {...}
 Function : returns TRUE if the single ID present in this linkset has a neighbor
            in the same database
 Returns  : boolean
 Args     : none
 Notes    : this checks cmd=ncheck (boolean for a neighbor in same database); no
            other checks performed at this time

=cut

sub has_neighbor {
    my $self = shift;
    if (exists $self->{'_hasneighbor'}) {
        return $self->{'_hasneighbor'} eq 'Y' ? 1 : 0;
    } else {
        return 0;
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
 Usage    : my @links = $linkset->get_LinkInfo
 Function : returns all LinkInfo objects
 Returns  : list of Bio::Tools::EUtilities::Link::LinkInfo
 Args     : 

=cut

sub get_LinkInfo {
    my $self = shift;
    return ref $self->{'_linkinfo'} ? @{ $self->{'_linkinfo'} } : return ();
}

=head2 rewind

 Title    : rewind
 Usage    : $info->rewind() # rewinds all (default)
            $info->rewind('links') # rewinds only links
 Function : 'rewinds' (resets) specified interators (all if no arg)
 Returns  : none
 Args     : [OPTIONAL] String: 
            'all'       - all iterators (default)
            'linkinfo' or 'linkinfos'  - LinkInfo objects only
            'urllinks'   - UrlLink objects only

=cut

{
    my %VALID_DATA = ('linkinfo' => 'linkinfo',
                      'linkinfos' => 'linkinfo',
                      'urllinks' => 'urllinks');
    
    sub rewind {
        my ($self, $arg) = @_;
        $arg ||= 'all';
        if (exists $VALID_DATA{$arg}) {
            delete $self->{'_'.$arg.'_it'};
        } elsif ($arg eq 'all') {
            delete $self->{'_'.$_.'_it'} for values %VALID_DATA;
        }
    }
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
    if (exists $data->{IdList}->{Id}) {
        @{$self->{'_submitted_ids'}} = @{$data->{IdList}->{Id}} ;
    }
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
        $urldata->{dbfrom} = $data->{DbFrom} if exists $data->{DbFrom};
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
        $linkinfo->{dbfrom} = $data->{DbFrom} if exists $data->{DbFrom};
        my $obj = Bio::Tools::EUtilities::Info::LinkInfo->new(-eutil => 'elink',
                                                             -datatype => 'linkinfo',
                                                             -verbose => $self->verbose
                                                             );
        $obj->_add_data($linkinfo);
        push @{$self->{'_linkinfo'}}, $obj;
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
    my $level = shift || 0;
    my $pad = 20 - $level;
    #        order     method                    name
    my %tags = (1 => ['get_databases'         => 'DB'],
                2 => ['get_ids'               => 'ID'],
                3 => ['get_link_names'        => 'Link Names'],
                5 => ['get_submitted_ids'     => 'Submitted IDs'],
                6 => ['has_scores'            => 'Scores?'],
                7 => ['has_linkout'           => 'LinkOut?'],
                8 => ['has_neighbor'          => 'DB Neighbors?'],
                9 => ['get_webenv'            => 'WebEnv'],
                10 => ['get_query_key'        => 'Key'],
                );
    my $string;
    for my $tag (sort {$a <=> $b} keys %tags) {
        my ($m, $nm) = (@{$tags{$tag}}[0..1]);
        # using this awkward little construct to deal with both lists and scalars
        my @content = grep {defined $_} $self->$m();
        next unless @content;
        $string .= $self->_text_wrap(
                 sprintf("%-*s%-*s:",$level, '',$pad, $nm,),
                 ' ' x ($pad).':',
                 join(', ',@content))."\n";
    }
    while (my $li = $self->next_LinkInfo) {
        $string .= $li->to_string(4);
    }
    while (my $ui = $self->next_UrlLink) {
        $string .= $ui->to_string(4);
    }
    if ($self->has_scores) {
        $string .= "Scores:\n";
        my %scores = $self->get_scores;
        $string .= sprintf("%-*s%-*s%s\n",
            $level + 4, '',
            $pad - 4, 'ID', 'Score'
            );        
        for my $id ($self->get_ids) {
            $string .= sprintf("%-*s%-*s%s\n",
                $level + 4, '',
                $pad - 4, $id, $scores{$id}
                );
        }
    }
    $string .= "\n";
    return $string;
}

1;

