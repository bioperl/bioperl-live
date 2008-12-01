# $Id$
#
# BioPerl module for Bio::Tools::EUtilities::Link::UrlLink
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

Bio::Tools::EUtilities::Link::UrlLink

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

package Bio::Tools::EUtilities::Link::UrlLink;

use base qw(Bio::Root::Root Bio::Tools::EUtilities::EUtilDataI);
use Data::Dumper;

=head2 get_dbfrom

 Title    : get_dbfrom
 Usage    : 
 Function : 
 Returns  : 
 Args     : 

=cut

sub get_dbfrom { return shift->{'_dbfrom'}; }

=head2 get_attribute

 Title    : get_attribute
 Usage    : 
 Function : 
 Returns  : 
 Args     : 

=cut

sub get_attribute { return shift->{'_attribute'}; }

=head2 get_icon_url

 Title    : get_icon_url
 Usage    : 
 Function : 
 Returns  : 
 Args     : 

=cut

sub get_iconurl { return shift->{'_iconurl'}; }

=head2 get_subject_type

 Title    : 
 Usage    : 
 Function : 
 Returns  : 
 Args     : 

=cut

sub get_subject_type { return shift->{'_subjecttype'}; }

=head2 get_url

 Title    : get_url
 Usage    : 
 Function : 
 Returns  : 
 Args     : 

=cut

sub get_url {
    my $self = shift;
    # fix Entrz LinkOut URLS without the full URL
    if ($self->{'_url'} && $self->{'_url'} =~ m{^/}) {
        $self->{'_url'} = 'http://www.ncbi.nih.gov'.$self->{'_url'};
    }
    return $self->{'_url'};
}

=head2 get_link_name

 Title    : get_link_name
 Usage    : 
 Function : 
 Returns  : 
 Args     : 

=cut

sub get_link_name { return shift->{'_linkname'};  }

=head2 get_provider_name

 Title    : get_provider_name
 Usage    : 
 Function : 
 Returns  : 
 Args     : 

=cut

sub get_provider_name { return shift->{'_provider_name'}; }

=head2 get_provider_abbr

 Title    : get_provider_abbr
 Usage    : 
 Function : 
 Returns  : 
 Args     : 

=cut

sub get_provider_abbr { return shift->{'_provider_nameabbr'}; }

=head2 get_provider_id

 Title    : get_provider_id
 Usage    : 
 Function : 
 Returns  : 
 Args     : 

=cut

sub get_provider_id { return shift->{'_provider_id'}[0]; }

=head2 get_provider_iconurl

 Title    : get_provider_iconurl
 Usage    : 
 Function : 
 Returns  : 
 Args     : 

=cut

sub get_provider_iconurl { return shift->{'_provider_iconurl'}; }

=head2 get_provider_url

 Title    : get_provider_url
 Usage    : 
 Function : 
 Returns  : 
 Args     : 

=cut

sub get_provider_url { return shift->{'_provider_url'}; }

# private method

sub _add_data {
    my ($self, $data) = @_;
    if (exists $data->{Provider}) {
        map {$self->{'_provider_'.lc $_} = $data->{Provider}->{$_};
            } keys %{$data->{Provider}};
        delete $data->{Provider};
    }
    map {$self->{'_'.lc $_} = $data->{$_} if $data->{$_}} keys %$data;
}

1;

__END__ 
