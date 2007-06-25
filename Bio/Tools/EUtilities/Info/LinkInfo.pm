# $Id$
#
# BioPerl module for Bio::Tools::EUtilities::Info::LinkInfo
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

Bio::Tools::EUtilities::Info::LinkInfo - class for storing einfo link data 

=head1 SYNOPSIS

    ## should not create instance directly; Bio::Tools::EUtilities does this ##

    # get a LinkInfo object using Bio:Tools::EUtilities    
    print "Link name: ",$link->get_link_name,"\n";
    print "Link name: ",$link->get_link_menu_name,"\n";
    print "Link desc: ",$link->get_link_description,"\n";
    print "DBFrom: ",$link->get_dbfrom,"\n"; # database linked from
    print "DBTo: ",$link->get_dbto,"\n"; # database linked to

=head1 DESCRIPTION

This class handles data output (XML) from einfo.

einfo is capable of returning two types of information: 1) a list of all
available databases (when called w/o parameters) and 2) information about a
specific database. The latter information includes the database description,
record count, and date/time stamp for the last update, among other things. It
also includes a list of fields (indices by which record data is stored which can
be used in queries) and links (crossrefs between related records in other
databases at NCBI). Data from the latter two are stored in two small subclasses
(Field and Link) which can be iterated through or retrieved all at once, as
demonstrated above. NOTE: Methods described for the Link and Field subclasses
are unique to those classes (as they retrieve data unique to those data types). 

Further documentation for Link and Field subclass methods is included below.

For more information on einfo see:

   http://eutils.ncbi.nlm.nih.gov/entrez/query/static/einfo_help.html

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other Bioperl
modules. Send your comments and suggestions preferably to one of the Bioperl
mailing lists. Your participation is much appreciated.

  bioperl-l@lists.open-bio.org               - General discussion
  http://www.bioperl.org/wiki/Mailing_lists  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track the bugs
and their resolution. Bug reports can be submitted via the web.

  http://bugzilla.open-bio.org/

=head1 AUTHOR 

Email cjfields at uiuc dot edu

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal
methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::Tools::EUtilities::Info::LinkInfo;

use base qw(Bio::Root::Root Bio::Tools::EUtilities::EUtilDataI);
use strict;
use warnings;

=head2 new

 Title    : new
 Note     : *** should not be called by end-users ***  
 Usage    : my $ct = Bio::Tools::EUtilities::Info::LinkInfo;
 Function : returns new LinkInfo instance
 Returns  : Bio::Tools::EUtilities::Info::LinkInfo instance
 Args     : none (all data added via _add_data, most methods are getters only)

=cut

sub new {
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(@args);
    my $eutil = $self->_rearrange([qw(EUTIL)], @args);
    $eutil ||= 'einfo';
    $self->eutil($eutil);
    $self->datatype('linkinfo');
    return $self;
}

=head2 get_dbto

 Title    : get_dbto
 Usage    : my $refd_db = $link->get_dbto;
 Function : returns database this link references (points to)
 Returns  : string
 Args     : none
 Note     : This is not the same as db()! (see DESCRIPTION for details)
 
=cut

sub get_dbto { return shift->{'_dbto'} }

=head2 get_dbfrom

 Title    : get_dbfrom
 Usage    : my $origdb = $link->get_dbfrom;
 Function : returns referring database
 Returns  : string
 Args     : none
 Note     : alias for get_db()

=cut

sub get_dbfrom { return shift->{'_dbfrom'} }

=head2 get_link_name

 Title    : get_link_name
 Usage    : $ln = $link->get_link_name;
 Function : returns raw link name (eutil-compatible)
 Returns  : string
 Args     : none

=cut

sub get_link_name {
    my $self = shift;
    if ($self->eutil eq 'elink') {
        return $self->{'_linkname'}
    } else {
        return $self->{'_name'}
    }
}

=head2 get_link_description

 Title    : get_link_description
 Usage    : $desc = $link->get_link_description;
 Function : returns the (more detailed) link description
 Returns  : string
 Args     : none

=cut

sub get_link_description { return shift->{'_description'} }

=head2 get_link_menu_name

 Title    : get_link_menu_name
 Usage    : my $mn = $link->get_link_menu_name;
 Function : returns formal menu name
 Returns  : string
 Args     : none

=cut

sub get_link_menu_name {
    my $self = shift;
    return $self->eutil eq 'elink' ? $self->{'_menutag'} : $self->{'_menu'};
}

=head2 get_priority

 Title    : get_priority
 Usage    : my $mn = $link->get_priority;
 Function : returns priority ranking
 Returns  : integer
 Args     : none
 Note     : only set when using elink and cmd set to 'acheck'
 
=cut

sub get_priority { return shift->{'_priority'} }

=head2 get_html_tag

 Title    : get_html_tag
 Usage    : my $tag = $link->get_html_tag;
 Function : returns HTML tag
 Returns  : string
 Args     : none
 Note     : only set when using elink and cmd set to 'acheck'
 
=cut

sub get_html_tag { return shift->{'_htmltag'} }

=head2 get_url

 Title    : get_url
 Usage    : my $url = $link->get_url;
 Function : returns URL string; note that the string isn't usable directly but
            has the ID replaced with the tag <@UID@>
 Returns  : string
 Args     : none
 Note     : only set when using elink and cmd set to 'acheck'
 
=cut

sub get_url { return shift->{'_url'} }

# private method

sub _add_data {
    my ($self, $simple) = @_;
    map { $self->{'_'.lc $_} = $simple->{$_} unless ref $simple->{$_}} keys %$simple;
}

1;

