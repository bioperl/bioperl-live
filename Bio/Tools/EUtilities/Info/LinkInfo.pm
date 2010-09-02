#
# BioPerl module for Bio::Tools::EUtilities::Info::LinkInfo
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

This class handles data output (XML) from both einfo and elink, and centers on
describing data that either describes how NCBI databases are linked together
via link names, or how databases are linked to outside databases (LinkOut).

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

=head2 Support 

Please direct usage questions or support issues to the mailing list:

I<bioperl-l@bioperl.org>

rather than to the module maintainer directly. Many experienced and 
reponsive experts will be able look at the problem and quickly 
address it. Please include a thorough description of the problem 
with code and data examples if at all possible.

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track the bugs
and their resolution. Bug reports can be submitted via the web.

  http://bugzilla.open-bio.org/

=head1 AUTHOR 

Email cjfields at bioperl dot org

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

=head2 get_database

 Title    : get_database
 Usage    : my $db = $info->get_database;
 Function : returns single database name (eutil-compatible).  This is the
            queried database. For elinks (which have 'db' and 'dbfrom')
            this is equivalent to db/dbto (use get_dbfrom() to for the latter)
 Returns  : string
 Args     : none

=cut

sub get_database {
    return shift->{'_dbto'};
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

=head2 get_dbfrom

 Title    : get_dbfrom
 Usage    : my $origdb = $link->get_dbfrom;
 Function : returns referring database
 Returns  : string
 Args     : none
 Note     :

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
    my %tags = (1 => ['get_link_name'         => 'Link Name'],
                2 => ['get_link_description'  => 'Description'],
                3 => ['get_dbfrom'            => 'DB From'],
                4 => ['get_dbto'              => 'DB To'],
                5 => ['get_link_menu_name'    => 'Menu Name'],
                6 => ['get_priority'          => 'Priority'],
                7 => ['get_html_tag'          => 'HTML Tag'],
                8 => ['get_url'               => 'URL'],
                );
    my $string = '';
    for my $tag (sort {$a <=> $b} keys %tags) {
        my ($m, $nm) = ($tags{$tag}->[0], $tags{$tag}->[1]);
        my $content = $self->$m();
        next unless $content;
        $string .= sprintf("%-*s%-*s%s\n",
            $level, '',
            $pad, $nm,
            $self->_text_wrap(':',
                 ' ' x ($pad).':',
                 $content ));
    }
    return $string;
}

1;

