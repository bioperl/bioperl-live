# $Id$
#
# BioPerl module for Bio::Tools::EUtilities::Info
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

Bio::Tools::EUtilities::Info - interface class for storing einfo data 

=head1 SYNOPSIS

  #### should not create instance directly; Bio::Tools::EUtilities does this ####
  
  my $info = Bio::Tools::EUtilities->new(-eutil => 'einfo',
                                         -file => 'einfo.xml');
  # can also use '-response' (for HTTP::Response objects) or '-fh' (for filehandles)
  
  # print available databases (if data is present)
  
  print join(', ',$info->get_available_databases),"\n";
  
  # get database info
  
  my $db = $info->get_database; # in case you forgot...
  my $desc = $info->get_description;
  my $nm = $info->get_menu_name;
  my $ct = $info->get_record_count;
  my $dt = $info->get_last_update;
  
  # EUtilDataI interface methods
  
  my $eutil = $info->eutil;
  my $type = $info->datatype;
  
  # iterate through Field and Link objects
  
  while (my $field = $info->next_Field) {
      print "Field code: ",$field->get_field_code,"\n";
      print "Field name: ",$field->get_field_name,"\n";
      print "Field desc: ",$field->get_field_description,"\n";
      print "DB  : ",$field->get_database,"\n";
      print "Term ct   : ",$field->get_term_count,"\n";
      for my $att (qw(is_date is_singletoken is_hierarchy is_hidden is_numerical)) {
          print "\tField $att\n" if $field->$att;
      }
  }
  
  my @fields = $info->get_Fields; # grab them all (useful for grep)
  
  while (my $link = $info->next_LinkInfo) {
      print "Link name: ",$link->get_link_name,"\n";
      print "Link desc: ",$link->get_link_description,"\n";
      print "DBFrom: ",$link->get_dbfrom,"\n"; # same as get_database()
      print "DBTo: ",$link->get_dbto,"\n"; # database linked to
  }
  
  my @links = $info->get_LinkInfo; # grab them all (useful for grep)
  
  $info->rewind(); # rewinds all iterators
  $info->rewind('links'); # rewinds Link iterator
  $info->rewind('fields'); # rewinds Field iterator

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

package Bio::Tools::EUtilities::Info;

use strict;
use warnings;
use base qw(Bio::Tools::EUtilities Bio::Tools::EUtilities::EUtilDataI);

use Bio::Tools::EUtilities::Info::LinkInfo;
use Bio::Tools::EUtilities::Info::FieldInfo;

=head2 new

 Title    : new
 Usage    : ***Should not be used directly; instance generated as follows***
            my $info = Bio::Tools::EUtilities->new(-eutil => 'einfo',
                                                -datatype => 'base');
 Function : create new Bio::Tools::EUtilities::Info data instance
 Returns  : new Bio::Tools::EUtilities::Info object based on type
 Args     : [REQUIRED] -datatype : string, 'info', 'field', 'link'

=cut

=head1 Bio::Tools::EUtilities::EUtilDataI methods

See Bio::Tools::EUtilities::EUtilDataI for details

=cut

=head2 add_data

=cut

sub _add_data {
    my ($self, $simple) = @_;
    if (exists $simple->{DbList} &&
        exists $simple->{DbList}->{DbName}) {
        $self->{'_available_databases'} = $simple->{DbList}->{DbName};
    }
    # start setting internal variables
    if (exists $simple->{DbInfo}) {
        for my $key (sort keys %{ $simple->{DbInfo} }) {
            my $data = 
            ($key eq 'FieldList') ? $simple->{DbInfo}->{$key}->{Field} :
            ($key eq 'LinkList' ) ? $simple->{DbInfo}->{$key}->{Link}  :
            $simple->{DbInfo}->{$key};
            if ($key eq 'FieldList' || $key eq 'LinkList') {
                for my $chunk (@{$data}) {
                    if (exists $simple->{DbInfo}->{DbName}) {
                        $chunk->{DbFrom} = $simple->{DbInfo}->{DbName};
                    }
                    my $type = ($key eq 'FieldList') ? 'FieldInfo' : 'LinkInfo';
                    my $obj = "Bio::Tools::EUtilities::Info::$type"->new(
                                           -eutil => 'einfo',
                                           -type => lc $type,
                                        -verbose => $self->verbose);
                    $obj->_add_data($chunk);
                    push @{ $self->{'_'.lc $type} }, $obj;
                }
            } else {
                $self->{'_'.lc $key} = $data;
            }
        }
    } else {
        map { $self->{'_'.lc $_} = $simple->{$_} unless ref $simple->{$_}} keys %$simple;
    }
    $self->{'_parsed'} = 1;
}

=head2 eutil

 Title    : eutil
 Usage    : $eutil->$foo->eutil
 Function : Get/Set eutil
 Returns  : string
 Args     : string (eutil)
 Throws   : on invalid eutil
 
=cut

=head2 datatype

 Title    : datatype
 Usage    : $type = $foo->datatype;
 Function : Get/Set data object type
 Returns  : string
 Args     : string

=cut

=head1 Bio::Tools::EUtilities::Info methods

=head2 get_available_databases

 Title    : get_available_databases
 Usage    : my @dbs = $info->get_available_databases
 Function : returns list of available eutil-compatible database names
 Returns  : Array of strings 
 Args     : none

=cut

sub get_available_databases {
    my $self = shift;
    $self->parse_data unless $self->data_parsed;
    ($self->{'_available_databases'}) ?
        return @{($self->{'_available_databases'})} :
        return ();
}

=head2 get_record_count

 Title    : get_record_count
 Usage    : my $ct = $eutil->get_record_count;
 Function : returns database record count
 Returns  : integer
 Args     : none

=cut

sub get_record_count {
    my $self = shift;
    $self->parse_data unless $self->data_parsed;
    return $self->{'_count'}
}

=head2 get_last_update

 Title    : get_last_update
 Usage    : my $time = $info->get_last_update;
 Function : returns string containing time/date stamp for last database update
 Returns  : integer
 Args     : none

=cut

sub get_last_update {
    my $self = shift;
    $self->parse_data unless $self->data_parsed;
    return $self->{'_lastupdate'}
}

=head2 get_menu_name

 Title    : get_menu_name
 Usage    : my $nm = $info->get_menu_name;
 Function : returns string of database menu name
 Returns  : string
 Args     : none
 
=cut

sub get_menu_name {
    my $self = shift;
    $self->parse_data unless $self->data_parsed;    
    exists $self->{'_menuname'} ? return $self->{'_menuname'} :
    exists $self->{'_menu'} ? return $self->{'_menu'} :
    return;
}

=head2 get_description

 Title    : get_description
 Usage    : my $desc = $info->get_description;
 Function : returns database description
 Returns  : string
 Args     : none

=cut

sub get_description {
    my $self = shift;
    $self->parse_data unless $self->data_parsed;        
    return $self->{'_description'};
}

=head2 get_db

 Title    : get_db
 Usage    : my $db = $info->get_db;
 Function : returns database name (eutil-compatible)
 Returns  : string
 Args     : none
 
=cut

sub get_db {
    my $self = shift;
    $self->parse_data unless $self->data_parsed;
    return $self->{'_dbname'};
}

=head2 get_database

 Title    : get_database
 Usage    : my $db = $eutil->get_database;
 Note     : alias for get_db()

=cut

*get_database = \&get_db;

=head2 next_FieldInfo

 Title    : next_FieldInfo
 Usage    : while (my $field = $info->next_FieldInfo) {...}
 Function : iterate through FieldInfo objects
 Returns  : Field object
 Args     : none
 Note     : uses callback() for filtering if defined for 'fields'
 
=cut

sub next_FieldInfo {
    my ($self, $cb) = @_;
    unless ($self->{'_fieldinfo_it'}) {
        $self->throw("Callback must be a code reference")
            if $cb && ref $cb ne 'CODE';
        my $fieldcount = $self->get_FieldInfo;
        my $current = 0;
        $self->{"_fieldinfo_it"} = sub {
            while ($current < $fieldcount) {
                if ($cb) {
                    $cb->($self->{'_fieldinfo'}->[$current++]) ?
                    return $self->{'_fieldinfo'}->[$current] :
                    next;
                } else {
                    return $self->{'_fieldinfo'}->[$current++]
                }
            }
        }
    }    
    $self->{'_fieldinfo_it'}->(); 
}

=head2 get_FieldInfo

 Title    : get_FieldInfo
 Usage    : my @fields = $info->get_FieldInfo;
 Function : returns list of FieldInfo objects
 Returns  : array (FieldInfo objects)
 Args     : none

=cut

sub get_FieldInfo {
    my $self = shift;
    $self->parse_data unless $self->data_parsed;        
    return ref $self->{'_fieldinfo'} ? @{ $self->{'_fieldinfo'} } : return 0;
}

=head2 next_LinkInfo

 Title    : next_LinkInfo
 Usage    : while (my $link = $info->next_LinkInfo) {...}
 Function : iterate through LinkInfo objects
 Returns  : LinkInfo object
 Args     : [OPTIONAL] callback; checks object and returns TRUE if wanted
 Note     : uses callback() for filtering if defined for 'links'
 
=cut

sub next_LinkInfo {
    my $self = shift;
    unless ($self->{'_linkinfo_it'}) {
        my $cb;
        $self->throw("Callback must be a code reference")
            if $cb && ref $cb ne 'CODE';
        my $linkcount = $self->get_LinkInfo;
        my $current = 0;
        $self->{"_linkinfo_it"} = sub {
            while ($current < $linkcount) {
                if ($cb) {
                    $cb->($self->{'_linkinfo'}->[$current++]) ?
                    return $self->{'_linkinfo'}->[$current] :
                    next;
                } else {
                    return $self->{'_linkinfo'}->[$current++]
                }
            }
        }
    }    
    $self->{'_linkinfo_it'}->();    
}

=head2 get_LinkInfo

 Title    : get_LinkInfo
 Usage    : my @links = $info->get_LinkInfo;
 Function : returns list of LinkInfo objects
 Returns  : array (LinkInfo objects)
 Args     : none

=cut

sub get_LinkInfo {
    my $self = shift;
    $self->parse_data unless $self->data_parsed;        
    return ref $self->{'_linkinfo'} ? @{ $self->{'_linkinfo'} } : return 0;
}

=head2 rewind

 Title    : rewind
 Usage    : $info->rewind() # rewinds all (default)
            $info->rewind('links') # rewinds only links
 Function : 'rewinds' (resets) specified interators (all if no arg)
 Returns  : none
 Args     : [OPTIONAL] String: 
            'all'    - all iterators (default)
            'links'  - link objects only
            'fields' - field objects only

=cut

{
    my %VALID_DATA = ('links' => 'linkinfo',
                      'fields' => 'fieldinfo');
    
    sub rewind {
        my ($self, $arg) = @_;
        $arg ||= 'all';
        if (exists $VALID_DATA{$arg}) {
            delete $self->{"_$VALID_DATA{$arg}_it"};
        } elsif ($arg eq 'all') {
            delete $self->{'_'.$_.'_it'} for values %VALID_DATA;
        }
    }
}

1;

