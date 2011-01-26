#
# BioPerl module for Bio::Tools::EUtilities::Info
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

Einfo is capable of returning two types of information:

=over 3

=item * A list of all available databases (when called w/o parameters)

=item * Information about a specific database.

=back

The latter information includes the database description, record count, and
date/time stamp for the last update, among other things. It also includes a list
of fields (indices by which record data is stored which can be used in queries)
and links (crossrefs between related records in other databases at NCBI). Data
from the latter two are stored in two small subclasses (FieldInfo and LinkInfo)
which can be iterated through or retrieved all at once, as demonstrated above.
NOTE: Methods described for the LinkInfo and FieldInfo subclasses are unique to
those classes (as they retrieve data unique to those data types).

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

package Bio::Tools::EUtilities::Info;

use strict;
use warnings;
use base qw(Bio::Tools::EUtilities Bio::Tools::EUtilities::EUtilDataI);

use Bio::Tools::EUtilities::Info::LinkInfo;
use Bio::Tools::EUtilities::Info::FieldInfo;

=head2 rewind

 Title    : rewind
 Usage    : $info->rewind() # rewinds all (default)
            $info->rewind('links') # rewinds only links
 Function : 'rewinds' (resets) specified interators (all if no arg)
 Returns  : none
 Args     : [OPTIONAL] String: 
            'all'    - all iterators (default)
            'linkinfo'  - LinkInfo objects only
            'fieldinfo' - FieldInfo objects only

=cut


# private EUtilDataI method

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
    if (my @dbs = $self->get_databases) {
        $string .= sprintf("%-20s:%s\n\n", 'DB',
            $self->_text_wrap('', ' 'x20 .':', join(', ',@dbs)));
    }
    while (my $fi = $self->next_FieldInfo) {
        $string .= $fi->to_string."\n";
    }
    while (my $li = $self->next_LinkInfo) {
        $string .= $li->to_string."\n";
    }
    return $string;
}

1;
