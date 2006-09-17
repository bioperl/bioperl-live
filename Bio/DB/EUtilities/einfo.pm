# $Id$
#
# BioPerl module for Bio::DB::EUtilities::einfo
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

Bio::DB::EUtilities::einfo - NCBI database information

=head1 SYNOPSIS

my $einfo = Bio::DB::EUtilities->new(
                                     -eutil    => 'einfo',
                                     -db       => 'pubmed'
                                      );

print $einfo->get_response->content;

=head1 DESCRIPTION

L<EInfo|Bio::DB::EUtilities::einfo> queries provide
information about NCBI databases.  At this time, data is postprocessed
for a complete list of Entrez databases (when '-C<db>' is not set) or for
specific database information, number of entries, date of the last update, or
Field or Link information.  Using the base URL with no
parameters returns a list of all Entrez databases.

=head2 Parameters

The following are a general list of parameters that can be used to take
advantage of EInfo.  Up-to-date help for EInfo is available at this URL
(the information below is a summary of the options found there):

  http://eutils.ncbi.nlm.nih.gov/entrez/query/static/einfo_help.html

=over 3

=item C<db>

Database parameter.  This is optional; not setting this will return a list of
all the available Entrez databases.  If a database is specified, returned XML
data will provide the following information : field names, index term counts,
last update, and available links for each Entrez database.

=back

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

package Bio::DB::EUtilities::einfo;
use strict;
use warnings;
use Bio::DB::EUtilities;
use XML::Simple;
use Data::Dumper;

use vars qw(@ISA $EUTIL);

@ISA = qw(Bio::DB::EUtilities);

BEGIN {
    $EUTIL = 'einfo';
}

sub _initialize {
    my ($self, @args ) = @_;
    $self->SUPER::_initialize(@args);
    # set by default
    $self->_eutil($EUTIL);
}

=head2 parse_response

 Title   : parse_response
 Usage   : $db->_parse_response($content)
 Function: parse out response for cookie
 Returns : empty
 Args    : none
 Throws  : 'unparseable output exception'

=cut

sub parse_response {
    my $self    = shift;
    my $response = shift if @_;
    if (!$response || !$response->isa("HTTP::Response")) {
        $self->throw("Need HTTP::Response object");
    }
    my $xs = XML::Simple->new();
    my $simple = $xs->XMLin($response->content,
                            forcearray => [qw(DbName Field Link)]);
    $self->debug("Response dumper:\n".Dumper($simple));
    # check for errors
    if ($simple->{ERROR}) {
        my $error = $simple->{ERROR} ? $simple->{ERROR} : 'No data returned';
        $self->throw("NCBI einfo nonrecoverable error: ".$error);
    }
    if (exists $simple->{DbList}->{DbName}) {
        $self->{'_einfo_dbname'} = $simple->{DbList}->{DbName};
        return;
    }
    # start setting internal variables
    for my $key (sort keys %{ $simple->{DbInfo} }) {
        my $data =
        ($key eq 'FieldList') ? $simple->{DbInfo}->{$key}->{Field} :
        ($key eq 'LinkList' ) ? $simple->{DbInfo}->{$key}->{Link}  :
        $simple->{DbInfo}->{$key};
        $self->_set_einfo_data($key, $data);
    }
}

sub einfo_dbs {
    my $self = shift;
    if (wantarray) {
        if( ref($self->{'_einfo_dbname'} ) =~ m{ARRAY}i  ) {
            return @{ $self->{'_einfo_dbname'} };
        }
    }
    return $self->{'_einfo_dbname'};
}

=head2 einfo_dbfield_info

 Title   : einfo_dbfield_info
 Usage   : @fields = $info->einfo_dbfield_info;
 Function: gets array of hashes with field information
 Returns : An array or array reference (based on wantarray) of hashes
           with information about each field 
 Args    : None (this is set using the _set_einfo_data method)

=cut

sub einfo_dbfield_info {
    my $self = shift;
    return @{ $self->{'_einfo_fieldlist'} } if wantarray;
    return $self->{'_einfo_fieldlist'};
}

=head2 einfo_dblink_info

 Title   : einfo_dblink_info
 Usage   : @links = $info->einfo_dblink_info;
 Function: gets array of hashes with link information
 Returns : An array or array reference (based on wantarray) of hashes
           with information about each link 
 Args    : None (this is set using the _set_einfo_data method)

=cut

sub einfo_dblink_info {
    my $self = shift;
    return @{ $self->{'_einfo_linklist'} } if wantarray;
    return $self->{'_einfo_linklist'};
}

=head2 einfo_db_lastupdate

 Title   : einfo_db_last_update
 Usage   : $date = $info->einfo_db_lastupdate;
 Function: returns last date database was updated
 Returns : String containing date
 Args    : None (this is set using the _set_einfo_data method)

=cut

sub einfo_db_lastupdate {
    my $self = shift;
    return $self->{'_einfo_lastupdate'};
}

=head2 einf_db_desc

 Title   : einfo_db_desc
 Usage   : $desc = $info->einfo_db_desc;
 Function: returns database description
 Returns : String containing descriptions
 Args    : None (this is set using the _set_einfo_data method)

=cut

sub einfo_db_desc {
    my $self = shift;
    return $self->{'_einfo_description'};
}

=head2 einfo_db_count

 Title   : einfo_db_count
 Usage   : $count = $info->einfo_db_count;
 Function: returns database record count
 Returns : Integer (number of database records)
 Args    : None (this is set using the _set_einfo_data method)

=cut

sub einfo_db_count {
    my $self = shift;
    return $self->{'_einfo_count'};
}

# no methods for MenuName

# set method

sub _set_einfo_data {
    my ($self, $key, $data) = @_;
    $self->throw('No data') if (!$data || !$key);
    my $info_key = '_einfo_'.lc($key);
    # fix to make a deep copy of data
    $self->{$info_key} = $data;
    return;
}    

1;
__END__