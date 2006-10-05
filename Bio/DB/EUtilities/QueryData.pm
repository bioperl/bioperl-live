# $Id$
#
# BioPerl module for Bio::DB::EUtilities::QueryData
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

Bio::DB::EUtilities::QueryData

=head1 SYNOPSIS

# use only in conjunction with Bio::DB::EUtilities::egquery

    my $eq = Bio::DB::EUtilities->new(
                                    -eutil  => 'egquery',
                                    -term   => 'Gallus gallus[organism]'
                                     );
    
    $eq->get_response;
    
    # get docsum objects
    while (my $qd  = $eq->next_query) {
        # do stuff here
    }

=head1 DESCRIPTION

This is a remedial object that acts as a container for QueryHit data from
Bio::DB::EUtilities::egquery.  It is in the very early stages of
development, so don't too be offended if the API changes.  It is possible
the various EUtilities container objects will be reorganized to have a
more consistent API; however, note that, due to the differences in the
actual data stored this may be next to impossible.

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

package Bio::DB::EUtilities::QueryData;
use strict;
use warnings;
use Data::Dumper;

use base qw(Bio::Root::Root);

sub new {
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(@args);
    $self->{'_querydata'} = [];
    return $self;
}

sub _add_data {
    my ($self, $qdata) = @_;
    if (!$qdata || ref($qdata) !~ /array/i) {
        $self->throw("Bad EGQuery data");
    }
    #$self->debug("Query data: ".Dumper($qdata));
    for my $item (@ {$qdata} ) {
        my %vals = %{ $item };
        push @{$self->{'_querydata'}}, \%vals;
    }
    return;
}

=head2 query_term

 Title   : query_term
 Usage   : $term = $qd->query_term();
 Function: get/set term used 
 Returns : string of term used
 Args    : OPTIONAL : string of term used

=cut

sub query_term {
    my $self = shift;
    return $self->{'_eq_term'} = shift if @_;
    return $self->{'_eq_term'};
}

=head2 get_all_dbs

 Title   : get_all_dbs
 Usage   : @dbs = $qd->get_all_dbs;
 Function: get array of databases
 Returns : array of database names
 Args    : none

=cut

sub get_all_dbs {
    my $self = shift;
    my @names = map {$_->{DbName}} @{ $self->{'_querydata'} };
    return @names;
}

=head2 get_querydata_by_db

 Title   : get_querydata_by_db
 Usage   : %data = $qd->get_querydata_by_db($db);
 Function: retrieve query data hash by database name
 Returns : hash containing all information for the query database
 Args    : REQUIRED: name of database 

=cut

sub get_querydata_by_db {
    my ($self, $db) = @_;
    $self->throw('Must supply database for get_querydata_by_db') if !$db;
    my ($data) = grep {$_->{DbName} eq $db} @{ $self->{'_querydata'} };
    return %{ $data } if $data;
    return;
}

=head2 get_Status_by_db

 Title   : get_Status_by_db
 Usage   : $st = $qd->get_Status_by_db($db);
 Function: retrieve the query status for a database
 Returns : string
 Args    : REQUIRED: name of database

=cut

sub get_Status_by_db {
    my ($self, $db) = @_;
    $self->throw('Must supply database for get_Status_by_db') if !$db;
    my ($data) = grep {$_->{DbName} eq $db} @{ $self->{'_querydata'} };
    return $data->{Status} if exists $data->{Status};
    return;
}

=head2 get_Count_by_db

 Title   : get_Count_by_db
 Usage   : $ct = $qd->get_Content_by_name($db);
 Function: retrieve the count (hits) for the database query
 Returns : Integer (number of hits)
 Args    : REQUIRED: name of database

=cut

sub get_Count_by_db {
    my ($self, $db) = @_;
    $self->throw('Must supply database for get_Count_by_db') if !$db;
    my ($data) = grep {$_->{DbName} eq $db} @{ $self->{'_querydata'} };
    return $data->{Count} if exists $data->{Count};
    return;
}

=head2 get_MenuName_by_db

 Title   : get_MenuName_by_db
 Usage   : $ct = $qd->get_MenuName_by_db($db);
 Function: retrieve the menu name tag for the database
 Returns : string
 Args    : REQUIRED: name of database

=cut

sub get_MenuName_by_db {
    my ($self, $db) = @_;
    $self->throw('Must supply database for get_Status_by_db') if !$db;
    my ($data) = grep {$_->{DbName} eq $db} @{ $self->{'_querydata'} };
    return $data->{MenuName} if exists $data->{MenuName};
    return;
}

# may add iterators; this will do for now

1;

__END__