# $Id$
#
# BioPerl module for Bio::DB::EUtilities::ElinkData
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

Bio::DB::EUtilities::ElinkData 

=head1 SYNOPSIS

*** Give standard usage here

=head1 DESCRIPTION

*** Describe the object here

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

package Bio::DB::EUtilities::ElinkData;
use strict;
use warnings;

use Bio::Root::Root;
use Data::Dumper;
use vars '@ISA';

@ISA = 'Bio::Root::Root';

sub new {
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(@args);
    my ($command) = $self->_rearrange([qw(COMMAND)], @args);
    $command    && $self->elink_command($command);
    $self->{'_dbindex'} = 0;
    $self->{'_scoreindex'} = 0;
    $self->{'_linkcount'} = 0;
    $self->{'_scoredb_key'} = '';
    $self->{'_databases'} = [];
    $self->{'_linksetdb'} = [];
    return $self;
}

# this should make a deep copy of the XML data for each ElinkData Linkset

sub _add_set {
    my ($self, $ls) = @_;
    if (!$ls) {
        $self->throw('No linkset data!');
    }
    # is there any data returned
    return 0 unless exists $ls->{LinkSetDb};
    my $dbfrom = $ls->{DbFrom};
    $self->elink_dbfrom($dbfrom);
    my $query_ids = $ls->{IdList}->{Id};
    if (!ref($query_ids)) {
        my $tempid = $query_ids;
        $query_ids = [$tempid];
    }
    $self->elink_queryids($query_ids);
    
    my $ct = 0;
    
    # each linkset database
    for my $ls_db (@{ $ls->{LinkSetDb} }) {
        $ct++;
        my $dbto = $ls_db->{DbTo} ;
        push @{ $self->{'_databases'}}, $dbto;
        my $linkname = $ls_db->{LinkName};
        if (exists $ls_db->{Info} || exists $ls->{ERROR} || !exists $ls_db->{Link}) {
            my $err_msg = $ls_db->{Info} || $ls->{ERROR} || 'No Links!';
            my $ids = (ref($query_ids) =~ /array/i) ?
                            join q(,), @{$query_ids}: $query_ids;
            $self->warn("ELink Error for $dbto and ids $ids: $err_msg");
            next;
        }
        my @ids;
        for my $id_ref (@{ $ls_db->{Link} } ) {
            my $id = $id_ref->{Id};
            my $score = exists $id_ref->{Score} ? $id_ref->{Score} : undef;
            push @ids, $id;
            # set up in case there are multiple databases that return scores
            if ($score) {
                $self->{'_scores'}->{$dbto}->{$id} = $score;
                if (!($self->{'_has_scores'})) {
                    push @{ $self->{'_has_scores'} }, $dbto;
                }
            }
        }
        my $linkset = {
                       'LinkName' => $linkname,
                       'DbTo'     => $dbto,
                       'Id'       => \@ids,
                      };
        $self->debug('Linkset:',Dumper($linkset));
        push @{ $self->{'_linksetdb'}}, $linkset;
    }
    return 1; # good linkset
}

=head2 elink_dbfrom

 Title   : elink_dbfrom
 Usage   : $dbfrom = $linkset->elink_dbfrom;
 Function: gets/sets dbfrom value
 Returns : originating database
 Args    : originating database

=cut

sub elink_dbfrom {
    my $self = shift;
    return $self->{'_elink_dbfrom'} = shift if @_;
    return $self->{'_elink_dbfrom'};
}

=head2 elink_queryids

 Title   : elink_queryids
 Usage   : @ids = $linkset->elink_queryids;
 Function: gets/sets original query ID values (ref to array)
 Returns : array or array ref of IDs (based on wantarray)
 Args    : array ref of IDs

=cut

sub elink_queryids {
    my $self = shift;
    return $self->{'_elink_queryids'} = shift if @_;
    return @{ $self->{'_elink_queryids'} } if wantarray;
    return $self->{'_elink_queryids'};
}

=head2 elink_command

 Title   : elink_command
 Usage   : $cmd = $linkset->elink_command;
 Function: gets/sets cmd used for elink query
 Returns : string (cmd parameter)
 Args    : string (cmd parameter)

=cut

sub elink_command {
    my $self = shift;
    return $self->{'_elink_command'} = shift if @_;
    return $self->{'_elink_command'};
}

=head2 get_LinkIds_by_db

 Title   : get_LinkIds_by_db
 Usage   : @ids = $linkset->get_LinkIds_by_db('protein');
 Function: retrieves primary ID list based on the database for the object
 Returns : array or array ref of IDs (based on wantarray)
 Args    : None

=cut

sub get_LinkIds_by_db {
    my $self = shift;
    my $db = shift if @_;
    $self->throw("Must use database to access IDs") if !$db;
    my $ct = scalar(@{ $self->{'_linksetdb'} });
    return [] if $ct == 0; # no linksets, blank anon array
    for my $linkset (@{ $self->{'_linksetdb'}}) {
        my $dbto = $linkset->{DbTo};
        if ($dbto eq $db) {
            return @{ $linkset->{Id} } if wantarray;
            return $linkset->{Id};
            
        }
    }
    $self->warn("Couldn't find ids for database $db");
}

=head2 next_linkdb

 Title   : next_linkdb
 Usage   : while (my $db = $linkset->next_linkdb) {
 Function: iterates through list of database names in internal queue
 Returns : String (name of database)
 Args    : None

=cut

sub next_linkdb {
    my $self = shift;
    my $index = $self->_next_db_index;
    return if ($index > scalar($self->{'_databases'}));
    return $self->{'_databases'}->[$index] ;
}

=head2 get_all_linkdbs

 Title   : get_all_linkdbs
 Usage   : @dbs = $linkset->get_all_linkdbs;
 Function: returns all database names which contain IDs
 Returns : array or array ref of databases (based on wantarray)
 Args    : None

=cut

sub get_all_linkdbs {
    my $self = shift;
    return @{ $self->{'_databases'} } if wantarray;
    return $self->{'_databases'};
}

=head2 next_scoredb

 Title   : next_scoredb
 Usage   : while (my $db = $linkset->next_scoredb) {
 Function: iterates through list of database with score values
 Returns : String (name of database)
 Args    : None

=cut

sub next_scoredb {
    my $self = shift;
    my $index = $self->_next_scoredb_index;
    return if ($index > scalar($self->{'_has_scores'}));
    my $db = $self->{'_has_scores'}->[$index];
    $self->set_scoredb($db);
    return $db;
}

=head2 get_all_scoredbs

 Title   : get_all_scoredbs
 Usage   : @dbs = $linkset->get_all_scoredbs;
 Function: returns database names which contain scores
 Returns : array or array ref of databases (based on wantarray)
 Args    : None

=cut

sub get_all_scoredbs {
    my $self = shift;
    return @{ $self->{'_has_scores'} } if wantarray;
    return $self->{'_has_scores'}->[0];
}

=head2 get_score

 Title   : get_score
 Usage   : $score = $linkset->get_score($id);
 Function: returns score value for ID
 Returns : score value
 Args    : ID
 Note    : if multiple databases are returned with scores (rare but possible),
         : you must set the default score database using set_scoredb.  If you
         : use next_scoredb to iterate through the databases, this is done for you

=cut

sub get_score {
    my $self = shift;
    my $id = shift if @_;
    if (!$self->get_all_scoredbs) {
        $self->warn("No scores!");
        return;
    }
    if (!$id) {
        $self->throw("Must use ID to access scores");
        return;
    }
    my $db = exists $self->{'_scoredb'} ? $self->{'_scoredb'} :
             $self->get_all_scoredbs;
    if ( exists $self->{'_scores'}->{$db}->{$id} ) {
        return $self->{'_scores'}->{$db}->{$id};
    }
}

=head2 get_score_hash

 Title   : get_score_hash
 Usage   : %scores = $linkset->get_score_hash($database);
 Function: returns ID(key)-score(value) hash based on database name
 Returns : score value
 Args    : OPTIONAL : database name.  If there is only one score hash, returns
         : that hash, otherwise throws an exception

=cut

sub get_score_hash {
    my $self = shift;
    $self->warn("No scores!") if !$self->has_scores;
    my $db = exists $self->{'_scoredb'} ? $self->{'_scoredb'} : $self->has_scores;
    if (exists $self->{'_scores'}->{$db}) {
        return %{ $self->{'_scores'}->{$db} };
    }
}

=head2 set_scoredb

 Title   : set_scoredb
 Usage   : $linkset->set_scoredb('protein');
 Function: sets the database to retrieve scores from
 Returns : None
 Args    : database name

=cut

sub set_scoredb {
    my ($self, $key) = shift;
    $self->{'_scoredb'} if $key;
}

=head2 rewind_linkdbs

 Title   : rewind_linkdbs
 Usage   : $linkset->rewind_linkdbs;
 Function: resets the iterator for next_database
 Returns : None
 Args    : None

=cut

sub rewind_linkdbs {
    my $self = shift;
    $self->{'_dbindex'} = 0;
}

=head2 rewind_scoredbs

 Title   : rewind_scoredbs
 Usage   : $linkset->rewind_scoredbs;
 Function: resets the iterator, current database for next_scoredb
 Returns : None
 Args    : None

=cut

sub rewind_scoredbs {
    my $self = shift;
    $self->{'_scoreindex'} = 0;
    $self->{'_scoredb'} = '';
}

# private methods

#iterator for full database list
sub _next_db_index {
    my $self = shift;
    return $self->{'_dbindex'}++;
}

#iterator for score database list
sub _next_scoredb_index {
    my $self = shift;
    return $self->{'_scoreindex'}++;
}

1;

__END__