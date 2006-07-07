# $Id$

# simple object to hold NCBI score information, links, booleans, and other info
# from elink queries; API to change dramatically

# this should hold all the linksets for one group of IDs and should
# accomodate all cmd types.

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
    $command    && $self->command($command);
    return $self;
}

# this makes a deep copy of the XML data

sub _add_set {
    my $self = shift;
    $self->throw('No linkset!') unless my $ls = shift;
    my $dbfrom = $ls->{DbFrom};
    $self->dbfrom($dbfrom);
    my $query_id = $ls->{IdList}->{Id};
    $self->query_id($query_id);
    for my $ls_db (@{ $ls->{LinkSetDb} }) {
        my $dbto = $ls_db->{DbTo};
        push @{ $self->{'_databases'}}, $dbto;
        my $linkname = $ls_db->{LinkName};
        if ( $ls_db->{Info} || $ls->{ERROR} || !($ls_db->{Link})) {
            my $err_msg = $ls_db->{Info} || $ls->{ERROR} || 'No Links!';
            my $ids = (ref($query_id) =~ /array/i) ?
                            join q(,), @{$query_id}: $query_id;
            $self->warn("ELink Error for $dbto and ids $ids: $err_msg");
            next;
        }
        my @ids;
        for my $id_ref (@{ $ls_db->{Link} } ) {
            my $id = $id_ref->{Id};
            my $score = $id_ref->{Score} ? $id_ref->{Score} : undef;
            push @ids, $id;
            if ($score) {
                $self->{'_scores'}->{$id} = $score;
                if (!($self->{'_has_scores'})) {
                    $self->{'_has_scores'} = $dbto;
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
}

sub dbfrom {
    my $self = shift;
    return $self->{'_dbfrom'} = shift if @_;
    return $self->{'_dbfrom'};
}

sub query_id {
    my $self = shift;
    return $self->{'_query_id'} = shift if @_;
    return $self->{'_query_id'};
}

sub command {
    my $self = shift;
    return $self->{'_command'} = shift if @_;
    return $self->{'_command'};
}

sub has_scores {
    my $self = shift;
    return $self->{'_has_scores'};
}

# remodel these to be something along the lines of next_cookie, but use closure

sub get_LinkIds_by_db {
    my $self = shift;
    my $db = shift if @_;
    $self->throw("Must use database to access IDs") if !$db;
    for my $linkset (@{ $self->{'_linksetdb'}}) {
        my $dbto = $linkset->{DbTo};
        if ($dbto == $db) {
            return @{ $linkset->{DbTo} } if wantarray;
            return $linkset->{DbTo};
        }
    }
    $self->warn("Couldn't find ids for database $db");
}

sub get_databases {
    my $self = shift;
    return @{ $self->{'_databases'} };
}

sub get_score {
    my $self = shift;
    my $id = shift if @_;
    $self->warn("No scores!") if !$self->has_scores;
    $self->warn("Must use ID to access scores") if !$id;
    if ($self->{'_scores'}->{$id} ) {
        return $self->{'_scores'}->{$id};
    }
}

sub get_score_hash {
    my $self = shift;
    $self->warn("No scores!") if !$self->has_scores;
    return %{ $self->{'_scores'} };
}

1;
__END__