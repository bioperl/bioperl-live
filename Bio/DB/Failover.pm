# POD documentation - main docs before the code

# $Id$


=head1 NAME

Bio::DB::Failover - A Bio::DB::RandomAccessI compliant class which wraps a priority list of DBs

=head1 SYNOPSIS

    $failover = Bio::DB::Failover->new();

    $failover->add_database($db);

    # fail over Bio::DB::RandomAccessI.pm

    # this will check each database in priority, returning when
    # the first one succeeds

    $seq = $failover->get_Seq_by_id($id);

=head1 DESCRIPTION

This module provides fail over access to a set of Bio::DB::RandomAccessI objects


=head1 CONTACT

Ewan Birney originally wrote this class.

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution. Bug reports can be submitted via email
or the web:

    bioperl-bugs@bio.perl.org                   
    http://bio.perl.org/bioperl-bugs/           

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...

package Bio::DB::Failover;

use vars qw(@ISA);
use strict;

use Bio::Root::Root;
use Bio::DB::RandomAccessI;
@ISA = qw(Bio::DB::RandomAccessI Bio::Root::Root);


sub new {
    my ($class,@args) = @_;

    my $self = Bio::Root::Root->new();
    bless $self,$class;

    $self->{'database'} = [];

    return $class;
}


sub add_database {
    my ($self,@db) = @_;


    foreach my $db ( @db ) {
	if( !ref $db || !$db->isa('Bio::DB::RandomAccessI') ) {
	    $self->throw("Database objects $db is a not a Bio::DB::RandomAccessI");
	    next;
	}

	push(@{$self->{'database'}},$db);
    }
    
}


sub get_Seq_by_id {
    my ($self,$id) = @_;

    if( !defined $id ) {
	$self->throw("no id is given!");
    }

    foreach my $db ( @{$self->{'database'}} ) {
	my $seq;
	eval {
	    $seq = $db->get_Seq_by_id($db);
	};
	if( defined $seq ) {
	    return $seq;
	}
    }

    return undef;
}



sub get_Seq_by_acc {
    my ($self,$id) = @_;

    if( !defined $id ) {
	$self->throw("no id is given!");
    }

    foreach my $db ( @{$self->{'database'}} ) {
	my $seq;
	eval {
	    $seq = $db->get_Seq_by_acc($db);
	};
	if( defined $seq ) {
	    return $seq;
	}
    }

    return undef;
}


## End of Package

1;

__END__

