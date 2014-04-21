
# POD documentation - main docs before the code

=head1 NAME

Bio::DB::Failover - A Bio::DB::RandomAccessI compliant class which
wraps a prioritized list of DBs

=head1 SYNOPSIS

    $failover = Bio::DB::Failover->new();

    $failover->add_database($db);

    # fail over Bio::DB::RandomAccessI.pm

    # this will check each database in priority, returning when
    # the first one succeeds

    $seq = $failover->get_Seq_by_id($id);

=head1 DESCRIPTION

This module provides fail over access to a set of Bio::DB::RandomAccessI
objects.

=head1 CONTACT

Ewan Birney E<lt>birney@ebi.ac.ukE<gt> originally wrote this class.

=head2 Support 

Please direct usage questions or support issues to the mailing list:

I<bioperl-l@bioperl.org>

rather than to the module maintainer directly. Many experienced and 
reponsive experts will be able look at the problem and quickly 
address it. Please include a thorough description of the problem 
with code and data examples if at all possible.

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution. Bug reports can be submitted via the
web:

  https://github.com/bioperl/bioperl-live/issues

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::DB::Failover;

use strict;

use base qw(Bio::Root::Root Bio::DB::RandomAccessI);

sub new {
    my ($class,@args) = @_;

    my $self = $class->SUPER::new(@args);

    $self->{'_database'} = [];
    return $self;
}

=head2 add_database

 Title   : add_database
 Usage   : add_database(%db)
 Function: Adds a database to the Failover object
 Returns : Count of number of databases
 Args    : Array of db resources
 Throws  : Not a RandomAccessI exception

=cut

sub add_database {
	my ($self,@db) = @_;
	for my $db ( @db ) {
		if ( !ref $db || !$db->isa('Bio::DB::RandomAccessI') ) {
			$self->throw("Database object $db is a not a Bio::DB::RandomAccessI");
			next;
		}

		push(@{$self->{'_database'}},$db);
	}
	scalar @{$self->{'_database'}};
}


=head2 get_Seq_by_id

 Title   : get_Seq_by_id
 Usage   : $seq = $db->get_Seq_by_id('ROA1_HUMAN')
 Function: Gets a Bio::Seq object by its name
 Returns : a Bio::Seq object
 Args    : the id (as a string) of a sequence
 Throws  : "no id" exception

=cut

sub get_Seq_by_id {
	my ($self,$id) = @_;

	if( !defined $id ) {
		$self->throw("no id is given!");
	}

	foreach my $db ( @{$self->{'_database'}} ) {
		my $seq;

		eval {
			$seq = $db->get_Seq_by_id($id);
		};
		$self->warn($@) if $@;
		if ( defined $seq ) {
			return $seq;
		} else {
			$self->warn("No sequence retrieved by database " . ref($db));
		}
	}

	return;
}

=head2 get_Seq_by_acc

 Title   : get_Seq_by_acc
 Usage   : $seq = $db->get_Seq_by_acc('X77802');
 Function: Gets a Bio::Seq object by accession number
 Returns : A Bio::Seq object
 Args    : accession number (as a string)
 Throws  : "no id" exception

=cut

sub get_Seq_by_acc {
	my ($self,$id) = @_;

	if( !defined $id ) {
		$self->throw("no id is given!");
	}

	foreach my $db ( @{$self->{'_database'}} ) {
		my $seq;
		eval {
			$seq = $db->get_Seq_by_acc($id);
		};
		$self->warn($@) if $@;
		if ( defined $seq ) {
			return $seq;
		} else {
			$self->warn("No sequence retrieved by database " . ref($db));
		}
	}
	return;
}

=head2 get_Seq_by_version

 Title   : get_Seq_by_version
 Usage   : $seq = $db->get_Seq_by_acc('X77802.2');
 Function: Gets a Bio::Seq object by versioned accession number
 Returns : A Bio::Seq object
 Args    : accession number (as a string)
 Throws  : "acc does not exist" exception

=cut

sub get_Seq_by_version {
	my ($self,$id) = @_;

	if( !defined $id ) {
		$self->throw("no acc is given!");
	}

	foreach my $db ( @{$self->{'_database'}} ) {
		my $seq;
		eval {
			$seq = $db->get_Seq_by_version($id);
		};
		$self->warn($@) if $@;
		if ( defined $seq ) {
			return $seq;
		} else {
			$self->warn("No sequence retrieved by database " . ref($db));
		}
	}
	return;
}

## End of Package

1;

__END__
