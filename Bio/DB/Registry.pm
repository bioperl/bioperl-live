# POD documentation - main docs before the code

# $Id$


=head1 NAME

Bio::DB::Registry - Access to the Open Bio Database Access registry scheme

=head1 SYNOPSIS

    $registry = Bio::DB::Registry->new();

    $db = $registry->get_database('embl');

    # $db is a Bio::DB::SeqI implementing class

=head1 DESCRIPTION

This module provides access to the Open Bio Database Access scheme,
which provides a cross language and cross platform specification of how
to get to databases. 


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

package Bio::DB::Registry;

use vars qw(@ISA);
use strict;

use Bio::Root::Root;
@ISA = qw(Bio::Root::Root);

my %implement = (
		 'bsane-corba'      => 'Bio::CorbaClient::SeqDB',
		 'index-berkeleydb' => 'XYZ',
		 'bioperldb' => 'Bio::DB::SQL::BioSeqDatabaseFetcher');

sub new {
    my ($class) = shift;

    my $self = Bio::Root::Root->new();
    bless $self,$class;

    # open files in order

    $self->_load_registry();

    return $self;
}


sub _load_registry {
    my ($self) = @_;

    my $home = (getpwuid($>))[7];
    
    if( -e "$home/.bioinformatics/seqdatabase.ini" ) {
	open(F,"$home/.bioinformatics/seqdatabase.ini");
    } elsif ( -e "/etc/bioinformatics/seqdatabase.ini" ) {
	open(F,"$home/.bioinformatics/seqdatabase.ini");
    } else {
	# waiting for information
	$self->throw("Oooops. We haven't implemented web fall back position yet");
    }

    while( <F> ) {
	if( /^#/ ) {
	    next;
	}
	if( /^\s/ ) {
	    next;
	}

	if( /\[(\w+)\]/ )  {
	    my $db;
	    $db = $1;
	    my $hash = {};
	    while( <F> ) {
		/^#/ && next;
		/^\s/ && last;
		my ($tag,$value) = split;
		$hash->{$tag} = $value;
	    }
	    if( !exists $self->{$db} ) {
		$self->{$db} = {};
		$self->{$db}->{'services'} = [];
	    }
	    
	    push(@{$self->{$db}->{'services'}},$hash);
	    next; # back to main loop
	}

	$self->warn("Uninterpretable line in registry, $_");
    }
}

	
sub get_database {
    my ($self,$dbname) = @_;

    if( !defined $dbname ) {
	$self->throw("must get_database with a database name");
    }

    if( !exists $self->{$dbname} ) {
	$self->throw("No database in with $dbname in registry");
    }

    if( exists $self->{$dbname}->{'active'} ) {
	return $self->{$dbname}->{'active'};
    }

    if( scalar(@{$self->{$dbname}->{'services'}}) > 1 ) {
	$self->throw("Sorry ... have not implemented multiple database failbacks yet. Complain to Ewan as he wrote this!");
    }

    my ($config) = @{$self->{$dbname}->{'services'}};
    my $class;
    unless ($class = $implement{$config->{'protocol'}}) {
	$self->throw("Registry does not support protocol ".$config->{'protocol'});
    }
    eval "require $class";
    $config->{'biodbname'}=$dbname;
    my $db = $class->new(%$config);

    $self->{$db}->{'active'} = $db;


    return $db;
}
    


## End of Package

1;

__END__

