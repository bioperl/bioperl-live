# POD documentation - main docs before the code

# $Id$


=head1 NAME

Bio::DB::Registry - Access to the Open Bio Database Access registry scheme

=head1 SYNOPSIS

    $registry = new Bio::DB::Registry();

    @available_services = $registry->services;

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
use Bio::DB::Failover;
use Bio::Root::HTTPget;

my %implement = (
		 'biocorba'      => 'Bio::CorbaClient::SeqDB',
		 'index-berkeleydb' => 'Bio::DB::Flat',
                 'index-flat'       => 'Bio::DB::Flat::OBDAIndex',
		 'biosql' => 'Bio::DB::SQL::BioDatabaseAdaptor',
		 'biofetch' => 'Bio::DB::BioFetch'
		 );

my $fallbackRegistryURL = 'http://www.open-bio.org/registry/seqdatabase.ini';


sub new {
    my ($class,@args) = shift;

    my $self = Bio::Root::Root->new();
    bless $self,$class;
    my ($verbose)= $self->_rearrange([qw(VERBOSE)],@args);
    
    # open files in order
    $self->_load_registry();
    $self->verbose($verbose);
    return $self;
}


sub _load_registry {
    my ($self) = @_;

    my $home = (getpwuid($>))[7];
    my $f;
    if( -e "$home/.bioinformatics/seqdatabase.ini" ) {
	open(F,"$home/.bioinformatics/seqdatabase.ini");
	$f = \*F;
    } elsif ( -e "/etc/bioinformatics/seqdatabase.ini" ) {
	open(F,"$home/.bioinformatics/seqdatabase.ini");
	$f = \*F;
    } else {
	# waiting for information
	$self->warn("No conf file found in ~/.bioinformatics/ \nor in /etc/.bioinformatics/ using web to get database registry from \n$fallbackRegistryURL\n");

	# Last gasp. Try to use HTTPget module to retrieve the registry from
        # the web...

	$f = Bio::Root::HTTPget::getFH($fallbackRegistryURL);

    }

    while( <$f> ) {
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
	    while( <$f> ) {
		chomp();
		/^#/ && next;
		/^$/ && last;
		my ($tag,$value) = split('=',$_);
		$value =~ s/\s//g;
		$tag =~ s/\s//g;
		$hash->{$tag} = $value;
	    }

	    if( !exists $self->{$db} ) {
		my $failover = Bio::DB::Failover->new();
		$self->{$db}=$failover;
	    }

	    my $class;
	    if (defined $implement{$hash->{'protocol'}}) {
		$class = $implement{$hash->{'protocol'}};
	    }
	    else {
		$self->warn("Registry does not support protocol ".$hash->{'protocol'});
		next;
	    }
	    eval "require $class";
	    
	    if ($@) {
		$self->verbose && $self->warn("Couldn't load $class");
		next;
	    }
	    
	    else {
		eval {
		    my $randi = $class->new_from_registry(%$hash);
		    $self->{$db}->add_database($randi);		};
		if ($@) {
		    $self->verbose && $self->warn("Couldn't call new_from_registry on [$class]\n$@");
		}
	    }
	    next; # back to main loop
	}

	$self->warn("Uninterpretable line in registry, $_");
    }
}


=head2 verbose

 Title   : verbose
 Usage   : 
 Function: get/set for verbose paramater
 Returns : integer
 Args    : none


=cut

sub verbose {
   my ($self,$value) = @_;

   if ( defined $value ) {
       $self->{'_verbose'} = $value;
   }
   return $self->{'_verbose'};
}

	
sub get_database {
    my ($self,$dbname) = @_;

    if( !defined $dbname ) {
	$self->throw("must get_database with a database name");
    }

    if( !exists $self->{$dbname} ) {
	$self->throw("No database in with $dbname in registry");
    }

    return $self->{$dbname};

}

## End of Package

1;

__END__

