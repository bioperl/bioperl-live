# POD documentation - main docs before the code

# $Id$


=head1 NAME

Bio::DB::Registry - Access to the Open Bio Database Access registry scheme

=head1 SYNOPSIS

    use Bio::DB::Registry();

    $registry = new Bio::DB::Registry();

    @available_services = $registry->services;

    $db = $registry->get_database('embl');
    # $db is a Bio::DB::SeqI implementing class

    $seq = $db->get_Seq_by_acc("J02231");

=head1 DESCRIPTION

This module provides access to the Open Bio Database Access scheme,
which provides a cross language and cross platform specification of how
to get to databases.

If the user or system administrator has not installed the default init file,
creating the first Registry object copies the default settings from the net.

=head1 CONTACT

Ewan Birney originally wrote this class.

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution. Bug reports can be submitted via email
or the web:

    bioperl-bugs@bio.perl.org
    http://bugzilla.bioperl.org/

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...

package Bio::DB::Registry;

use vars qw(@ISA $OBDA_SPEC_VERSION $OBDA_SEARCH_PATH);
use strict;

use Bio::Root::Root;
@ISA = qw(Bio::Root::Root);
use Bio::DB::Failover;
use Bio::Root::HTTPget;

BEGIN {
    $OBDA_SPEC_VERSION = 1.0;
    if (defined $ENV{OBDA_SEARCH_PATH}) {
        $OBDA_SEARCH_PATH = $ENV{OBDA_SEARCH_PATH} || '';

    }
}

my %implement = (
		 'biocorba'         => 'Bio::CorbaClient::SeqDB',
		 'flat'             => 'Bio::DB::Flat',
		 'biosql'           => 'Bio::DB::BioSQL::BioDatabaseAdaptor',
		 'biofetch'         => 'Bio::DB::BioFetch'
		 );

my $fallbackRegistryURL = 'http://www.open-bio.org/registry/seqdatabase.ini';


sub new {
    my ($class,@args) = shift;
    my $self = $class->SUPER::new(@args);

    # open files in order
    $self->{'_dbs'} = {};
    $self->_load_registry();
    return $self;
}


sub _load_registry {
    my ($self) = @_;

    my $home = (getpwuid($>))[7];
    my $f;

    if ( $OBDA_SEARCH_PATH ) {
        foreach ( split /\+/, $OBDA_SEARCH_PATH ) {
            next unless -e $_;
            open(F,"$OBDA_SEARCH_PATH/seqdatabase.ini");
            $f = \*F;
            last;
        }
    }
    elsif( -e "$home/.bioinformatics/seqdatabase.ini" ) {
	open(F,"$home/.bioinformatics/seqdatabase.ini");
	$f = \*F;
    } elsif ( -e "/etc/bioinformatics/seqdatabase.ini" ) {
	open(F,"/etc/bioinformatics/seqdatabase.ini");
	$f = \*F;
    } else {
	# waiting for information
	$self->warn("No conf file found in ~/.bioinformatics/ \nor in /etc/.bioinformatics/.\n".
                    "Using web to get database registry from \n$fallbackRegistryURL");

	# Last gasp. Try to use HTTPget module to retrieve the registry from
        # the web...

	$f = Bio::Root::HTTPget::getFH($fallbackRegistryURL);

        # store the default registry file
        mkdir "$home/.bioinformatics" unless -e "$home/.bioinformatics";
	open(F,">$home/.bioinformatics/seqdatabase.ini");
        print F while (<$f>);
        close F;

	$self->warn("Stored the default registry configuration into:\n".
                    "  $home/.bioinformatics/seqdatabase.ini");

	open(F,"$home/.bioinformatics/seqdatabase.ini");
	$f = \*F;

    }

    while( <$f> ) {
	/^VERSION=([\d\.]+)/;
        $self->throw("Do not know about this version [$1] > $OBDA_SPEC_VERSION")
            if $1 > $OBDA_SPEC_VERSION or !$1;
        last;
    }

    while( <$f> ) {
      if( /^#/ ) {
	  next;
       }
	if( /^\s/ ) {
	  next;
	}

      if( /\[(\w+)\]/ )  {
	my $db = $1;
	my $hash = {};
	while( <$f> ) {
	  chomp();
	  /^#/ && next;
	    /^$/ && last;
	  my ($tag,$value) = split('=',$_);
	  $value =~ s/\s//g;
	  $tag =~ s/\s//g;
	  $hash->{"\L$tag"} = lc $value;
	}

	if( !exists $self->{'_dbs'}->{$db} ) {
	  my $failover = Bio::DB::Failover->new();
	  $self->{'_dbs'}->{$db}=$failover;
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
	    $self->{'_dbs'}->{$db}->add_database($randi); };
	  if ($@) {
	    $self->warn("Couldn't call new_from_registry on [$class]\n$@");
	  }
	}
	next; # back to main loop
      }
      $self->warn("Uninterpretable line in registry, $_");
    }
}

=head2 get_database

 Title   : get_database
 Usage   : my $db = $registry->get_database($dbname);
 Function: Retrieve a Database object which implements Bio::DB::SeqI interface
 Returns : Bio::DB::SeqI object
 Args    : string describing the name of the database

=cut

sub get_database {
    my ($self,$dbname) = @_;

    $dbname = lc $dbname;
    if( !defined $dbname ) {
	$self->warn("must get_database with a database name");
	return undef;
    }
    if( !exists $self->{'_dbs'}->{$dbname} ) {
	$self->warn("No database in with $dbname in registry");
	return undef;
    }
    return $self->{'_dbs'}->{$dbname};
}

=head2 services

 Title   : services
 Usage   : my @available = $registry->services();
 Function: returns list of possible services 
 Returns : list of strings
 Args    : none


=cut

sub services{ 
    my ($self) = @_;
    return () unless ( defined $self->{'_dbs'} &&
		       ref( $self->{'_dbs'} ) =~ /HASH/i);
    return keys %{$self->{'_dbs'}};
}


## End of Package

1;

__END__

