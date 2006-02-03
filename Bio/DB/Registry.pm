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

This module provides access to the Open Bio Database Access (OBDA)
scheme, which provides a single cross-language and cross-platform 
specification of how to get to databases. These databases may be 
accessible through the Web, they may be BioSQL databases, or
they may be local, indexed flatfile databases.

If the user or system administrator has not installed the default init 
file, seqdatabase.ini, in /etc/bioinformatics or ${HOME}/.bioinformatics 
then creating the first Registry object copies the default settings from 
the net. The Registry object will attempt to store these settings in
${HOME}/.bioinformatics/seqdatabase.ini.

Users can specify one or more custom locations for the init file by 
setting $OBDA_SEARCH_PATH to those directories, where multiple 
directories should be separated by ';'.

Please see the OBDA Access HOWTO for more information
(http://bioperl.open-bio.org/wiki/HOWTO:OBDA).

=head1 CONTACT

Ewan Birney originally wrote this class.

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution. Bug reports can be submitted via the
web:

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

my %implement = ('biocorba'         => 'Bio::CorbaClient::SeqDB',
		 'flat'             => 'Bio::DB::Flat',
		 'biosql'           => 'Bio::DB::BioSQL::BioDatabaseAdaptor',
		 'biofetch'         => 'Bio::DB::BioFetch' );

my $fallbackRegistryURL = 'http://www.open-bio.org/registry/seqdatabase.ini';

sub new {
    my ($class,@args) = shift;
    my $self = $class->SUPER::new(@args);
    # open files in order
    $self->{'_dbs'} = {};
    $self->_load_registry();
    return $self;
}

=head2 _load_registry

 Title   : _load_registry
 Usage   :
 Function: Looks for seqdatabase.ini files in the expected locations and
           in the directories specified by $OBDA_SEARCH_PATH. If no files
           are found it downloads a default file from www.open-bio.org
 Returns : nothing
 Args    : none

=cut

sub _load_registry {
   my ($self) = @_;
   my $home = "";
	$home = $ENV{"HOME"} if defined $ENV{"HOME"};
	eval {$home = (getpwuid($>))[7];} unless $home;
	if ($@) {
		warn "This Perl doesn't implement function getpwuid(). Skipping...\n"
	}
	my @ini_files = _get_ini_files($home);

   unless (@ini_files) {
      $self->warn("No seqdatabase.ini file found in ~/.bioinformatics/\nnor in /etc/bioinformatics/ nor in directory specified by\n$OBDA_SEARCH_PATH. Using web to get database registry from\n$fallbackRegistryURL");

      # Last gasp. Try to use HTTPget module to retrieve the registry from
      # the web...
      my $f = Bio::Root::HTTPget::getFH($fallbackRegistryURL);

      # store the default registry file
      mkdir "$home/.bioinformatics" unless -e "$home/.bioinformatics";
      open(F,">$home/.bioinformatics/seqdatabase.ini");
      print F while (<$f>);
      close F;

      $self->warn("Stored the default registry configuration in\n" .
		  "$home/.bioinformatics/seqdatabase.ini");

      push @ini_files,"$home/.bioinformatics/seqdatabase.ini";
   }

   my ($db,$hash) = ();
   foreach my $file (@ini_files) {
      open FH,"$file";
      while( <FH> ) {
			if (/^VERSION=([\d\.]+)/) {
				if ($1 > $OBDA_SPEC_VERSION or !$1) {
					$self->throw("Do not know about this version [$1] > $OBDA_SPEC_VERSION");
					last;
				}
				next;
         }
			next if( /^#/ );
			next if( /^\s/ );
			if ( /^\[(\w+)\]/ ) {
				$db = $1;
				next;
			}
			my ($tag,$value) = split('=',$_);
			$value =~ s/\s//g;
			$tag =~ s/\s//g;
			$hash->{$db}->{"\L$tag"} = lc $value;
      }
   }

   foreach my $db( keys %{$hash} ) {
      if ( !exists $self->{'_dbs'}->{$db} ) {
			my $failover = Bio::DB::Failover->new();
			$self->{'_dbs'}->{$db} = $failover;
      }
      my $class;
      if (defined $implement{$hash->{$db}->{'protocol'}}) {
			$class = $implement{$hash->{$db}->{'protocol'}};
      } else {
			$self->warn("Registry does not support protocol " .
							$hash->{$db}->{'protocol'});
			next;
      }
      eval "require $class";
      if ($@) {
			$self->verbose && $self->warn("Couldn't load $class");
			next;
      } else {
	  eval {
	      my $randi = $class->new_from_registry( %{$hash->{$db}} );
	      $self->{'_dbs'}->{$db}->add_database($randi); };
	  if ($@) {
	      $self->warn("Couldn't call new_from_registry on [$class]\n$@");
	  }
      }
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
		$self->warn("No database with name $dbname in Registry");
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

sub services {
    my ($self) = @_;
    return () unless ( defined $self->{'_dbs'} &&
		       ref( $self->{'_dbs'} ) =~ /HASH/i);
    return keys %{$self->{'_dbs'}};
}

=head2 _get_ini_files

 Title   : _get_ini_files
 Usage   :
 Function: To find all the seqdatabase.ini files
 Returns : list of seqdatabase.ini paths
 Args    : $home

=cut

sub _get_ini_files {
   my $home = shift;
   my @ini_files = ();
   if ( $OBDA_SEARCH_PATH ) {
      foreach my $dir ( split /;/,$OBDA_SEARCH_PATH ) {
			my $file = $dir . "/" . "seqdatabase.ini";
			next unless -e $file;
			push @ini_files,$file;
      }
   }
   push @ini_files,"$home/.bioinformatics/seqdatabase.ini" 
     if ( $home && -e "$home/.bioinformatics/seqdatabase.ini" );
   push @ini_files,"/etc/bioinformatics/seqdatabase.ini"
     if ( $home && -e "/etc/bioinformatics/seqdatabase.ini" );
   @ini_files;
}

## End of Package

1;

__END__

