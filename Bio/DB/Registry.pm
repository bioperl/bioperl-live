#
# POD documentation - main docs before the code

=head1 NAME

Bio::DB::Registry - Access to the Open Bio Database Access registry scheme

=head1 SYNOPSIS

    use Bio::DB::Registry();

    $registry = Bio::DB::Registry->new();

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
the www.open-bio.org. The Registry object will attempt to store these 
settings in a new file, ${HOME}/.bioinformatics/seqdatabase.ini.

Users can specify one or more custom locations for the init file by 
setting $OBDA_SEARCH_PATH to those directories, where multiple 
directories should be separated by ';'.

Please see the OBDA Access HOWTO for more information
(L<http://bioperl.open-bio.org/wiki/HOWTO:OBDA>).

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

package Bio::DB::Registry;

use vars qw($OBDA_SPEC_VERSION $OBDA_SEARCH_PATH
			   $HOME $PRIVATE_DIR $PUBLIC_DIR $REGISTRY 
			   $FALLBACK_REGISTRY);
use strict;

use Bio::DB::Failover;
use Bio::Root::HTTPget;
use base qw(Bio::Root::Root);

BEGIN {
   $OBDA_SPEC_VERSION = 1.0;
	$HOME = $ENV{HOME} if (defined $ENV{HOME});
	if (defined $ENV{OBDA_SEARCH_PATH}) {
		$OBDA_SEARCH_PATH = $ENV{OBDA_SEARCH_PATH} || '';
   }
}

my %implement = ('flat'     => 'Bio::DB::Flat',
					   'biosql'   => 'Bio::DB::BioSQL::OBDA',
					   'biofetch' => 'Bio::DB::BioFetch'
					   # 'biocorba' => 'Bio::CorbaClient::SeqDB',
					   );

$FALLBACK_REGISTRY = 'http://www.open-bio.org/registry/seqdatabase.ini';
$PRIVATE_DIR = '.bioinformatics';
$PUBLIC_DIR = '/etc/bioinformatics';
$REGISTRY = 'seqdatabase.ini';

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
           are found download a default file from www.open-bio.org
 Returns : nothing
 Args    : none

=cut

sub _load_registry {
   my $self = shift;
	eval { $HOME = (getpwuid($>))[7]; } unless $HOME;
	if ($@) {
		$self->warn("This Perl doesn't implement function getpwuid(), no \$HOME");
	}
	my @ini_files = $self->_get_ini_files();

	@ini_files = $self->_make_private_registry() unless (@ini_files);

   my ($db,$hash) = ();
   for my $file (@ini_files) {
      open my $FH, '<', $file or $self->throw("Could not read file '$file': $!");
      while( <$FH> ) {
			if (/^VERSION=([\d\.]+)/) {
				if ($1 > $OBDA_SPEC_VERSION or !$1) {
					$self->throw("Do not know about this version [$1] > $OBDA_SPEC_VERSION");
					last;
				}
				next;
         }
			next if( /^#/ );
			next if( /^\s/ );
			if ( /^\[(\S+)\]/ ) {
				$db = $1;
				next;
			}
			my ($tag,$value) = split('=',$_);
			$value =~ s/\s//g;
			$tag =~ s/\s//g;
			$hash->{$db}->{"\L$tag"} = $value;
      }
   }

   for my $db ( keys %{$hash} ) {
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
			$self->warn("Couldn't load $class");
			next;
      } else {
			eval {
				my $randi = $class->new_from_registry( %{$hash->{$db}} );
				$self->{'_dbs'}->{$db}->add_database($randi); 
			};
			if ($@) {
				$self->warn("Couldn't call new_from_registry() on [$class]\n$@");
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
		return;
	}
	if( !exists $self->{'_dbs'}->{$dbname} ) {
		$self->warn("No database with name $dbname in Registry");
		return;
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
 Usage   : my @files = $self->_get_ini_files
 Function: To find all the seqdatabase.ini files
 Returns : list of seqdatabase.ini paths
 Args    : None

=cut

sub _get_ini_files {
	my $self = shift;
   my @ini_files = ();
   if ( $OBDA_SEARCH_PATH ) {
      foreach my $dir ( split /;/, $OBDA_SEARCH_PATH ) {
			my $file = $dir . "/" . $REGISTRY;
			next unless -e $file;
			push @ini_files,$file;
      }
   }
   push @ini_files,"$HOME/$PRIVATE_DIR/$REGISTRY" 
     if ( $HOME && -e "$HOME/$PRIVATE_DIR/$REGISTRY" );
   push @ini_files, "$PUBLIC_DIR/$REGISTRY"
     if ( -e "$PUBLIC_DIR/$REGISTRY" );
   @ini_files;
}

=head2 _make_private_registry

 Title   : _make_private_registry
 Usage   :
 Function: Make private registry in file in $HOME
 Returns : Path to private registry file
 Args    : None

=cut

sub _make_private_registry {
	my $self = shift;
   my @ini_file;

	my $nor_in = $OBDA_SEARCH_PATH ? 
	  "nor in directory specified by\n$OBDA_SEARCH_PATH" : 
	  "and environment variable OBDA_SEARCH_PATH wasn't set";

	$self->warn("No $REGISTRY file found in $HOME/$PRIVATE_DIR/\n" . 
					"nor in $PUBLIC_DIR $nor_in.\n" .
					"Using web to get registry from\n$FALLBACK_REGISTRY");

	# Last gasp. Try to use HTTPget module to retrieve the registry from
	# the web...
	my $f = Bio::Root::HTTPget::getFH($FALLBACK_REGISTRY);

	# store the default registry file
	eval {
		mkdir "$HOME/$PRIVATE_DIR" unless -e "$HOME/$PRIVATE_DIR";
	};
	$self->throw("Could not make directory $HOME/$PRIVATE_DIR, " .
					 "no $REGISTRY file available") if $@;

	open my $F, '>', "$HOME/$PRIVATE_DIR/$REGISTRY"
	  or $self->throw("Could not write file '$HOME/$PRIVATE_DIR/$REGISTRY': $!");
	print $F while (<$F>);
	close $F;

	$self->warn("Stored $REGISTRY file in $HOME/$PRIVATE_DIR");

	push @ini_file,"$HOME/$PRIVATE_DIR/$REGISTRY";
	@ini_file;
}

1;

__END__
