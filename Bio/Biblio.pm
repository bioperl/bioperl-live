# $Id$
#
# BioPerl module Bio::Biblio
#
# Cared for by Martin Senger <senger@ebi.ac.uk>
# For copyright and disclaimer see below.
#

# POD documentation - main docs before the code

=head1 NAME

Bio::Biblio - A Bibliographic Query Service module

=head1 SYNOPSIS

  use Bio::Biblio;
  my $biblio = new Bio::Biblio;

  print $biblio->find ('perl')->get_count . "\n";

  my $collection = $biblio->find ('brazma', 'authors');
  while ( $collection->has_next ) {
      print $collection->get_next;
  }

Here are some one-liners:

  perl -MBio::Biblio -e 'print new Bio::Biblio->get_by_id ("94033980")' 
  perl -MBio::Biblio \
       -e 'print join ("\n", @{ new Bio::Biblio->find ("brazma")->get_all_ids })' 
  perl -MBio::Biblio \
       -e 'print new Bio::Biblio->find ("Java")->find ("perl")->get_count'

The C<new> method can get parameters, for example:

  my $biblio = Bio::Biblio 
    (-access          => 'soap',
     -location        => 'http://industry.ebi.ac.uk/soap/openBQS',
     -destroy_on_exit => '0');

=head1 DESCRIPTION

This is a class whose instances can access bibliographic
repositories. It allows to query a bibliographic database (such as
MEDLINE) and then to retrieve resulting citations from it. The
citations are returned in an XML format which is native to the
repository but there are also supporting modules for converting them
into Perl objects.

The detailed descriptions of all query and retrieval methods are in
L<Bio::DB::BiblioI> (an interface). All those methods should be
called on instances of this (Bio::Biblio) module.

The module complies (with some simplifications) with the specification
described in the B<OpenBQS> project. Its home page is at
I<http://industry.ebi.ac.uk/openBQS>. There are also links to
available servers providing access to the bibliographic repositories
(namely to I<MEDLINE>).

The module also gives an access to a set of controlled vocabularies
and their values. It allows to introspect bibliographic repositories
and to find what citation resource types (such as journal and book
articles, patents or technical reports) are provided, and what
attributes they have, eventually what attribute values are allowed.

=head1 OVERVIEW OF CLASSES AND PACKAGES

=over

=item B<Bio::Biblio>

This is the main class to be used by the end users. It
loads a real implementation for a particular access protocol according
to the argument I<-access>. At the time of writing this documentation
there is only one available access module implementing all query and
retrieval methods:

   -access => soap

This module implements all methods defined in the interface
I<Bio::DB::BiblioI> (see L<Bio::DB::BiblioI>) by delegating
calls to a loaded low-level module (e.g. see
L<Bio::DB::Biblio::soap>).

Note that there is also another module (and perhaps more) which does
not use SOAP protocol and do not implement all query methods -
nevertheless it has retrieval methods and it can be used in the same
way:

   -access => biofetch


=item Bio::DB::BiblioI

This is an interface defining all methods that can be called on
I<Bio::Biblio> instances.

=item Bio::DB::Biblio::soap

This is a real implementation of all methods defined in
Bio::DB::BiblioI using SOAP protocol (calling a WebService
based on SOAP). This class should not be instantiated directly (use
I<Bio::Biblio> instead). See L<Bio::DB::BiblioI> for details.

=item Bio::Biblio::IO

This module instantiates and uses a converter of the citations read by
any of the access methods mentioned above. See L<Bio::Biblio::IO> for
details.

=item Bio::Biblio::IO::medlinexml and Bio::Biblio::IO::medline2ref

A converter of MEDLINE citations in XML into Perl objects.

=item Bio::Biblio::IO::pubmedxml and Bio::Biblio::IO::pubmed2ref

A converter of PUBMED citations in XML into Perl objects.

=back

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

  bioperl-l@bioperl.org              - General discussion
  http://bioperl.org/MailList.shtml  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution. Bug reports can be submitted via
email or the web:

  bioperl-bugs@bioperl.org
  http://bugzilla.bioperl.org/

=head1 AUTHOR

Martin Senger (senger@ebi.ac.uk)

=head1 COPYRIGHT

Copyright (c) 2002 European Bioinformatics Institute. All Rights Reserved.

This module is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.

=head1 DISCLAIMER

This software is provided "as is" without warranty of any kind.

=head1 SEE ALSO

=over

=item *

OpenBQS home page: http://industry.ebi.ac.uk/openBQS

=item *

Comments to the Perl client: http://industry.ebi.ac.uk/openBQS/Client_perl.html

=back

=head1 APPENDIX

The main documentation details are to be found in
L<Bio::DB::BiblioI>.

Here is the rest of the object methods.  Internal methods are preceded
with an underscore _.

=cut


# Let the code begin...


package Bio::Biblio;
use vars qw(@ISA $VERSION $Revision);
use strict;

use Bio::Root::Root;
use Bio::DB::BiblioI;

@ISA = qw(Bio::Root::Root Bio::DB::BiblioI);


BEGIN { 
    $VERSION = do { my @r = (q$Revision$ =~ /\d+/g); sprintf "%d.%-02d", @r };
    $Revision = q$Id$;
}

# -----------------------------------------------------------------------------

=head2 new

 Usage   : my $obj = new Bio::Biblio (@args);
 Returns : Bio::Biblio object on success, or undef on failure
 Args    : This module recognizes and uses:

             -access => 'soap'
               It indicates what lower-level module to load.
               Default is 'soap'.

             -location => 'http://...'
                It says where to find a bibliographic query service.
                The format and contents of this argument is dependent
                on the '-access' argument.

                For 'soap' access it is a URL of a WebService.
                Default is http://industry.ebi.ac.uk/soap/openBQS

             -data_source synonym for -access # NG 03-06-09
               Only works for 'ncbi_eutils'

           Other arguments can be given here but they are
           recognized by the lower-level module
           (e.g. see Bio::DB::Biblio::soap).

It builds, populates and returns a new I<Bio::Biblio> object. This is
how it is seen from the outside. But in fact, it builds, populates and
returns a more specific lower-level object, for example
I<Bio::DB::Biblio::soap> object - which one it is depends on the
parameter I<-access>.

The real initialization is done in the method I<_initialize> of the
lower-level object.

This method can also be used for I<cloning> an existing object and
changing or adding new attributes to it in the same time. This is,
however, not particulary useful for the casual users of this module,
because the query methods (see L<Bio::DB::BiblioI>) themselves
already return cloned objects with more refined query
collections. Anyway this is how the cloning can be done:

  use Bio::Biblio;
  my $biblio = new Bio::Biblio;

  # this will create a new object which will NOT send a 'destroy'
  # message to the remote server when its life ends
  my $clone = $biblio->new (-destroy-on-exit => '0'); 

=cut

sub new {
    my ($caller,@args) = @_;
    my $class = ref($caller) || $caller;
  
    # if $caller is an object, or if it is an underlying
    # 'real-work-doing' class (e.g. Bio::DB::Biblio::soap) then
    # we want to call SUPER to create and bless an object

    if ($class =~ /Bio::DB::Biblio::(\S+)/) {
	my ($self) = $class->SUPER::new (@args);

	# now the $self is an empty object - we will populate it from
	# the $caller - if $caller is an object

	if (ref ($caller)) {
	    %{ $self } = %{ $caller };
	}

	# and finally add values from '@args' into the newly created
	# object (the values will overwrite the values copied above)

	$self->_initialize (@args);
	return $self;

    # this is called only the first time when somebody calls: 'new
    # Bio::Biblio (...)', and it actually loads a 'real-work-doing'
    # module and call this new() method again (unless the loaded
    # module has its own new() method)

    } else { 
	my %param = @args;
	@param { map { lc $_ } keys %param } = values %param; # lowercase keys
	
	# NG 03-06-15. Changed to handle ncbi_eutils parameter combinations
	my $access;
	$access= 'ncbi_eutils' if
	  (lc $param {'-access'} eq 'ncbi_eutils' || lc $param {'-format'} eq 'ncbi_eutils') && 
	    !($param {'-file'} || $param {'-filehandle'});
	$access or $access='ncbi_eutils_file' if 
	  lc $param {'-format'} eq 'ncbi_eutils' && ($param {'-file'} || $param {'-filehandle'});
	$access or $access='alzforum_file' if 
	  lc $param {'-format'} eq 'alzforum' && ($param {'-file'} || $param {'-filehandle'});
	$access or $access =
	  $param {'-access'} || 
	    $class->_guess_access ( $param {'-location'} ) ||
	      'soap';
	$access = "\L$access";	# normalize capitalization to lower case

	# load module with the real implementation - as defined in $access
	return undef unless (&_load_access_module ($access));

	# this will call this same method new() - but rather its the
	# upper (object) branche
	return "Bio::DB::Biblio::$access"->new (@args);
    }
}

# NG 03-06-15
# Code adapted from Bio::Biblio::IO::medlinexml by Martin Senger
sub newFh {
  my $class = shift;
  my $self = $class->new(@_);
  $self->throw("Sorry, newFh only defined for ncbi_eutils or alzforum access module") unless $self->isa('Bio::DB::Biblio::ncbi_eutils') || $self->isa('Bio::DB::Biblio::alzforum');
  return $self->fh;
}
#sub fh {
#  my $self = shift;
#  my $class = ref($self) || $self;
#  my $s = Symbol::gensym;
#  tie $$s,$class,$self;
#  return $s;
#}
#sub TIEHANDLE {
#  my ($class,$val) = @_;
#  return bless {'biblio' => $val}, $class;
#}
#sub READLINE {$_[0]->throw("Not implemented");}

# -----------------------------------------------------------------------------

=head2 _load_access_module

 Usage   : $class->_load_access_module ($access)
 Returns : 1 on success, undef on failure
 Args    : 'access' should contain the last part of the
           name of a module who does the real implementation

It does (in run-time) a similar thing as

   require Bio::DB::Biblio::$access

It prints an error on STDERR if it fails to find and load the module
(for example, because of the compilation errors in the module).

=cut

sub _load_access_module {
  my ($access) = @_;
  my ($module, $load, $m);

  $module = "_<Bio/DB/Biblio/$access.pm";
  $load = "Bio/DB/Biblio/$access.pm";

  return 1 if $main::{$module};
  eval {
    require $load;
  };

  if ( $@ ) {
    Bio::Root::Root->throw (<<END);
$load: $access cannot be found or loaded
Exception $@
For more information about the Biblio system please see the Bio::Biblio docs.
END
  ;
    return;
  }
  return 1;
}

# -----------------------------------------------------------------------------

=head2 _guess_access

 Usage   : $class->_guess_access ($location)
 Returns : string with a guessed access protocol (e.g. 'soap')
 Args    : 'location' defines where to find a bibliographic service
           in a protocol-dependent manner (e.g. for SOAP it is
           a URL of a bibliographic WebService)

It makes an expert guess what kind of access/transport protocol should
be used based on the I<location> of the service (e.g. if the
I<location> looks like an IOR then the access protocol is probably
CORBA).

=cut

# this is kept here for the future when more access protocols
# (e.g. CORBA) may be available for accessing bibliographic query
# services

sub _guess_access {
#   my ($class, $location) = @_;
   return 'soap';
}

=head2 VERSION and Revision

 Usage   : print $Bio::Biblio::VERSION;
           print $Bio::Biblio::Revision;

=cut

1;
__END__
