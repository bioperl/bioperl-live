#
# BioPerl module Bio::DB::Biblio::soap.pm
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Martin Senger <senger@ebi.ac.uk>
# For copyright and disclaimer see below.

# POD documentation - main docs before the code

=head1 NAME

Bio::DB::Biblio::soap - A SOAP-based access to a bibliographic query service

=head1 SYNOPSIS

Do not use this object directly, it is recommended to access it and use
it through the I<Bio::Biblio> module:

  use Bio::Biblio;
  my $biblio = Bio::Biblio->new (-access => 'soap');

=head1 DESCRIPTION

This object contains the real implementation of a Bibliographic Query
Service as defined in L<Bio::DB::BiblioI> - using a SOAP protocol
to access a WebService (a remote server) that represents a
bibliographic repository.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

  bioperl-l@bioperl.org                  - General discussion
  http://bioperl.org/wiki/Mailing_lists  - About the mailing lists

=head2 Support 

Please direct usage questions or support issues to the mailing list:

I<bioperl-l@bioperl.org>

rather than to the module maintainer directly. Many experienced and 
reponsive experts will be able look at the problem and quickly 
address it. Please include a thorough description of the problem 
with code and data examples if at all possible.

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution. Bug reports can be submitted via the
web:

  https://redmine.open-bio.org/projects/bioperl/

=head1 AUTHOR

Martin Senger (martin.senger@gmail.com)

=head1 COPYRIGHT

Copyright (c) 2002 European Bioinformatics Institute. All Rights Reserved.

This module is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.

=head1 DISCLAIMER

This software is provided "as is" without warranty of any kind.

=head1 BUGS AND LIMITATIONS

=over

=item *

Methods returning a boolean value (I<has_next>, I<exists> and
I<contains>) can be used only with SOAP::Lite version 0.52 and newer
(probably due to a bug in the older SOAP::Lite).

=item *

It does not use WSDL.

=item *

More testing and debugging needed to ensure that returned citations
are properly transferred even if they contain foreign characters.

=back

=head1 APPENDIX

The main documentation details are to be found in
L<Bio::DB::BiblioI>.

Here is the rest of the object methods.  Internal methods are preceded
with an underscore _.

=cut


# Let the code begin...


package Bio::DB::Biblio::soap;
use vars qw($DEFAULT_SERVICE $DEFAULT_NAMESPACE);
use strict;

use SOAP::Lite
    on_fault => sub {
	my $soap = shift;
	my $res = shift;
	my $msg =
	    ref $res ? "--- SOAP FAULT ---\n" . $res->faultcode . " " . $res->faultstring
		     : "--- TRANSPORT ERROR ---\n" . $soap->transport->status . "\n$res\n";
        Bio::DB::Biblio::soap->throw ( -text => $msg );
    }
;

use base qw(Bio::Biblio);

BEGIN {
    # where to go...
    $DEFAULT_SERVICE = 'http://www.ebi.ac.uk/openbqs/services/MedlineSRS';

    # ...and what to find there
    
    ## TODO: This namespace is no longer valid (check for deprecation or update)
    $DEFAULT_NAMESPACE = 'http://industry.ebi.ac.uk/openBQS';
}

# -----------------------------------------------------------------------------

=head2 _initialize

 Usage   : my $obj = Bio::Biblio->new(-access => 'soap' ...);
           (_initialize is internally called from this constructor)
 Returns : nothing interesting
 Args    : This module recognises and uses following arguments:

             -namespace => 'urn'
               The namespace used by the WebService that is being
               accessed. It is a string which guarantees its world-wide
               uniqueness - therefore it often has a style of a URL -
               but it does not mean that such pseudo-URL really exists.

               ## TODO: This namespace is no longer valid (check for deprecation
               ## or update)
               
               Default is 'http://industry.ebi.ac.uk/openBQS'.

             -destroy_on_exit => '0'
                Default value is '1' which means that all Bio::Biblio
                objects - when being finalised - will send a request
                to the remote WebService to forget the query collections
                they represent.

                If you change it to '0' make sure that you know the
                query collection identification - otherwise you will
                not be able to re-established connection with it.
                This can be done by calling method get_collection_id.

              -collection_id => '...'
                It defines what query collection will this object work
                with. Use this argument when you know a collection ID
                of an existing query collection and when you wish to
                re-established connection with it.

                By default, the collection IDs are set automatically
                by the query methods - they return Bio::Biblio objects
                already having a collection ID.

                A missing or undefined collection ID means that the
                object represents the whole bibliographic repository
                (which again means that some methods, like get_all,
                will be probably refused).

              -soap => a SOAP::Lite object
                Usually all Bio::Biblio objects share an instance of
                the underlying SOAP::Lite module. But you are free
                to have more - perhaps with different characteristics.

                See the code for attributes of the default SOAP::Lite
                object.

              -httpproxy => 'http://server:port'
                 In addition to the 'location' parameter, you may need
                 to specify also a location/URL of a HTTP proxy server
                 (if your site requires one).

	   Additionally, the main module Bio::Biblio recognises
	   also:
             -access => '...'
             -location => '...'

It populates calling object with the given arguments, and then - for
some attributes and only if they are not yet populated - it assigns
some default values.

This is an actual new() method (except for the real object creation
and its blessing which is done in the parent class Bio::Root::Root in
method _create_object).

Note that this method is called always as an I<object> method (never as
a I<class> method) - and that the object who calls this method may
already be partly initiated (from Bio::Biblio::new method); so if you
need to do some tricks with the 'class invocation' you need to change
Bio::Biblio::new method, not this one.

=cut

sub _initialize {
    my ($self, @args) = @_;
    
    # make a hashtable from @args
    my %param = @args;
    @param { map { lc $_ } keys %param } = values %param; # lowercase keys

    # copy all @args into this object (overwriting what may already be
    # there) - changing '-key' into '_key'
    my $new_key;
    foreach my $key (keys %param) {
	($new_key = $key) =~ s/^-/_/;
	$self->{ $new_key } = $param { $key };
    }

    # finally add default values for those keys who have default value
    # and who are not yet in the object
    $self->{'_location'} = $DEFAULT_SERVICE unless $self->{'_location'};
    $self->{'_namespace'} = $DEFAULT_NAMESPACE unless $self->{'_namespace'};
    $self->{'_destroy_on_exit'} = 1 unless defined $self->{'_destroy_on_exit'};
    unless ($self->{'_soap'}) {
	if (defined $self->{'_httpproxy'}) {
	    $self->{'_soap'} = SOAP::Lite
	                          -> uri ($self->{'_namespace'})
		                  -> proxy ($self->{'_location'},
				            proxy => ['http' => $self->{'_httpproxy'}]);
	} else {
	    $self->{'_soap'} = SOAP::Lite
	                          -> uri ($self->{'_namespace'})
				  -> proxy ($self->{'_location'});
	}
#	$self->{'_soap'}->soapversion (1.2);
    }
}

# -----------------------------------------------------------------------------

#
# objects representing query collections are being destroyed if they
# have attribute '_destroy_on_exit' set to true - which is a default
# value
#
sub DESTROY {
    my $self = shift;
    my $soap = $self->{'_soap'};
    my $destroy = $self->{'_destroy_on_exit'};
    return unless $destroy;
    my $collection_id = $self->{'_collection_id'};
    return unless $collection_id;

    # ignore all errors here
    eval {
	$soap->destroy (SOAP::Data->type (string => $collection_id));
    }
}

#
# some methods must be called with an argument containing a collection
# ID; here we return a proper error message explaining it
#
sub _no_id_msg {
    my $self = shift;
    my $package = ref $self;
    my $method = (caller(1))[3];
    my $strip_method = $method;
    $strip_method =~ s/^$package\:\://;

    return <<"END_OF_MSG";
Method '$method' works only if its object has a query collection ID.
Perhaps you need to use:
\tBio::Biblio->new(-collection_id => '1234567')->$strip_method;
or to obtain a collection ID indirectly from a query method:
\tBio::Biblio->new->find ('keyword')->$strip_method;
END_OF_MSG
}
    
#
# some methods do not work with older SOAP::Lite version; here we
#return message explaining it
#
sub _old_version_msg {
    my $self = shift;
    my $method = (caller(1))[3];

    return <<"END_OF_MSG";
Method '$method' works only with SOAP::Lite
version 0.52 and newer (the problem is with returning a boolean value from the server).
END_OF_MSG
}

#
# some controlled vocabulary methods needs two parameters; here we
# return message explaining it
#
sub _two_params_msg {
    my $self = shift;
    my $method = (caller(1))[3];

    return <<"END_OF_MSG";
Method '$method' expects two parameters: vocabulary name and a value.
END_OF_MSG
}

#
# some controlled vocabulary methods needs a vocabulary name; here we
# return message explaining it
#
sub _missing_name_msg {
    my $self = shift;
    my $method = (caller(1))[3];

    return <<"END_OF_MSG";
Method '$method' expects vocabulary name as parameter.
END_OF_MSG
}

# 
# return a copy of a given array, with all its elements replaced
# with the SOAP-Data objects defining elements type as 'string'
#
sub _as_strings {
    my ($ref_input_array) = @_;
    my (@result) = map { SOAP::Data->new (type => 'string', value => $_) } @$ref_input_array;
    return \@result;
}
    
# ---------------------------------------------------------------------
#
#   Here are the methods implementing Bio::DB::BiblioI interface
#   (documentation is in Bio::DB::BiblioI)
#
# ---------------------------------------------------------------------

sub get_collection_id {
   my ($self) = @_;
   $self->{'_collection_id'};
}

sub get_count {
   my ($self) = @_;
   my $soap = $self->{'_soap'};
   my ($collection_id) = $self->{'_collection_id'};
   if ($collection_id) {
       $soap->getBibRefCountOfCollection (SOAP::Data->type (string => $collection_id))->result;
   } else {
       $soap->getBibRefCount->result;
   }
}

# try: 12368254 (it's a Bioperl article)
sub get_by_id {
   my ($self, $citation_id) = @_;
   $self->throw ("Citation ID is expected as a parameter of method 'get_by_id'.")
       unless $citation_id;
   my $soap = $self->{'_soap'};
   $soap->getById (SOAP::Data->type (string => $citation_id))->result;
}

sub find {
   my ($self, $keywords, $attrs) = @_;
   my (@keywords, @attrs);

   # $keywords can be a comma-delimited scalar or a reference to an array
   if ($keywords) {
       my $ref = ref $keywords;
       @keywords = split (/,/, $keywords) unless $ref;
       @keywords = @$keywords if $ref =~ /ARRAY/;
   }
   $self->throw ("No keywords given in 'find' method.\n")
       unless (@keywords);

   # ...and the same with $attrs
   if ($attrs) {
       my $ref = ref $attrs;
       @attrs = split (/,/, $attrs) unless $ref;
       @attrs = @$attrs if $ref =~ /ARRAY/;
   }

   my $soap = $self->{'_soap'};
   my $collection_id = $self->{'_collection_id'};
   my $new_id;
   if ($collection_id) {
       if (@attrs) {
	   $new_id = $soap->reFindInAttrs (SOAP::Data->name ('arg0')->type (string => $collection_id),
				           SOAP::Data->name ('arg1')->value (&_as_strings (\@keywords)),
				           SOAP::Data->name ('arg2')->value (&_as_strings (\@attrs)))
	       ->result;
       } else {
	   $new_id = $soap->reFind (SOAP::Data->name ('arg0')->type (string => $collection_id),
				    SOAP::Data->name ('arg1')->value (&_as_strings (\@keywords)))
	       ->result;
       }
   } else {
       if (@attrs) {
	   $new_id = $soap->findInAttrs (SOAP::Data->name ('arg0')->value (&_as_strings (\@keywords)),
				         SOAP::Data->name ('arg1')->value (&_as_strings (\@attrs)))
	       ->result;
       } else {
	   $new_id = $soap->find (SOAP::Data->name ('arg0')->value (&_as_strings (\@keywords)))
	       ->result;
       }
   }

   # clone itself but change the collection ID to a new one
   return $self->new (-collection_id        => $new_id,
		      -parent_collection_id => $collection_id);
}

sub get_all_ids {
   my ($self) = @_;
   my $soap = $self->{'_soap'};
   my ($collection_id) = $self->{'_collection_id'};
   $self->throw ($self->_no_id_msg) unless $collection_id;
   $soap->getAllIDs (SOAP::Data->type (string => $collection_id))->result;
}

sub get_all {
   my ($self) = @_;
   my $soap = $self->{'_soap'};
   my ($collection_id) = $self->{'_collection_id'};
   $self->throw ($self->_no_id_msg) unless $collection_id;
   $soap->getAllBibRefs (SOAP::Data->type (string => $collection_id))->result;
}

sub has_next {
   my ($self) = @_;
   my $soap = $self->{'_soap'};
   my ($collection_id) = $self->{'_collection_id'};
   $self->throw ($self->_no_id_msg) unless $collection_id;
   $self->throw ($self->_old_version_msg) if $SOAP::Lite::VERSION lt '0.52';
   $soap->hasNext (SOAP::Data->type (string => $collection_id))->result;
}

sub get_next {
   my ($self) = @_;
   my $soap = $self->{'_soap'};
   my ($collection_id) = $self->{'_collection_id'};
   $self->throw ($self->_no_id_msg) unless $collection_id;
   $soap->getNext (SOAP::Data->type (string => $collection_id))->result;
}

sub get_more {
   my ($self, $how_many) = @_;
   my $soap = $self->{'_soap'};
   my $collection_id = $self->{'_collection_id'};
   $self->throw ($self->_no_id_msg) unless $collection_id;

   unless (defined ($how_many) and $how_many =~ /^\d+$/) {
       $self->warn ("Method 'get_more' expects a numeric argument. Changing to 1.\n");
       $how_many = 1;
   }
   unless ($how_many > 0) {
       $self->warn ("Method 'get_more' expects a positive argument. Changing to 1.\n");
       $how_many = 1;
   }

   my $ra = $soap->getMore (SOAP::Data->type (string => $collection_id),
			    SOAP::Data->type (int    => $how_many))->result;
   $self->{'_collection_id'} = shift @{ $ra };
   $ra;
}

sub reset_retrieval {
   my ($self) = @_;
   my $soap = $self->{'_soap'};
   my ($collection_id) = $self->{'_collection_id'};
   $self->throw ($self->_no_id_msg) unless $collection_id;
   $self->{'_collection_id'} = $soap->resetRetrieval (SOAP::Data->type (string => $collection_id))->result;
}

sub exists {
   my ($self) = @_;
   my $soap = $self->{'_soap'};
   my ($collection_id) = $self->{'_collection_id'};
   $self->throw ($self->_no_id_msg) unless $collection_id;
   $self->throw ($self->_old_version_msg) if $SOAP::Lite::VERSION lt '0.52';
   $soap->exists (SOAP::Data->type (string => $collection_id))->result;
}

sub destroy {
   my ($self) = @_;
   my $soap = $self->{'_soap'};
   my ($collection_id) = $self->{'_collection_id'};
   $self->throw ($self->_no_id_msg) unless $collection_id;
   $soap->destroy (SOAP::Data->type (string => $collection_id));
}

sub get_vocabulary_names {
   my ($self) = @_;
   my $soap = $self->{'_soap'};
   $soap->getAllVocabularyNames->result;
}

sub contains {
   my ($self, $vocabulary_name, $value) = @_;
   my $soap = $self->{'_soap'};
   $self->throw ($self->_old_version_msg) if $SOAP::Lite::VERSION lt '0.52';
   $self->throw ($self->_two_params_msg)
       unless defined $vocabulary_name and defined $value;
   $soap->contains (SOAP::Data->type (string => $vocabulary_name),
		    SOAP::Data->type (string => $value))->result;
}

sub get_entry_description {
   my ($self, $vocabulary_name, $value) = @_;
   my $soap = $self->{'_soap'};
   $self->throw ($self->_two_params_msg)
       unless defined $vocabulary_name and defined $value;
   $soap->getEntryDescription (SOAP::Data->type (string => $vocabulary_name),
			         SOAP::Data->type (string => $value))->result;
}

sub get_all_values {
   my ($self, $vocabulary_name) = @_;
   my $soap = $self->{'_soap'};
   $self->throw ($self->_missing_name_msg)
       unless defined $vocabulary_name;
   $soap->getAllValues (SOAP::Data->type (string => $vocabulary_name))->result;
}

sub get_all_entries {
   my ($self, $vocabulary_name) = @_;
   my $soap = $self->{'_soap'};
   $self->throw ($self->_missing_name_msg)
       unless defined $vocabulary_name;
   $soap->getAllEntries (SOAP::Data->type (string => $vocabulary_name))->result;
}

=head2 VERSION and Revision

 Usage   : print $Bio::DB::Biblio::soap::VERSION;
           print $Bio::DB::Biblio::soap::Revision;

=cut

=head2 Defaults

 Usage   : print $Bio::DB::Biblio::soap::DEFAULT_SERVICE;
           print $Bio::DB::Biblio::soap::DEFAULT_NAMESPACE;

=cut

1;
__END__
