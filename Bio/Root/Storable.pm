package Bio::Root::Storable;
use strict;
use Bio::Root::IO;
use Data::Dumper qw( Dumper );
use File::Spec;
use base qw(Bio::Root::Root);

# ABSTRACT: object serialisation methods
# AUTHOR:   Will Spooner <whs@sanger.ac.uk>
# OWNER:    Will Spooner
# LICENSE:  Perl_5

=head1 SYNOPSIS

  my $storable = Bio::Root::Storable->new();

  # Store/retrieve using class retriever
  my $token     = $storable->store();
  my $storable2 = Bio::Root::Storable->retrieve( $token );

  # Store/retrieve using object retriever
  my $storable2 = $storable->new_retrievable();
  $storable2->retrieve();


=head1 DESCRIPTION

Generic module that allows objects to be safely stored/retrieved from
disk.  Can be inhereted by any BioPerl object. As it will not usually
be the first class in the inheretence list, _initialise_storable()
should be called during object instantiation.

Object storage is recursive; If the object being stored contains other
storable objects, these will be stored separately, and replaced by a
skeleton object in the parent heirarchy. When the parent is later
retrieved, its children remain in the skeleton state until explicitly
retrieved by the parent. This lazy-retrieve approach has obvious
memory efficiency benefits for certain applications.


By default, objects are stored in binary format (using the Perl
Storable module). Earlier versions of Perl5 do not include Storable as
a core module. If this is the case, ASCII object storage (using the
Perl Data::Dumper module) is used instead.

ASCII storage can be enabled by default by setting the value of
$Bio::Root::Storable::BINARY to false.

=cut

use vars qw( $BINARY );

BEGIN{
    if( eval "require Storable" ){
        Storable->import( 'freeze', 'thaw' );
        $BINARY = 1;
    }
}

#----------------------------------------------------------------------

=head2 new

  Arg [1]   : -workdir  => filesystem path,
              -template => tmpfile template,
              -suffix   => tmpfile suffix,
  Function  : Builds a new Bio::Root::Storable inhereting object
  Returntype: Bio::Root::Storable inhereting object
  Exceptions:
  Caller    :
  Example   : $storable = Bio::Root::Storable->new()

=cut

sub new {
    my ($caller, @args) = @_;
    my $self = $caller->SUPER::new(@args);
    $self->_initialise_storable;
    return $self;
}

#----------------------------------------------------------------------

=head2 _initialise_storable

  Arg [1]   : See 'new' method
  Function  : Initialises storable-specific attributes
  Returntype: boolean
  Exceptions:
  Caller    :
  Example   :

=cut

sub _initialise_storable {
    my $self = shift;
    my( $workdir, $template, $suffix ) =
        $self->_rearrange([qw(WORKDIR TEMPLATE SUFFIX)], @_ );
    $workdir  && $self->workdir ( $workdir );
    $template && $self->template( $template );
    $suffix   && $self->suffix  ( $suffix   );
    return 1;
}



#----------------------------------------------------------------------

=head2 statefile

  Arg [1]   : string (optional)
  Function  : Accessor for the file to write state into.
              Should not normaly use as a setter - let Root::IO
              do this for you.
  Returntype: string
  Exceptions:
  Caller    : Bio::Root::Storable->store
  Example   : my $statefile = $obj->statefile();

=cut

sub statefile{
    my $key = '_statefile';
    my $self  = shift;

    if( @_ ){ $self->{$key} = shift }

    if( ! $self->{$key} ){ # Create a new statefile
        my $workdir  = $self->workdir;
        my $template = $self->template;
        my $suffix   = $self->suffix;

        # TODO: add cleanup and unlink methods. For now, we'll keep the
        # statefile hanging around.
        my @args = ( CLEANUP=>0, UNLINK=>0 );
        if( $template ){ push( @args, 'TEMPLATE' => $template )};
        if( $workdir  ){ push( @args, 'DIR'      => $workdir  )};
        if( $suffix   ){ push( @args, 'SUFFIX'   => $suffix   )};
        my( $fh, $file ) = Bio::Root::IO->new->tempfile( @args );
        # If filehandle is not stored, don't leave it open
        $fh->close;

        $self->{$key} = $file;
    }

    return $self->{$key};
}

#----------------------------------------------------------------------

=head2 workdir

  Arg [1]   : string (optional) (TODO - convert to array for x-platform)
  Function  : Accessor for the statefile directory. Defaults to File::Spec->tmpdir
  Returntype: string
  Exceptions:
  Caller    :
  Example   : $obj->workdir('/tmp/foo');

=cut

sub workdir {
    my $key = '_workdir';
    my $self = shift;
    if( @_ ){
        my $caller = join( ', ', (caller(0))[1..2] );
        $self->{$key} && $self->debug("Overwriting workdir: probably bad!");
        $self->{$key} = shift
    }
    #$self->{$key} ||= $Bio::Root::IO::TEMPDIR;
    $self->{$key} ||= File::Spec->tmpdir();
    return $self->{$key};
}

#----------------------------------------------------------------------

=head2 template

  Arg [1]   : string (optional)
  Function  : Accessor for the statefile template. Defaults to XXXXXXXX
  Returntype: string
  Exceptions:
  Caller    :
  Example   : $obj->workdir('RES_XXXXXXXX');

=cut

sub template {
    my $key = '_template';
    my $self = shift;
    if( @_ ){ $self->{$key} = shift }
    $self->{$key} ||= 'XXXXXXXX';
    return $self->{$key};
}

#----------------------------------------------------------------------

=head2 suffix

  Arg [1]   : string (optional)
  Function  : Accessor for the statefile template.
  Returntype: string
  Exceptions:
  Caller    :
  Example   : $obj->suffix('.state');

=cut

sub suffix {
    my $key = '_suffix';
    my $self = shift;
    if( @_ ){ $self->{$key} = shift }
    return $self->{$key};
}

#----------------------------------------------------------------------

=head2 new_retrievable

  Arg [1]   : Same as for 'new'
  Function  : Similar to store, except returns a 'skeleton' of the calling
              object, rather than the statefile.
              The skeleton can be repopulated by calling 'retrieve'. This
              will be a clone of the original object.
  Returntype: Bio::Root::Storable inhereting object
  Exceptions:
  Caller    :
  Example   : my $skel = $obj->new_retrievable(); # skeleton
              $skel->retrieve();                  # clone

=cut

sub new_retrievable{
    my $self = shift;
    my @args = @_;

    $self->_initialise_storable( @args );

    if( $self->retrievable ){ return $self->clone } # Clone retrievable
    return bless( { _statefile   => $self->store(@args),
                    _workdir     => $self->workdir,
                    _suffix      => $self->suffix,
                    _template    => $self->template,
                    _retrievable => 1 },
                 ref( $self ) );
}

#----------------------------------------------------------------------

=head2 retrievable

  Arg [1]   : none
  Function  : Reports whether the object is in 'skeleton' state, and the
              'retrieve' method can be called.
  Returntype: boolean
  Exceptions:
  Caller    :
  Example   : if( $obj->retrievable ){ $obj->retrieve }

=cut

sub retrievable {
    my $self = shift;
    if( @_ ){ $self->{_retrievable} = shift }
    return $self->{_retrievable};
}

#----------------------------------------------------------------------

=head2 token

  Arg [1]   : None
  Function  : Accessor for token attribute
  Returntype: string. Whatever retrieve needs to retrieve.
              This base implementation returns the statefile
  Exceptions:
  Caller    :
  Example   : my $token = $obj->token();

=cut

sub token{
    my $self = shift;
    return $self->statefile;
}


#----------------------------------------------------------------------

=head2 store

  Arg [1]   : none
  Function  : Saves a serialised representation of the object structure
              to disk. Returns the name of the file that the object was
              saved to.
  Returntype: string

  Exceptions:
  Caller    :
  Example   : my $token = $obj->store();

=cut

sub store{
    my $self = shift;
    my $statefile = $self->statefile;
    my $store_obj = $self->serialise;
    my $io = Bio::Root::IO->new( ">$statefile" );
    $io->_print( $store_obj );
    $self->debug( "STORING $self to $statefile\n" );
    # If filehandle is not stored, don't leave it open
    $io->close;
    return $statefile;
}

#----------------------------------------------------------------------

=head2 serialise

  Arg [1]   : none
  Function  : Prepares the the serialised representation of the object.
              Object attribute names starting with '__' are skipped.
              This is useful for those that do not serialise too well
              (e.g. filehandles).
              Attributes are examined for other storable objects. If these
              are found they are serialised separately using 'new_retrievable'
  Returntype: string
  Exceptions:
  Caller    :
  Example   : my $serialised = $obj->serialise();

=cut

sub serialise{
    my $self = shift;

    # Create a new object of same class that is going to be serialised
    my $store_obj = bless( {}, ref( $self ) );

    my %retargs = ( -workdir =>$self->workdir,
                    -suffix  =>$self->suffix,
                    -template=>$self->template );
    # Assume that other storable bio objects held by this object are
    # only 1-deep.

    foreach my $key( keys( %$self ) ){
        if( $key =~ /^__/ ){ next } # Ignore keys starting with '__'
        my $value = $self->{$key};

        # Scalar value
        if( ! ref( $value ) ){
            $store_obj->{$key} = $value;
        }

        # Bio::Root::Storable obj: save placeholder
        elsif( ref($value) =~ /^Bio::/ and $value->isa('Bio::Root::Storable') ){
            # Bio::Root::Storable
            $store_obj->{$key} = $value->new_retrievable( %retargs );
            next;
        }

        # Arrayref value. Look for Bio::Root::Storable objs
        elsif( ref( $value ) eq 'ARRAY' ){
            my @ary;
            foreach my $val( @$value ){
                if( ref($val) =~ /^Bio::/ and $val->isa('Bio::Root::Storable') ){
                    push(  @ary, $val->new_retrievable( %retargs ) );
                }
                else{ push(  @ary, $val ) }
            }
            $store_obj->{$key} = \@ary;
        }

        # Hashref value. Look for Bio::Root::Storable objs
        elsif( ref( $value ) eq 'HASH' ){
            my %hash;
            foreach my $k2( keys %$value ){
                my $val = $value->{$k2};
                if( ref($val) =~ /^Bio::/ and $val->isa('Bio::Root::Storable') ){
                    $hash{$k2} = $val->new_retrievable( %retargs );
                }
                else{ $hash{$k2} = $val }
            }
            $store_obj->{$key} = \%hash;
        }

        # Unknown, just add to the store object regardless
        else{ $store_obj->{$key} = $value }
    }
    $store_obj->retrievable(0); # Once deserialised, obj not retrievable
    return $self->_freeze( $store_obj );
}


#----------------------------------------------------------------------

=head2 retrieve

  Arg [1]   : string; filesystem location of the state file to be retrieved
  Function  : Retrieves a stored object from disk.
              Note that the retrieved object will be blessed into its original
              class, and not the
  Returntype: Bio::Root::Storable inhereting object
  Exceptions:
  Caller    :
  Example   : my $obj = Bio::Root::Storable->retrieve( $token );

=cut

sub retrieve{
    my( $caller, $statefile ) = @_;

    my $self = {};
    my $class = ref( $caller ) || $caller;

    # Is this a call on a retrievable object?
    if (    ref( $caller )
        and $caller->retrievable
        ){
        $self = $caller;
        $statefile = $self->statefile;
    }
    bless( $self, $class );

    # Recover serialised object
    if( ! -f $statefile ){
        $self->throw( "Token $statefile is not found" );
    }
    my $io = Bio::Root::IO->new( $statefile );
    local $/ = undef;
    my $state_str = $io->_readline('-raw'=>1);
    # If filehandle is not stored, don't leave it open
    $io->close;

    # Dynamic-load modules required by stored object
    my $stored_obj;
    my $success;
    for( my $i=0; $i<10; $i++ ){
        eval{ $stored_obj = $self->_thaw( $state_str ) };
        if( ! $@ ){
            $success = 1;
            last;
        }
        my $package;
        if( $@ =~ /Cannot restore overloading(.*)/i ){
            my $postmatch = $1; #'
            if( $postmatch =~ /\(package +([\w\:]+)\)/ ) {
                $package = $1;
            }
        }
        if( $package ){
            eval "require $package";
            $self->throw($@) if $@;
        }
        else{ $self->throw($@) }
    }
    if( ! $success ){ $self->throw("maximum number of requires exceeded" ) }

    if( ! ref( $stored_obj ) ){
        $self->throw( "Token $statefile returned no data" );
    }
    map{ $self->{$_} = $stored_obj->{$_} } keys %$stored_obj; # Copy hasheys
    $self->retrievable(0);

    # Maintain class of stored obj
    return $self;
}

#----------------------------------------------------------------------


=head2 clone

  Arg [1]   : none
  Function  : Returns a clone of the calling object
  Returntype: Bio::Root::Storable inhereting object
  Exceptions:
  Caller    :
  Example   : my $clone = $obj->clone();

=cut

sub clone {
    my $self = shift;
    my $frozen = $self->_freeze( $self );
    return $self->_thaw( $frozen );
}



#----------------------------------------------------------------------

=head2 remove

  Arg [1]   : none
  Function  : Clears the stored object from disk
  Returntype: boolean
  Exceptions:
  Caller    :
  Example   : $obj->remove();

=cut

sub remove {
    my $self = shift;
    if( -e $self->statefile ){
        unlink( $self->statefile );
    }
    return 1;
}

#----------------------------------------------------------------------

=head2 _freeze

  Arg [1]   : variable
  Function  : Converts whatever is in the the arg into a string.
              Uses either Storable::freeze or Data::Dumper::Dump
              depending on the value of $Bio::Root::BINARY
  Returntype:
  Exceptions:
  Caller    :
  Example   :

=cut

sub _freeze {
    my $self = shift;
    my $data = shift;
    if( $BINARY ){
        return freeze( $data );
    }
    else{
        $Data::Dumper::Purity = 1;
        return Data::Dumper->Dump( [\$data],["*code"] );
    }
}

#----------------------------------------------------------------------

=head2 _thaw

  Arg [1]   : string
  Function  : Converts the string into a perl 'whatever'.
              Uses either Storable::thaw or eval depending on the
              value of $Bio::Root::BINARY.
              Note; the string arg should have been created with
              the _freeze method, or strange things may occur!
  Returntype: variable
  Exceptions:
  Caller    :
  Example   :

=cut

sub _thaw {
    my $self = shift;
    my $data = shift;
    if( $BINARY ){
        return thaw( $data )
    }
    else{
        my $code;
        $code = eval( $data ) ;
        if($@) {
            $self->throw( "eval: $@" );
        }
        ref( $code ) eq 'REF'
            or $self->throw( "Serialised string was not a scalar ref" );
        return $$code;
    }
}

1;
