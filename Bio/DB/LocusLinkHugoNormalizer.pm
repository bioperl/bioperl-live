package Bio::DB::LocusLinkHugoNormalizer;

# $Id$

=head1 NAME

Bio::SeqFeature::LocusLinkHugoNormalizer - A name normalizer for
L<Bio::SeqFeatureI> objects that uses a LocusLink database backend to
normalize SeqFeature display_names() to the Hugo name.
L<LocusLinkHugoNormalizer> can take a SeqFeatureI object and set its
display_name to the Hugo name for that feature.

=head1 SYNOPSIS

  my $normalizer =
    Bio::DB::LocusLinkHugoNormalizer->new( @args );
  my $iterator = $seq_feature_collection->features( '-iterator' => 1 );
  while( $iterator->has_more_features() ) {
    my $seq_feature = $iterator->next_feature();
    $normalizer->normalize( $seq_feature );
    print "The feature's normalized name is ".$seq_feature->display_name();
  }

=head1 DESCRIPTION

  A name normalizer for L<Bio::SeqFeatureI> objects that uses a
  LocusLink database backend to normalize SeqFeature display_names() to
  the Hugo name.  L<LocusLinkHugoNormalizer> can take a SeqFeatureI
  object and set its display_name to the Hugo name for that feature.  A
  name normalizer for L<Bio::SeqFeatureI> objects.  Of course, to do
  this, the SeqFeatureI object must already have its unique_id() or
  display_name() set to some value known by the normalizer to be an
  alias for some standard name.

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

Paul Edlefsen E<lt>paul@systemsbiology.orgE<gt>.

Copyright (c) 2003 Institute for Systems Biology

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.  See DISCLAIMER.txt for
disclaimers of warranty.

=head1 CONTRIBUTORS

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

# Let the code begin...

use strict;
use vars qw( $VERSION @ISA );

$VERSION = '0.01';
use Bio::Root::Root;
use Bio::SeqFeature::NameNormalizerI;
@ISA = qw( Bio::Root::Root Bio::SeqFeature::NameNormalizerI );

use Bio::DB::GFF::Util::Rearrange; # for rearrange()
use DBI;

## TODO: REMOVE?
use constant DEFAULT_DB_SERVER   => 'mssql';
use constant DEFAULT_DB_DATABASE => 'JDRFdev';

use constant DEFAULT_DB_DSN	 => ("DBI:Sybase:server=".DEFAULT_DB_SERVER.";database=".DEFAULT_DB_DATABASE);
use constant DEFAULT_DB_USER     => 'seqDataLoader';
use constant DEFAULT_DB_PASS     => 'DL444';


=head2 new

 Title   : new
 Usage   : $normalizer = Bio::DB::LocusLinkHugoNormalizer->new( @args )
 Function: Instantiates a new name normalizer that connects to a LocusLink db.
 Returns : a new Bio::DB::LocusLinkHugoNormalizer
 Args    : -dsn => The name of the Sybase LocusLink Hugo database to connect to
           -username => The username to use to connect to the database
           -password => The password to use to connect to the database
 Status  : Public

=cut

sub new {
  my $class = shift;
  my $self = $class->SUPER::new( @_ );

  my ( $dsn, $username, $password, $other ) =
    rearrange([
               'dsn',
               [qw(USERNAME USER UNAME)],
               [qw(PASSWORD PASS PW)],
              ], @_
             );

  $dsn      |= DEFAULT_DB_DSN;
  $username |= DEFAULT_DB_USER;
  $password |= DEFAULT_DB_PASS;

  my @args;
  push( @args, $dsn );
  push( @args, $username );
  push( @args, $password );
  my $db_handle = DBI->connect( @args )
    or $class->throw( "new(): Failed to connect to $dsn: ".$DBI::errstr );

  $self->db_handle( $db_handle );

  return $self;
} # new(..)

sub db_handle {
  my $self = shift;
  my ( $new_handle ) = @_;
  my $current_value = $self->{ '_db_handle' };
  if( defined $new_handle ) {
    $self->{ '_db_handle' } = $new_handle;
  }
  return $current_value;
} # db_handle(..)

=head2 normalize

 Title   : normalize
 Usage   : $normalizer->normalize( $seq_feature );
 Function: Changes the given feature\'s display_name to the Hugo
           standard if a standard name can be found that corresponds
           to the feature\'s unique_id() or display_name() or string
           value, tested in that order.
 Returns : the given L<Bio::SeqFeatureI>
 Args    : a L<Bio::SeqFeatureI> object
 Status  : Public

=cut

sub normalize {
  my $self = shift;
  my ( $seq_feature ) = @_;

  ## TODO: REMOVE
  #warn "LocusLinkHugoNormalizer->normalize( $seq_feature )";

  my $dbh = $self->db_handle();

  my @name_possibilities = 
    ( ( $seq_feature->can( 'unique_id' ) ? $seq_feature->unique_id() : undef ),
      ( $seq_feature->can( 'display_name' ) ? $seq_feature->display_name() : undef ),
      "$seq_feature" );
  foreach my $name ( @name_possibilities ) {

    next unless defined( $name );

    $name = uc( $name );

    ## TODO: REMOVE
    #warn "Hey, got current name '$name'.";

    if( my $hugo_name = $self->{ '_hugo_names' }->{ $name } ) {
      ## TODO: REMOVE.  This is a hack.
      $seq_feature->add_tag_value( 'alias', $name );
      $seq_feature->display_name( $hugo_name );
      return $seq_feature;
    }

    my $q_name = $dbh->quote( $name );

    my $sth		= $dbh->prepare("
    	SELECT	        gnh.hugo_name
    	FROM            genename_hugoII gnh
    	WHERE	        gnh.field_name     = $q_name
    ");
    ####
    #    SELECT	        gi.hugo_name
    #    FROM            gene_information gi LEFT OUTER JOIN gene_refseq gr
    #       ON             ( gi.locus_link_id = gr.locus_link_id )
    #    WHERE	        gr.refseq_id     = $q_name
    ####
    #    SELECT	        gi.hugo_name
    #    FROM            gene_information gi LEFT OUTER JOIN gene_refseq gr
    #       ON             ( gi.locus_link_id = gr.locus_link_id )
    #    WHERE	        gr.refseq_id     = $q_name
    #                 OR gi.mim_id        = $q_name
    #                 OR gi.gene_cards_id = $q_name
    $sth->execute() || warn "Got SQL problem: ".$sth->errstr();

    ## TODO: Deal with the fact that there could be multiple of these.
    my $found_it;
    while (my ($hugo_name) = $sth->fetchrow_array() ) {
      if( defined $found_it ) {
        ## TODO: REMOVE
        warn "****!!!!***** Yo yo yo, there's multiple hugo names for this feature: $seq_feature, with this name: $name.  Earlier we got $found_it, and now we have $hugo_name.";
      }
      #warn "Dude! got $hugo_name";
      $found_it = $hugo_name;
      ## Save it for later.
      $self->{ '_hugo_names' }->{ $name } = $hugo_name;
      ## TODO: REMOVE.  This is a hack.
      $seq_feature->add_tag_value( 'alias', $name );
      $seq_feature->display_name( $hugo_name );


      ############ TODO: REMOVE.  TESTING. ###################
      #my $q_hugo_name = $dbh->quote( $hugo_name );
      #
      #my $sth		= $dbh->prepare("
      #    SELECT	        srh.location_string
      #    FROM            single_region_hugo srh
      #    WHERE	        srh.hugo_name     = $q_hugo_name
      #");
      #$sth->execute() || warn "Got SQL problem: ".$sth->errstr();
      #
      ### TODO: Deal with the fact that there could be multiple of these.
      #while (my ($location) = $sth->fetchrow_array() ) {
      #  warn "Dude! got $location";
      #}
    }

  }

  return $seq_feature;
} # normalize()

=head2 locate

 Title   : locate
 Usage   : $normalizer->locate( $feature_string );
 Function: If the given string corresponds to a known location, return that location.
 Returns : the location corresponding to the given string.
 Args    : A string like 'ccr5' or 'NM_000567'
 Status  : Public

=cut

sub locate {
  my $self = shift;
  my ( $feature_string ) = @_;

  $feature_string = uc( $feature_string );

  my $dbh = $self->db_handle();
  my $hugo_name;
  $hugo_name = $self->{ '_hugo_names' }->{ $feature_string };
  unless( defined $hugo_name ) {
    ## TODO: REMOVE
    #warn "Looking up feature string '$feature_string'";

    my $q_feature_string = $dbh->quote( $feature_string );
    my $sth		= $dbh->prepare("
    	SELECT	        gnh.hugo_name
    	FROM            genename_hugoII gnh
    	WHERE	        gnh.field_name     = $q_feature_string
    ");
    $sth->execute() || warn "Got SQL problem: ".$sth->errstr();
    ## TODO: Deal with the fact that there could be multiple of these.
    my $found_it;
    while( my ( $a_hugo_name ) = $sth->fetchrow_array() ) {
      if( defined $found_it ) {
        ## TODO: REMOVE
        warn "****!!!!***** Yo yo yo, there's multiple hugo names for this string: $feature_string.  Earlier we got $found_it, and now we have $a_hugo_name.";
      } else {
        ## TODO: REMOVE
        #warn "Found hugo name $a_hugo_name";
      }
      $hugo_name = $found_it = $a_hugo_name;
      ## Save it for later.
      $self->{ '_hugo_names' }->{ $feature_string } = $a_hugo_name;
    }
  }
  ## TODO: REMOVE?
  unless( defined $hugo_name ) {
    return undef;
  }
  ## TODO: REMOVE
  #warn "Looking up hugo name $hugo_name";
  if( defined( $hugo_name ) ) {
    my $q_hugo_name = $dbh->quote( uc( $hugo_name ) );
    
    my $sth		= $dbh->prepare("
        SELECT	        srh.location_string
        FROM            single_region_hugo srh
        WHERE	        srh.hugo_name     = $q_hugo_name
    ");
    $sth->execute() || warn "Got SQL problem: ".$sth->errstr();
    
    ## TODO: Deal with the fact that there could be multiple of these.
    my $location;
    while( my ( $a_location ) = $sth->fetchrow_array() ) {
      if( defined $location ) {
        ## TODO: REMOVE
        warn "****!!!!***** Yo yo yo, there's multiple locations for this string: $feature_string.  Earlier we got $location, and now we have $a_location.";
      }
      $location = $a_location;
    }
    if( defined $location ) {
      return $location;
    }
  }

  ## If we got to here then there's a hugo name for the
  ## feature_string, but there's no corresponding location.
  return 'none';
} # locate()

sub DESTROY {
  my $self = shift;
  ## We have to remember to kill the handle.
  my $db_handle = $self->db_handle();
  if( defined $db_handle ) {
    $db_handle->disconnect();
  }
} # DESTROY()

1;

__END__
