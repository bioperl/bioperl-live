# $Id$
#
# BioPerl module Bio::DB::Biblio::biofetch.pm
#
# Cared for by Heikki Lehvaslaiho <heikki@ebi.ac.uk>
# For copyright and disclaimer see below.

# POD documentation - main docs before the code

=head1 NAME

Bio::DB::Biblio::biofetch - A BioFetch-based access to a bibliographic 
  citation retrieval

=head1 SYNOPSIS

Do not use this object directly, only access it through the
I<Bio::Biblio> module:

  use Bio::Biblio;
  my $biblio = new Bio::Biblio (-access => 'biofetch');
  my $ref = $biblio->get_by_id('20063307'));

=head1 DESCRIPTION

This class uses BioFetch protocol based service to retrieve Medline
references by their ID.

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
  http://bioperl.org/bioperl-bugs/

=head1 AUTHOR

Heikki Lehvaslaiho (heikki@ebi.ac.uk)

=head1 COPYRIGHT

Copyright (c) 2002 European Bioinformatics Institute. All Rights Reserved.

This module is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.

=head1 DISCLAIMER

This software is provided "as is" without warranty of any kind.

=head1 BUGS AND LIMITATIONS

=over 1

=item *

Only method get_by_id() is supported.

=back

=head1 APPENDIX

The main documentation details are to be found in
L<Bio::DB::BiblioI>.

Here is the rest of the object methods.  Internal methods are preceded
with an underscore _.

=cut


# Let the code begin...


package Bio::DB::Biblio::biofetch;
use vars qw(@ISA $VERSION %HOSTS  %FORMATMAP  $DEFAULTFORMAT 
	    $Revision $DEFAULT_SERVICE $DEFAULT_NAMESPACE);
use strict;

use Bio::Biblio;
use Bio::DB::DBFetch;
use Bio::Biblio::IO;

@ISA = qw( Bio::DB::DBFetch Bio::Biblio);

BEGIN {

    # set the version for version checking
    $VERSION = do { my @r = (q$Revision$ =~ /\d+/g); sprintf "%d.%-02d", @r };
    $Revision = q$Id$;

    # you can add your own here theoretically.
    %HOSTS = (
	       'dbfetch' => {
		   baseurl => 'http://%s/cgi-bin/dbfetch?db=medline&style=raw',
		   hosts   => {
		       'ebi'  => 'www.ebi.ac.uk'
		       }
	       }
	      );
    %FORMATMAP = ( 'default' => 'medlinexml'
		   );
    $DEFAULTFORMAT = 'default';

    $DEFAULT_SERVICE = 'http://www.ebi.ac.uk/cgi-bin/dbfetch';

}


sub new {
    my ($class, @args ) = @_;
    my $self = $class->SUPER::new(@args);

    $self->{ '_hosts' } = {};
    $self->{ '_formatmap' } = {};

    $self->hosts(\%HOSTS);
    $self->formatmap(\%FORMATMAP);
    $self->{'_default_format'} = $DEFAULTFORMAT;

    return $self;
}

=head2 get_by_id

 Title   : get_by_id
 Usage   : $entry = $db->get__by_id('20063307')
 Function: Gets a Bio::Biblio::RefI object by its name
 Returns : a Bio::Biblio::Medline object
 Args    : the id (as a string) of the reference

=cut

sub get_by_id {
    my ($self,$id) = @_;
    my $io = $self->get_Stream_by_id([$id]);
    $self->throw("id does not exist") if( !defined $io ) ;
    return $io->next_bibref();
}



=head2 get_seq_stream

 Title   : get_seq_stream
 Usage   : my $seqio = $self->get_seq_sream(%qualifiers)
 Function: builds a url and queries a web db
 Returns : a Bio::SeqIO stream capable of producing sequence
 Args    : %qualifiers = a hash qualifiers that the implementing class 
           will process to make a url suitable for web querying 

=cut

sub get_seq_stream {
    my ($self, %qualifiers) = @_;
    my ($rformat, $ioformat) = $self->request_format();
    my $seen = 0;
    foreach my $key ( keys %qualifiers ) {
	if( $key =~ /format/i ) {
	    $rformat = $qualifiers{$key};
	    $seen = 1;
	}
    }
    $qualifiers{'-format'} = $rformat if( !$seen);
    ($rformat, $ioformat) = $self->request_format($rformat);
    
    my $request = $self->get_request(%qualifiers);
    my ($stream,$resp);
    if( $self->retrieval_type =~ /temp/i ) {
	my $dir = $self->io()->tempdir( CLEANUP => 1);
	my ( $fh, $tmpfile) = $self->io()->tempfile( DIR => $dir );
	close $fh;
	my ($resp) = $self->_request($request, $tmpfile);		
	if( ! -e $tmpfile || -z $tmpfile || ! $resp->is_success() ) {
            $self->throw("WebDBSeqI Error - check query sequences!\n");
	}
	$self->postprocess_data('type' => 'file',
				'location' => $tmpfile);	
	# this may get reset when requesting batch mode
	($rformat,$ioformat) = $self->request_format();
	if( $self->verbose > 0 ) {
	    open(ERR, "<$tmpfile");
	    while(<ERR>) { $self->debug($_);}
	} 
	$stream = new Bio::Biblio::IO('-format' => $ioformat,
				      '-file'   => $tmpfile);
    } elsif( $self->retrieval_type =~ /io_string/i ) {
	my ($resp) = $self->_request($request);
        my $content = $resp->content_ref;
	$self->debug( "content is $$content\n");
	if( ! $resp->is_success() || length(${$resp->content_ref()}) == 0 ) {
	    $self->throw("WebDBSeqI Error - check query sequences!\n");	
        }  
	($rformat,$ioformat) = $self->request_format();
	$self->postprocess_data('type'=> 'string',
				'location' => $content);
        print STDERR "str is $$content\n" if ( $self->verbose > 0);
	$stream = new Bio::Biblio::IO('-format' => $ioformat,
				      '-data'   => $$content);
    } else { 
	$self->throw("retrieval type " . $self->retrieval_type . 
		     " unsupported\n");
    }
    return $stream;
}

=head2 VERSION and Revision

 Usage   : print $Bio::DB::Biblio::biofetch::VERSION;
           print $Bio::DB::Biblio::biofetch::Revision;

=cut

=head2 Defaults

 Usage   : print $Bio::DB::Biblio::biofetch::DEFAULT_SERVICE;

=cut

1;
__END__
