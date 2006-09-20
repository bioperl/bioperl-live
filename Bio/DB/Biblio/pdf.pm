# $Id$
#
# BioPerl module Bio::DB::Biblio::pdf.pm
#
# Cared for by Allen Day <allenday@ucla.edu>
# For copyright and disclaimer see below.

# POD documentation - main docs before the code

=head1 NAME

Bio::DB::Biblio::pdf - Fetch PDF for a PubMed ID

=head1 SYNOPSIS

Do not use this object directly, it is recommended to access it and use
it through the I<Bio::Biblio> module:

  use Bio::Biblio;
  my $biblio = new Bio::Biblio (-access => 'pdf');

=head1 DESCRIPTION

This object contains the real implementation of a Bibliographic Query
Service as defined in L<Bio::DB::BiblioI>.

L<Bio::DB::BiblioI> is not implemented as documented in the interface,
particularly the find() method, which is not compatible with PubMed's
query language.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

  bioperl-l@bioperl.org                  - General discussion
  http://bioperl.org/wiki/Mailing_lists  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution. Bug reports can be submitted via the
web:

  http://bugzilla.open-bio.org/

=head1 AUTHOR

Allen Day E<lt>allenday@ucla.eduE<gt>

=head1 COPYRIGHT

Copyright (c) 2004 Allen Day, University of California, Los Angeles.

This module is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.

=head1 DISCLAIMER

This software is provided "as is" without warranty of any kind.

=head1 BUGS AND LIMITATIONS

=over

=item *

Of course, you'll need access to the sites hosting the PDFs to download
them.

= item *

If you're having problems retrieving PDF from a site you have access to,
you might try adjusting the max_depth() attribute.  It is default set to 3,
and affects how many links deep will be recursively followed in page
fetches to try to find your PDF.

=back

=head1 SEE ALSO

 Pub Med Help
 http://eutils.ncbi.nlm.nih.gov/entrez/query/static/help/pmhelp.html

 Entrez Utilities
 http://eutils.ncbi.nlm.nih.gov/entrez/query/static/eutils_help.html

=head1 APPENDIX

The main documentation details are to be found in
L<Bio::DB::BiblioI>.

Here is the rest of the object methods.  Interface methods first,
followed by internal methods.

=cut

# Let the code begin...


package Bio::DB::Biblio::pdf;
use vars qw($DEFAULT_URN);
use strict;

use Data::Dumper;
use WWW::Mechanize;
use base qw(Bio::Biblio Bio::DB::BiblioI);

use constant DEBUG         => 1;
use constant NCBI_BASE     => 'http://www.ncbi.nlm.nih.gov';
use constant ABSTRACT_BASE => NCBI_BASE . '/entrez/query.fcgi?cmd=Retrieve&db=pubmed&dopt=Abstract&list_uids=';

my %visit = (); #for spidering

# -----------------------------------------------------------------------------

=head2 _initialize

 Usage   : my $obj = new Bio::Biblio (-access => 'pdf' ...);
           (_initialize is internally called from this constructor)
 Returns : 1 on success
 Args    : none

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

    #pdf doesn't need this code, but it doesn't hurt to leave it here... -ad

    # make a hashtable from @args
    my %param = @args;
    @param { map { lc $_ } keys %param } = values %param; # lowercase keys

    $self->max_depth(3);
    $self->depth(0);
    $self->ua( WWW::Mechanize->new());
    $self->ua->agent_alias('Linux Mozilla');

    #now override with passed hash

    # copy all @args into this object (overwriting what may already be
    # there) - changing '-key' into '_key'
    my $new_key;
    foreach my $key (keys %param) {
      ($new_key = $key) =~ s/^-/_/;
      $self->{ $new_key } = $param { $key };
    }

    #AOK
    return 1;
}

=head2 get_next

  Title   : get_next
  Usage   : $xml = $biblio->get_next();
  Function: return next record as xml
  Returns : an xml string
  Args    : none


=cut

sub get_next {
  my $self = shift;

  return $self->pdf();

  return;
}

=head2 find

  Title   : find
  Usage   : $biblio = $biblio->find(1234);
  Function: perform a PubMed query by PubMed ID
  Returns : a reference to the object on which the method was called
  Args    : a PubMed ID

=cut

sub find {
  my ($self,$id) = @_;

  $self->{pdf} = undef;
  $self->_process_pubmed_html($id);
}

=head2 exists

  Title   : exists
  Usage   : do not use
  Function: no-op.  this is here only for interface compatibility
  Returns : undef
  Args    : none


=cut

sub exists {
  return;

}

=head2 destroy

  Title   : destroy
  Usage   : do not use
  Function: no-op.  this is here only for interface compatibility
  Returns : undef
  Args    : none


=cut

sub destroy {
  return;

}

=head2 get_vocabulary_names

  Title   : get_vocabulary_names
  Usage   : do not use
  Function: no-op.  this is here only for interface compatibility
  Returns : empty arrayref
  Args    : none


=cut

sub get_vocabulary_names {
  return [];
}

=head2 contains

  Title   : contains
  Usage   : do not use
  Function: no-op.  this is here only for interface compatibility
  Returns : undef
  Args    : none


=cut

sub contains {
  return;
}

=head2 get_entry_description

  Title   : get_entry_description
  Usage   : do not use
  Function: no-op.  this is here only for interface compatibility
  Returns : undef
  Args    : none


=cut

sub get_entry_description {
  return;
}

=head2 get_all_values

  Title   : get_all_values
  Usage   : do not use
  Function: no-op.  this is here only for interface compatibility
  Returns : undef
  Args    : none


=cut

sub get_all_values {
  return;
}

=head2 get_all_entries

  Title   : get_all_entries
  Usage   : do not use
  Function: no-op.  this is here only for interface compatibility
  Returns : undef
  Args    : none


=cut

sub get_all_entries {
  return;
}

=head1 Internal methods unrelated to Bio::DB::BiblioI

=cut

=head2 depth()

 Usage   : $obj->depth($newval)
 Function: track link recursion depth
 Example : 
 Returns : value of depth (a scalar)
 Args    : on set, new value (a scalar or undef, optional)


=cut

sub depth {
  my($self,$val) = @_;
  $self->{'depth'} = $val if defined($val);
  return $self->{'depth'};
}

=head2 max_depth()

 Usage   : $obj->max_depth($newval)
 Function: how far should link recursion go?
 Example : 
 Returns : value of max_depth (a scalar)
 Args    : on set, new value (a scalar or undef, optional)


=cut

sub max_depth {
  my($self,$val) = @_;
  $self->{'max_depth'} = $val if defined($val);
  return $self->{'max_depth'};
}


=head2 ua()

 Usage   : $obj->ua($newval)
 Function: holds an LWP::UserAgent instance
 Example : 
 Returns : value of ua (a scalar)
 Args    : on set, new value (a scalar or undef, optional)


=cut

sub ua {
  my($self,$val) = @_;
  $self->{'ua'} = $val if defined($val);
  return $self->{'ua'};
}

sub _process_pubmed_html {
  my ($self,$id) = @_;

  $self->ua->get( ABSTRACT_BASE . $id );
  my $page = $self->ua->content();

  #here is the treasure
  $page =~ m|<!---- Pager -- \(page header\) -- end ------>.+?<SPAN><a href="(.+?)" onClick="window.open|s;

  if( ! defined($1) ) {
    return;
  }

  $self->ua->follow_link( url => $1 );

  #uncomment this to do site crawl -- old style
  #$self->_crawl();

  my $pdf_url = $self->guess_pdf_url($self->ua->uri);
  $self->throw( "didn't recognize pattern in '".$self->ua->uri."', please patch module" ) unless $pdf_url;
  $self->ua->get( $pdf_url );
  my $content = $self->ua->content();
  $self->pdf( $content );
}

sub guess_pdf_url {
  my($self,$url) = @_;

  #cancer research
  if( $url =~ m!^(.+?)/cgi/content/full/(\d+)/(\d+)/(\d+)/?$! ) {
    return qq($1/cgi/reprint/$2/$3/$4.pdf);
  }
  #nature
  #http://www.nature.com/cgi-taf/DynaPage.taf?file=/onc/journal/v18/n27/abs/1202776a.html&dynoptions=doi1108513968
  #http://www.nature.com/cgi-taf/DynaPage.taf?file=/onc/journal/v18/n27/full/1202776a.html&filetype=pdf
  elsif( $url =~ m!^(.+?cgi-taf/DynaPage.taf.+?)/journal/(.+?)/abs/(.+?\.html)! ) {
    return qq($1/journal/$2/full/$3\&filetype=pdf);
  }
  #science direct
  #these pages contain some unpredictable md5 bullshit, so we need to parse the page
  #http://www.sciencedirect.com/science?_ob=ArticleURL&_udi=B6VPM-480CTTS-5&_coverDate=04%2F30%2F2003&_alid=247076467&_rdoc=1&_fmt=&_orig=search&_qd=1&_cdi=6210&_sort=d&view=c&_acct=C000059605&_version=1&_urlVersion=0&_userid=4423&md5=8054dea49e32e98a6b30b206ea47fbfe
  #http://www.sciencedirect.com/science?_ob=MImg&_imagekey=B6VPM-480CTTS-5-5&_cdi=6210&_user=4423&_orig=search&_coverDate=04%2F30%2F2003&_qd=1&_sk=999779997&view=c&wchp=dGLbVtz-zSkzV&md5=5b04979d84dab066be5cde52fd2affa7&ie=/sdarticle.pdf
  elsif( $url =~ m!^(.+?science\?_ob=)ArticleURL(.+?)$! ) {
    my $link = $self->ua->find_link( text_regex => qr/PDF \(.+?\)/s );
    return unless $link;
    return $link->url_abs();
  }
  #genome biology
  #http://genomebiology.com/2003/4/7/R43
  #http://genomebiology.com/content/pdf/gb-2003-4-7-r43.pdf
  elsif( $url =~ m!^(.+?genomebiology.com)/(\d+)/(\d+)/(\d+)/(.+?)/?$! ) {
    my $file = lc(sprintf("gb-%d-%d-%d-%s.pdf",$2,$3,$4,$5));
    return qq($1/content/pdf/$file);
  }
  #wiley interscience
  #http://www3.interscience.wiley.com/cgi-bin/abstract/91013753/ABSTRACT
  #http://www3.interscience.wiley.com/cgi-bin/fulltext/91013753/PDFSTART
  #http://download.interscience.wiley.com/cgi-bin/fulltext?ID=96515300&PLACEBO=IE.pdf&mode=pdf
  elsif( $url =~ m!^(.+?/cgi-bin)/abstract/(\d+?)/ABSTRACT$! ) {
    $self->ua->get( qq($1/fulltext/$2/PDFSTART) );
    my $link = $self->ua->find_link( url_regex => qr/fulltext/ );
    return unless $link;
    return $link->url_abs();
  }
  #nar, bioinformatics
  #http://nar.oupjournals.org/cgi/content/full/32/suppl_1/D258
  #http://nar.oupjournals.org/cgi/reprint/32/suppl_1/D258.pdf
  elsif( $url =~ m!^(.+?oupjournals.org/cgi)/reprint/(.+?)$! ) {
    return qq($1/reprint/$2.pdf);
  }
  elsif( $url =~ m!^(.+?oupjournals.org/cgi)/content/full/(.+?)$! ) {
    return qq($1/reprint/$2.pdf);
  }

  #plos
  #http://biology.plosjournals.org/plosonline/?request=get-document&doi=10.1371/journal.pbio.0020009
  #http://www.plosbiology.org/archive/1545-7885/2/1/pdf/10.1371_journal.pbio.0020009-S.pdf
  elsif( $url =~ m!^http://[^.]+?\.plos! ) {
    my $link = $self->ua->find_link( text_regex => qr/^Screen/s );
    return unless $link;
    return $link->url_abs();
  }
  #bmc bioinformatics
  #http://www.biomedcentral.com/1471-2105/2/7
  #http://www.biomedcentral.com/content/pdf/1471-2105-2-7.pdf
  elsif( $url =~ m!^(.+?biomedcentral.+?)/(\d+\-\d+)/(\d+)/(\d+)/?$! ) {
    my $file = lc(sprintf("%s-%d-%d.pdf",$2,$3,$4));
    return qq($1/content/pdf/$file);
  }


warn $url;
  return;
}


# sub _crawl {
#   my( $self ) = @_;

#   return if $self->depth() == $self->max_depth();
#   return if $self->pdf();

#   $self->depth( $self->depth + 1 );

#   #try to find "PDF" link first
#   my ( $link ) = $self->ua->find_link( text_regex => qr/PDF|View article/ );
#   if ( $link ) {
#     $self->_fetch_pdf( $link );
#   }
#   else {
#     foreach my $link ( $self->ua->find_all_links ) {
#       $self->_fetch_pdf( $link );
#     }
#   }

#   $self->depth( $self->depth - 1 );

#   return;
# }

# sub _fetch_pdf {
#   my $self = shift;
#   my $link = shift;

#   return if $visit{ $link->url_abs };
#   $visit{ $link->url_abs }++;

#   $self->ua->get( $link->url_abs );
#   print "[" . $self->depth() . "] fetching: " . $link->url_abs . " " . $self->ua->ct() . "\n" if DEBUG;

#   #test for a likely string "href", because some misconfigured webservers will send pdf
#   #as text/html
#   if ( $self->ua->ct() eq 'application/pdf' or
#        ( $self->ua->ct() =~ /text/ and $self->ua->content !~ /href|src/is )
#      ) {
#     print "*****FOUND IT (" . $link->url_abs . ") *****\n" if DEBUG;

#     $self->ua->get( $link->url_abs );
#     my $content = $self->ua()->content();
#     $self->pdf( $content );
#   }
#   else {
#     $self->_crawl();
#   }
# }

=head2 pdf()

 Usage   : $obj->pdf($newval)
 Function: holds pdf data
 Example : 
 Returns : value of pdf (a scalar)
 Args    : on set, new value (a scalar or undef, optional)


=cut

sub pdf {
  my($self,$val) = @_;
  $self->{'pdf'} = $val if defined($val);
  return $self->{'pdf'};
}


1;
