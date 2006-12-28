# $Id$
#
# BioPerl module for Bio::SearchIO::blastxml
#
# Cared for by Jason Stajich <jason@bioperl.org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::SearchIO::blastxml - A SearchIO implementation of NCBI Blast XML parsing. 

=head1 SYNOPSIS

    use Bio::SearchIO;
    my $searchin = new Bio::SearchIO(-format => 'blastxml',
				     -file   => 't/data/plague_yeast.bls.xml');
    while( my $result = $searchin->next_result ) {
    }

    # one can also request that the parser NOT keep the XML data in memory
    # by using the tempfile initialization flag.
    my $searchin = new Bio::SearchIO(-tempfile => 1,
				     -format => 'blastxml',
				     -file   => 't/data/plague_yeast.bls.xml');
    while( my $result = $searchin->next_result ) {
    }

=head1 DESCRIPTION

This object implements a NCBI Blast XML parser.  It requires XML::SAX; it is
also recommended (for faster parsing) that XML::SAX::ExpatXS be installed and
set as the default parser in ParserDetails.ini.  This file is located in the
SAX subdirectory of XML in your local perl library (normally in the 'site'
directory).  Currently, XML::SAX::Expat will NOT work as expected if set as
default; you must have local copies of the NCBI DTDs if using XML::SAX::Expat.

There is one additional initialization flag from the SearchIO defaults
- that is the -tempfile flag.  If specified as true, then the parser
will write out each report to a temporary filehandle rather than
holding the entire report as a string in memory.  The reason this is
done in the first place is NCBI reports have an uncessary E<lt>?xml
version="1.0"?E<gt> at the beginning of each report and RPS-BLAST reports
have an additional unecessary RPS-BLAST tag at the top of each report.
So we currently have implemented the work around by preparsing the
file (yes it makes the process slower, but it works).

=head1 DEPENDENCIES

In addition to parts of the Bio:: hierarchy, this module uses:

 XML::SAX

It is also recommended that XML::SAX::ExpatXS be installed and made the default
XML::SAX parser using , along with the
Expat library () for faster parsing.  XML::SAX::Expat is not recommended; 
XML::SAX::ExpatXS is considered the current replacement for XML::SAX:Expat
and is actively being considered to replace XML::SAX::Expat.  XML::SAX::Expat
will work, but only if you have local copies of the NCBI BLAST DTDs. This is
due to issues with NCBI's BLAST XML format.  The DTDs and the web address to
obtain them are:

  NCBI_BlastOutput.dtd	    
  NCBI_BlastOutput.mod.dtd

  http://www.ncbi.nlm.nih.gov/data_specs/dtd/

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

=head1 AUTHOR - Jason Stajich

Email jason-at-bioperl.org

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::SearchIO::blastxml;
use strict;
# Object preamble - inherits from Bio::Root::Root

use base qw(Bio::SearchIO);
use Bio::Root::Root;
use XML::SAX;
use IO::File;
use Bio::SearchIO::XML::BlastHandler;

our $DTD = 'ftp://ftp.ncbi.nlm.nih.gov/blast/documents/NCBI_BlastOutput.dtd';

our $DEBUG;

# mapping of NCBI Blast terms to Bioperl hash keys

=head2 new

 Title   : new
 Usage   : my $searchio = new Bio::SearchIO(-format => 'blastxml',
					    -file   => 'filename',
					    -tempfile => 1);
 Function: Initializes the object - this is chained through new in SearchIO
 Returns : Bio::SearchIO::blastxml object
 Args    : One additional argument from the format and file/fh parameters.
           -tempfile    => boolean.  Defaults to false.  Write out XML data
                           to a temporary filehandle to send to PerlSAX parser.
                    
=cut

=head2 _initialize

 Title   : _initialize
 Usage   : private
 Function: Initializes the object - this is chained through new in SearchIO

=cut

sub _initialize{
    my ($self,@args) = @_;   
    $self->SUPER::_initialize(@args);
    my ($usetempfile) = $self->_rearrange([qw(TEMPFILE)],@args);
    defined $usetempfile && $self->use_tempfile($usetempfile);
    
    # uncomment only for testing XML::SAX backend parsers
    #$XML::SAX::ParserPackage = 'XML::SAX::PurePerl';
    
    # BlastHandler does the heavy lifting 
    my $xmlhandler = Bio::SearchIO::XML::BlastHandler->new(-verbose => $self->verbose);
    
    # Pass the SearchIO eventhandler to the XML handler
    # The XML handler does the heavy work, passes data to object handler
    $xmlhandler->_eventHandler($self->_eventHandler());
    
    # start up the parser factory
    my $parserfactory = XML::SAX::ParserFactory->parser(
        Handler => $xmlhandler);
    
    if (ref($parserfactory) eq 'XML::SAX::Expat') {
        $self->throw('XML::SAX::Expat not supported as it is no '.
                     'longer maintained.  Please use any other XML::SAX '.
                     'backend (such as XML::SAX::ExpatXS or XML::LibXML)');
    } elsif (ref($parserfactory) eq 'XML::SAX::PurePerl' && $self->verbose > -1) {
        $self->warn("XML::SAX::PurePerl installed as default XML::SAX parser.\n".
                     "This works but has a small bug which breaks ".
                     "with character encoding (Bug 2159). \n".
                     "We recommend using a different ".
                     "backend (such as XML::SAX::ExpatXS or XML::LibXML)");
    }
    
    $self->{'_xmlparser'} = $parserfactory;
    $self->{'_result_cache'} = [];
    eval {  require Time::HiRes };	
    if( $@ ) { $DEBUG = 0; }
    $DEBUG = 1 if( ! defined $DEBUG && ($self->verbose > 0));
}

=head2 next_result

 Title   : next_result
 Usage   : my $hit = $searchio->next_result;
 Function: Returns the next Result from a search
 Returns : Bio::Search::Result::ResultI object
 Args    : none

=cut

sub next_result {
    my ($self) = @_;
    
    my $result;
    
    if ($result = shift @{$self->{'_result_cache'} }) {
        return $result;
    }
    
    local $/ = "\n";
    local $_;

    my $data = '';
    my $firstline = 1;
    my ($tfh);
    if( $self->use_tempfile ) {
	$tfh = IO::File->new_tmpfile or $self->throw("Unable to open temp file: $!");	
	$tfh->autoflush(1);
    }
   
    my ($sawxmlheader,$okaytoprocess,$sawdoctype);
    while( defined( $_ = $self->_readline) ) {
	if( /^RPS-BLAST/i ) {
	    $self->{'_type'} = 'RPS-BLAST';
	    next;
	}
	if( /^<\?xml version/ ) {
	    if( ! $firstline ) {
		$self->_pushback($_);
		last;
	    }
	    $sawxmlheader = 1;
	}
	# for the non xml version prefixed in each section
	if( /DOCTYPE/ ) { #|| /<BlastOutput>/
	    if(  $sawdoctype ) {
		if( ! $sawxmlheader ) { 
		    $self->_pushback("<?xml version=\"1.0\"?>\n");
		}
		$self->_pushback($_);
		last;
	    }
	    $sawdoctype = 1;
	    unless( $sawxmlheader ) {
		$self->debug( "matched here\n");
		$self->_pushback("<?xml version=\"1.0\"?>\n");
		$self->_pushback($_);
		next;
	    }
	}
	$okaytoprocess = 1;
	if( defined $tfh ) {
	    print $tfh $_;
	} else {
	    $data .= $_;
	}
	$firstline = 0;
    }
    return unless( $okaytoprocess);
    
    my %parser_args;
    if( defined $tfh ) {
	seek($tfh,0,0);
	%parser_args = ('Source' => { 'ByteStream' => $tfh });
    } else {
	%parser_args = ('Source' => { 'String' => $data });
    }

    my $starttime;
    if(  $DEBUG ) {  $starttime = [ Time::HiRes::gettimeofday() ]; }
    
    eval {
	$self->{'_result_cache'} = $self->{'_xmlparser'}->parse(%parser_args);
        $self->{'_result_count'} += scalar(@{ $self->{'_result_cache'} });
        # remove result refs from handler
        $self->{'_xmlparser'}->get_handler->reset_results;
    };
    
    if( $@ ) {
	$self->warn("error in parsing a report:\n $@");
	$result = undef;
    }    
    if( $DEBUG ) {
	$self->debug( sprintf("parsing took %f seconds\n", Time::HiRes::tv_interval($starttime)));
    }
    # parsing magic here - but we call event handlers rather than 
    # instantiating things 
    return shift @{ $self->{'_result_cache'} };
}

=head2 use_tempfile

 Title   : use_tempfile
 Usage   : $obj->use_tempfile($newval)
 Function: Get/Set boolean flag on whether or not use a tempfile
 Example : 
 Returns : value of use_tempfile
 Args    : newvalue (optional)


=cut

sub use_tempfile{
   my ($self,$value) = @_;
   if( defined $value) {
      $self->{'_use_tempfile'} = $value;
    }
    return $self->{'_use_tempfile'};
}

sub result_count {
    my $self = shift;
    return $self->{'_result_count'};
}

sub no_preparse {
    my $self = shift;
    return $self->{'_result_count'} = shift if @_;
    return $self->{'_result_count'};
}

sub saxparser {
    my $self = shift;
    return ref($self->{'_xmlparser'});
}

1;
