#
# BioPerl module for Bio::SearchIO::blastxml
#
# Please direct questions and support issues to <bioperl-l@bioperl.org>
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
    my $searchin = Bio::SearchIO->new(-format => 'blastxml',
				     -file   => 't/data/plague_yeast.bls.xml');

    while( my $result = $searchin->next_result ) {
        ....
    }

    # one can also request that the parser NOT keep the XML data in memory
    # by using the tempfile initialization flag.

    $searchin = Bio::SearchIO->new(-tempfile => 1,
				     -format => 'blastxml',
				     -file   => 't/data/plague_yeast.bls.xml');

    while( my $result = $searchin->next_result ) {
       ....
    }

    # PSI-BLAST parsing (default is normal BLAST)
    $searchin = Bio::SearchIO->new(
                     -format => 'blastxml',
                     -blasttype => 'psiblast',
				     -file   => 't/data/plague_yeast.bls.xml');

    while( my $result = $searchin->next_result ) {
       ....
    }

=head1 DESCRIPTION

This object implements a NCBI Blast XML parser.  It requires XML::SAX; it is
also recommended (for faster parsing) that XML::SAX::ExpatXS or XML::LibXML
be installed.  Either 'XML::SAX::ExpatXS' or 'XML::LibXML::SAX::Parser' should
be set as the default parser in ParserDetails.ini.  This file is located in the
SAX subdirectory of XML in your local perl library (normally in the 'site'
directory).

Two different XML handlers currently exist to deal with logical differences
between how normal BLAST reports and PSI-BLAST reports are logically parsed into
BioPerl objects; this is explicitly settable using the B<-blasttype> parameter.
The default is for parsing a normal BLAST report ('blast'), but if one is
expecting PSI-BLAST report parsing, -blasttype B<must> be set explicitly to
'psiblast'. This is due to a lack of any information in the XML output which
tells the parser the report is derived from a PSI-BLAST run vs. a normal BLAST
run.

There is one additional initialization flag from the SearchIO defaults. That is
the B<-tempfile> flag. If specified as true, then the parser will write out each
report to a temporary filehandle rather than holding the entire report as a
string in memory. The reason this is done in the first place is NCBI reports
have an uncessary E<lt>?xml version="1.0"?E<gt> at the beginning of each report
and RPS-BLAST reports have an additional unnecessary RPS-BLAST tag at the top of
each report. So we currently have implemented the work around by preparsing the
file (yes it makes the process slower, but it works). We are open to suggestions
on how to optimize this in the future.

=head1 DEPENDENCIES

In addition to parts of the Bio:: hierarchy, this module uses:

 XML::SAX

It is also recommended that XML::SAX::ExpatXS be installed and made the default
XML::SAX parser using , along with the Expat library () for faster parsing.
XML::SAX::Expat is not recommended; XML::SAX::ExpatXS is considered the current
replacement for XML::SAX:Expat and is actively being considered to replace
XML::SAX::Expat. XML::SAX::Expat will work, but only if you have local copies of
the NCBI BLAST DTDs. This is due to issues with NCBI's BLAST XML format. The
DTDs and the web address to obtain them are:

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

  https://github.com/bioperl/bioperl-live/issues

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
use Bio::SearchIO::IteratedSearchResultEventBuilder;

our $DEBUG;

my %VALID_TYPE = (
    'BLAST'      => 'Bio::SearchIO::XML::BlastHandler',
    'PSIBLAST'   => 'Bio::SearchIO::XML::PsiBlastHandler',
    'PSI-BLAST'  => 'Bio::SearchIO::XML::PsiBlastHandler'
    );

# mapping of NCBI Blast terms to Bioperl hash keys

=head2 new

 Title   : new
 Usage   : my $searchio = Bio::SearchIO->new(-format => 'blastxml',
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
    my ($usetempfile, $blasttype,$xmlcompact) = $self->_rearrange([qw(
                                            TEMPFILE
                                            BLASTTYPE
                                            XMLCOMPACT)],@args);
    $blasttype ||= 'BLAST';
    $self->{_xml_compact} = $xmlcompact || 0;
    $self->blasttype(uc $blasttype);
    defined $usetempfile && $self->use_tempfile($usetempfile);
    $self->{_result_count} = 0;
    eval {  require Time::HiRes };
    if( $@ ) { $DEBUG = 0; }
    $DEBUG = 1 if( ! defined $DEBUG && ($self->verbose > 0));
}

sub attach_EventHandler {
    my ($self,$handler) = @_;

    $self->SUPER::attach_EventHandler($handler);

	# Make sure if there is an XML parser present already, the internal Handler
	# is set
	if (exists $self->{'_xmlparser'}) {
		$self->{'_xmlparser'}->get_handler->eventHandler($handler);
	}

    # Optimization: caching the EventHandler since it is used a lot
    # during the parse.

    $self->{'_handler_cache'} = $handler;
    return;
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

    my ($tfh);

    # XMLCOMPACT
    # WU-BLAST has an XML_COMPACT option which needs to be preprocessed before
    # passing on to the parser.
    if ($self->{_xml_compact}) {
        $self->debug("XMLCOMPACT mode\n");
        my ($tfh2, $filename) = IO::File->new_tmpfile or $self->throw("Unable to open temp file: $!");
        $tfh2->autoflush(1);
        my $fh = $self->_fh;
        while (my $line = <$fh>) {
            $line =~ s/></>\n</g;
            print $tfh2 $line;
        }
        seek($tfh2,0,0);
        close $fh;
        # redirect self's IO to use new tempfile
        $self->_fh($tfh2);
    }

    if( $self->use_tempfile ) {
        $tfh = IO::File->new_tmpfile or $self->throw("Unable to open temp file: $!");
        $tfh->autoflush(1);
    }

    my $okaytoprocess = ($self->blasttype =~ /PSI/) ? $self->_chunk_psiblast($tfh) :
        $self->_chunk_normalblast($tfh);

    return unless( $okaytoprocess);

    my %parser_args;
    if( defined $tfh ) {
	seek($tfh,0,0);
	%parser_args = ('Source' => { 'ByteStream' => $tfh });
    } else {
	%parser_args = ('Source' => { 'String' => $self->{'_blastdata'} });
    }

    my $starttime;
    if(  $DEBUG ) {  $starttime = [ Time::HiRes::gettimeofday() ]; }

    eval {
	$result = $self->{'_xmlparser'}->parse(%parser_args);
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
    if (defined $result) {
        # result count is handled here, as the BLASTXML reports are
        # broken up into smaller easier to digest bits
        $self->{_result_count}++;
        return $result;
    } else {
        return;
    }
}

=head2 result_count

 Title   : result_count
 Usage   : $num = $stream->result_count;
 Function: Gets the number of Blast results that have been successfully parsed
           at the point of the method call.  This is not the total # of results
           in the file.
 Returns : integer
 Args    : none
 Throws  : none

=cut

sub result_count {
    my $self = shift;
    return $self->{_result_count};
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

=head2 blasttype

 Title   : blasttype
 Usage   : $obj->blasttype($newtype)
 Function: Get/Set BLAST report type.
 Returns : BLAST report type
 Args    : case-insensitive string of types BLAST or PSIBLAST (default: BLAST)
 Note    : this is used to determine how reports are 'chunked' (in cases
           where multiple queries are submitted) and which XML handler
           to use when parsing the report(s)

=cut

sub blasttype{
    my ($self,$value) = @_;
    if ($value) {
        $self->throw("$value is not a supported BLAST type") unless exists $VALID_TYPE{$value};
        my $ok;
        eval {
            $ok = $self->_load_module($VALID_TYPE{$value});
        };
        if ($@) {
            print STDERR <<END;
$self: data module $VALID_TYPE{$value} cannot be found
Exception $@
For more information about the Bio::SearchIO::blastxml system please see the Bio::SearchIO::blastxml.
END
            return unless $ok;
        }
        # BlastHandler does the heavy lifting
        my $xmlhandler = $VALID_TYPE{$value}->new(-verbose => $self->verbose);

        # The XML handler does the heavy work, passes data to object handler
        if ($value =~ /^PSI/) {
            my $handler = Bio::SearchIO::IteratedSearchResultEventBuilder->new();
            $self->{'_handler'} = $handler; # cache
        }
        $xmlhandler->eventHandler($self->_eventHandler());

        # start up the parser factory
        my $parserfactory = XML::SAX::ParserFactory->parser(
            Handler => $xmlhandler);
        $self->{'_xmlparser'} = $parserfactory;
        $self->saxparser(ref($parserfactory));

        $self->{'_blasttype'} = $value;
    }
    return $self->{'_blasttype'};
}

sub saxparser {
    my $self = shift;
    return ref($self->{'_xmlparser'});
}

sub _chunk_normalblast {
    my ($self, $tfh) = @_;

    local $/ = "\n";
    local $_;
    $self->{'_blastdata'} = '';

    my ($sawxmlheader, $okaytoprocess);

    my $mode = 'header';

    my $tail = << 'XML_END';
  </BlastOutput_iterations>
</BlastOutput>
XML_END

    # no buffering needed (famous last words...)
    my $fh = $self->_fh;

    #chop up XML into edible bits for the parser
    while( defined( my $line = <$fh>) ) {
        next if $line =~ /^\s*$/;
        next if $line =~ m{^\s*</BlastOutput_iterations>}xmso || $line =~ m{^</BlastOutput>}xmso;
        if( $line =~ m{^RPS-BLAST}i ) {
            $self->{'_type'} = 'RPS-BLAST';
            next;
        } elsif ($line =~ m{^<\?xml\sversion="1.0"}xms) {# <?xml version="1.0"?> & <?xml version="1.0" encoding="UTF-8"?>
            delete $self->{'_header'} if exists $self->{'_header'};
            $sawxmlheader++;
            $mode = 'header';
        } elsif ($line =~ m{^\s*<Iteration>}xmso) {
            if (!$sawxmlheader) {
                if (defined $tfh) {
                    print $tfh $self->{'_header'}
                } else {
                    $self->{'_blastdata'} .= $self->{'_header'};
                }
            }
            $mode = 'iteration';
        } elsif ($line =~ m{^\s*</Iteration>}xmso) {
            if (defined $tfh) {
                print $tfh $line.$tail;
            } else {
                $self->{'_blastdata'} .= $line.$tail;
            }
            $okaytoprocess++;
            last;
        }
        if (defined $tfh) {
            print $tfh $line;
        } else {
            $self->{'_blastdata'} .= $line;
        }
        $self->{"_$mode"} .= $line if $mode eq 'header';
    }
    return $okaytoprocess;
}

sub _chunk_psiblast {
    my ($self, $tfh) = @_;

    local $/ = "\n";
    local $_;
    $self->{'_blastdata'} = '';

    my ($sawxmlheader, $okaytoprocess);

    # no buffering needed (famous last words...)
    my $fh = $self->_fh;

    #chop up XML into edible bits for the parser
    while( defined( my $line = <$fh>) ) {
        if (defined $tfh) {
            print $tfh $line;
        } else {
            $self->{'_blastdata'} .= $line;
        }
        #$self->{"_$mode"} .= $line;
        if ($line =~ m{^</BlastOutput>}xmso) {
            $okaytoprocess++;
            last;
        }
    }
    #$self->debug($self->{'_blastdata'}."\n");
    return $okaytoprocess;
}

1;
