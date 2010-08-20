# $Id$
#
# BioPerl module for Bio::DB::AlignIO::Prosite
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Jun Yin <jun dot yin at ucd dot ie>
#
# Copyright Jun Yin
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code
# 
# RESTful service-based extension of GenericWebDBI interface 

=head1 NAME

Bio::DB::AlignIO::Prosite - webagent which interacts with and retrieves alignment 
sequences from Prosite

=head1 SYNOPSIS

	use Bio::DB::Align;
	use Bio::DB::Align::Prosite;
	use Bio::SimpleAlign; 
  
	my $dbobj=Bio::DB::Align::Prosite->new(); #create a db object
	my $dbobj2=Bio::DB::Align->new(-db=>"Prosite"); #create a db object
	                                           #the same with above

	#retrieve a Bio::SimpleAlign object
	my $aln=$dbobj->get_Aln_by_acc("PS51092"); 
	
	#do something here with the align object
	print $aln->length,"\n";
	print $aln->num_sequences,"\n";	
	
	foreach my $seq ($aln->next_Seq) {
		#do something with $seq
	}
	
	#or parameter based calling
	my $aln2=$dbobj->get_Aln_by_acc(-accession=>"PS00023",
	                                -format=>"clustalw");
	print $aln2->accession,"\n";

=head1 DESCRIPTION

This package uses the RESTful service provided by Prosite. It retrieves
online alignment sequences and save it into Bio::SimpleAlign object.
Bio::DB::Align::Prosite only supports alignment sequence retrieval using
protein domain accession number, e.g. "RF00360".

=head1 TODO

=over 3

=item * Finish documentation

HOWTOs (both standard and Cookbook).

=item * Cookbook tests

Set up dev-only tests for Cookbook examples to make sure they are consistently
updated.

=back

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the 
evolution of this and other Bioperl modules. Send
your comments and suggestions preferably to one
of the Bioperl mailing lists. Your participation
is much appreciated.

  bioperl-l@lists.open-bio.org               - General discussion
  http://www.bioperl.org/wiki/Mailing_lists  - About the mailing lists

=head2 Support 

Please direct usage questions or support issues to the mailing list:

I<bioperl-l@bioperl.org>

rather than to the module maintainer directly. Many experienced and 
reponsive experts will be able look at the problem and quickly 
address it. Please include a thorough description of the problem 
with code and data examples if at all possible.

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to
help us keep track the bugs and their resolution.
Bug reports can be submitted via the web.

  http://bugzilla.open-bio.org/

=head1 AUTHOR 

Email jun dot yin at ucd dot ie

=head1 APPENDIX

The rest of the documentation details each of the
object methods. Internal methods are usually
preceded with a _

=cut

# Let the code begin...

package Bio::DB::Align::Prosite;
use strict;
use warnings;
use LWP::UserAgent;
use HTTP::Request;
use Bio::AlignIO;
use Bio::Root::IO;
use vars qw(%FORMATS %ALNTYPE $HOSTBASE $TIMEOUT);

use base qw(Bio::Root::Root Bio::Root::IO Bio::DB::GenericWebAgent);

BEGIN {
	$HOSTBASE = 'http://expasy.org';
	%FORMATS=qw(fasta 1 clustalw 1); #supported formats in Prosite
	$TIMEOUT=1000;
}

=head2 new

 Title   : new
 Usage   : $dbobj = Bio::DB::Align->new(-db=>"Prosite");
             Or, it can be called through specific package
             $dbobj = Bio::DB::Align::Prosite->new();
 Function: Returns a Bio::DB::Align::Prosite stream
 Returns : Bio::DB::Align::Prosite object
 Args    : 
 Note    : 
=cut


sub new {
	my($class,@args) = @_;
	my $self = $class->SUPER::new(@args);
   my $ua = new LWP::UserAgent(env_proxy => 1);
   #$ua->agent(ref($self) ."/$MODVERSION");
   $self->ua($ua);  
   $self->ua->timeout($TIMEOUT);
   $self->{'_authentication'} = [];	
	return $self;	
}


=head2 get_Aln_by_id

 Title   : get_Aln_by_id
 Usage   : $aln = $dbobj->get_Aln_by_id($id)
 Function: Gets a Bio::SimpleAlign object by id
 Returns : a Bio::SimpleAlign object
 Args    : -id  the id as a string
          -format     the output format from Prosite. This will decide which
                      package to use in the Bio::AlignIO
                      possible options can be fasta (default) or clustalw
 Note    : 
 Throws  : "Bio::Root::NotImplemented" exception
=cut

sub get_Aln_by_id {
	my ($self,@args)=@_;
	$self->throw_not_implemented();
}


=head2 get_Aln_by_acc

 Title   : get_Aln_by_acc
 Usage   : $aln = $dbobj->get_Aln_by_acc("PS51092");
 Function: Gets a Bio::SimpleAlign object by accession numbers
 Returns : a Bio::SimpleAlign object
 Args    : -accession  the accession number as a string
           -format     the output format from Prosite. This will decide which
                        package to use in the Bio::AlignIO
                        possible options can be fasta (default) or clustalw
 Note    : 
 Throws  : "Bio::DB::Align::Prosite Request Error" exception
=cut

sub get_Aln_by_acc {
	my ($self,@args)=@_;
	
	my ($acc,$format)=$self->_checkparameter(@args);
	
	my $CGI_location= '/cgi-bin/aligner';
	my %params;
	if($format eq "fasta") {
		%params= (
			"psa"=>$acc,
			"format"=>$format,
		);
	}
	else {
		#default format is clustal
		$params{"psa"}=$acc;
	}

	my $url = URI->new($HOSTBASE . $CGI_location);
	$url->query_form(%params);
	
	#save the retrieved file into a tempfile
	my $dir = $self->tempdir( CLEANUP => 1);
	my ( $fh, $tmpfile) = $self->tempfile( DIR => $dir );
	close $fh;
	
	my $request = HTTP::Request->new(GET => $url);
	#my $request = $self->ua->get($url,content_file=>$tmpfile);
	
	$request->proxy_authorization_basic($self->authentication)
	  if ( $self->authentication);
	$self->debug("request is ". $request->as_string(). "\n");
	
	my $respond = $self->ua->request($request,$tmpfile);
	
	if( $respond->is_error  ) {
		$self->throw("Bio::DB::Align::Prosite Request Error:\n".$request->as_string);
	}
	
	my $alnobj=Bio::AlignIO->new(-format=>$format,-file=>$tmpfile);
	my $aln=$alnobj->next_aln;
	
	#record the accession and id
	unless ($aln->accession) {
		$aln->accession($acc);
	}
	unless ($aln->source) {
		$aln->source("Prosite");
	}	
	
	return $aln;
}

=head2 id2acc

 Title   : id2acc
 Usage   : $acc = $db->id2acc($id)
 Function: Convert id to accession
 Returns : Accession
 Args    : ID (as a string)
 Throws  : "Bio::Root::NotImplemented" exception
=cut

sub id2acc {
	my ($self,@args)=@_;
	$self->throw_not_implemented();
}

=head2 acc2id

 Title   : acc2id
 Usage   : $id = $db->acc2id($acc)
 Function: Convert Accession to ID
 Returns : ID
 Args    : Accession (as a string)
 Throws  : "Bio::Root::NotImplemented" exception
=cut

sub acc2id {
	my ($self,@args)=@_;
	$self->throw_not_implemented();
}


sub _checkparameter {
	my ($self,@args)=@_;
	my ($acc,$format)=$self->_rearrange([qw(ACCESSION FORMAT)],@args);

	#check format
	if($format) {
		$format=lc $format;
		unless(defined $FORMATS{$format}) {
			$self->throw("Only [".join(" ",sort keys %FORMATS)."] are supported by Prosite, not [".$format."]");
		}
	}
	else {
		$format="fasta"; #default format
	}
		
	return ($acc,$format);
}

1;


                  