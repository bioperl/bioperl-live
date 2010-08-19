# $Id$
#
# BioPerl module for Bio::DB::AlignIO::ProtClustDB
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

Bio::DB::AlignIO::ProtClustDB - webagent which interacts with and retrieves alignment 
sequences from ProtClustDB (Entrez Protein Clusters Database)

=head1 SYNOPSIS

  # ...To be added!

=head1 DESCRIPTION

	# ...To be added!

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

package Bio::DB::Align::ProtClustDB;
use strict;
use warnings;
use LWP::UserAgent;
use HTTP::Request;
use Bio::AlignIO;
use Bio::SimpleAlign;
use Bio::Root::IO;
use Bio::DB::EUtilities;
use vars qw(%FORMATS %ALNTYPE $HOSTBASE);

use base qw(Bio::Root::Root Bio::Root::IO Bio::DB::GenericWebAgent);

BEGIN {
	$HOSTBASE = 'http://www.ncbi.nlm.nih.gov';
	#http://www.ncbi.nlm.nih.gov/sutils/prkview.cgi?result=align&cluster=PRK00286&image=download
	%FORMATS=qw(fasta 1 stockholm 1 pfam 1 fastau 1); #supported formats in pfam
	%ALNTYPE=qw(seed 1 full 1); #supported alignment types
}


=head2 new

 Title   : new
 Usage   : $dbobj = Bio::DB::Align->new(-db=>"ProtClustDB",
                                        -email=>'my@foo.bar');
             Or, it can be called through specific package
             $dbobj = Bio::DB::Align::ProtClustDB->new(-email=>'my@foo.bar');
 Function: Returns a Bio::DB::Align::ProtClustDB stream
 Returns : Bio::DB::Align::ProtClustDB object
 Args    : -db     database used
           -email  email address requested per NCBI policy
 Note    : 
=cut


sub new {
	my($class,@args) = @_;
	my $self = $class->SUPER::new(@args);
	my ($email)=$self->_rearrange(["EMAIL"],@args);
	if($email) {
		$self->email($email);
	}
	else {
		$self->warn('The -email parameter is now required for Bio::DB::Align::ProtClustDB, per NCBI E-utilities policy');
	}
   my $ua = new LWP::UserAgent(env_proxy => 1);
   #$ua->agent(ref($self) ."/$MODVERSION");
   $self->ua($ua);  
   $self->{'_authentication'} = [];	
	return $self;	
}


=head2 get_Aln_by_id

 Title   : get_Aln_by_id
 Usage   : $aln = $db->get_Aln_by_id(2725839)
 Function: Gets a Bio::SimpleAlign object by id
 Returns : a Bio::SimpleAlign object
 Args    : -id  the ID as a string
 Note    : 
 Throws  : "Bio::DB::Align::ProtClustDB Request Error" exception
=cut

sub get_Aln_by_id {
	my ($self,@args)=@_;
	my ($id)=$self->_rearrange(["ID"],@args);
	
	my $acc=$self->id2acc($id);
	
	my $aln=$self->_get_request($acc,$id);
	
	$aln->accession($acc);
	$aln->id($id);	
	
	return $aln;
}

=head2 get_Aln_by_acc

 Title   : get_Aln_by_acc
 Usage   : $seq = $db->get_Aln_by_acc("CLSN2725839");
 Function: Gets a Bio::SimpleAlign object by accession numbers
 Returns : a Bio::SimpleAlign object
 Args    : -accession  the accession number as a string
 Note    : 
 Throws  : "Bio::DB::Align::ProtClustDB Request Error" exception
=cut

sub get_Aln_by_acc {
	my ($self,@args)=@_;
	
	my ($acc)=$self->_rearrange(["ACCESSION"],@args);
	
	my $id=$self->acc2id($acc);
	
	my $aln=$self->_get_request($acc,$id);
	
	return $aln;
}

=head2 id2acc

 Title   : id2acc
 Usage   : $acc = $db->id2acc($id)
 Function: Convert ID to Accession
 Returns : Accession
 Args    : Protein cluster ID (as a string)
 Throws  : "Converting ID failed" exception
=cut

sub id2acc {
	my ($self,@args)=@_;
	my $id=shift @args;
	my $acc;
	my $factory = Bio::DB::EUtilities->new(
                        -eutil => 'esummary',
                        -db => 'proteinclusters',
                        -email => $self->email,
                        -id=>[$id],
                       );
	while (my $ds = $factory->next_DocSum) {
	   while (my $item = $ds->next_Item('flattened'))  {
			if ($item->get_name eq "ACCN") {
				if($item->get_content) {
					$acc=$item->get_content;
				}
				last;
			}
		}
	}
	if($acc) {
		return $acc;
	}
	else {
		$self->throw("Converting ID [$id] to Accession failed");
	}
}

=head2 acc2id

 Title   : acc2id
 Usage   : $id = $db->acc2id($acc)
 Function: Convert Accession to ID
 Returns : ID
 Args    : Protein cluster accession (as a string)
 Throws  : "Converting Accession failed" exception
=cut

sub acc2id {
	my ($self,@args)=@_;
	my $acc=shift @args;
	my $factory = Bio::DB::EUtilities->new(
                        -eutil => 'esearch',
                        -db => 'proteinclusters',
                        -email => $self->email,
                        -term=>$acc,
                       );
   my @ids=$factory->get_ids;
	if(@ids) {
		return shift @ids;
	}
	else {
		$self->throw("Converting Accession [$acc] to ID failed");
	}
}

=head2 email

 Title   : email
 Usage   : $dbobj->email('my@foo.bar')
 Function: Assign email address as requested by NCBI
 Returns : Assigned email address
 Args    : Your email address
 Throws  : 
=cut

sub email {
	my ($self,@args)=@_;
	if(@args) {
		$self->{_email}=shift @args;
	}
	return $self->{_email};
}


sub _get_request {
	#get the request in Protein clusters database
	my ($self,$acc,$id)=@_;
	
	my $CGI_location= '/sutils/prkview.cgi';
	my %params= (
		"cluster"=>$acc,
		"clusterid"=>$id,
		"image"=>"download",
		"result"=>"align",
	);

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
		$self->throw("Bio::DB::Align::ProtClustDB Request Error:\n".$request->as_string);
	}
	
	my $alnobj;
	$alnobj=Bio::AlignIO->new(-format=>"fasta",-file=>$tmpfile);
	my $aln=$alnobj->next_aln;


	$aln->accession($acc);
	$aln->id($id);
	unless ($aln->source) {
		$aln->source("ProtClustDB");
	}
	
	return $aln;
}



sub _checkparameter {
	my ($self,@args)=@_;
	return $self->throw_not_implemented();
}

1;


                  