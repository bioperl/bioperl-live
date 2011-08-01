# $Id$
#
# BioPerl module for Bio::DB::Align::Pfam
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

Bio::DB::Align::Pfam - webagent which interacts with and retrieves alignment 
sequences from Pfam

=head1 SYNOPSIS

	use Bio::DB::Align;
	use Bio::DB::Align::Pfam;
	use Bio::SimpleAlign; 
  
	my $dbobj=Bio::DB::Align::Pfam->new(); #create a db object
	my $dbobj2=Bio::DB::Align->new(-db=>"Pfam"); #create a db object
	                                           #the same with above

	#retrieve a Bio::SimpleAlign object
	my $aln=$dbobj->get_Aln_by_id("Piwi"); 
	
	#do something here with the align object
	print $aln->length,"\n";
	print $aln->num_sequences,"\n";	
	
	foreach my $seq ($aln->next_Seq) {
		#do something with $seq
	}
	
	#or parameter based calling
	my $aln2=$dbobj->get_Aln_by_acc(-accession=>"PF02171",
	                                -alignment=>"full",
	                                -format=>"stockholm",
	                                -order=>"a",
	                                -case=>"u",
	                                -gap=>"dots");
	print $aln2->id,"\n";
	print $aln2->accession,"\n";
	
	#id accession conversion
	my $acc=$dbobj->id2acc("Piwi");
	print $acc,"\n";
	my $id=$dbobj->acc2id("PF02171");
	print $id,"\n";


=head1 DESCRIPTION

This package uses the RESTful service provided by Pfam. It retrieves
online alignment sequences and save it into Bio::SimpleAlign object.

The retrieval can be based on protein domain id, e.g. "Piwi", or 
accession number, e.g. "PF02171".

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

package Bio::DB::Align::Pfam;
use strict;
use warnings;
use LWP::UserAgent;
use HTTP::Request;
use Bio::AlignIO;
use Bio::Root::IO;
use vars qw(%FORMATS %ALNTYPE $HOSTBASE);

use base qw(Bio::Root::Root Bio::Root::IO Bio::DB::GenericWebAgent);

BEGIN {
	$HOSTBASE = 'http://pfam.sanger.ac.uk';
	%FORMATS=qw(fasta 1 stockholm 1 selex 1 msf 1); #supported formats in pfam
	%ALNTYPE=qw(seed 1 full 1 ncbi 1 metagenomics 1); #supported alignment types
}


=head2 new

 Title   : new
 Usage   : $dbobj = Bio::DB::Align->new(-db=>"Pfam");
             Or, it can be called through specific package
             $dbobj = Bio::DB::Align::Pfam->new();
 Function: Returns a Bio::DB::Align::Pfam stream
 Returns : Bio::DB::Align::Pfam object
 Args    : 
 Note    : 
=cut


sub new {
	my($class,@args) = @_;
	my $self = $class->SUPER::new(@args);
   my $ua = new LWP::UserAgent(env_proxy => 1);
   #$ua->agent(ref($self) ."/$MODVERSION");
   $self->ua($ua);  
   $self->{'_authentication'} = [];	
	return $self;	
}


=head2 get_Aln_by_id

 Title   : get_Aln_by_id
 Usage   : $aln = $dbobj->get_Aln_by_id('Piwi')
 Function: This method uses Pfam id conversion service id2acc to convert 
           id to accession. Then, it gets a Bio::SimpleAlign object 
	          using get_Aln_by_acc
 Returns : a Bio::SimpleAlign object
 Args    : -id  the id as a string
          -alignment  Seed(default), Full, NCBI or metagenomics
          -format     the output format from Pfam. This will decide which
                      package to use in the Bio::AlignIO
                      possible options can be fasta (default), stockholm, 
                      selex and MSF
          -order      t (default)   Order by tree 
                      a             Order alphabetically
          -case       i (default)   Inserts lower case  
                      a             All upper case
          -gap        dashes (default) "-" as gap char
                      dots           "." as gap char
                      mixed         "-" and "." mixed 
                                    (not recommended, this may cause bug 
                                    in BioPerl)
                      none          Unaligned
 
 
 Note    : 
 Throws  : "Bio::DB::Align::Pfam Request Error" exception
=cut

sub get_Aln_by_id {
	my ($self,@args)=@_;
	my ($id,$alignment,$format, $order, $case, $gap)=$self->_rearrange([qw(ID ALIGNMENT FORMAT ORDER CASE GAP)],@args);
	#id 2 accession convertion
	my $acc=$self->id2acc($id);
	
	#give new -accession argument
	return $self->get_Aln_by_acc($acc,$alignment,$format, $order, $case, $gap);
}

=head2 get_Aln_by_acc

 Title   : get_Aln_by_acc
 Usage   : $aln = $dbobj->get_Aln_by_acc("PF02171");
 Function: Gets a Bio::SimpleAlign object by accession numbers
  Returns : a Bio::SimpleAlign object
 Args    : -accession  the accession number as a string
           -alignment  Seed(default), Full, NCBI or metagenomics
           -format     the output format from Pfam. This will decide 
                       which package to use in the Bio::AlignIO
                       possible options can be fasta (default), 
                       stockholm, selex and MSF
           -order      t (default)   Order by tree 
                       a             Order alphabetically
           -case       l (default)   Inserts lower case  
                       u             All upper case
           -gap        dashes (default) "-" as gap char
                       dots           "." as gap char
                       mixed         "-" and "." mixed 
                                     (not recommended, this may cause 
                                     bug in BioPerl)
                       none          Unaligned
  
  
  Note    : 
  Throws  : "Bio::DB::Align::Pfam Request Error" exception
=cut

sub get_Aln_by_acc {
	my ($self,@args)=@_;
	
	my ($acc,$alignment,$format, $order, $case, $gap)=$self->_checkparameter(@args);
	
	my $CGI_location= '/family/alignment/download/format';
	
	my %params= (
		"acc"=>$acc,
		"alnType"=>$alignment,
		"format"=>$format,
		"order"=>$order,
		"case"=>$case,
		"gaps"=>$gap,
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
		$self->throw("Bio::DB::Align::Pfam Request Error:\n".$request->as_string);
	}
	
	my $alnobj=Bio::AlignIO->new(-format=>$format,-file=>$tmpfile);
	my $aln=$alnobj->next_aln;
	
	#gap char
	if($gap eq "dashes") {
		$aln->gap_char("-");
	}
	elsif($gap eq "dots") {
		$aln->gap_char(".");
	}
	elsif($gap eq "none") {
		$aln->gap_char("-");
	}
	
	#record the accession and id
	unless($aln->accession) {
		$aln->accession($acc);
	}
	unless($aln->id && $aln->id ne "NoName") {
		$aln->id($self->acc2id($acc));
	}
	unless($aln->source) {
		$aln->source("Pfam");
	}
	
	return $aln;
}

=head2 id2acc

 Title   : id2acc
 Usage   : $acc = $dbobj->id2acc($id)
 Function: Convert ID to Accession
 Returns : Accession
 Args    : ID (as a string)
 Throws  : "Converting ID failed" exception
=cut

sub id2acc {
	my ($self,@args)=@_;
	my $id=shift @args;
	
	my $CGI_location= '/family/acc/';
	
	my $url = URI->new($HOSTBASE . $CGI_location. $id);
	
	my $request = $self->ua->get($url);
	
	if($request->is_success && $request->content ne "No such family") {
		return $request->content;
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
 Args    : Accession (as a string)
 Throws  : "Converting Accession failed" exception
=cut

sub acc2id {
	my ($self,@args)=@_;
	my $acc=shift @args;
	
	my $CGI_location= '/family/id/';
	
	my $url = URI->new($HOSTBASE . $CGI_location. $acc);
	
	my $request = $self->ua->get($url);
	
	if($request->is_success && $request->content ne "No such family") {
		return $request->content;
	}
	else {
		$self->throw("Converting Accession [$acc] to ID failed");
	}
}


sub _checkparameter {
	my ($self,@args)=@_;
	my ($acc,$alignment,$format, $order, $case, $gap)=$self->_rearrange([qw(ACCESSION ALIGNMENT FORMAT ORDER CASE GAP)],@args);
	#check alntype
	if($alignment) {
		$alignment=lc $alignment;
		unless(defined $ALNTYPE{$alignment}) {
			$self->throw("Only [".join(" ",sort keys %ALNTYPE)."] are supported by Pfam, not [".$alignment."]");
		}
	}else {
		$alignment="seed";
	}
	#check format
	if($format) {
		$format=lc $format;
		unless(defined $FORMATS{$format}) {
			$self->throw("Only [".join(" ",sort keys %FORMATS)."] are supported by Pfam, not [".$format."]");
		}
	}
	else {
		$format="fasta"; #default format
	}
	#check order
	if($order&&$order=~/^a/i) {
		$order="a";
	}
	else {
		$order="t"; #default value
	}
	#check case
	if($case && $case=~/^u/i) {
		$case="u";
	}
	else {
		$case="l";
	}
	#check gaps
	if($gap) {
		if($gap=~/dashes/i||$gap eq "-") {
			$gap="dashes";
		}
		elsif($gap=~/dots/i||$gap eq ".") {
			$gap="dots";
		}
		elsif($gap=~/mixed/i) {
			$gap="default";
		}
		elsif ($gap=~/none/i) {
			$gap="none";
		}
		else {
			$gap="dashes";
		}
	}
	else {
		$gap="dashes";
	}
		
	return ($acc,$alignment,$format, $order, $case, $gap);
}

1;
