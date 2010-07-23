# $Id$
#
# BioPerl module for Bio::DB::AlignIO::Pfam
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

Bio::DB::AlignIO::Pfam - webagent which interacts with and retrieves alignment 
sequences from Pfam

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

package Bio::DB::Align::Pfam;
use strict;
use warnings;
use LWP::UserAgent;
use HTTP::Request;
use Bio::AlignIO;
use Bio::Root::IO;
use vars qw(%FORMATS %ALNTYPE $HOSTBASE $CGILOCATION $CGILOCATION_CONV);

use base qw(Bio::Root::Root Bio::Root::IO Bio::DB::GenericWebAgent);

BEGIN {
	$HOSTBASE = 'http://pfam.sanger.ac.uk';
	$CGILOCATION= '/family/alignment/download/format';
	$CGILOCATION_CONV= '/family/acc/';
	%FORMATS=qw(fasta 1 stockholm 1 selex 1 msf 1); #supported formats in pfam
	%ALNTYPE=qw(seed 1 full 1 ncbi 1 metagenomics 1); #supported alignment types
}

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
	Usage   : $aln = $db->get_Aln_by_id('Piwi')
	Function: This method uses Pfam id conversion service to convert id to accession
	          Then, it gets a Bio::SimpleAlign object using get_Aln_by_acc
	Returns : a Bio::SimpleAlign object
	Args    : -id  the id as a string
	         -alignment  Seed(default), Full, NCBI or metagenomics
	         -format     the output format from Pfam. This will decide which
	                     package to use in the Bio::AlignIO
	                     possible options can be fasta (default), stockholm, selex and MSF
	         -order      t (default)   Order by tree 
	                     a             Order alphabetically
	         -case       i (default)   Inserts lower case  
	                     a             All upper case
	         -gap        dashes (default) "-" as gap char
	                     dots           "." as gap char
	                     mixed         "-" and "." mixed 
	                                   (not recommended, this may cause bug in BioPerl)
	                     none          Unaligned
	
	
	Note    : 
	Throws  : "Bio::DB::Align::Pfam Request Error" exception
=cut

sub get_Aln_by_id {
	my ($self,@args)=@_;
	my $id=$self->_rearrange(["ID"],@args);
	
	#id 2 accession convertion
	my $acc=$self->id2acc($id);
	
	#give new -accession argument
	push @args,"-accession",$acc;
	return $self->get_Aln_by_acc(@args);
}


=head2 id2acc

 Title   : id2acc
 Usage   : $acc = $db->id2acc('Piwi')
 Function: Convert id to accession
 Returns : Accession
 Args    : the id (as a string) of a sequence for the alignment
 Throws  : "Bio::DB::Align::Pfam Request Error" exception
=cut

sub id2acc {
	my ($self,@args)=@_;
	my $id=shift @args;
	
	my $url = URI->new($HOSTBASE . $CGILOCATION_CONV . $id);
	
	my $request = $self->ua->get($url);
	
	if($request->is_success) {
		return $request->content;
	}
	else {
		$self->throw("Bio::DB::Align::Pfam Request Error:\n",$request->to_string);
	}
}

=head2 get_Aln_by_acc

  Title   : get_Aln_by_acc
  Usage   : $seq = $db->get_Aln_by_acc($acc);
  Function: Gets a Bio::SimpleAlign object by accession numbers
  Returns : a Bio::SimpleAlign object
  Args    : -accession  the accession number as a string
            -alignment  Seed(default), Full, NCBI or metagenomics
            -format     the output format from Pfam. This will decide which
                        package to use in the Bio::AlignIO
                        possible options can be fasta (default), stockholm, selex and MSF
            -order      t (default)   Order by tree 
                        a             Order alphabetically
            -case       i (default)   Inserts lower case  
                        a             All upper case
            -gap        dashes (default) "-" as gap char
                        dots           "." as gap char
                        mixed         "-" and "." mixed 
                                      (not recommended, this may cause bug in BioPerl)
                        none          Unaligned
  
  
  Note    : 
  Throws  : "Bio::DB::Align::Pfam Request Error" exception
=cut

sub get_Aln_by_acc {
	my ($self,@args)=@_;
	
	my ($acc,$alignment,$format, $order, $case, $gap)=$self->_checkparameter(@args);
	
	my %params= (
		"acc"=>$acc,
		"alnType"=>$alignment,
		"format"=>$format,
		"order"=>$order,
		"case"=>$case,
		"gaps"=>$gap,
		);

	my $url = URI->new($HOSTBASE . $CGILOCATION);
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
	
	return $aln;
}


sub _checkparameter {
	my ($self,@args)=@_;
	my ($acc,$alignment,$format, $order, $case, $gap)=$self->_rearrange([qw(ACCESSION ALIGNMENT FORMAT ORDER CASE GAP)],@args);
	#check alntype
	if($alignment) {
		$alignment=lc $alignment;
		unless(defined $ALNTYPE{$alignment}) {
			$self->throw("Only ".join(" ",sort keys %ALNTYPE)." are supported by Pfam, not [".$alignment."]");
		}
	}else {
		$alignment="seed";
	}
	#check format
	if($format) {
		$format=lc $format;
		unless(defined $FORMATS{$format}) {
			$self->throw("Only ".join(" ",sort keys %FORMATS)." are supported by Pfam, not [".$format."]");
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
	if($case && $case=~/^a/i) {
		$case="a";
	}
	else {
		$case="i";
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


                  