# $Id$
#
# BioPerl module for Bio::DB::AlignIO::Rfam
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

Bio::DB::AlignIO::Rfam - webagent which interacts with and retrieves alignment 
sequences from Rfam

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

package Bio::DB::Align::Rfam;
use strict;
use warnings;
use LWP::UserAgent;
use HTTP::Request;
use Bio::AlignIO;
use Bio::SimpleAlign;
use Bio::Root::IO;
use vars qw(%FORMATS %ALNTYPE $HOSTBASE);

use base qw(Bio::Root::Root Bio::Root::IO Bio::DB::GenericWebAgent);

BEGIN {
	$HOSTBASE = 'http://rfam.sanger.ac.uk';
	%FORMATS=qw(fasta 1 stockholm 1 pfam 1 fastau 1); #supported formats in pfam
	%ALNTYPE=qw(seed 1 full 1); #supported alignment types
}


=head2 new

 Title   : new
 Usage   : $dbobj = Bio::DB::Align->new(-db=>"Rfam");
             Or, it can be called through specific package
             $dbobj = Bio::DB::Align::Rfam->new();
 Function: Returns a Bio::DB::Align::Rfam stream
 Returns : Bio::DB::Align::Rfam object
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
 Usage   : $aln = $db->get_Aln_by_id($id)
 Function: Gets a Bio::SimpleAlign object by id
 Returns : a Bio::SimpleAlign object
 Args    : 
 Note    : 
 Throws  : "Bio::Root::NotImplemented" exception
=cut

sub get_Aln_by_id {
	my ($self,@args)=@_;
	$self->throw_not_implemented();
}


=head2 get_Aln_by_acc

 Title   : get_Aln_by_acc
 Usage   : $seq = $db->get_Aln_by_acc("RF00360");
 Function: Gets a Bio::SimpleAlign object by accession numbers
 Returns : a Bio::SimpleAlign object
 Args    : -accession  the accession number as a string
           -alignment  Seed(default) or Full
           -format     the output format from Rfam. This will decide 
                       which package to use in the Bio::AlignIO
                       possible options can be fasta (default), 
                       stockholm or pfam
           -nselabel   0 (default)   Label by species name
                       1             Label by name/start-end
  	        -gap        dashes (default) "-" as gap char
	                    dots           "." as gap char
	                    none          Unaligned
  
 Note    : 
 Throws  : "Bio::DB::Align::Rfam Request Error" exception
=cut

sub get_Aln_by_acc {
	my ($self,@args)=@_;
	
	my ($acc,$alignment,$format, $nselabel,$gap)=$self->_checkparameter(@args);
	
	my $CGI_location= '/family/alignment/download/format';
	
	my %params= (
		"acc"=>$acc,
		"alnType"=>$alignment,
		"format"=>$format,
		"nseLabels"=>$nselabel,
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
		$self->throw("Bio::DB::Align::Rfam Request Error:\n".$request->as_string);
	}
	
	
	my $alnobj;
	unless($format eq "fastau") {
		$alnobj=Bio::AlignIO->new(-format=>$format,-file=>$tmpfile);
	}
	else {
		$alnobj=Bio::AlignIO->new(-format=>"fasta",-file=>$tmpfile);
	}
	my $aln=$alnobj->next_aln;
	
	#gap char
	if($gap eq "dashes" ) {
		$aln->gap_char("-");
		$aln->map_chars('\.','-');#the default gap_char in Rfam is ".", need to be changed to "-" upon request
	}
	elsif($gap eq "dots") {
		if($format eq "fastau") {
			$aln->map_chars('-','\.');#The panding gap char can only be "-", thus need to be substituted to "."
		}
		$aln->gap_char(".");
	}
	elsif($gap eq "none") {
		$aln->gap_char("-");
	}
	
	#record the accession and id
	unless ($aln->accession) {
		$aln->accession($acc);
	}
	unless($aln->source) {
		$aln->source("Rfam");
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
	my ($acc,$alignment,$format, $nselabel, $gap)=$self->_rearrange([qw(ACCESSION ALIGNMENT FORMAT NSELABEL GAP)],@args);
	#check alntype
	if($alignment) {
		$alignment=lc $alignment;
		unless(defined $ALNTYPE{$alignment}) {
			$self->throw("Only [".join(" ",sort keys %ALNTYPE)."] are supported by Rfam, not [".$alignment."]");
		}
	}
	else {
		$alignment="seed"; #default
	}
	#check format
	if($format) {
		$format=lc $format;
		unless(defined $FORMATS{$format}) {
			$self->throw("Only [".join(" ",sort keys %FORMATS)."] are supported by Rfam, not [".$format."]");
		}
	}
	else {
		$format="fasta"; #default format
	}
	#check nselabel
	if($nselabel) {
		unless($nselabel==0 || $nselabel==1) {
			$self->throw("Only nseLabel [0/1] are supported by Rfam, not [".$nselabel."]");
		}
	}
	else {
		$nselabel=0; #default value
	}
	#check gaps
	if($gap) {
		if($gap=~/dashes/i||$gap eq "-") {
			$gap="dashes";
		}
		elsif($gap=~/dots/i||$gap eq ".") {
			$gap="dots";
		}
		elsif ($gap=~/none/i) {
			$gap="none"; #ungapped
			$format="fastau";
		}
		else {
			$gap="dashes";
		}
	}
	else {
		$gap="dashes";
	}
		
	return ($acc,$alignment,$format, $nselabel, $gap);
}

1;


                  