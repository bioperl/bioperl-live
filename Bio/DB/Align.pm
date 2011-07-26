# $Id$
#
# BioPerl module for Bio::DB::Align
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
# Interface from Bio::DB::Align packages

=head1 NAME

Bio::DB::Align - interface fro webagent packages retrieves alignment sequences 
from online databases, e.g. Pfam

=head1 SYNOPSIS
	
	#Take Bio::DB::Align::Pfam as example
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

This package is an interface to Bio::DB::Align packages. See also
Bio::DB::Align::Pfam
Bio::DB::Align::Rfam
Bio::DB::Align::Prosite
Bio::DB::Align::ProtClustDB

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

package Bio::DB::Align;
use strict;
use warnings;
use LWP::UserAgent;
use HTTP::Request;
use Bio::AlignIO;
use Bio::Root::IO;

use base qw(Bio::Root::Root Bio::Root::IO Bio::DB::GenericWebAgent);


=head2 new

 Title   : new
 Usage   : $dbobj = Bio::DB::Align->new(-db=>"Pfam");
             Or, it can be called through specific package
             $dbobj = Bio::DB::Align::Pfam->new();
 Function: Returns a db stream
 Returns : Bio::DB::Align initialised with the appropriate db object
 Args    : -db     database query(case sensitive)
                   Currently Bio::DB::Align supports Pfam, Rfam, 
                   Prosite and ProtClustDB
           -email  email address requested per NCBI policy 
                   (only for Bio::DB::Align::ProtClustDB)
 Note    : 
=cut

sub new {
	my($class,@args) = @_;
	my %param = @args;
	if(defined $param{"-db"}) {
		my $db=$param{"-db"};
		#transform the db name to the right format
		#$db=lc $db;
		#substr($db,0,1)=uc substr($db,0,1);
		
		my $module = "Bio::DB::Align::" . $db;
		my $ok;
	
		eval {
		   $ok = $class->_load_module($module);
		};
		if ( $@ ) {
			$class->throw("$module not supported in BioPerl");
		}
		else {
			return $module->new(@args);
		}
	}
	else {
		$class->throw("-db parameter must be defined");
	}
}

=head2 get_Aln_by_id

 Title   : get_Aln_by_id
 Usage   : $aln = $dbobj->get_Aln_by_id($id)
 Function: Gets a Bio::SimpleAlign object by id
 Returns : a Bio::SimpleAlign object
 Args    : -id  the id as a string
 Note    : 
=cut

sub get_Aln_by_id {
    my ($self,$aln) = @_;
    $self->throw("Sorry, you cannot read from a generic Bio::DB::Align object.");
}


=head2 get_Aln_by_acc

 Title   : get_Aln_by_acc
 Usage   : $seq = $dbobj->get_Aln_by_acc($acc);
 Function: Gets a Bio::SimpleAlign object by accession numbers
 Returns : a Bio::SimpleAlign object
 Args    : -accession  the accession number as a string
 Note    : 
=cut

sub get_Aln_by_acc {
    my ($self,$aln) = @_;
    $self->throw("Sorry, you cannot read from a generic Bio::DB::Align object.");
}

=head2 id2acc

 Title   : id2acc
 Usage   : $acc = $db->id2acc($id)
 Function: Convert id to accession
 Returns : Accession
 Args    : ID (as a string)
 Note    : 
=cut

sub id2acc {
	my ($self,@args)=@_;
	$self->throw("Sorry, you cannot read from a generic Bio::DB::Align object.");
}

=head2 acc2id

 Title   : acc2id
 Usage   : $id = $db->acc2id($acc)
 Function: Convert Accession to ID
 Returns : ID
 Args    : Accession (as a string)
 Note    : 
=cut

sub acc2id {
	my ($self,@args)=@_;
	$self->throw("Sorry, you cannot read from a generic Bio::DB::Align object.");
}


1;


                  