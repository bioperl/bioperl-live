# $Id$
#
# BioPerl module for Bio::Cluster::UniGene.pm
#
# Cared for by Andrew Macgregor <andrew@anatomy.otago.ac.nz>
#
# Copyright Andrew Macgregor, Jo-Ann Stanton, David Green
# Molecular Embryology Group, Anatomy & Structural Biology, University of Otago
# http://anatomy.otago.ac.nz/meg
#
# You may distribute this module under the same terms as perl itself
#
# _history
# April 17, 2002 - Initial implementation by Andrew Macgregor
# POD documentation - main docs before the code

=head1 NAME

Bio::Cluster::UniGene - UniGene object

=head1 SYNOPSIS

	use Bio::Cluster::UniGene;
	use Bio::ClusterIO;

	$stream  = Bio::ClusterIO->new('-file' => "Hs.data", 
                                       '-format' => "unigene");
	# note: we quote -format to keep older perl's from complaining.

	while ( my $in = $stream->next_unigene() ) {
		print $in->unigene_id() . "\n";
		while ( my $sequence = $in->next_seq() ) {
			print $sequence->accession_number() . "\n";
		}

=head1 DESCRIPTION

This UniGene object is returned by ClusterIO and contains all the
data associated with one UniGene record.

Available methods (see below for details):

new() - standard new call

unigene_id() - set/get unigene_id

title() - set/get title (description)

gene() - set/get gene

cytoband() - set/get cytoband

locuslink() - set/get locuslink

gnm_terminus() - set/get gnm_terminus

scount() - set/get scount

express() - set/get express, currently takes/returns a reference to an
array of expressed tissues

next_express() - returns the next tissue expression from the expressed
tissue array

chromosome() - set/get chromosome, currently takes/returns a reference
to an array of chromosome lines

next_chromosome() - returns the next chromosome line from the array of
chromosome lines

sts() - set/get sts, currently takes/returns a reference to an array
of sts lines

next_sts() - returns the next sts line from the array of sts lines

txmap() - set/get txmap, currently takes/returns a reference to an
array of txmap lines

next_txmap() - returns the next txmap line from the array of txmap
lines

protsim() - set/get protsim, currently takes/returns a reference to an
array of protsim lines

next_protsim() - returns the next protsim line from the array of
protsim lines

sequence() - set/get sequence, currently takes/returns a reference to
an array of references to seq info

next_seq() - returns a Seq object that currently only contains an
accession number


=head1 FEEDBACK


=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists. Your participation is much appreciated.

  bioperl-l@bioperl.org              - General discussion
  http://bio.perl.org/MailList.html  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via email
or the web:

  bioperl-bugs@bioperl.org
  http://bio.perl.org/bioperl-bugs/

=head1 AUTHOR - Andrew Macgregor

Email andrew@anatomy.otago.ac.nz


=head1 APPENDIX


The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a "_".

=cut

# Let the code begin...


package Bio::Cluster::UniGene;
use vars qw(@ISA $VERSION);
use strict;


use Bio::Root::Root;
use Bio::Seq;
use Bio::Seq::RichSeq;
use Bio::Factory::SequenceStreamI;
use Bio::Seq::SeqFactory;
$VERSION = '1.0';
@ISA = qw(Bio::Root::Root Bio::Factory::SequenceStreamI);


=head2 new

 Title   : new
 Usage   : used by ClusterIO
 Returns : a new Bio::Cluster::Unigene object

=cut

sub new {
    # standard new call..
    my($caller,@args) = @_;
    my $self = $caller->SUPER::new(@args);
    my ($seqfact) = $self->_rearrange([qw(SEQFACTORY)], @args);
    $self->{'_alphabet'} = 'dna';
    if( ! defined $seqfact ) {
	$seqfact = new Bio::Seq::SeqFactory
	    (-verbose => $self->verbose(), 
	     -type => 'Bio::Seq::RichSeq');
    }
    $self->sequence_factory($seqfact);
    return $self;
}


=head2 unigene_id

 Title   : unigene_id
 Usage   : unigene_id();
 Function: Returns the unigene_id associated with the object.
 Example : $id = $unigene->unigene_id or $unigene->unigene_id($id)
 Returns : A string
 Args    : None or an id


=cut

sub unigene_id {
   my ($obj,$value) = @_;

   if( defined $value) {
      $obj->{'unigene_id'} = $value;
    }
   if( ! exists $obj->{'unigene_id'} ) {
       return "$obj";
   }
   return $obj->{'unigene_id'};
}



=head2 title

 Title   : title
 Usage   : title();
 Function: Returns the title associated with the object.
 Example : $title = $unigene->title or $unigene->title($title)
 Returns : A string
 Args    : None or a title


=cut

sub title {
   my ($obj,$value) = @_;

   if( defined $value) {
      $obj->{'title'} = $value;
    }
   if( ! exists $obj->{'title'} ) {
       return "$obj";
   }
   return $obj->{'title'};
}


=head2 gene

 Title   : gene
 Usage   : gene();
 Function: Returns the gene associated with the object.
 Example : $gene = $unigene->gene or $unigene->gene($gene)
 Returns : A string
 Args    : None or a gene


=cut

sub gene {
   my ($obj,$value) = @_;

   if( defined $value) {
      $obj->{'gene'} = $value;
    }
   if( ! exists $obj->{'gene'} ) {
       return "$obj";
   }
   return $obj->{'gene'};
}


=head2 cytoband

 Title   : cytoband
 Usage   : cytoband();
 Function: Returns the cytoband associated with the object.
 Example : $cytoband = $unigene->cytoband or $unigene->cytoband($cytoband)
 Returns : A string
 Args    : None or a cytoband


=cut

sub cytoband {
   my ($obj,$value) = @_;

   if( defined $value) {
      $obj->{'cytoband'} = $value;
    }
   if( ! exists $obj->{'cytoband'} ) {
       return "$obj";
   }
   return $obj->{'cytoband'};
}


=head2 locuslink

 Title   : locuslink
 Usage   : locuslink();
 Function: Returns or stores a reference to an array containing locuslink data.
 Returns : An array reference
 Args    : None or an array reference

=cut

sub locuslink {
   my ($obj,$value) = @_;

   if( defined $value) {
      $obj->{'locuslink'} = $value;
    }
   if( ! exists $obj->{'locuslink'} ) {
       return "$obj";
   }
   return $obj->{'locuslink'};
}


=head2 next_locuslink

 Title   : next_locuslink
 Usage   : next_locuslink();
 Function: Returns the next locuslink from an array referred 
           to using $obj->{'locuslink'}
 Example : 	while ( my $locuslink = $in->next_locuslink() ) {
				print "$locuslink\n";
			}
 Returns : String
 Args    : None

=cut

sub next_locuslink {
	my ($obj) = @_;
	shift @{$obj->{'locuslink'}};
}


=head2 gnm_terminus

 Title   : gnm_terminus
 Usage   : gnm_terminus();
 Function: Returns the gnm_terminus associated with the object.
 Example : $gnm_terminus = $unigene->gnm_terminus or 
           $unigene->gnm_terminus($gnm_terminus)
 Returns : A string
 Args    : None or a gnm_terminus

=cut

sub gnm_terminus {
   my ($obj,$value) = @_;

   if( defined $value) {
      $obj->{'gnm_terminus'} = $value;
    }
   if( ! exists $obj->{'gnm_terminus'} ) {
       return "$obj";
   }
   return $obj->{'gnm_terminus'};
}

=head2 scount

 Title   : scount
 Usage   : scount();
 Function: Returns the scount associated with the object.
 Example : $scount = $unigene->scount or $unigene->scount($scount)
 Returns : A string
 Args    : None or a scount

=cut

sub scount {
   my ($obj,$value) = @_;

   if( defined $value) {
      $obj->{'scount'} = $value;
    }
   if( ! exists $obj->{'scount'} ) {
       return "$obj";
   }
   return $obj->{'scount'};
}



=head2 express

 Title   : express
 Usage   : express();
 Function: Returns or stores a reference to an array containing 
           tissue expression data
 	   This should really only be used by ClusterIO, not directly
 Returns : An array reference
 Args    : None or an array reference

=cut

sub express {
   my ($obj,$value) = @_;

   if( defined $value) {
      $obj->{'express'} = $value;
    }
   if( ! exists $obj->{'express'} ) {
       return "$obj";
   }
   return $obj->{'express'};
}


=head2 next_express

 Title   : next_express
 Usage   : next_express();
 Function: Returns the next tissue from an array referred 
           to using $obj->{'express'}
 Example : 	while ( my $express = $in->next_express() ) {
				print "$express\n";
			}
 Returns : String
 Args    : None

=cut

sub next_express {
	my ($obj) = @_;
	shift @{$obj->{'express'}};
}


=head2 chromosome

 Title   : chromosome
 Usage   : chromosome();
 Function: Returns or stores a reference to an array containing chromosome lines
 		   This should really only be used by ClusterIO, not directly
 Returns : An array reference
 Args    : None or an array reference

=cut

sub chromosome {
   my ($obj,$value) = @_;

   if( defined $value) {
      $obj->{'chromosome'} = $value;
    }
   if( ! exists $obj->{'chromosome'} ) {
       return "$obj";
   }
   return $obj->{'chromosome'};
}


=head2 next_chromosome

 Title   : next_chromosome
 Usage   : next_chromosome();
 Function: Returns the next chromosome line from an array referred to using $obj->{'chromosome'}
 Example : 	while ( my $chromosome = $in->next_chromosome() ) {
				print "$chromosome\n";
			}
 Returns : String
 Args    : None

=cut

sub next_chromosome {
	my ($obj) = @_;
	shift @{$obj->{'chromosome'}};
}



=head2 sts

 Title   : sts
 Usage   : sts();
 Function: Returns or stores a reference to an array containing sts lines
 	   This should really only be used by ClusterIO, not directly
 Returns : An array reference
 Args    : None or an array reference

=cut

sub sts {
   my ($obj,$value) = @_;

   if( defined $value) {
      $obj->{'sts'} = $value;
    }
   if( ! exists $obj->{'sts'} ) {
       return "$obj";
   }
   return $obj->{'sts'};
}


=head2 next_sts

 Title   : next_sts
 Usage   : next_sts();
 Function: Returns the next sts line from an array referred 
           to using $obj->{'sts'}
 Example : 	while ( my $sts = $in->next_sts() ) {
				print "$sts\n";
			}
 Returns : String
 Args    : None

=cut

sub next_sts {
	my ($obj) = @_;
	shift @{$obj->{'sts'}};
}


=head2 txmap

 Title   : txmap
 Usage   : txmap();
 Function: Returns or stores a reference to an array containing txmap lines
	   This should really only be used by ClusterIO, not directly
 Returns : An array reference
 Args    : None or an array reference

=cut

sub txmap {
   my ($obj,$value) = @_;

   if( defined $value) {
      $obj->{'txmap'} = $value;
    }
   if( ! exists $obj->{'txmap'} ) {
       return "$obj";
   }
   return $obj->{'txmap'};
}


=head2 next_txmap

 Title   : next_txmap
 Usage   : next_txmap();
 Function: Returns the next txmap line from an array 
           referred to using $obj->{'txmap'}
 Example : 	while ( my $tsmap = $in->next_txmap() ) {
				print "$txmap\n";
			}
 Returns : String
 Args    : None

=cut

sub next_txmap {
	my ($obj) = @_;
	shift @{$obj->{'txmap'}};
}


=head2 protsim

 Title   : protsim
 Usage   : protsim();
 Function: Returns or stores a reference to an array containing protsim lines
	   This should really only be used by ClusterIO, not directly
 Returns : An array reference
 Args    : None or an array reference

=cut

sub protsim {
   my ($obj,$value) = @_;

   if( defined $value) {
      $obj->{'protsim'} = $value;
    }
   if( ! exists $obj->{'protsim'} ) {
       return "$obj";
   }
   return $obj->{'protsim'};
}


=head2 next_protsim

 Title   : next_protsim
 Usage   : next_protsim();
 Function: Returns the next protsim line from an array referred 
           to using $obj->{'protsim'}
 Example : 	while ( my $protsim = $in->next_protsim() ) {
				print "$protsim\n";
			}
 Returns : String
 Args    : None

=cut

sub next_protsim {
	my ($obj) = @_;
	shift @{$obj->{'protsim'}};
}



=head2 sequence

 Title   : sequence
 Usage   : sequence();
 Function: Returns or stores a reference to an array containing sequence data
 	   This should really only be used by ClusterIO, not directly
 Returns : An array reference
 Args    : None or an array reference

=cut

sub sequence {
   my ($obj,$value) = @_;

   if( defined $value) {
      $obj->{'sequence'} = $value;
    }
   if( ! exists $obj->{'sequence'} ) {
       return "$obj";
   }
   return $obj->{'sequence'};
}


=head2 next_seq

 Title   : next_seq
 Usage   : next_seq();
 Function: Returns the next seq as a Seq object as defined by 
           $seq->sequence_factory(), 
           at present an empty Bio::Seq::RichSeq object with 
           just the accession_number() and pid() set
 Example :  while ( my $sequence = $in->next_seq() ) {
             print $sequence->accession_number() . "\n";
	    }
 Returns : Bio::PrimarySeqI object
 Args    : None

=cut

sub next_seq {
    my ($obj) = @_;
    return unless (my $seq = shift @{$obj->{'sequence'}});
    my $seqobj = $obj->sequence_factory->create_sequence
	( -accession_number => $seq->{acc},
	  -pid => $seq->{pid},
	  -id => $seq->{acc},
	  -desc => join(' ', map { uc($_) ."=". $seq->{$_}} sort keys %{$seq} ));
    return $seqobj;
}

=head2 sequence_factory

 Title   : sequence_factory
 Usage   : $seqio->sequence_factory($seqfactory)
 Function: Get/Set the Bio::Factory::SequenceFactoryI
 Returns : Bio::Factory::SequenceFactoryI
 Args    : [optional] Bio::Factory::SequenceFactoryI


=cut

sub sequence_factory {
    my ($self,$obj) = @_;   
    if( defined $obj ) {
	if( ! ref($obj) || ! $obj->isa('Bio::Factory::SequenceFactoryI') ) {
	    $self->throw("Must provide a valid Bio::Factory::SequenceFactoryI object to ".ref($self)." sequence_factory()");
	}
	$self->{'_seqfactory'} = $obj;
    }
    $self->{'_seqfactory'};
}

# keep AUTOLOAD happy
sub DESTROY {
    my ($self) = @_;
}

1;
