# $Id$
# BioPerl module for Bio::ClusterIO::unigene
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

Bio::ClusterIO::unigene - UniGene input stream

=head1 SYNOPSIS

Do not use this module directly.  Use it via the Bio::ClusterIO class.

=head1 DESCRIPTION

This object reads from Unigene *.data files downloaded from ftp://ftp.ncbi.nih.gov/repository/UniGene/.
It doesn't download and decompress the file, you have to do that yourself.


=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists.  Your participation is much appreciated.

  bioperl-l@bioperl.org            - General discussion
  http://bioperl.org/MailList.shtml - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
 the bugs and their resolution.
 Bug reports can be submitted via email or the web:

  bioperl-bugs@bio.perl.org
  http://bio.perl.org/bioperl-bugs/

=head1 AUTHORS - Andrew Macgregor

Email: andrew@anatomy.otago.ac.nz


=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

#'
# Let the code begin...

package Bio::ClusterIO::unigene;
use vars qw(@ISA);
use strict;

use Bio::ClusterIO;
use Bio::Cluster::Unigene;

@ISA = qw(Bio::ClusterIO);

=head2 next_cluster

 Title   : next_cluster
 Usage   : $unigene = $stream->next_cluster()
 Function: returns the next unigene in the stream
 Returns : Bio::Cluster::UniGene object
 Args    : NONE

=cut

sub next_cluster {
	my( $self) = @_;
	local $/ = "//";
	return unless my $entry = $self->_readline;
	
	my $UGobj = Bio::Cluster::UniGene->new();

# set up the variables we'll need
my (%unigene,@express,@locuslink,@chromosome,@sts,@txmap,@protsim,@sequence);

# set up the regexes
my $data = qr/(?:.+)/;
my $num  = qr/(?:\d+)/;
my $organism = qr/(?:At|Bt|Dr|Hs|Hv|Mm|Os|Rn|Ta|Xl|Zm)/;

my %line_is = (
			ID				=> 	qr/ID ($organism\.$num)/,
			TITLE			=>	qr/TITLE ($data)/,
			GENE			=>	qr/GENE ($data)/,
			CYTOBAND		=>	qr/CYTOBAND ($data)/,
			MGI				=>	qr/MGI ($data)/,
			LOCUSLINK		=>	qr/LOCUSLINK ($data)/,
			EXPRESS			=>	qr/EXPRESS ($data)/,
			GNM_TERMINUS	=>	qr/GNM_TERMINUS ($data)/,
			CHROMOSOME		=>	qr/CHROMOSOME ($data)/,
			STS				=>	qr/STS ($data)/,
			TXMAP			=>	qr/TXMAP ($data)/,
			PROTSIM			=>	qr/PROTSIM ($data)/,
			SCOUNT			=>	qr/SCOUNT ($num)/,
			SEQUENCE		=>	qr/SEQUENCE ($data)/,
			ACC				=>	qr/ACC=($data)/,
			NID				=>	qr/NID=($data)/,
			PID				=>	qr/PID=($data)/,
			CLONE			=>	qr/CLONE=($data)/,
			END				=>	qr/END=($data)/,
			LID				=>	qr/LID=($data)/,
			MGC				=>	qr/MGC=($data)/,
			DELIMITER		=>	qr/^\/\//
);


foreach (values %line_is) {
	$_ =~ s/\s+/\\s+/g;
	$_ = qr/$_/x;
}



# run each line in an entry against the regexes
	foreach my $line (split /\n/, $entry) {
		if ($line =~ /$line_is{ID}/gcx) {
				$unigene{ID} = $1;
		}
		elsif ($line =~ /$line_is{TITLE}/gcx) {
				$unigene{TITLE} = $1;
		}
		elsif ($line =~ /$line_is{GENE}/gcx) {
				$unigene{GENE} = $1;
		}
		elsif ($line =~ /$line_is{CYTOBAND}/gcx) {
				$unigene{CYTOBAND} = $1;
		}
		elsif ($line =~ /$line_is{MGI}/gcx) {
				$unigene{MGI} = $1;
		}
		elsif ($line =~ /$line_is{LOCUSLINK}/gcx) {
				@locuslink = split /;/, $1;
		}
		elsif ($line =~ /$line_is{EXPRESS}/gcx) {
				my $express = $1;
				$express =~ s/^;//;	# remove initial semicolon if present
				@express = split /;/, $express;
		}
		elsif ($line =~ /$line_is{GNM_TERMINUS}/gcx) {
				$unigene{GNM_TERMINUS} = $1;
		}
		elsif ($line =~ /$line_is{CHROMOSOME}/gcx) {
				push @chromosome, $1;
		}
		elsif ($line =~ /$line_is{TXMAP}/gcx) {
				push @txmap, $1;
		}
		elsif ($line =~ /$line_is{STS}/gcx) {
				push @sts, $1;
		}
		elsif ($line =~ /$line_is{PROTSIM}/gcx) {
				push @protsim, $1;
		}
		elsif ($line =~ /$line_is{SCOUNT}/gcx) {
				$unigene{SCOUNT} = $1;
		}
		elsif ($line =~ /$line_is{SEQUENCE}/gcx) {		# parse into each sequence line
				my $seq = {};
				$seq->{unigene_id} = $unigene{ID}; # add unigene id to each seq
				my @items = split /;/,$1;
				foreach (@items) {
					if (/$line_is{ACC}/gcx) {
						$seq->{acc} = $1;
					}
					elsif (/$line_is{NID}/gcx) {
						$seq->{nid} = $1;
					}
					elsif (/$line_is{PID}/gcx) {
						$seq->{pid} = $1;
					}
					elsif (/$line_is{CLONE}/gcx) {
						$seq->{clone} = $1;
					}
					elsif (/$line_is{END}/gcx) {
						$seq->{end} = $1;
					}
					elsif (/$line_is{LID}/gcx) {
						$seq->{lid} = $1;
					}
					elsif (/$line_is{MGC}/gcx) {
						$seq->{mgc} = $1;
					}
				}
				push @sequence, $seq;			
		}
		elsif ($line =~ /$line_is{DELIMITER}/gcx) {		# at the end of the record, add data to the object
					$UGobj->unigene_id($unigene{ID});
					$UGobj->title($unigene{TITLE});
					if ( defined ($unigene{GENE}) ) { $UGobj->gene($unigene{GENE}) };
					if ( defined ($unigene{CYTOBAND}) ) { $UGobj->cytoband($unigene{CYTOBAND}) };
					if ( defined ($unigene{MGI}) ) { $UGobj->mgi($unigene{MGI}) };
					$UGobj->locuslink(\@locuslink);
					$UGobj->express(\@express);
					if ( defined ($unigene{GNM_TERMINUS}) ) { $UGobj->gnm_terminus($unigene{GNM_TERMINUS}) };
					$UGobj->chromosome(\@chromosome);
					$UGobj->sts(\@sts);
					$UGobj->txmap(\@txmap);
					$UGobj->protsim(\@protsim);
					$UGobj->scount($unigene{SCOUNT});
					$UGobj->sequence(\@sequence);
					
					
					
		}
	}
	return $UGobj;
}

1;

