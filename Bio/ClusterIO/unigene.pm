# BioPerl module for Bio::ClusterIO::unigene
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Andrew Macgregor <andrew at cbbc.murdoch.edu.au>
#
# Copyright Andrew Macgregor, Jo-Ann Stanton, David Green
# Molecular Embryology Group, Anatomy & Structural Biology, University of Otago
# http://meg.otago.ac.nz
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

This object reads from Unigene *.data files downloaded from
ftp://ftp.ncbi.nih.gov/repository/UniGene/.  It does not download and
decompress the file, you have to do that yourself.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists.  Your participation is much appreciated.

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
the bugs and their resolution.  Bug reports can be submitted via the
web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHORS - Andrew Macgregor

Email: andrew at cbbc.murdoch.edu.au


=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

#'
# Let the code begin...

package Bio::ClusterIO::unigene;
use strict;

use Bio::Cluster::UniGene;
use Bio::Cluster::ClusterFactory;

use base qw(Bio::ClusterIO);

my %line_is = (
		ID			=>	q/ID\s+(\w{2,3}\.\d+)/,
		TITLE			=>	q/TITLE\s+(\S.*)/,
		GENE			=>	q/GENE\s+(\S.*)/,
		CYTOBAND		=>	q/CYTOBAND\s+(\S.*)/,
		MGI			=>	q/MGI\s+(\S.*)/,
		LOCUSLINK		=>	q/LOCUSLINK\s+(\S.*)/,
		HOMOL		=>	q/HOMOL\s+(\S.*)/,
		EXPRESS			=>	q/EXPRESS\s+(\S.*)/,
		RESTR_EXPR		=>	q/RESTR_EXPR\s+(\S.*)/,
		GNM_TERMINUS		=>	q/GNM_TERMINUS\s+(\S.*)/,
		CHROMOSOME		=>	q/CHROMOSOME\s+(\S.*)/,
		STS			=>	q/STS\s+(\S.*)/,
		TXMAP			=>	q/TXMAP\s+(\S.*)/,
		PROTSIM			=>	q/PROTSIM\s+(\S.*)/,
		SCOUNT			=>	q/SCOUNT\s+(\S.*)/,
		SEQUENCE		=>	q/SEQUENCE\s+(\S.*)/,
		ACC			=>	q/ACC=(\w+)(\.\d+)?/,
		NID			=>	q/NID=\s*(\S.*)/,
		PID			=>	q/PID=\s*(\S.*)/,
		CLONE			=>	q/CLONE=\s*(\S.*)/,
		END			=>	q/END=\s*(\S.*)/,
		LID			=>	q/LID=\s*(\S.*)/,
		MGC			=>	q/MGC=\s*(\S.*)/,
		SEQTYPE		=>	q/SEQTYPE=\s*(\S.*)/,
		TRACE			=>	q/TRACE=\s*(\S.*)/,
		PERIPHERAL		=>	q/PERIPHERAL=\s*(\S.*)/,
		DELIMITER		=>	q{^//},
);

# we set the right factory here
sub _initialize {
	my($self, @args) = @_;

	$self->SUPER::_initialize(@args);
	if(! $self->cluster_factory()) {
	$self->cluster_factory(Bio::Cluster::ClusterFactory->new(
						-type => 'Bio::Cluster::UniGene'));
	}
}

=head2 next_cluster

 Title	 : next_cluster
 Usage	 : $unigene = $stream->next_cluster()
 Function: returns the next unigene in the stream
 Returns : Bio::Cluster::UniGene object
 Args	 : NONE

=cut

sub next_cluster {
	my( $self) = @_;
	local $/ = "\n//";
	return unless my $entry = $self->_readline;
	
# set up the variables we'll need
	my (%unigene,@express,@locuslink,@chromosome,
		@sts,@txmap,@protsim,@sequence);
	my $UGobj;
	
# set up the regexes

# add whitespace parsing and precompile regexes
#foreach (values %line_is) {
#	$_ =~ s/\s+/\\s+/g;
#	print STDERR "Regex is $_\n";
#	#$_ = qr/$_/x;
#}

#$line_is{'TITLE'} = qq/TITLE\\s+(\\S.+)/;

# run each line in an entry against the regexes
	foreach my $line (split /\n/, $entry) {
	  #print STDERR "Wanting to match $line\n";
		if ($line =~ /$line_is{ID}/gcx) {
			$unigene{ID} = $1;
		}
		elsif ($line =~ /$line_is{TITLE}/gcx ) {
		  #print STDERR "MATCHED with [$1]\n";
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
		elsif ($line =~ /$line_is{HOMOL}/gcx) {
			$unigene{HOMOL} = $1;
		}
		elsif ($line =~ /$line_is{EXPRESS}/gcx) {
			my $express = $1;
			# remove initial semicolon if present
			$express =~ s/^;//; 
			@express = split /\s*;/, $express;
		}
		elsif ($line =~ /$line_is{RESTR_EXPR}/gcx) {
			$unigene{RESTR_EXPR} = $1;
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
		elsif ($line =~ /$line_is{SEQUENCE}/gcx) { 
			# parse into each sequence line
			my $seq = {};
			# add unigene id to each seq
			#$seq->{unigene_id} = $unigene{ID}; 
			my @items = split(/;/, $1);
			foreach (@items) {
                            if (/$line_is{ACC}/gcx) {
                                $seq->{acc} = $1;
                                # remove leading dot if version pattern matched
                                $seq->{version} = substr($2,1) if defined $2;
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
                            elsif (/$line_is{SEQTYPE}/gcx) {
                                $seq->{seqtype} = $1;
                            }
                            elsif (/$line_is{TRACE}/gcx) {
                                $seq->{trace} = $1;
                            }
                            elsif (/$line_is{PERIPHERAL}/gcx) {
                                $seq->{peripheral} = $1;
                            }
			}
			push @sequence, $seq;			
		}
		elsif ($line =~ /$line_is{DELIMITER}/gcx) {
			# at the end of the record, add data to the object
			$UGobj = $self->cluster_factory->create_object(
				  -display_id  => $unigene{ID},
				  -description => $unigene{TITLE},
				  -size		   => $unigene{SCOUNT},
				  -members	   => \@sequence);
			$UGobj->gene($unigene{GENE}) if defined ($unigene{GENE});
			$UGobj->cytoband($unigene{CYTOBAND}) if defined($unigene{CYTOBAND});
			$UGobj->mgi($unigene{MGI}) if defined ($unigene{MGI});
			$UGobj->locuslink(\@locuslink);
			$UGobj->homol($unigene{HOMOL}) if defined ($unigene{HOMOL});
			$UGobj->express(\@express);
			$UGobj->restr_expr($unigene{RESTR_EXPR}) if defined ($unigene{RESTR_EXPR});
			$UGobj->gnm_terminus($unigene{GNM_TERMINUS}) if defined ($unigene{GNM_TERMINUS});
			$UGobj->chromosome(\@chromosome);
			$UGobj->sts(\@sts);
			$UGobj->txmap(\@txmap);
			$UGobj->protsim(\@protsim);
		}
	}
	return $UGobj;
}

1;

