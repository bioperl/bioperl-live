# $Id$
# BioPerl module for Bio::SeqIO::tigr
#
# Cared for by Josh Lauricha (laurichj@bioinfo.ucr.edu)
#
# Copyright Josh Lauricha
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::SeqIO::tigr - TIGR XML sequence input/output stream

=head1 SYNOPSIS

Do not use this module directly.  Use it via the Bio::SeqIO class.

=head1 DESCRIPTION

This object can transform Bio::Seq objects to and from efa flat
file databases.

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

=head1 AUTHORS - Josh Lauricha

Email: laurichj@bioinfo.ucr.edu


=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

# TODO:
#  - Clean up code
#  - Find and fix bugs ;)

# Let the code begin...
package Bio::SeqIO::tigr;
use vars qw(@ISA);
use strict;

use Bio::SeqIO;
use Bio::Seq::RichSeq;
use Bio::Species;
use Bio::Annotation::Comment;
use Bio::SeqFeature::Generic;
use Bio::Seq::SeqFactory;
use Bio::Seq::RichSeq;
use Data::Dumper;

@ISA = qw(Bio::SeqIO);

sub _initialize
{
	my($self, @args) = @_;

	$self->SUPER::_initialize(@args);
	($self->{'tu'}, $self->{'model'})
	 = $self->_rearrange([qw(TU MODEL)], @args);
	$self->{'tu'}      = 1 if !defined($self->{'tu'});
	$self->{'model'}   = 1 if !defined($self->{'model'});
	$self->{'eof'}     = 0;
	$self->sequence_factory(new Bio::Seq::SeqFactory(
			-type => 'Bio::Seq::RichSeq'));

	$self->_process();
	my $list=$self->{'ASSEMBLY'}->{'GENE_LIST'}->{'PROTEIN_CODING'}->{'TU'};
}

=head2 next_seq

 Title   : next_seq
 Usage   : $seq = $stream->next_seq()
 Function: returns the next sequence in the stream
 Returns : Bio::Seq object
 Args    : NONE

=cut

sub next_seq()
{
	my ($self) = @_;
	my $list=$self->{'ASSEMBLY'}->{'GENE_LIST'}->{'PROTEIN_CODING'}->{'TU'};

	return 0 if $self->{'eof'} == 1;

	if(defined(@{$self->{'seq'}}) and @{$self->{'seq'}} > 0) {
		return shift @{$self->{'seq'}};
	}

	while($self->{'eof'} == 0) {
		until((defined(@$list) and @$list > 0) or $self->{'eof'} == 1) {
			$self->_process();
		}

		while(defined(@$list) and @$list > 0 and
			(!defined(@{$self->{'seq'}}) or @{$self->{'seq'}} == 0)
		) {
			$self->_process_seq();
		}

		if(defined(@{$self->{'seq'}}) and @{$self->{'seq'}} > 0) {
			return shift @{$self->{'seq'}};
		}
	}

	return undef;
}

sub _get_sequence
{
	my($self, $start, $end) = @_;
	my $dir   = ($start < $end)? 1 : -1;
	my $len   = ($end - $start) * $dir + 1;
	my $seqstr;


	if($dir == -1) {
		$seqstr = reverse(substr(
			$self->{'ASSEMBLY'}->{'ASSEMBLY_SEQUENCE'},
			$end - 1, $len));
		$seqstr =~ tr/acgtrymkswhbvdnxACGTRYMKSWHBVDNX/tgcayrkmswdvbhnxTGCAYRKMSWDVBHNX/;
	} else {
		$seqstr = substr(
			$self->{'ASSEMBLY'}->{'ASSEMBLY_SEQUENCE'},
			$start - 1, $len);
	}
	return $seqstr;
}

sub _process_seq
{
	my ($self) = @_;
	my $list=$self->{'ASSEMBLY'}->{'GENE_LIST'}->{'PROTEIN_CODING'}->{'TU'};
	my $kw = join(' ', @{$self->{'ASSEMBLY'}->{'HEADER'}->{'KEYWORDS'}})
		if defined($self->{'ASSEMBLY'}->{'HEADER'}->{'KEYWORDS'});

	my $tu = shift(@$list);
	push(@{$self->{'seq'}}, $self->_process_seq_tu($tu))
		if defined($self->{'tu'});
}

sub _process_seq_tu()
{
	my ($self, $tu) = @_;
	my $start  = $tu->{'COORDSET'}->{'END5'};
	my $end    = $tu->{'COORDSET'}->{'END3'};
	my $seqstr = $self->_get_sequence($start, $end);
	
	my $id = $tu->{GENE_INFO}->{PUB_LOCUS};
	if(!defined($id) or $id =~ /^\s*$/o) {
		$id = $tu->{FEAT_NAME};
	}

	my $seq = $self->sequence_factory()->create(
		-id => $id,
		-seq => $seqstr,
		-accession_number => $id
	);

	$seq->primary_id($id);
	$seq->add_date($tu->{'DATE'}) if defined($tu->{'DATE'});
	$seq->add_date($tu->{'GENE_INFO'}->{'DATE'});
	$seq->add_secondary_accession($tu->{'FEAT_NAME'}) unless $tu->{FEAT_NAME} eq $id;
	$seq->add_secondary_accession($tu->{'GENE_INFO'}->{'LOCUS'});
	$seq->add_secondary_accession($tu->{'GENE_INFO'}->{'ALT_LOCUS'})
		if defined($tu->{'GENE_INFO'}->{'ALT_LOCUS'});

	$seq->description(join(' | ',
		'genomic nucleotide',
		$tu->{'GENE_INFO'}->{'COM_NAME'}
	));

	my $strand = 1;
	$strand = -1 if $start > $end;
	
	my $source = new Bio::SeqFeature::Generic(
		-primary => 'source',
		-start   => 1,
		-end     => $seq->length(),
		-tag     => {
			date    =>  $tu->{DATE},
			end5    =>  $tu->{COORDSET}->{END5},
			end3    =>  $tu->{COORDSET}->{END3},
			strand  => ($tu->{COORDSET}->{END5} < $tu->{COORDSET}->{END3} ? "positive" : "negative"),
		}
	);

	$self->_add_source($source);
	$seq->add_SeqFeature($source);

	foreach my $model (@{$tu->{'MODEL'}}) {
		my @exons;
		my @cds;
		my $mstart = $strand*($model->{'COORDSET'}->{'END5'} - $start) + 1;
		my $mend   = $strand*($model->{'COORDSET'}->{'END3'} - $start) + 1;

		my $locus = $model->{PUB_LOCUS};
		if(!defined($locus) or $locus =~ /^\s*$/) {
			$locus = $model->{FEAT_NAME};
		}

		$seq->add_secondary_accession($locus);
		$seq->add_secondary_accession($model->{FEAT_NAME}) unless $locus eq $model->{FEAT_NAME};

		my $mfeat = new Bio::SeqFeature::Generic(
			-primary     => 'MODEL',
			-start       => $mstart,
			-end         => $mend,
			-tag => {
				locus => $locus,
				feat_name => $model->{'FEAT_NAME'}
			}
		);
		$seq->add_SeqFeature($mfeat);

		foreach my $exon (@{$model->{'EXON'}}) {
			my $estart = $strand * ($exon->{'COORDSET'}->{'END5'} - $start) + 1;
			my $eend   = $strand * ($exon->{'COORDSET'}->{'END3'} - $start) + 1;

			push @exons, {
				"start" => $estart,
				"end"   => $eend,
				"date"  => $exon->{'DATE'},
				"locus" => $exon->{'FEAT_NAME'},
			};


			if(defined($exon->{'CDS'})) {
				my $c = $exon->{'CDS'};
				my $cstart = $strand * ($c->{'COORDSET'}->{'END5'} - $start) + 1;
				my $cend   = $strand * ($c->{'COORDSET'}->{'END3'} - $start) + 1;

				push @cds, {
					"start" => $cstart,
					"end"   => $cend,
					"date"  => $c->{'DATE'},
					"locus" => $c->{'FEAT_NAME'},
				};
			}

			if(defined($exon->{'UTRS'})) {
				my $u = $exon->{'UTRS'};

				while( defined @{$u->{'LEFT_UTR'}} and 0 < @{$u->{'LEFT_UTR'}}) {
					my $lutr = shift @{$u->{'LEFT_UTR'}};
					my $us = $strand * ($lutr->{'END5'} - $start) + 1;
					my $ue = $strand * ($lutr->{'END3'} - $start) + 1;

					$seq->add_SeqFeature(new Bio::SeqFeature::Generic(
						-primary => 'LEFT_UTR',
						-start   => $us,
						-end     => $ue,
						-tag     => {
							locus => $model->{'PUB_LOCUS'} || $model->{'FEAT_NAME'}
						}
					));
				}

				while( defined @{$u->{'RIGHT_UTR'}} and 0 < @{$u->{'RIGHT_UTR'}}) {
					my $rutr = shift @{$u->{'RIGHT_UTR'}};
					my $us = $strand * ($rutr->{'END5'} - $start) + 1;
					my $ue = $strand * ($rutr->{'END3'} - $start) + 1;

					$seq->add_SeqFeature(new Bio::SeqFeature::Generic(
						-primary => 'RIGHT_UTR',
						-start   => $us,
						-end     => $ue,
						-tag => {
							locus => $model->{'PUB_LOCUS'} || $model->{'FEAT_NAME'}
						}
					));
				}

				while( defined @{$u->{'EXTENDED_UTR'}} and 0 < @{$u->{'EXTENDED_UTR'}}) {
					my $eutr = shift @{$u->{'EXTENDED_UTR'}};
					my $us = $strand * ($eutr->{'END5'} - $start) + 1;
					my $ue = $strand * ($eutr->{'END3'} - $start) + 1;

					$seq->add_SeqFeature(new Bio::SeqFeature::Generic(
						-primary => 'EXTENDED_UTR',
						-start   => $us,
						-end     => $ue,
						-tag => {
							locus => $model->{'PUB_LOCUS'} || $model->{'FEAT_NAME'}
						}
					));
				}
			}
		}

		my $loc = new Bio::Location::Split();
		foreach my $e (@exons) {
			$loc->add_sub_Location(new Bio::Location::Simple(
					-start => $e->{'start'},
					-end   => $e->{'end'}
			));
		}

		my $efeat = new Bio::SeqFeature::Generic(
			-primary      => 'EXON',
			-location     => $loc,
			-tag => {
				locus => $model->{'PUB_LOCUS'} || $model->{'FEAT_NAME'}
			}
		);
		$seq->add_SeqFeature($efeat);

		if(scalar(@cds) > 0) {
			$loc = new Bio::Location::Split();
			foreach my $c (@cds) {
				$loc->add_sub_Location(new Bio::Location::Simple(
						-start => $c->{'start'},
						-end   => $c->{'end'}
				));
			}

			my $cfeat = new Bio::SeqFeature::Generic(
				-primary     => 'CDS',
				-location    => $loc,
				-tag => {
					locus => $model->{'PUB_LOCUS'} || $model->{'FEAT_NAME'}
				}
			);

			my $trans;
			eval {
				$trans = new Bio::PrimarySeq(
					-id => $seq->primary_id(),
					-seq => $seq->subseq($loc)
				);
				$cfeat->add_tag_value('translation',
					$trans->translate(undef, undef, undef, undef, 1, 0)->seq()
				);
			};
			if($@) {
				$self->warn("Unable to translate protein. Probably a psuedo-protein.");
			} else {
				$seq->add_SeqFeature($cfeat);
			}
		}
	}
	
	my $kw = $self->{'ASSEMBLY'}->{'HEADER'}->{'KEYWORDS'};
	$seq->keywords(join(' ', @$kw)) if defined(@$kw);
	$self->_add_species($seq);
	
	while(my $comment = shift @{$tu->{'GENE_INFO'}->{'COMMENTS'}}) {
		my $com = new Bio::Annotation::Comment('-text' => $comment);
		$seq->annotation()->add_Annotation('comment', $com);
	}
	
	return $seq;
}


sub _add_species
{
	my($self, $seq) = @_;
	my $o = $self->{'ASSEMBLY'}->{'HEADER'}->{'ORGANISM'};
	my $lineage = $self->{'ASSEMBLY'}->{'HEADER'}->{'LINEAGE'};
	my $species = new Bio::Species(-classification => $lineage);
	my ($genus, $spec, @sub) = split(/ /, $o);

	$species->genus($genus);
	$species->species($spec);
	$species->sub_species(join(' ', @sub));

	$seq->species($species);
}

sub _add_source
{
	my($self, $source) = @_;
	$source->add_tag_value("chromosome", $self->{ASSEMBLY}->{CHROMOSOME});
	$source->add_tag_value("clone", $self->{ASSEMBLY}->{HEADER}->{CLONE_NAME});
}

sub _process
{
	my($self) = @_;
	my $line;
	my $tu = undef;

	return if $self->{'eof'} == 1;

	do {
		$line = $self->_readline();
		if($line =~ /<\?xml\s+version\s+=\s+"\d+\.\d+"\?>/o) {
			# do nothing
		} elsif ($line =~ /<!DOCTYPE (\w+) SYSTEM "[\w\.]+">/o) {
			$self->throw("DOCTYPE of $1, not TIGR!")
				if $1 ne "TIGR" ;
		} elsif ($line =~ /<TIGR>/o) {
			$self->_pushback($line);
			$self->_process_tigr();
		} elsif ($line =~ /<ASSEMBLY.*?>/o) {
			$self->_pushback($line);
			$self->_process_assembly();
		} elsif ($line =~ /<\/TIGR>/o) {
			$self->{'eof'}     = 1;
			return;
		} else {
			$self->throw("Unknown or Invalid process directive:",
				join('', ($line =~ /^\s*(<[^>]+>)/o)));
		}
		$tu = $self->{'ASSEMBLY'}->{'GENE_LIST'}->{'PROTEIN_CODING'};
	} until((defined($tu->{'TU'}) and @{$tu->{'TU'}} != 0)
			or !defined($line));
}

sub _process_tigr
{
	my($self) = @_;
	my $line;

	$line = $self->_readline();
	if($line !~ /<TIGR>/o) {
		$self->throw("Bio::SeqIO::tigr::_process_tigr called but no ",
		             "<TIGR> found in stream");
	}

	$line = $self->_readline();
	if($line =~ /<PSEUDOCHROMOSOME>/o) {
		$self->_pushback($line);
		$self->_process_pseudochromosome();
	} elsif ($line =~ /<ASSEMBLY.*?>/o) {
		$self->_pushback($line);
		$self->{'ASSEMBLY'} = $self->_process_assembly();
	}
}

sub _process_pseudochromosome
{
	my($self) = @_;
	my $line;

	$line = $self->_readline();
	return if $line !~ /<PSEUDOCHROMOSOME>/o;

	$line = $self->_readline();

	if($line =~ /<SCAFFOLD>/o) {
		$self->_pushback($line);
		$self->_process_scaffold();
		$line = $self->_readline();
	} else {
		$self->warn( "No Scaffold found in <PSUEDOCHROMOSOME> this " .
		             "is a violation of the TIGR dtd, but we ignore " .
			     "it so we are ignoring the error\n"
		);
	}

	if($line =~ /<ASSEMBLY.*>/o) {
		$self->_pushback($line);
		$self->{'ASSEMBLY'} = $self->_process_assembly();
		$line = $self->_readline();
	} else {
		$self->throw("Missing required ASSEMBLY in <PSEUDOCHROMOSOME>");
	}

	if($line =~ /<\/PSEUDOCHROMOSOME>/) {
		return;
	}

	$self->throw("Reached end of _process_psuedochromosome");
}

sub _process_assembly
{
	my($self) = @_;
	my $line;
	my $assembly;

	
	$line = $self->_readline();
	if($line !~ /<ASSEMBLY([^>]*)>/o) {
		$self->throw("Bio::SeqIO::tigr::_process_assembly called ",
		             "but no <ASSEMBLY> found in stream");
	}


	my %attribs = ($1 =~ /(\w+)\s*=\s+"(.*?)"/og);
	foreach my $key (keys(%attribs)) {
		$assembly->{$key} =  $attribs{$key};
	}

	$line = $self->_readline();
	my($attr, $val); 
	if(($attr, $val) = ($line =~ /<ASMBL_ID([^>]*)>([^<]*)<\/ASMBL_ID>/o)) {
		%attribs = ($attr =~ /(\w+)\s*=\s+"(.*?)"/og);
		foreach my $key (keys(%attribs)) {
			$assembly->{$key} =  $attribs{$key};
		}
		$assembly->{'ASMBL_ID'} = $val;
		$line = $self->_readtag();
	} else {
		$self->throw("Required <ASMBL_ID> missing");
	}

	if($line =~ /<COORDSET>/o) {
		$self->_pushback($line);
		$assembly->{'COORDSET'} = $self->_process_coordset();
		$line = $self->_readline();
	} else {
		$self->throw("Required <COORDSET> missing");
	}

	if($line =~ /<HEADER>/o) {
		$self->_pushback($line);
		$assembly->{'HEADER'} = $self->_process_header();
		$line = $self->_readline();
	} else {
		$self->throw("Required <HEADER> missing");
	}

	if($line =~ /<TILING_PATH>/o) {
		$self->_pushback($line);
		$assembly->{'TILING_PATH'} = $self->_process_tiling_path();
		$line = $self->_readline();
	}

	if($line =~ /<GENE_LIST>/o) {
		$self->_pushback($line);
		$assembly->{'GENE_LIST'} = $self->_process_gene_list();
		$line = $self->_readline();
	} else {
		$self->throw("Required <GENE_LIST> missing");
	}

	if($line =~ /<MISC_INFO>/o) {
		$self->_pushback($line);
		$assembly->{'MISC_INFO'} = $self->_process_misc_info();
		$line = $self->_readline();
	}

	if($line =~ /<REPEAT_LIST>/o) {
		$self->_pushback($line);
		$assembly->{'REPEAT_LIST'} = $self->_process_repeat_list();
		$line = $self->_readline();
	}

	if($line =~ /<ASSEMBLY_SEQUENCE>/o) {
		$self->_pushback($line);
		$assembly->{'ASSEMBLY_SEQUENCE'}=$self->_process_assembly_seq();
		$line = $self->_readline();
	} else {
		$self->throw("Required <ASSEMBLY_SEQUENCE> missing");
	}

	if($line =~ /<\/ASSEMBLY>/o) {
		return $assembly;
	}

	$self->throw("Reached the end of <ASSEMBLY>");
}

sub _process_assembly_seq()
{
	my ($self) = @_;
	my $line;
	
	$line = $self->_readline();
	if($line !~ /<ASSEMBLY_SEQUENCE>/o) {
		$self->throw("Bio::SeqIO::tigr::_process_assembly_seq called ".
			     "with no <ASSEMBLY_SEQUENCE> in the stream");
	}

	$line = $self->_readline();
	if($line =~ /^(.+)<\/ASSEMBLY_SEQUENCE>/o) {
		return $1;
	}

	$self->throw("Reached end of _proces_assembly");
}

sub _process_coordset($)
{
	my ($self) = @_;
	my $line;
	my $h;

	$line = $self->_readline();
	if($line =~ /<COORDSET>/o) {
		$self->_pushback($line);
		$line = $self->_readtag();
		($h->{'END5'}, $h->{'END3'}) = ($line =~ /<COORDSET>\s*<END5>\s*(\d+)\s*<\/END5>\s*<END3>\s*(\d+)\s*<\/END3>/os);
		if(!defined($h->{'END5'}) or !defined($h->{'END3'})) {
			$self->throw("Invalid <COORDSEt>");
		}
		return $h;
	} else {
		$self->throw("Bio::SeqIO::tigr::_process_coordset() called ",
		             "but no <COORDSET> found in stream");
	}
}

sub _process_header
{
	my ($self) = @_;
	my $header;
	my $line = $self->_readline();

	if($line !~ /<HEADER>/o) {
		$self->throw("Bio::SeqIO::tigr::_process_header called ",
		             "but no <HEADER> found in stream");
	}

	$line = $self->_readtag();
	if($line =~ /<CLONE_NAME>([^>]+)<\/CLONE_NAME>/o) {
		$header->{'CLONE_NAME'} = $1;
		$line = $self->_readtag();
	} else {
		$self->throw("Required <CLONE_NAME> missing");
	}

	if($line =~ /<SEQ_LAST_TOUCHED>/o) {
		($header->{'SEQ_LAST_TOUCHED'}) =
			($line =~ /<DATE>([^<]+)<\/DATE>/o);
		$line = $self->_readtag();
	} else {
		$self->throw("Reqired <SEQ_LAST_TOUCHED> missing");
	}

	if($line =~ /<GB_ACCESSION>([^<]*)<\/GB_ACCESSION>/o) {
		$header->{'GB_ACCESSION'} = $1;
		$line = $self->_readtag();
	} else {
		$self->throw("Required <GB_ACCESSION> missing");
	}

	if($line =~ /<ORGANISM>([^<]*)<\/ORGANISM>/o) {
		$header->{'ORGANISM'} = $1;
		$line = $self->_readtag();
	} else {
		$self->throw("Required <ORGANISM> missing");
	}

	if($line =~ /<LINEAGE>([^<]*)<\/LINEAGE>/o) {
		@{$header->{'LINEAGE'}}
			= reverse(split(/\s*;\s*/o, $1));
		$line = $self->_readtag();
	} else {
		$self->throw("Required <LINEAGE> missing");
	}

	if($line =~ /<SEQ_GROUP>([^<]*)<\/SEQ_GROUP>/o) {
		$header->{'SEQ_GROUP'} = $1;
		$line = $self->_readtag();
	} else {
		$self->throw("Required <SEQ_GROUP> missing");
	}

	while($line =~ /<KEYWORDS>[^<]*<\/KEYWORDS>/o) {
		push(@{$header->{'KEYWORDS'}}, $1);
		$line = $self->_readtag();
	}

	while($line =~ /<GB_DESCRIPTION>([^<]+)<\/GB_DESCRIPTION>/o) {
		push(@{$header->{'GB_DESCRIPTION'}},$1);
		$line = $self->_readtag();
	}

	while($line =~ /<GB_COMMENT>([^<]+)<\/GB_COMMENT>/o) {
		push(@{$header->{'GB_COMMENT'}}, $1);
		$line = $self->_readtag();
	}

	if(my %h = ($line =~ /<AUTHOR_LIST(?:\s*(\w+)\s*=\s*"([^"]+)"\s*)*>/o)) {
		$header->{'AUTHOR_LIST'}=$h{'CONTACT'};
		while($line !~ /<\/AUTHOR_LIST>/o) {
			$self->_readtag();
		}
		$line = $self->_readline();
	} else {
		$self->throw("Required <AUTHOR_LIST> missing");
	}

	if($line =~ /<\/HEADER>/o) {
		return $header;
	}

	$self->throw("Reached end of header\n");
}

sub _process_gene_list
{
	my($self) = @_;
	my $gene;
	my $line;

	$line = $self->_readline();
	if($line !~ /<GENE_LIST>/o) {
		$self->throw("Bio::SeqIO::tigr::_process_gene_list called ",
		             "but no <GENE_LIST> in the stream");
	}

	$line = $self->_readline();
	if($line =~ /<PROTEIN_CODING>/o) {
		$self->_pushback($line);
		$gene->{'PROTEIN_CODING'} = $self->_process_protein_coding();
		$line = $self->_readline();
	} else {
		$self->throw("Required <PROTEIN_CODING> missing");
	}

	if($line =~ /<RNA_GENES>/o) {
		$self->_pushback($line);
		$gene->{'RNA_GENES'} = $self->_process_rna_genes();
		$line = $self->_readline();
	} else {
		$self->throw("Required <RNA_GENES> missing");
	}

	if($line =~ /<\/GENE_LIST>/o) {
		return $gene;
	}

	$self->throw("Reached end of _process_gene_list");
}

sub _process_protein_coding
{
	my ($self) = @_;
	my $prot;
	my $line = $self->_readline();

	if($line !~ /<PROTEIN_CODING>/o) {
		$self->throw("Bio::SeqIO::tigr::_process_protein_coding called"
		             . "but no <GENE_LIST> in the stream");
	}

	$line = $self->_readline();
	while($line =~ /<TU>/o) {
		$self->_pushback($line);
		push(@{$prot->{'TU'}}, $self->_process_tu());
		$line = $self->_readline();
	}

	if($line =~ /<\/PROTEIN_CODING>/o) {
		return $prot;
	}

	$self->throw("Reached end of _process_protein_coding");
}


sub _process_rna_genes
{
	my ($self) = @_;
	my $line = $self->_readline();

	if($line =~ /<RNA_GENES>/o) {
		while($line !~ /<\/RNA_GENES>/o) {
			$line = $self->_readline();
		}
	} else {
		$self->throw("Bio::SeqIO::tigr::_process_rna_genes called ",
		             "but no <RNA_GENES> in the stream");
	}
}

sub _process_misc_info
{
	my ($self) = @_;
	my $line = $self->_readline();

	if($line =~ /<MISC_INFO>/o) {
		while($line !~ /<\/MISC_INFO>/o) {
			$line = $self->_readline();
		}
	} else {
		$self->throw("Bio::SeqIO::tigr::_process_misc_info called ",
		             "but no <MISC_INFO> in the stream");
	}
}

sub _process_repeat_list
{
	my ($self) = @_;
	my $line = $self->_readline();

	if($line =~ /<REPEAT_LIST>/o) {
		while($line !~ /<\/REPEAT_LIST>/o) {
			$line = $self->_readline();
		}
	} else {
		$self->throw("Bio::SeqIO::tigr::_process_repeat_list called ",
		             "but no <MISC_INFO> in the stream");
	}
}

sub _process_tiling_path
{
	my($self) = @_;
	my $line = 4self->_readine();


	if($line =~ /<TILING_PATH>/o) {
		while($line !~ /<\/TILING_PATH>/o) {
			$line = $self->_readline();
		}
	} else {
		$self->throw("Bio::SeqIO::tigr::_process_repeat_list called ",
		             "but no <MISC_INFO> in the stream");
	}
}

sub _process_scaffold
{
	my ($self) = @_;
	my $line;

	# for now we just skip them
	$line = $self->_readline();
	return if $line !~ /<SCAFFOLD>/o;
	do {
		$line = $self->_readline();
	} while(defined($line) && $line !~ /<\/SCAFFOLD>/o);
}

sub _process_tu
{
	my($self) = @_;
	my $line = $self->_readline();
	my $tu;

	if($line !~ /<TU>/o) {
		$self->throw("Process_tu called when no <TU> tag");
	}

	$line = $self->_readtag();
	if ($line =~ /<FEAT_NAME>([\w\.]+)<\/FEAT_NAME>/o) {
		$tu->{'FEAT_NAME'} = $1;
		$line = $self->_readtag();
	} else {
		$self->throw("Invalid Feat_Name");
	}

	while($line =~ /<GENE_SYNONYM>/o) {
		$line = $self->_readtag();
	}
	
	while($line =~ /<CHROMO_LINK>([\w\.]+)<\/CHROMO_LINK>/o) {
		 push @{$tu->{'CHROMO_LINK'}}, $1;
		 $line = $self->_readtag();
	}

	if ($line =~ /<DATE>([^>]+)<\/DATE>/o) {
		$tu->{'DATE'} = $1;
		$line = $self->_readline();
	} else {
		$self->throw("Invalid Date");
	}

	if ($line =~ /<GENE_INFO>/o) {
		$self->_pushback($line);
		$tu->{'GENE_INFO'} = $self->_process_gene_info();
		$line = $self->_readline();
	} else {
		$self->throw("Invalid Gene_Info");
	}

	if($line =~ /<COORDSET>/o) {
		$self->_pushback($line);
		$tu->{'COORDSET'} = $self->_process_coordset();
		$line = $self->_readline();
	} else {
		$self->throw("Invalid Coordset");
	}

	if($line =~ /<MODEL[^>]*>/o) {
		do {
			$self->_pushback($line);
			push(@{$tu->{'MODEL'}},
				$self->_process_model());
			$line = $self->_readline();
		} while($line =~ /<MODEL[^>]*>/o);
		$self->_pushback($line);
		$line = $self->_readtag();
	} else {
		$self->throw("Expected <MODEL> not found");
	}

	if($line =~ /<TRANSCRIPT_SEQUENCE>/o) {
		$line = $self->_readtag();
	}

	if($line =~ /<GENE_EVIDENCE>/o) {
		$line = $self->_readtag();
	}

	while($line =~ /<URL[^>]*>[^<]*<\/URL>/o) {
		$line = $self->_readtag();
	}

	if($line =~ /<\/TU>/o) {
		return $tu;
	} else {
		$self->throw("Expected </TU> not found: $line");
	}
}

sub _process_gene_info
{
	my($self) = @_;
	my $line = $self->_readline();
	my $geneinfo;

	$self->throw("Invalid Gene Info: $line") if $line !~ /<GENE_INFO>/o;
	$line = $self->_readline();

	if($line =~ /<LOCUS>([^>]*)<\/LOCUS>/o) {
		$geneinfo->{'LOCUS'} = $1;
		$line = $self->_readline();
	} else {
		$self->throw("Invalid Locus");
	}

	if($line =~ /<ALT_LOCUS>([^>]+)<\/ALT_LOCUS>/o) {
		$geneinfo->{'ALT_LOCUS'} = $1;
		$line = $self->_readline();
	}

	if($line =~ /<PUB_LOCUS>([^>]+)<\/PUB_LOCUS>/o) {
		$geneinfo->{'PUB_LOCUS'} = $1;
		$line = $self->_readtag();
	} else {
#		$self->throw("Invalid Pub_Locus");
	}

	if($line =~ /<COM_NAME\s+CURATED="(\d+)">([^>]+)<\/COM_NAME>/o) {
		$geneinfo->{'CURATED'} = $1;
		$geneinfo->{'COM_NAME'} = $2;
		$line = $self->_readtag();
	} else {
		$self->throw("invalid com_name");
	}

	while($line =~ /<COMMENT>([^<]+)<\/COMMENT>/o) {
		push(@{$geneinfo->{'COMMENTS'}}, $1);
		$line = $self->_readtag();
	}

	while($line =~ /<PUB_COMMENT>([^<]+)<\/PUB_COMMENT>/o) {
		push(@{$geneinfo->{'COMMENTS'}}, $1);
		$line = $self->_readtag();
	}

	if($line =~ /<EC_NUM>([\w\-\\\.]+)<\/EC_NUM>/o) {
		$geneinfo->{'EC_NUM'} = $1;
		$line = $self->_readtag();
	}

	if($line =~ /<GENE_SYM>([^<]+)<\/GENE_SYM>/o) {
		$geneinfo->{'GENE_SYM'} = $1;
		$line = $self->_readtag();
	}

	if($line =~ /<IS_PSEUDOGENE>([^>]+)<\/IS_PSEUDOGENE>/o) {
		$geneinfo->{'IS_PSEUDOGENE'} = $1;
		$line = $self->_readtag();
	} else {
		$self->throw("invalid is_pseudogene: $line");
	}

	if($line =~ /<FUNCT_ANNOT_EVIDENCE/o) {
		$line = $self->_readtag();
	}

	if($line =~ /<DATE>([^>]+)<\/DATE>/o) {
		$geneinfo->{'DATE'} = $1;
		$line = $self->_readtag();
	}

	if($line =~ /<GENE_ONTOLOGY>/o) {
		until($line =~ /<\/GENE_ONTOLOGY>/o) {
			$line = $self->_readline();
		}
		$line = $self->_readline();
	}
	
	if($line =~ /<\/GENE_INFO/o) {
		return $geneinfo;
	}

	$self->throw("unexpected end of gene_info");
}

sub _process_model
{
	my($self) = @_;
	my $line;
	my $model;

	$line = $self->_readline();
	if($line !~ /<MODEL ([^>]+)>/o) {
		$self->throw("Invalid Model: $line")
	}
	my %attribs = ($1 =~ /(\w+)\s*=\s*"([^"]*)"/og);
	$model->{'CURATED'} = $attribs{'CURATED'};
	$line = $self->_readline();

	if($line =~ /<FEAT_NAME>([^>]+)<\/FEAT_NAME>/o) {
		$model->{'FEAT_NAME'} = $1;
		$line = $self->_readline();
	} else {
		$self->throw("Invalid Feature Name: $line");
	}

	if($line =~ /<PUB_LOCUS>([^>]+)<\/PUB_LOCUS>/o) {
		$model->{'PUB_LOCUS'} = $1;
		$line = $self->_readline();
	} else {
#		$self->throw("Invalid Pub_Locus: $line");
	}

	if($line =~ /<CDNA_SUPPORT>/o) {
		$self->_pushback($line);
		$self->_process_cdna_support();
		$line = $self->_readline();
	}

	while($line =~ /<CHROMO_LINK>([^>]+)<\/CHROMO_LINK>/o) {
		push @{$model->{'CHROMO_LINK'}}, $1;
		$line = $self->_readline();
	} 

	if($line =~ /<DATE>([^>]+)<\/DATE>/o) {
		$model->{'DATE'} = $1;
		$line = $self->_readline();
	} else {
		$self->throw("Invalid Date: $line");
	}

	if($line =~ /<COORDSET>/o) {
		$self->_pushback($line);
		$model->{'COORDSET'} = $self->_process_coordset();
		$line = $self->_readline();
	} else {
		$self->throw("Invalid Coordset: $line");
	}

	if($line =~ /<EXON>/o) {
		do {
			$self->_pushback($line);
			push(@{$model->{'EXON'}}, $self->_process_exon());
			$line = $self->_readline();
		} while($line =~ /<EXON>/o);
	} else {
		$self->throw("Required <EXON> missing");
	}
	
	until($line =~ /<\/MODEL>/o) {
		$line = $self->_readline();
	}

	return $model;
}


sub _process_cdna_support
{
	my($self) = @_;
	my $line = $self->_readline();

	if($line !~ /<CDNA_SUPPORT>/o) {
		$self->throw("Bio::SeqIO::tigr::_process_cdna_support called ",
		             "but no <CDNA_SUPPORT> in the stream");
	}

	# TODO Add CDNA Support
	do {
		$line = $self->_readline();
	} while($line !~ /<\/CDNA_SUPPORT>/o);
}


sub _process_exon
{
	my($self) = @_;
	my $line = $self->_readline();
	my $exon;

	if($line !~ /<EXON>/o) {
		$self->throw("Bio::SeqIO::tigr::_process_exon called ",
		             "but no <EXON> in the stream");
	}

	$line = $self->_readtag();
	if($line =~ /<FEAT_NAME>([^<]+)<\/FEAT_NAME>/o) {
		$exon->{'FEAT_NAME'} = $1;
		$line = $self->_readtag();
	} else {
		$self->throw("Required <FEAT_NAME> missing");
	}

	if($line =~ /<DATE>([^<]+)<\/DATE>/o) {
		$exon->{'DATE'} = $1;
		$line = $self->_readtag();
	} else {
		$self->throw("Required <DATE> missing");
	}

	if($line =~ /<COORDSET>/o) {
		$self->_pushback($line);
		$exon->{'COORDSET'} = $self->_process_coordset();
		$line = $self->_readline();
	} else {
		$self->throw("Required <COORDSET> missing");
	}

	if($line =~ /<CDS>/o) {
		$self->_pushback($line);
		$exon->{'CDS'} = $self->_process_cds();
		$line = $self->_readline();
	}

	if($line =~ /<UTRS>/o) {
		$self->_pushback($line);
		$exon->{'UTRS'} = $self->_process_utrs();
		$line = $self->_readline();
	}

	if($line =~ /<\/EXON>/o) {
		return $exon;
	}

	$self->throw("Reached End of Bio::SeqIO::tigr::_process_exon");
}

sub _process_cds
{
	my($self) = @_;
	my $line = $self->_readline();
	my $cds;

	if($line !~ /<CDS>/o) {
		$self->throw("Bio::SeqIO::tigr::_process_cda_support called ",
		             "but no <CDS> in the stream");
	}
	
	$line = $self->_readtag();
	if($line =~ /<FEAT_NAME>([^<]+)<\/FEAT_NAME>/o) {
		$cds->{'FEAT_NAME'} = $1;
		$line = $self->_readtag();
	} else {
		$self->throw("Required <FEAT_NAME> missing");
	}

	if($line =~ /<DATE>([^<]+)<\/DATE>/o) {
		$cds->{'DATE'} = $1;
		$line = $self->_readtag();
	} else {
		$self->throw("Required <DATE> missing");
	}

	if($line =~ /<COORDSET>/o) {
		$self->_pushback($line);
		$cds->{'COORDSET'} = $self->_process_coordset();
		$line = $self->_readtag();
	} else {
		$self->throw("Required <COORDSET> missing");
	}

	if($line =~ /<\/CDS>/o) {
		return $cds;
	}

	$self->throw("Reached onf of Bio::SeqIO::tigr::_process_cds");
}

sub _process_utrs
{
	my($self) = @_;
	my $line = $self->_readline();
	my $utrs;

	if($line !~ /<UTRS/o) {
		$self->throw("Bio::SeqIO::tigr::_process_utrs called but no ",
		             "<UTRS> found in stream");
	}

	$line = $self->_readline();
	while($line !~ /<\/UTRS>/o) {
		$self->_pushback($line);
		if($line =~ /<LEFT_UTR>/o) {
			push(@{$utrs->{'LEFT_UTR'}}, $self->_process_left_utr());
		} elsif ($line =~ /<RIGHT_UTR>/o) {
			push(@{$utrs->{'RIGHT_UTR'}}, $self->_process_right_utr());
		} elsif ($line =~ /<EXTENDED_UTR>/o) {
			push(@{$utrs->{'EXTENDED_UTR'}}, $self->_process_ext_utr());
		} else {
			$self->throw("Unexpected tag");
		}
	
		$line = $self->_readline();
	}

	if($line =~ /<\/UTRS>/o) {
		return $utrs;
	}
	$self->throw("Reached end of Bio::SeqIO::tigr::_process_utrs");
}

sub _process_left_utr
{
	my($self) = @_;
	my $line = $self->_readline();
	my $coordset;

	if($line !~ /<LEFT_UTR>/o) {
		$self->throw("Bio::SeqIO::tigr::_process_left_utr called but ",
		             "no <LEFT_UTR> found in stream");
	}

	$line = $self->_readtag();
	if($line =~ /<COORDSET>/o) {
		$self->_pushback($line);
		$coordset = $self->_process_coordset();
	} else {
		$self->throw("Required <COORDSET> missing");
	}

	$line = $self->_readline();
	if($line =~ /<\/LEFT_UTR>/o) {
		return $coordset;
	}
	$self->throw("Reached end of Bio::SeqIO::tigr::_process_left_utr");
}

sub _process_right_utr
{
	my($self) = @_;
	my $line = $self->_readline();
	my $coordset;

	if($line !~ /<RIGHT_UTR>/o) {
		$self->throw("Bio::SeqIO::tigr::_process_right_utr called but ",
		             "no <RIGHT_UTR> found in stream");
	}

	$line = $self->_readtag();
	if($line =~ /<COORDSET>/o) {
		$self->_pushback($line);
		$coordset = $self->_process_coordset();
	} else {
		$self->throw("Required <COORDSET> missing");
	}

	$line = $self->_readline();
	if($line =~ /<\/RIGHT_UTR>/o) {
		return $coordset;
	}
	$self->throw("Reached end of Bio::SeqIO::tigr::_process_right_utr");
}

sub _process_ext_utr
{
	my($self) = @_;
	my $line = $self->_readline();
	my $coordset;

	if($line !~ /<EXTENDED_UTR>/o) {
		$self->throw("Bio::SeqIO::tigr::_process_ext_utr called but ",
		             "no <EXTENDED_UTR> found in stream");
	}

	$line = $self->_readtag();
	if($line =~ /<COORDSET>/o) {
		$self->_pushback($line);
		$coordset = $self->_process_coordset();
	} else {
		$self->throw("Required <COORDSET> missing");
	}

	$line = $self->_readline();
	if($line =~ /<\/EXTENDED_UTR>/o) {
		return $coordset;
	}
	$self->throw("Reached end of Bio::SeqIO::tigr::_process_ext_utr");
}

sub _readtag
{
	my($self) = @_;
	my $line = $self->_readline();
	chomp($line);

	my $tag;
	if(($tag) = ($line =~ /^[^<]*<\/(\w+)/o)) {
		$self->_pushback($1) if $line =~ /<\/$tag>(.+)$/;
		return "</$tag>";
	}
 
	until(($tag) = ($line =~ /<(\w+)[^>]*>/o)) {
		$line = $self->_readline();
		chomp $line;
	}

	until($line =~ /<\/$tag>/) {
		$line .= $self->_readline();
	}

	if(my ($val) = ($line =~ /(<$tag.*>.*?<\/$tag>)/s)) {
		if($line =~ /<\/$tag>\s*(\w+[\s\w]*?)\s*$/s) {
			$self->_pushback($1)
		}
		return $val;
	}
	$self->throw("summerror");
}

sub _readline
{
	my($self) = @_;
	my $line;
	do {
		$line = $self->SUPER::_readline();
	} while(defined($line) and $line =~ /^\s*$/o);

	return $line;
}

sub throw
{
	my($self, @s) = @_;
	my $string = "[$.]" . join('', @s);
	$self->SUPER::throw($string);
}
