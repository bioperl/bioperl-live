# $Id$
#
# BioPerl module for Bio::Structure::IO::pdb
#
# Cared for by Kris Boulez <kris.boulez@algonomics.com>
#
# Copyright 2001 Kris Boulez
#
# Framework is a copy of Bio::SeqIO::embl.pm
# 
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Structure::IO::embl - PDB input/output stream

=head1 SYNOPSIS

It is probably best not to use this object directly, but
rather go through the Bio::Structure::IO handler system. Go:

    $stream = Bio::Structure::IO->new(-file => $filename, -format => 'PDB');

    while ( (my $structure = $stream->next_structure()) ) {
	# do something with $structure
    }

=head1 DESCRIPTION

This object can transform Bio::Structure objects to and from PDB flat
file databases (at this time only reading is supported, most of the
framework for writing is in place).

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this
and other Bioperl modules. Send your comments and suggestions preferably
 to one of the Bioperl mailing lists.
Your participation is much appreciated.

  bioperl-l@bioperl.org                 - General discussion
  http://www.bioperl.org/MailList.shtml - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
 the bugs and their resolution.
 Bug reports can be submitted via email or the web:

  bioperl-bugs@bio.perl.org
  http://bio.perl.org/bioperl-bugs/

=head1 AUTHOR - Kris Boulez

Email kris.boulez@algonomics.com

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Structure::IO::pdb;
use vars qw(@ISA);
use strict;
use Bio::Structure::IO;
use Bio::Structure::Entry;
#use Bio::Structure::Model;
#use Bio::Structure::Chain;
#use Bio::Structure::Residue;
use Bio::Structure::Atom;
use Bio::SeqFeature::Generic;

@ISA = qw(Bio::Structure::IO);

sub _initialize {
  my($self,@args) = @_;

  $self->SUPER::_initialize(@args);  

  my ($noheader, $noatom) =
  	$self->_rearrange([qw(
			NOHEADER
			NOATOM
		)],
		@args);
  $noheader && $self->_noheader($noheader);
  $noatom   && $self->_noatom($noatom);
}


=head2 next_structure;

 Title   : next_structure
 Usage   : $seq = $stream->next_seq()
 Function: returns the next structure in the stream
 Returns : Bio::Seq object
 Args    :


=cut

sub next_structure {
   my ($self,@args) = @_;
   my ($line);
   my ($obslte, $title, $caveat, $compnd, $source, $keywds,
	$expdta, $author, %revdat, $revdat, $sprsde, $jrnl, %remark, $dbref,
	$seqadv, $seqres, $modres, $het, $hetnam, $hetsyn, $formul, $helix, 
	$sheet, $turn, $ssbond, $link, $hydbnd, $sltbrg, $cispep,
	$site, $cryst1, $tvect,);
   my $struc = Bio::Structure::Entry->new(-id => 'created from pdb.pm');
   my $all_headers = ( !$self->_noheader );  # we'll parse all headers and store as annotation
   my %header;  # stores all header RECORDs an is stored as annotations when ATOM is reached


   $line = $self->_readline;   # This needs to be before the first eof() test

   if( !defined $line ) {
       return undef; # no throws - end of file
   }

   if( $line =~ /^\s+$/ ) {
       while( defined ($line = $self->_readline) ) {
	   $line =~/\S/ && last;
       }
   }   
   if( !defined $line ) {
       return undef; # end of file
   }
   $line =~ /^HEADER\s+\S+/ || $self->throw("PDB stream with no HEADER. Not pdb in my book");
   my($class, $depdate, $idcode) = unpack "x10 A40 A9 x3 A4", $line;
   $idcode =~ s/^\s*(\S+)\s*$/$1/;
   $struc->id($idcode);
$self->debug("PBD c $class d $depdate id $idcode\n"); # XXX KB

   my $buffer = $line;
 
   
   BEFORE_COORDINATES :
   until( !defined $buffer ) {
       $_ = $buffer;

       # Exit at start of coordinate section
       last if /^(MODEL|ATOM|HETATM)/;

       # OBSLTE line(s)
       if (/^OBSLTE / && $all_headers) {
		$obslte = $self->_read_PDB_singlecontline("OBSLTE","12-70",\$buffer);
		$header{'obslte'} = $obslte;
       }

       # TITLE line(s)
       if (/^TITLE / && $all_headers) {
		$title = $self->_read_PDB_singlecontline("TITLE","11-70",\$buffer);
		$header{'title'} = $title;
       }

       # CAVEAT line(s)
       if (/^CAVEAT / && $all_headers) {
		$caveat = $self->_read_PDB_singlecontline("CAVEAT","12-70",\$buffer);
		$header{'caveat'} = $caveat;
       }

       # COMPND line(s)
       if (/^COMPND / && $all_headers) {
		$compnd = $self->_read_PDB_singlecontline("COMPND","11-70",\$buffer);
		$header{'compnd'} = $compnd;
$self->debug("get COMPND $compnd\n");
       }

	# SOURCE line(s)
	if (/^SOURCE / && $all_headers) {
		$source = $self->_read_PDB_singlecontline("SOURCE","11-70",\$buffer);
		$header{'source'} = $source;
	}

	# KEYWDS line(s)
	if (/^KEYWDS / && $all_headers) {
		$keywds = $self->_read_PDB_singlecontline("KEYWDS","11-70",\$buffer);
		$header{'keywds'} = $keywds;
	}

	# EXPDTA line(s)
	if (/^EXPDTA / && $all_headers) {
		$expdta = $self->_read_PDB_singlecontline("EXPDTA","11-70",\$buffer);
		$header{'expdta'} = $expdta;
	}
	
	# AUTHOR line(s)
	if (/^AUTHOR / && $all_headers) {
		$author = $self->_read_PDB_singlecontline("AUTHOR","11-70",\$buffer);
		$header{'author'} = $author;
	}

	# REVDAT line(s)
	#  a bit more elaborate as we also store the modification number
	if (/^REVDAT / && $all_headers) {
		##my($modnum,$rol) = unpack "x7 A3 x3 A53", $_;
		##$modnum =~ s/\s+//; # remove  spaces
		##$revdat{$modnum} .= $rol;
		my ($rol) = unpack "x7 a59", $_;
		$revdat .= $rol;
		$header{'revdat'} = $revdat;
	}

	# SPRSDE line(s)
	if (/^SPRSDE / && $all_headers) {
		$sprsde = $self->_read_PDB_singlecontline("SPRSDE","12-70",\$buffer);
		$header{'sprsde'} = $sprsde;
	}

	# jRNL line(s)
	if (/^JRNL / && $all_headers) {
		$jrnl = $self->_read_PDB_jrnl(\$buffer);
		$struc->annotation->add_Annotation('reference',$jrnl);
	}

	# REMARK line(s)
	#  we only parse the "REMARK   1" lines (additional references)
	#  thre rest is stored in %remark (indexed on remarkNum) (pack does space-padding)
	if (/^REMARK\s+(\d+)\s*/ && $all_headers) {
		my $remark_num = $1;
		if ($remark_num == 1) {
			my @refs = $self->_read_PDB_remark_1(\$buffer);
			# How can we find the primary reference when writing (JRNL record) XXX KB
			foreach my $ref (@refs) {
				$struc->annotation->add_Annotation('reference', $ref);
			}
		} else { # other remarks, we store literlly at the moment
			my ($rol) = unpack "x11 a59", $_;
			$remark{$remark_num} .= $rol;
		}
	} # REMARK

	# DBREF line(s)
	#  references to sequences in other databases
	#  we store as 'dblink' annotations and whole line as simple annotation (round-trip)
	if (/^DBREF / && $all_headers) {
		my ($rol) = unpack "x7 a61", $_;
		$dbref .= $rol;
		my ($db, $acc) = unpack "x26 a6 x1 a8", $_;
		$db =~ s/\s*$//;
		$acc =~ s/\s*$//;
		my $link = Bio::Annotation::DBLink->new;
		$link->database($db);
		$link->primary_id($acc);
		$struc->annotation->add_Annotation('dblink', $link);
	} # DBREF

	# SEQADV line(s)
	if (/^SEQADV / && $all_headers) {
		my ($rol) = unpack "x7 a63", $_;
		$self->_concatenate_lines($seqadv, $rol);
		$header{'seqadv'} = $seqadv;
	} # SEQADV

	# SEQRES line(s)
	#  this is (I think) the sequence of macromolecule that was analysed
	#  this will be returned when doing $struc->seq
	if (/^SEQRES / && $all_headers) {
		my ($rol) = unpack "x7 a63", $_;
		$seqres .= $rol;
		$header{'seqres'} = $seqres;
	} # SEQRES
	
	# MODRES line(s)
	if (/^MODRES / && $all_headers) {
		my ($rol) = unpack "x7 a63", $_;
		$modres .= $rol;
		$header{'modres'} = $modres;
	} # MODRES

	# HET line(s)
	if (/^HET / && $all_headers) {
		my ($rol) = unpack "x7 a63", $_;
		$het .= $rol;
		$header{'het'} = $het;
	} # HET

	# HETNAM line(s)
	if (/^HETNAM / && $all_headers) {
		my ($rol) = unpack "x8 a62", $_;
		$hetnam .= $rol;
		$header{'hetnam'} = $hetnam;
	} # HETNAM

	# HETSYN line(s)
	if (/^HETSYN / && $all_headers) {
		my ($rol) = unpack "x8 a62", $_;
		$hetsyn .= $rol;
		$header{'hetsyn'} = $hetsyn;
	} # HETSYN

	# FORMUL line(s)
	if (/^FORMUL / && $all_headers) {
		my ($rol) = unpack "x8 a62", $_;
		$formul .= $rol;
		$header{'formul'} = $formul;
	} # FORMUL
	
	# HELIX line(s)
	#  store as specific object ??
	if (/^HELIX / && $all_headers) {
		my ($rol) = unpack "x7 a69", $_;
		$helix .= $rol;
		$header{'helix'} = $helix;
	} # HELIX
	
	# SHEET line(s)
	#  store as specific object ??
	if (/^SHEET / && $all_headers) {
		my ($rol) = unpack "x7 a63", $_;
		$sheet .= $rol;
		$header{'sheet'} = $sheet;
	} # SHEET

	# TURN line(s)
	#  store as specific object ??
	if (/^TURN / && $all_headers) {
		my ($rol) = unpack "x7 a63", $_;
		$turn .= $rol;
		$header{'turn'} = $turn;
	} # TURN

	# SSBOND line(s)
	#  store in connection-like object (see parsing of CONECT record)
	if (/^SSBOND / && $all_headers) {
		my ($rol) = unpack "x7 a65", $_;
		$ssbond .= $rol;
		$header{'ssbond'} = $ssbond;
	} # SSBOND

	# LINK
	#  store like SSBOND ?
	if (/^LINK / && $all_headers) {
		my ($rol) = unpack "x12 a60", $_;
		$link .= $rol;
		$header{'link'} = $link;
	} # LINK

	# HYDBND
	#  store like SSBOND
	if (/^HYDBND / && $all_headers) {
		my ($rol) = unpack "x12 a60", $_;
		$hydbnd .= $rol;
		$header{'hydbnd'} = $hydbnd;
	} # HYDBND

	# SLTBRG
	#  store like SSBOND ?
	if (/^SLTBRG / && $all_headers) {
		my ($rol) = unpack "x12 a60",$_;
		$sltbrg .= $rol;
		$header{'sltbrg'} = $sltbrg;
	} # SLTBRG

	# CISPEP
	#   store like SSBOND ?
	if (/^CISPEP / && $all_headers) {
		my ($rol) = unpack "x7 a52", $_;
		$cispep .= $rol;
		$header{'cispep'} = $cispep;
	}

	# SITE line(s)
	if (/^SITE / && $all_headers) {
		my ($rol) = unpack "x7 a54", $_;
		$site .= $rol;
		$header{'site'} = $site;
	} # SITE

	# CRYST1 line
	#  store in some crystallographic subobject ?
	if (/^CRYST1/ && $all_headers) {
		my ($rol) = unpack "x6 a64", $_;
		$cryst1 .= $rol;
		$header{'cryst1'} = $cryst1;
	} # CRYST1

	# ORIGXn line(s) (n=1,2,3)
	if (/^(ORIGX\d) / && $all_headers) {
		my $origxn = lc($1);
		my ($rol) = unpack "x10 a45", $_;
		$header{$origxn} .= $rol;
	} # ORIGXn

	# SCALEn line(s) (n=1,2,3)
	if (/^(SCALE\d) / && $all_headers) {
		my $scalen = lc($1);
		my ($rol) = unpack "x10 a45", $_;
		$header{$scalen} .= $rol;
	} # SCALEn
	
	# MTRIXn line(s) (n=1,2,3)
	if (/^(MTRIX\d) / && $all_headers) {
		my $mtrixn = lc($1);
		my ($rol) = unpack "x7 a53", $_;
		$header{$mtrixn} .= $rol;
	} # MTRIXn

	# TVECT line(s)
	if (/^TVECT / && $all_headers) {
		my ($rol) = unpack "x7 a63", $_;
		$tvect .= $rol;
		$header{'tvect'} = $tvect;
	}

	# Get next line.
	$buffer = $self->_readline;
   }
   
   # store %header entries a annotations
   if (%header) {
	for my $record (keys %header) {
		my $sim = Bio::Annotation::SimpleValue->new();
		$sim->value($header{$record});
		$struc->annotation->add_Annotation($record, $sim);
	}
   }
   # store %remark entries as annotations
   if (%remark) {
	my $sim = Bio::Annotation::SimpleValue->new();
	for my $remark_num (keys %remark) {
		$sim->value($remark{$remark_num});
		$struc->annotation->add_Annotation("remark_$remark_num", $sim);
	}
   }
	
   # Coordinate section, the real meat
   #
   #  $_ contains a line beginning with (ATOM|MODEL)

   $buffer = $_;


   if (defined($buffer) && $buffer =~ /^(ATOM |MODEL |HETATM)/ ) {  # can you have an entry without ATOM ?
	until( !defined ($buffer) ) {				 #  (yes : 1a7z )
		   # read in one model at a time
		   my $model = $self->_read_PDB_coordinate_section(\$buffer, $struc);
		   # add this to $struc
		   $struc->add_model($struc, $model);

		   if ($buffer !~ /^MODEL /) { # if we get here we have multiple MODELs
			   last;
		   }
	}	
   } 
   else {
	   $self->throw("Could not find a coordinate section in this record\n");
   }
      

   until( !defined $buffer ) {
	$_ = $buffer;
	
   	# CONNECT records
	if (/^CONECT/) {
		# do not differentiate between different type of connect (column dependant)
		my $conect_unpack = "x6 a5 a5 a5 a5 a5 a5 a5 a5 a5 a5 a5";
		my (@conect) = unpack $conect_unpack, $_;
		for my $k (0 .. $#conect) {
			$conect[$k] =~ s/\s//g;
		}
		my $source = shift @conect;
		for my $conect (@conect) {
			next unless ($conect =~ /^\d+$/);
			$struc->conect($source, $conect);
		}
	}

	# MASTER record
	if (/^MASTER /) {
		# the numbers in here a checksums, we should use them :)
		my ($rol) = unpack "x10 a60", $_;
		$struc->master($rol);
	}

	if (/^END/) {
		# this it the end ...
	}
   
   	$buffer = $self->_readline;
   }
   	

   return $struc;
}

=head2 write_structure

 Title   : write_structure
 Usage   : $stream->write_seq($seq)
 Function: writes the $seq object (must be seq) to the stream
 Returns : 1 for success and 0 for error
 Args    : Bio::Seq


=cut

sub write_structure {
	my ($self) = @_;

	$self->throw("write_structure is not yet implemented, start holding your breath\n");
}


=head2 _filehandle

 Title   : _filehandle
 Usage   : $obj->_filehandle($newval)
 Function: 
 Example : 
 Returns : value of _filehandle
 Args    : newvalue (optional)


=cut

sub _filehandle{
   my ($obj,$value) = @_;
   if( defined $value) {
      $obj->{'_filehandle'} = $value;
    }
    return $obj->{'_filehandle'};

}

=head2 _noatom

 Title   : _noatom
 Usage   : $obj->_noatom($newval)
 Function: 
 Example : 
 Returns : value of _noatom 
 Args    : newvalue (optional)


=cut

sub _noatom{
   my ($obj,$value) = @_;
   if( defined $value) {
      $obj->{'_noatom'} = $value;
    }
    return $obj->{'_noatom'};

}

=head2 _noheader

 Title   : _noheader
 Usage   : $obj->_noheader($newval)
 Function: 
 Example : 
 Returns : value of _noheader 
 Args    : newvalue (optional)


=cut

sub _noheader{
   my ($obj,$value) = @_;
   if( defined $value) {
      $obj->{'_noheader'} = $value;
    }
    return $obj->{'_noheader'};

}

=head2 _read_PDB_singlecontline

 Title   : _read_PDB_singlecontline
 Usage   : $obj->_read_PDB_singlecontline($record, $fromto, $buffer))
 Function: read single continued record from PDB
 Returns : concatenated record entry (between $fromto columns)
 Args    : record, colunm delimiters, buffer

=cut

sub _read_PDB_singlecontline {
	my ($self, $record, $fromto, $buffer) = @_;
	my $concat_line;

	my ($begin, $end) = (split (/-/, $fromto));
	my $unpack_string = "x8 a2 ";
	if($begin == 12) { # one additional space
		$unpack_string .= "x1 a59";
	} else {
		$unpack_string .= "a60";
	}
	$_ = $$buffer;
	while (defined( $_ ||= $self->_readline ) ) {
		if ( /^$record/ ) {
			my($cont, $rol) = unpack $unpack_string, $_;
			if($cont =~ /\d$/ && $begin == 11) { # continuation line
			     			# and text normally at pos 11
		       		$rol =~ s/^\s//; # strip leading space
			}
			## no space (store litteraly) $concat_line .= $rol . " ";
			$concat_line .= $rol; 
		} else {
			last;
		}
	
		$_ = undef;
	}
		$concat_line =~ s/\s$//;  # remove trailing space
	$$buffer = $_;
	
	return $concat_line;
}


=head2 _read_PDB_jrnl

 Title   : _read_PDB_jrnl
 Usage   : $obj->_read_PDB_jrnl($\buffer))
 Function: read jrnl record from PDB
 Returns : Bio::Annotation::Reference object
 Args    : 

=cut

sub _read_PDB_jrnl {
	my ($self, $buffer) = @_;
	
	$_ = $$buffer;
	my ($auth, $titl,$edit,$ref,$publ,$refn);
	while (defined( $_ ||= $self->_readline )) {
		if (/^JRNL /) {
			# this code belgons in a seperate method (shared with
			# remark 1 parsing)
			my ($rec, $subr, $cont, $rol) = unpack "a6 x6 a4 a2 x1 a51", $_;
			$auth = $self->_concatenate_lines($auth,$rol) if ($subr eq "AUTH");
			$titl = $self->_concatenate_lines($titl,$rol) if ($subr eq "TITL");
			$edit = $self->_concatenate_lines($edit,$rol) if ($subr eq "EDIT");
			$ref  = $self->_concatenate_lines($ref ,$rol) if ($subr eq "REF");
			$publ = $self->_concatenate_lines($publ,$rol) if ($subr eq "PUBL");
			$refn = $self->_concatenate_lines($refn,$rol) if ($subr eq "REFN");
		} else {
			last;
		}

		$_ = undef; # trigger reading of next line
	} # while

	$$buffer = $_;
	my $jrnl_ref = Bio::Annotation::Reference->new;

	$jrnl_ref->authors($auth);
	$jrnl_ref->title($titl);
	$jrnl_ref->location($ref);
	$jrnl_ref->publisher($publ);
	$jrnl_ref->editors($edit);
	$jrnl_ref->encoded_ref($refn);
	# XXX KB waht to do with $publ (publisher), $edit (editor) and $refn (ASTM code) ?
	
	return $jrnl_ref;
} # sub _read_PDB_jrnl


=head2 _read_PDB_remark_1

 Title   : _read_PDB_remark_1
 Usage   : $obj->_read_PDB_remark_1($\buffer))
 Function: read "remark 1"  record from PDB
 Returns : array of Bio::Annotation::Reference objects
 Args    : 

=cut

sub _read_PDB_remark_1 {
	my ($self, $buffer) = @_;
	
	$_ = $$buffer;
	my ($auth, $titl,$edit,$ref,$publ,$refn);
	my @refs;

	while (defined( $_ ||= $self->_readline )) {
		if (/^REMARK   1 /) {
			if (/^REMARK   1\s+(\d+)\s*/) {
				my $refnum = $1;
				if ($refnum != 1) { # this is first line of a reference
					my $rref = Bio::Annotation::Reference->new;
					$rref->authors($auth);
					$rref->title($titl);
					$rref->location($ref);
					$rref->publisher($publ);
					$rref->editors($edit);
					$rref->encoded_ref($refn);
					push @refs, $rref;
				}
			} else {
				# this code belgons in a seperate method (shared with
				# remark 1 parsing)
				my ($rec, $subr, $cont, $rol) = unpack "A6 x6 A4 A2 x1 A51", $_;
				$auth = $self->_concatenate_lines($auth,$rol) if ($subr eq "AUTH");
				$titl = $self->_concatenate_lines($titl,$rol) if ($subr eq "TITL");
				$edit = $self->_concatenate_lines($edit,$rol) if ($subr eq "EDIT");
				$ref  = $self->_concatenate_lines($ref ,$rol) if ($subr eq "REF");
				$publ = $self->_concatenate_lines($publ,$rol) if ($subr eq "PUBL");
				$refn = $self->_concatenate_lines($refn,$rol) if ($subr eq "REFN");
			}
		} else {
			# create last reference
                        my $rref = Bio::Annotation::Reference->new;
		        $rref->authors($auth);
		        $rref->title($titl);
		        $rref->location($ref);
			$rref->publisher($publ);
			$rref->editors($edit);
			$rref->encoded_ref($refn);
			push @refs, $rref;
			last;
		}

		$_ = undef; # trigger reading of next line
	} # while

	$$buffer = $_;
	
	return @refs;
} # sub _read_PDB_jrnl


=head2 _read_PDB_coordinate_section

 Title   : _read_PDB_coordinate_section
 Usage   : $obj->_read_PDB_coordinate_section($\buffer))
 Function: read one model from a PDB
 Returns : Bio::Structure::Model object
 Args    : 

=cut

sub _read_PDB_coordinate_section {
	my ($self, $buffer, $struc) = @_;
	my ($model_num, $chain_name, $residue_name, $atom_name);  # to keep track of state
	$model_num = "";
	$chain_name = "";
	$residue_name = "";
	$atom_name = "";

	my $atom_unpack =   "x6 a5 x1 a4 a1 a3 x1 a1 a4 a1 x3 a8 a8 a8 a6 a6 x6 a4 a2 a2";
	my $anisou_unpack = "x6 a5 x1 a4 a1 a3 x1 a1 a4 a1 x1 a7 a7 a7 a7 a7 a7 a4 a2 a2";

	my $model = Bio::Structure::Model->new;
	$model->id('default');
	my $noatom = $self->_noatom;
	my ($chain, $residue, $atom, $old);

	$_ = $$buffer;
	while (defined( $_ ||= $self->_readline )) {
		# start of a new model
		if (/^MODEL\s+(\d+)/) {
			$model_num = $1;
$self->debug("_read_PDB_coor: parsing model $model_num\n");
			$model->id($model_num);
			if (/^MODEL\s+\d+\s+\S+/) { # old format (pre 2.1)
				$old = 1;
			}
		} 

		# ATOM lines, if first set chain
		if (/^(ATOM |HETATM|SIGATM)/) {
			my @line_elements = unpack $atom_unpack, $_;
			for my $k (0 .. $#line_elements) {
				$line_elements[$k] =~ s/^\s+//; # remove leading space
				$line_elements[$k] =~ s/\s+$//; # remove trailing space
				$line_elements[$k] = undef if ($line_elements[$k] =~ /^\s*$/);
			}
			my ($serial, $atomname, $altloc, $resname, $chainID, $resseq, $icode, $x, $y, $z, 
				$occupancy, $tempfactor, $segID, $element, $charge) = @line_elements;
			$chainID = 'default' if ( !defined $chainID ); 
			if ($chainID ne $chain_name) { # new chain
				$chain = Bio::Structure::Chain->new;
				$struc->add_chain($model,$chain);
				$chain->id($chainID);
				$chain_name = $chainID;
			}
			my $res_name_num = $resname."-".$resseq;
			if ($res_name_num ne $residue_name) { # new residue
				$residue = Bio::Structure::Residue->new;
				$struc->add_residue($chain,$residue);
				$residue->id($res_name_num);
				$residue_name = $res_name_num;
				$atom_name = ""; # only needed inside a residue 
			}
			# get out of here if we don't want the atom objects
			if ($noatom) {
				$_ = undef;
				next;
			}
			# alternative location: only take first one
			if ( $altloc && ($altloc =~ /\S+/) && ($atomname eq $atom_name) ) {
				$_ = undef; # trigger reading next line
				next;
			}
			if (/^(ATOM |HETATM)/) { # ATOM  / HETATM
				$atom_name = $atomname;
				$atom = Bio::Structure::Atom->new;
				$struc->add_atom($residue,$atom);
				$atom->id($atomname);
				$atom->serial($serial);
				$atom->icode($icode);
				$atom->x($x);
				$atom->y($y);
				$atom->z($z);
				$atom->occupancy($occupancy);
				$atom->tempfactor($tempfactor);
				# ? segment ID ? (deprecated)
				if (! $old ) {
					$atom->element($element);
					$atom->charge($charge);
				}
			}
			else {  # SIGATM
				my $sigx = $x;
				my $sigy = $y;
				my $sigz = $z;
				my $sigocc = $occupancy;
				my $sigtemp = $tempfactor;
				if ($atom_name ne $atomname) {  # something wrong with PDB file
					$self->throw("A SIGATM record should have the same $atomname as the previous record $atom_name\n");
				}
				$atom->sigx($sigx);
				$atom->sigy($sigy);
				$atom->sigz($sigz);
				$atom->sigocc($sigocc);
				$atom->sigtemp($sigtemp);

			} 
		} # ATOM|HETARM|SIGATM

		# ANISOU | SIGUIJ  lines
		if (/^(ANISOU|SIGUIJ)/) {
			if ($noatom) {
				$_ = undef;
				next;
			}
			my @line_elements = unpack $anisou_unpack, $_;
			for my $k (0 .. $#line_elements) {
				$line_elements[$k] =~ s/^\s+//; # remove leading space
				$line_elements[$k] =~ s/\s+$//; # remove trailing space
				$line_elements[$k] = undef if ($line_elements[$k] =~ /^\s*$/);
			}
			my ($serial, $atomname, $altloc, $resname, $chainID, $resseq, $icode, 
				$u11,$u22, $u33, $u12, $u13, $u23, $segID, $element, $charge) = @line_elements;
$self->debug("read_PDB_coor: parsing ANISOU record: $serial $atomname\n");
			if ( $altloc && ($altloc =~ /\S+/) && ($atomname eq $atom_name) ) {
				$_ = undef;
				next;
			}
			if (/^ANISOU/) {
				if ($atom_name ne $atomname) {  # something wrong with PDB file
					$self->throw("A ANISOU record should have the same $atomname as the previous record $atom_name\n");
				}
				$atom->aniso("u11",$u11);
				$atom->aniso("u22",$u22);
				$atom->aniso("u33",$u33);
				$atom->aniso("u12",$u12);
				$atom->aniso("u13",$u13);
				$atom->aniso("u23",$u23);
			} 
			else { # SIGUIJ
				if ($atom_name ne $atomname) {  # something wrong with PDB file
					$self->throw("A SIGUIJ record should have the same $atomname as the previous record $atom_name\n");
				}
				# could use different variable names, but hey ...
				$atom->aniso("sigu11",$u11);
				$atom->aniso("sigu22",$u22);
				$atom->aniso("sigu33",$u33);
				$atom->aniso("sigu12",$u12);
				$atom->aniso("sigu13",$u13);
				$atom->aniso("sigu23",$u23);
			}
		} # ANISOU | SIGUIJ

		if (/^TER /) {
			$_ = undef;
			next;
		}

		if (/^ENDMDL/) {
			$_ = $self->_readline;
			last;
		}

		if (/^(CONECT|MASTER)/) { # get out of here
			# current line is OK
			last;
		}
		$_ = undef;

	} # while 

	$$buffer = $_;

	return $model;
} # _read_PDB_coordinate_section

1;
