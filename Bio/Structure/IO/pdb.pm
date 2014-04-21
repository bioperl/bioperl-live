#
# BioPerl module for Bio::Structure::IO::pdb
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Kris Boulez <kris.boulez@algonomics.com>
#
# Copyright 2001, 2002 Kris Boulez
#
# Framework is a copy of Bio::SeqIO::embl.pm
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Structure::IO::pdb - PDB input/output stream

=head1 SYNOPSIS

It is probably best not to use this object directly, but
rather go through the Bio::Structure::IO handler system. Go:

    $stream = Bio::Structure::IO->new(-file => $filename,
                                      -format => 'PDB');

    while (my $structure = $stream->next_structure) {
	    # do something with $structure
    }

=head1 DESCRIPTION

This object can transform Bio::Structure objects to and from PDB flat
file databases. The working is similar to that of the Bio::SeqIO handlers.

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

=head1 AUTHOR - Kris Boulez

Email kris.boulez@algonomics.com

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Structure::IO::pdb;
use strict;
use Bio::Structure::Entry;
#use Bio::Structure::Model;
#use Bio::Structure::Chain;
#use Bio::Structure::Residue;
use Bio::Structure::Atom;
use Bio::SeqFeature::Generic;
use Bio::Annotation::Reference;

use base qw(Bio::Structure::IO);

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
 Usage   : $struc = $stream->next_structure()
 Function: returns the next structure in the stream
 Returns : Bio::Structure object
 Args    :

=cut

sub next_structure {
   my ($self,@args) = @_;
   my ($line);
   my ($obslte, $title, $caveat, $compnd, $source, $keywds,
	$expdta, $author, %revdat, $revdat, $sprsde, $jrnl, %remark, $dbref,
	$turn, $ssbond, $link, $hydbnd, $sltbrg, $cispep,
	$site, $cryst1, $tvect,);
   my $struc = Bio::Structure::Entry->new(-id => 'created from pdb.pm');
   my $all_headers = ( !$self->_noheader );  # we'll parse all headers and store as annotation
   my %header;  # stores all header RECORDs an is stored as annotations when ATOM is reached


   $line = $self->_readline;   # This needs to be before the first eof() test

   if( !defined $line ) {
       return; # no throws - end of file
   }

   if( $line =~ /^\s+$/ ) {
       while( defined ($line = $self->_readline) ) {
	   $line =~/\S/ && last;
       }
   }
   if( !defined $line ) {
       return; # end of file
   }
   $line =~ /^HEADER\s+\S+/ || $self->throw("PDB stream with no HEADER. Not pdb in my book");
   my($header_line) = unpack "x10 a56", $line;
   $header{'header'} = $header_line;
   my($class, $depdate, $idcode) = unpack "x10 a40 a9 x3 a4", $line;
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
		$header{'jrnl'} = 1; # when writing out, we need a way to check there was a JRNL record (not mandatory)
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
			# $_ still holds the REMARK_1 line, $buffer now contains the first non
			#  REMARK_1 line. We need to parse it in this pass (so no else block)
			$_ = $buffer;
		}
		# for the moment I don't see a better solution (other then using goto)
		if (/^REMARK\s+(\d+)\s*/) {
			my $r_num = $1;
			if ($r_num != 1) { # other remarks, we store literlly at the moment
				my ($rol) = unpack "x11 a59", $_;
				$remark{$r_num} .= $rol;
			}
		}
	} # REMARK

	# DBREF line(s)
	#  references to sequences in other databases
	#  we store as 'dblink' annotations and whole line as simple annotation (round-trip)
	if (/^DBREF / && $all_headers) {
		my ($rol) = unpack "x7 a61", $_;
		$dbref .= $rol;
		$header{'dbref'} = $dbref;
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
		$header{'seqadv'} .= $rol;
	} # SEQADV

	# SEQRES line(s)
	#  this is (I think) the sequence of macromolecule that was analysed
	#  this will be returned when doing $struc->seq
	if (/^SEQRES / && $all_headers) {
		my ($rol) = unpack "x8 a62", $_;
		$header{'seqres'} .= $rol;
	} # SEQRES

	# MODRES line(s)
	if (/^MODRES / && $all_headers) {
		my ($rol) = unpack "x7 a63", $_;
		$header{'modres'} .= $rol;
	} # MODRES

	# HET line(s)
	if (/^HET / && $all_headers) {
		my ($rol) = unpack "x7 a63", $_;
		$header{'het'} .= $rol;
	} # HET

	# HETNAM line(s)
	if (/^HETNAM / && $all_headers) {
		my ($rol) = unpack "x8 a62", $_;
		$header{'hetnam'} .= $rol;
	} # HETNAM

	# HETSYN line(s)
	if (/^HETSYN / && $all_headers) {
		my ($rol) = unpack "x8 a62", $_;
		$header{'hetsyn'} .= $rol;
	} # HETSYN

	# FORMUL line(s)
	if (/^FORMUL / && $all_headers) {
		my ($rol) = unpack "x8 a62", $_;
		$header{'formul'} .= $rol;
	} # FORMUL

	# HELIX line(s)
	#  store as specific object ??
	if (/^HELIX / && $all_headers) {
		my ($rol) = unpack "x7 a69", $_;
		$header{'helix'} .= $rol;
	} # HELIX

	# SHEET line(s)
	#  store as specific object ??
	if (/^SHEET / && $all_headers) {
		my ($rol) = unpack "x7 a63", $_;
		$header{'sheet'} .= $rol;
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
	for my $remark_num (keys %remark) {
		my $sim = Bio::Annotation::SimpleValue->new();
		$sim->value($remark{$remark_num});
		$struc->annotation->add_Annotation("remark_$remark_num", $sim);
	}
   }

   # Coordinate section, the real meat
   #
   #  $_ contains a line beginning with (ATOM|MODEL)

   $buffer = $_;


   if (defined($buffer) && $buffer =~ /^(ATOM |MODEL |HETATM)/ ) {  # can you have an entry without ATOM ?
	while( defined ($buffer) ) {				 #  (yes : 1a7z )
		   # read in one model at a time
		   my $model = $self->_read_PDB_coordinate_section(\$buffer, $struc);
		   # add this to $struc
		   $struc->add_model($struc, $model);

		   if ($buffer && $buffer !~ /^MODEL /) { # if we get here we have multiple MODELs
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
		my $type;
		for my $k (0 .. 9) {
			next unless ($conect[$k] =~ /^\d+$/);
			# 0..3 		bond
			if( $k <= 3 ) {
				$type = "bond";
			}
			# 4..5,7..8 	hydrogen bonded
			elsif( ($k >= 4 && $k <= 5) || ($k >= 7 && $k <= 8) ) {
				$type = "hydrogen";
			}
			# 6, 9		salt bridged
			elsif( $k == 6 || $k == 9 ) {
				$type = "saltbridged";
			} else {
				$self->throw("k has impossible value ($k), check brain");
			}
			$struc->conect($source, $conect[$k], $type);
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
 Usage   : $stream->write_structure($struc)
 Function: writes the $struc object (must be a Bio::Structure) to the stream
 Returns : 1 for success and 0 for error
 Args    : Bio::Structure object

=cut

sub write_structure {
	my ($self, $struc) = @_;
	if( !defined $struc ) {
		$self->throw("Attempting to write with no structure!");
	}

	if( ! ref $struc || ! $struc->isa('Bio::Structure::StructureI') ) {
		$self->throw(" $struc is not a StructureI compliant module.");
	}
	my ($ann, $string, $output_string, $key);
	# HEADER
	($ann) = $struc->annotation->get_Annotations("header");
	if (defined $ann) {
		$string = $ann->as_text;
		$string =~ s/^Value: //;
		$output_string = pack ("A10 A56", "HEADER", $string);
	} else {	# not read in via read_structure, create HEADER line
		my $id = $struc->id;
		if (!$id) {
			$id = "UNK1";
		}
		if (length($id) > 4) {
			$id = substr($id,0,4);
		}
		my $classification = "DEFAULT CLASSIFICATION";
		my $dep_date       = "24-JAN-70";
		$output_string = pack ("A10 A40 A12 A4", "HEADER", $classification, $dep_date, $id);
	}
	$output_string .= " " x (80 - length($output_string) );
	$self->_print("$output_string\n");

	my (%header);
	for  $key ($struc->annotation->get_all_annotation_keys) {
		$header{$key} = 1;;
	}

	exists $header{'obslte'} && $self->_write_PDB_simple_record(-name => "OBSLTE  ", -cont => "9-10",
					-annotation => $struc->annotation->get_Annotations("obslte"), -rol => "11-70");

	exists $header{'title'} && $self->_write_PDB_simple_record(-name => "TITLE   ", -cont => "9-10",
					-annotation => $struc->annotation->get_Annotations("title"), -rol => "11-70");

	exists $header{'caveat'} && $self->_write_PDB_simple_record(-name => "CAVEAT  ", -cont => "9-10",
					-annotation => $struc->annotation->get_Annotations("caveat"), -rol => "12-70");

	exists $header{'compnd'} && $self->_write_PDB_simple_record(-name => "COMPND  ", -cont => "9-10",
					-annotation => $struc->annotation->get_Annotations("compnd"), -rol => "11-70");

	exists $header{'source'} && $self->_write_PDB_simple_record(-name => "SOURCE  ", -cont => "9-10",
					-annotation => $struc->annotation->get_Annotations("source"), -rol => "11-70");

	exists $header{'keywds'} && $self->_write_PDB_simple_record(-name => "KEYWDS  ", -cont => "9-10",
					-annotation => $struc->annotation->get_Annotations("keywds"), -rol => "11-70");

	exists $header{'expdta'} && $self->_write_PDB_simple_record(-name => "EXPDTA  ", -cont => "9-10",
					-annotation => $struc->annotation->get_Annotations("expdta"), -rol => "11-70");

	exists $header{'author'} && $self->_write_PDB_simple_record(-name => "AUTHOR  ", -cont => "9-10",
					-annotation => $struc->annotation->get_Annotations("author"), -rol => "11-70");

	exists $header{'revdat'} && $self->_write_PDB_simple_record(-name => "REVDAT ",
					-annotation => $struc->annotation->get_Annotations("revdat"), -rol => "8-66");

	exists $header{'sprsde'} && $self->_write_PDB_simple_record(-name => "SPRSDE  ", -cont => "9-10",
					-annotation => $struc->annotation->get_Annotations("sprsde"), -rol => "12-70");

	# JRNL en REMARK 1
	my ($jrnl_done, $remark_1_counter);
	if ( !exists $header{'jrnl'} ) {
		$jrnl_done = 1;
	}
	foreach my $ref ($struc->annotation->get_Annotations('reference') ) {
		if( !$jrnl_done ) { # JRNL record
			$ref->authors && $self->_write_PDB_simple_record(-name => "JRNL        AUTH",
					-cont => "17-18", -rol => "20-70", -string => $ref->authors );
			$ref->title && $self->_write_PDB_simple_record(-name => "JRNL        TITL",
					-cont => "17-18", -rol => "20-70", -string => $ref->title );
			$ref->editors && $self->_write_PDB_simple_record(-name => "JRNL        EDIT",
					-cont => "17-18", -rol => "20-70", -string => $ref->editors );
			$ref->location && $self->_write_PDB_simple_record(-name => "JRNL        REF ",
					-cont => "17-18", -rol => "20-70", -string => $ref->location );
			$ref->editors && $self->_write_PDB_simple_record(-name => "JRNL        EDIT",
					-cont => "17-18", -rol => "20-70", -string => $ref->editors );
			$ref->encoded_ref && $self->_write_PDB_simple_record(-name => "JRNL        REFN",
					-cont => "17-18", -rol => "20-70", -string => $ref->encoded_ref );
			$jrnl_done = 1;
		} else { # REMARK 1
			if (!$remark_1_counter) { # header line
				my $remark_1_header_line = "REMARK   1" . " " x 70;
				$self->_print("$remark_1_header_line\n");
				$remark_1_counter = 1;
			}
			# per reference header
			my $rem_line = "REMARK   1 REFERENCE " . $remark_1_counter;
			$rem_line .= " " x (80 - length($rem_line) );
			$self->_print($rem_line,"\n");
			$ref->authors && $self->_write_PDB_simple_record(-name => "REMARK   1  AUTH",
					-cont => "17-18", -rol => "20-70", -string => $ref->authors );
			$ref->title && $self->_write_PDB_simple_record(-name => "REMARK   1  TITL",
					-cont => "17-18", -rol => "20-70", -string => $ref->title );
			$ref->editors && $self->_write_PDB_simple_record(-name => "REMARK   1  EDIT",
					-cont => "17-18", -rol => "20-70", -string => $ref->editors );
			$ref->location && $self->_write_PDB_simple_record(-name => "REMARK   1  REF ",
					-cont => "17-18", -rol => "20-70", -string => $ref->location );
			$ref->editors && $self->_write_PDB_simple_record(-name => "REMARK   1  EDIT",
					-cont => "17-18", -rol => "20-70", -string => $ref->editors );
			$ref->encoded_ref && $self->_write_PDB_simple_record(-name => "REMARK   1  REFN",
					-cont => "17-18", -rol => "20-70", -string => $ref->encoded_ref );
			$remark_1_counter++;
		}
	}
	if (! defined $remark_1_counter ) { 	# no remark 1 record written yet
		my $remark_1_header_line = "REMARK   1" . " " x 70;
		$self->_print("$remark_1_header_line\n");  # write dummy  (we need this line)
	}

	# REMARK's  (not 1 at the moment, references)
	my (%remarks, $remark_num);
	for  $key (keys %header) {
		next unless ($key =~ /^remark_(\d+)$/);
		next if ($1 == 1);
		$remarks{$1} = 1;
	}
	for $remark_num (sort {$a <=> $b} keys %remarks) {
		$self->_write_PDB_remark_record($struc, $remark_num);
	}

	exists $header{'dbref'} && $self->_write_PDB_simple_record(-name =>  "DBREF  ",
					-annotation => $struc->annotation->get_Annotations("dbref"), -rol => "8-68");
	exists $header{'seqadv'} && $self->_write_PDB_simple_record(-name => "SEQADV ",
					-annotation => $struc->annotation->get_Annotations("seqadv"), -rol => "8-70");
	exists $header{'seqres'} && $self->_write_PDB_simple_record(-name => "SEQRES  ",
					-annotation => $struc->annotation->get_Annotations("seqres"), -rol => "9-70");
	exists $header{'modres'} && $self->_write_PDB_simple_record(-name => "MODRES ",
					-annotation => $struc->annotation->get_Annotations("modres"), -rol => "8-70");
	exists $header{'het'} && $self->_write_PDB_simple_record(-name =>    "HET    ",
					-annotation => $struc->annotation->get_Annotations("het"), -rol => "8-70");
	exists $header{'hetnam'} && $self->_write_PDB_simple_record(-name => "HETNAM  ",
					-annotation => $struc->annotation->get_Annotations("hetnam"), -rol => "9-70");
	exists $header{'hetsyn'} && $self->_write_PDB_simple_record(-name => "HETSYN  ",
					-annotation => $struc->annotation->get_Annotations("hetsyn"), -rol => "9-70");
	exists $header{'formul'} && $self->_write_PDB_simple_record(-name => "FORMUL  ",
					-annotation => $struc->annotation->get_Annotations("formul"), -rol => "9-70");
	exists $header{'helix'} && $self->_write_PDB_simple_record(-name =>  "HELIX  ",
					-annotation => $struc->annotation->get_Annotations("helix"), -rol => "8-76");
	exists $header{'sheet'} && $self->_write_PDB_simple_record(-name =>  "SHEET  ",
					-annotation => $struc->annotation->get_Annotations("sheet"), -rol => "8-70");
	exists $header{'turn'} && $self->_write_PDB_simple_record(-name =>   "TURN   ",
					-annotation => $struc->annotation->get_Annotations("turn"), -rol => "8-70");
	exists $header{'ssbond'} && $self->_write_PDB_simple_record(-name => "SSBOND ",
					-annotation => $struc->annotation->get_Annotations("ssbond"), -rol => "8-72");
	exists $header{'link'} && $self->_write_PDB_simple_record(-name =>   "LINK        ",
					-annotation => $struc->annotation->get_Annotations("link"), -rol => "13-72");
	exists $header{'hydbnd'} && $self->_write_PDB_simple_record(-name => "HYDBND      ",
					-annotation => $struc->annotation->get_Annotations("hydbnd"), -rol => "13-72");
	exists $header{'sltbrg'} && $self->_write_PDB_simple_record(-name => "SLTBRG      ",
					-annotation => $struc->annotation->get_Annotations("sltbrg"), -rol => "13-72");
	exists $header{'cispep'} && $self->_write_PDB_simple_record(-name => "CISPEP ",
					-annotation => $struc->annotation->get_Annotations("cispep"), -rol => "8-59");
	exists $header{'site'} && $self->_write_PDB_simple_record(-name =>   "SITE   ",
					-annotation => $struc->annotation->get_Annotations("site"), -rol => "8-61");
	exists $header{'cryst1'} && $self->_write_PDB_simple_record(-name => "CRYST1",
					-annotation => $struc->annotation->get_Annotations("cryst1"), -rol => "7-70");
	for my $k (1..3) {
		my $origxn = "origx".$k;
		my $ORIGXN = uc($origxn)."    ";
		exists $header{$origxn} && $self->_write_PDB_simple_record(-name => $ORIGXN,
			-annotation => $struc->annotation->get_Annotations($origxn), -rol => "11-55");
	}
	for my $k (1..3) {
		my $scalen = "scale".$k;
		my $SCALEN = uc($scalen)."    ";
		exists $header{$scalen} && $self->_write_PDB_simple_record(-name => $SCALEN,
			-annotation => $struc->annotation->get_Annotations($scalen), -rol => "11-55");
	}
	for my $k (1..3) {
		my $mtrixn = "mtrix".$k;
		my $MTRIXN = uc($mtrixn)." ";
		exists $header{$mtrixn} && $self->_write_PDB_simple_record(-name => $MTRIXN,
			-annotation => $struc->annotation->get_Annotations($mtrixn), -rol => "8-60");
	}
	exists $header{'tvect'} && $self->_write_PDB_simple_record(-name => "TVECT  ",
					-annotation => $struc->annotation->get_Annotations("tvect"), -rol => "8-70");

	# write out coordinate section
	#
	my %het_res;  # hetero residues
	$het_res{'HOH'} = 1;  # water is default
	if (exists $header{'het'}) {
		my ($het_line) = ($struc->annotation->get_Annotations("het"))[0]->as_text;
		$het_line =~ s/^Value: //;
		for ( my $k = 0; $k <= length $het_line ; $k += 63) {
			my $l = substr $het_line, $k, 63;
			$l =~ s/^\s*(\S+)\s+.*$/$1/;
			$het_res{$l} = 1;
		}
	}
	for my $model ($struc->get_models) {
		# more then one model ?
		if ($struc->get_models > 1) {
			my $model_line = sprintf("MODEL     %4d", $model->id);
			$model_line .= " " x (80 - length($model_line) );
			$self->_print($model_line, "\n");
		}
		for my $chain ($struc->get_chains($model)) {
			my ($residue, $atom, $resname, $resnum, $atom_line, $atom_serial, $atom_icode, $chain_id);
			my ($prev_resname, $prev_resnum, $prev_atomicode); # need these for TER record
			my $last_record = ""; # Used to spot an ATOM -> HETATM change within a chain
			$chain_id = $chain->id;
			if ( $chain_id eq "default" ) {
				$chain_id = " ";
			}
			$self->debug("model_id: $model->id chain_id: $chain_id\n");
			for $residue ($struc->get_residues($chain)) {
				($resname, $resnum) = split /-/, $residue->id;
				for $atom ($struc->get_atoms($residue)) {
					if ($het_res{$resname}) {  # HETATM
						if ( $resname ne "HOH" && $last_record eq "ATOM  " ) {
							# going from ATOM -> HETATM, we have to write TER
							my $ter_line = "TER   ";
							$ter_line .= sprintf("%5d", $atom_serial + 1);
							$ter_line .= "      ";
							$ter_line .= sprintf("%3s ", $prev_resname);
							$ter_line .= $chain_id;
							$ter_line .= sprintf("%4d", $prev_resnum);
							$ter_line .= $atom_icode ? $prev_atomicode : " "; # 27
							$ter_line .= " " x (80 - length $ter_line);  # extend to 80 chars
							$self->_print($ter_line,"\n");
						}
						$atom_line = "HETATM";
					} else {
						$atom_line = "ATOM  ";
					}
					$last_record = $atom_line;
					$atom_line .= sprintf("%5d ", $atom->serial);
					$atom_serial = $atom->serial; # we need it for TER record
					$atom_icode = $atom->icode;
					# remember some stuff if next iteration needs writing TER
					$prev_resname = $resname;
					$prev_resnum  = $resnum;
					$prev_atomicode = $atom_icode;
					# getting the name of the atom correct is subtrivial
					my $atom_id = $atom->id;
					# is pdb_atomname set, then use this (most probably set when
					# reading in the PDB record)
					my $pdb_atomname = $atom->pdb_atomname;
					if( defined $pdb_atomname ) {
						$atom_line .= sprintf("%-4s", $pdb_atomname);
					} else {
						# start (educated) guessing
						my $element = $atom->element;
						if( defined $element && $element ne "H") {
							# element should be at first two positions (right justified)
							# ie. Calcium should be "CA  "
							#     C alpha should be " CA "
							if( length($element) == 2 ) {
								$atom_line .= sprintf("%-4s", $atom->id);
							} else {
								$atom_line .= sprintf(" %-3s", $atom->id);
							}
						} else { # old behaviour do a best guess
							if ($atom->id =~ /^\dH/) { # H: four positions, left justified
								$atom_line .= sprintf("%-4s", $atom->id);
							} elsif (length($atom_id) == 4) {
								if ($atom_id =~ /^(H\d\d)(\d)$/) {  # turn H123 into 3H12
									$atom_line .= $2.$1;
								} else {	# no more guesses, no more alternatives
									$atom_line .= $atom_id;
								}
							} else { # if we get here and it is not correct let me know
								$atom_line .= sprintf(" %-3s", $atom->id);
							}
						}
					}
					# we don't do alternate location at this moment
					$atom_line .= " "; 				# 17
					$atom_line .= sprintf("%3s",$resname);		# 18-20
					$atom_line .= " ".$chain_id; 			# 21, 22
					$atom_line .= sprintf("%4d", $resnum); 		# 23-26
					$atom_line .= $atom->icode ? $atom->icode : " "; # 27
					$atom_line .= "   ";				# 28-30
					$atom_line .= sprintf("%8.3f", $atom->x);	# 31-38
					$atom_line .= sprintf("%8.3f", $atom->y);	# 39-46
					$atom_line .= sprintf("%8.3f", $atom->z);	# 47-54
					$atom_line .= sprintf("%6.2f", $atom->occupancy); # 55-60
					$atom_line .= sprintf("%6.2f", $atom->tempfactor); # 61-66
					$atom_line .= "      ";				# 67-72
					$atom_line .= $atom->segID ? 			# segID 73-76
							sprintf("%-4s",  $atom->segID) :
							"    ";
					$atom_line .= $atom->element ?
							sprintf("%2s", $atom->element) :
							"  ";
					$atom_line .= $atom->charge ?
							sprintf("%2s", $atom->charge) :
							"  ";

					$self->_print($atom_line,"\n");
				}
			}
			# write out TER record
			if ( $resname ne "HOH" ) {
				my $ter_line = "TER   ";
				$ter_line .= sprintf("%5d", $atom_serial + 1);
				$ter_line .= "      ";
				$ter_line .= sprintf("%3s ", $resname);
				$ter_line .= $chain_id;
				$ter_line .= sprintf("%4d", $resnum);
				$ter_line .= $atom_icode ? $atom_icode : " "; # 27
				$ter_line .= " " x (80 - length $ter_line);  # extend to 80 chars
				$self->_print($ter_line,"\n");
			}
		}
		if ($struc->get_models > 1) { # we need ENDMDL
			my $endmdl_line = "ENDMDL" . " " x 74;
			$self->_print($endmdl_line, "\n");
		}
	} # for my $model

	# CONECT
	my @sources = $struc->get_all_conect_source;
	my ($conect_line,@conect, @bond, @hydbond, @saltbridge, $to, $type);
	for my $source (@sources) {
		# get all conect's
		my @conect = $struc->conect($source);
		# classify
		for my $con (@conect) {
			($to, $type) = split /_/, $con;
			if($type eq "bond") {
				push @bond, $to;
			} elsif($type eq "hydrogenbonded") {
				push @hydbond, $to;
			} elsif($type eq "saltbridged") {
				push @saltbridge, $to;
			} else {
				$self->throw("type $type is unknown for conect");
			}
		}
		# and write out CONECT lines as long as there is something
		# in one of the arrays
		while ( @bond || @hydbond ||  @saltbridge) {
			my ($b, $hb, $sb);
			$conect_line = "CONECT". sprintf("%5d", $source);
			for my $k (0..3) {
				$b = shift @bond;
				$conect_line .= $b ? sprintf("%5d", $b) : "    ";
			}
			for my $k (4..5) {
				$hb = shift @hydbond;
				$conect_line .= $hb ? sprintf("%5d", $hb) : "    ";
			}
			$sb = shift @saltbridge;
			$conect_line .= $sb ? sprintf("%5d", $sb) : "    ";
			for my $k (7..8) {
				$hb = shift @hydbond;
				$conect_line .= $hb ? sprintf("%5d", $hb) : "    ";
			}
			$sb = shift @saltbridge;
			$conect_line .= $sb ? sprintf("%5d", $sb) : "    ";

			$conect_line .= " " x (80 - length($conect_line) );
			$self->_print($conect_line, "\n");
		}
	}

	# MASTER line contains checksums, we should calculate them of course :)
	my $master_line = "MASTER    " . $struc->master;
	$master_line .= " " x (80 - length($master_line) );
	$self->_print($master_line, "\n");

	my $end_line = "END" . " " x 77;
	$self->_print($end_line,"\n");

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
	my ($auth, $titl,$edit,$ref,$publ,$refn, $pmid, $doi);
	while (defined( $_ ||= $self->_readline )) {
		if (/^JRNL /) {
			# this code belgons in a seperate method (shared with
			# remark 1 parsing)
			my ($rec, $subr, $cont, $rol) = unpack "A6 x6 A4 A2 x1 A51", $_;
			$auth = $self->_concatenate_lines($auth,$rol) if ($subr eq "AUTH");
			$titl = $self->_concatenate_lines($titl,$rol) if ($subr eq "TITL");
			$edit = $self->_concatenate_lines($edit,$rol) if ($subr eq "EDIT");
			$ref  = $self->_concatenate_lines($ref ,$rol) if ($subr eq "REF");
			$publ = $self->_concatenate_lines($publ,$rol) if ($subr eq "PUBL");
			$refn = $self->_concatenate_lines($refn,$rol) if ($subr eq "REFN");
			$pmid = $self->_concatenate_lines($pmid,$rol) if ($subr eq "PMID");
			$doi = $self->_concatenate_lines($doi,$rol) if ($subr eq "DOI");
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
	$jrnl_ref->pubmed($pmid);
	$jrnl_ref->doi($doi);

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
	my ($auth, $titl,$edit,$ref,$publ,$refn,$refnum,$pmid, $doi);
	my @refs;

	while (defined( $_ ||= $self->_readline )) {
		if (/^REMARK   1 /) {
			if (/^REMARK   1\s+REFERENCE\s+(\d+)\s*/) {
				$refnum = $1;
				if ($refnum != 1) { # this is first line of a reference
					my $rref = Bio::Annotation::Reference->new;
					$rref->authors($auth);
					$rref->title($titl);
					$rref->location($ref);
					$rref->publisher($publ);
					$rref->editors($edit);
					$rref->encoded_ref($refn);
					$rref->pubmed($pmid);
					$rref->doi($doi);
					$auth = $titl = $edit = $ref = $publ = $refn = undef;
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
				$pmid = $self->_concatenate_lines($pmid,$rol) if ($subr eq "PMID");
				$doi = $self->_concatenate_lines($doi,$rol) if ($subr eq "DOI");
			}
		} else {
			# have we seen any reference at all (could be single REMARK  1 line
			if ( ! defined ($refnum) ) {
				last; # get out of while()
			}

			# create last reference
                        my $rref = Bio::Annotation::Reference->new;
		        $rref->authors($auth);
		        $rref->title($titl);
		        $rref->location($ref);
			$rref->publisher($publ);
			$rref->editors($edit);
			$rref->encoded_ref($refn);
			$rref->pubmed($pmid);
			$rref->doi($doi);
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
	my (%_ch_in_model);  # which chains are already in this model

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
		# old hier ook setten XXX
		# ATOM lines, if first set chain
		if (/^(ATOM |HETATM|SIGATM)/) {
			my @line_elements = unpack $atom_unpack, $_;
			my $pdb_atomname = $line_elements[1]; # need to get this before removing spaces
			for my $k (0 .. $#line_elements) {
				$line_elements[$k] =~ s/^\s+//; # remove leading space
				$line_elements[$k] =~ s/\s+$//; # remove trailing space
				$line_elements[$k] = undef if ($line_elements[$k] =~ /^\s*$/);
			}
			my ($serial, $atomname, $altloc, $resname, $chainID, $resseq, $icode, $x, $y, $z,
				$occupancy, $tempfactor, $segID, $element, $charge) = @line_elements;
			$chainID = 'default' if ( !defined $chainID );
			if ($chainID ne $chain_name) { # possibly a new chain
				# fix for bug #1187
				#  we can have ATOM/HETATM of an already defined chain (A B A B)
				#  e.g. 1abm

				if (exists $_ch_in_model{$chainID} ) { # we have already seen this chain in this model
					$chain = $_ch_in_model{$chainID};
				} else {  # we create a new chain
					$chain = Bio::Structure::Chain->new;
					$struc->add_chain($model,$chain);
					$chain->id($chainID);
					$_ch_in_model{$chainID} = $chain;
				}
				$chain_name = $chain->id;
			}
			#my $res_name_num = $resname."-".$resseq;
			my $res_name_num = $resname."-".$resseq;
			$res_name_num .= '.'.$icode if $icode;
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
				$atom->pdb_atomname($pdb_atomname); # store away PDB atomname for writing out
				$atom->serial($serial);
				$atom->icode($icode);
				$atom->x($x);
				$atom->y($y);
				$atom->z($z);
				$atom->occupancy($occupancy);
				$atom->tempfactor($tempfactor);
				$atom->segID($segID); # deprecated but used by people
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


sub _write_PDB_simple_record {
	my ($self, @args) = @_;
	my ($name, $cont , $annotation, $rol, $string) =
		$self->_rearrange([qw(
				NAME
				CONT
				ANNOTATION
				ROL
				STRING
			)],
			@args);
	if (defined $string && defined $annotation) {
		$self->throw("you can only supply one of -annoation or -string");
	}
	my ($output_string, $ann_string, $t_string);
	my ($rol_begin, $rol_end) = $rol =~ /^(\d+)-(\d+)$/;
	my $rol_length = $rol_end - $rol_begin +1;
	if ($string) {
		if (length $string > $rol_length) {
			# we might need to split $string in multiple lines
			while (length $string > $rol_length) {
				# other option might be to go for a bunch of substr's
				my @c = split//,$string;
				my $t = $rol_length; # index into @c
				while ($c[$t] ne " ") { # find first space, going backwards
$self->debug("c[t]: $c[$t] $t\n");
					$t--;
					if ($t == 0) { $self->throw("Found no space for $string\n"); }
				}
$self->debug("t: $t rol_length: $rol_length\n");
				$ann_string .= substr($string, 0, $t);
$self->debug("ann_string: $ann_string\n");
				$ann_string .= " " x ($rol_length - $t );
				$string = substr($string, $t+1);
				$string =~ s/^\s+//;
$self->debug("ann_string: $ann_string~~\nstring: $string~~\n");
			}
			$ann_string .= $string;
		} else {
			$ann_string = $string;
		}
	} else {
		$ann_string = $annotation->as_text;
		$ann_string =~ s/^Value: //;
	}
	# ann_string contains the thing to write out, writing out happens below
	my $ann_length = length $ann_string;

$self->debug("ann_string: $ann_string\n");
	if ($cont) {
		my ($c_begin, $c_end) = $cont =~ /^(\d+)-(\d+)$/;
		if ( $ann_length > $rol_length ) { # we need to continuation lines
			my $first_line = 1;
			my $cont_number = 2;
			my $out_line;
			my $num_pos = $rol_length;
			my $i = 0;
			while( $i < $ann_length ) {
				$t_string = substr($ann_string, $i, $num_pos);
$self->debug("t_string: $t_string~~$i $num_pos\n");
				if ($first_line) {
					$out_line = $name . " " x ($rol_begin - $c_begin) . $t_string;
					$out_line .= " " x (80 - length($out_line) ) . "\n";
					$first_line = 0;
					$output_string = $out_line;
					$i += $num_pos;	# first do counter
					if ($rol_begin - $c_end == 1) { # next line one character less
						$num_pos--;
					}
				} else {
					$out_line = $name . sprintf("%2d",$cont_number);
					# a space after continuation number
					if ($rol_begin - $c_end == 1) {  # one space after cont number
						$out_line .= " ";
						$out_line .=  $t_string;
					} else {
						$out_line .= " " x ($rol_begin - $c_end - 1) . $t_string;
					}
					$out_line .= " " x (80 -length($out_line) ) . "\n";
					$cont_number++;
					$output_string .= $out_line;
					$i += $num_pos;
				}
			}
		} else { # no continuation
			my $spaces = $rol_begin - $c_begin; # number of spaces need to insert
			$output_string = $name . " " x $spaces . $ann_string;
			$output_string .= " " x (80 - length($output_string) );
		}
	} else { # no contintuation lines
		if ($ann_length < $rol_length) {
			$output_string = $name . $ann_string;
			$output_string .= " " x (80 - length($output_string) );
		} else {
			for (my $i = 0; $i < $ann_length; $i += $rol_length) {
				my $out_line;
				$t_string = substr($ann_string, $i, $rol_length);
				$out_line = $name . $t_string;
				$out_line .= " " x (80 -length($out_line) ) . "\n";
				$output_string .= $out_line;
			}
		}
	}
	$output_string =~ s/\n$//;  # remove trailing newline
	$self->_print("$output_string\n");

}

sub _write_PDB_remark_record {
	my ($self, $struc, $remark_num) = @_;
	my ($ann) = $struc->annotation->get_Annotations("remark_$remark_num");
	my $name = sprintf("REMARK %3d ",$remark_num);
	$self->_write_PDB_simple_record(-name => $name, -annotation => $ann, -rol => "12-70");
}

1;
