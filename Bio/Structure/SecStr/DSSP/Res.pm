# $Header$
#
# bioperl module for Bio::Structure::SecStr::DSSP::Res.pm
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Ed Green <ed@compbio.berkeley.edu>
#
# Copyright Univ. of California
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=head1 NAME

Bio::Structure::SecStr::DSSP::Res - Module for parsing/accessing dssp output

=head1 SYNOPSIS

  my $dssp_obj = Bio::Structure::SecStr::DSSP::Res->new('-file'=>'filename.dssp');

  # or

  my $dssp_obj = Bio::Structure::SecStr::DSSP::Res->new('-fh'=>\*STDOUT);

  # get DSSP defined Secondary Structure for residue 20
  $sec_str = $dssp_obj->resSecStr( 20 );

  # get dssp defined sec. structure summary for PDB residue  # 10 of chain A

  $sec_str = $dssp_obj->resSecStrSum( '10:A' );

=head1 DESCRIPTION

DSSP::Res is a module for objectifying DSSP output.  Methods are then
available for extracting all the information within the output file
and convenient subsets of it.
The principal purpose of DSSP is to determine secondary structural
elements of a given structure.

    ( Dictionary of protein secondary structure: pattern recognition
      of hydrogen-bonded and geometrical features.
      Biopolymers. 1983 Dec;22(12):2577-637. )

The DSSP program is available from:
  http://www.cmbi.kun.nl/swift/dssp

This information is available on a per residue basis ( see resSecStr
and resSecStrSum methods ) or on a per chain basis ( see secBounds
method ).

resSecStr() & secBounds() return one of the following:
    'H' = alpha helix
    'B' = residue in isolated beta-bridge
    'E' = extended strand, participates in beta ladder
    'G' = 3-helix (3/10 helix)
    'I' = 5 helix (pi helix)
    'T' = hydrogen bonded turn
    'S' = bend
    ''  = no assignment

A more general classification is returned using the resSecStrSum()
method.  The purpose of this is to have a method for DSSP and STRIDE
derived output whose range is the same.
Its output is one of the following:

    'H' = helix         ( => 'H', 'G', or 'I' from above )
    'B' = beta          ( => 'B' or 'E' from above )
    'T' = turn          ( => 'T' or 'S' from above )
    ' ' = no assignment ( => ' ' from above )

The methods are roughly divided into 3 sections:
1.  Global features of this structure (PDB ID, total surface area,
    etc.).  These methods do not require an argument.
2.  Residue specific features ( amino acid, secondary structure,
    solvent exposed surface area, etc. ).  These methods do require an
    argument.  The argument is supposed to uniquely identify a
    residue described within the structure.  It can be of any of the
    following forms:
    ('#A:B') or ( #, 'A', 'B' )
      || |
      || - Chain ID (blank for single chain)
      |--- Insertion code for this residue.  Blank for most residues.
      |--- Numeric portion of residue ID.

    (#)
     |
     --- Numeric portion of residue ID.  If there is only one chain and
         it has no ID AND there is no residue with an insertion code at this
         number, then this can uniquely specify a residue.

    ('#:C') or ( #, 'C' )
      | |
      | -Chain ID
      ---Numeric portion of residue ID.

  If a residue is incompletely specified then the first residue that
  fits the arguments is returned.  For example, if 19 is the argument
  and there are three chains, A, B, and C with a residue whose number
  is 19, then 19:A will be returned (assuming its listed first).

  Since neither DSSP nor STRIDE correctly handle alt-loc codes, they
  are not supported by these modules.

3.  Value-added methods.  Return values are not verbatem strings
    parsed from DSSP or STRIDE output.


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

=head1 AUTHOR - Ed Green

Email ed@compbio.berkeley.edu


=head1 APPENDIX

The rest of the documentation details each method.
Internal methods are preceded with a _

=cut

package Bio::Structure::SecStr::DSSP::Res;
use strict;
use Bio::Root::IO;
use Bio::PrimarySeq;

use base qw(Bio::Root::Root);

# Would be a class variable if Perl had them

               #attribute        begin col        # columns
our %lookUp = ( 'pdb_resnum'     => [  5,           5 ],
		'insertionco'    => [  10,          1 ],
		'pdb_chain'      => [  11,          1 ],
		
		'amino_acid'     => [  13,          1 ],
		'term_sig'       => [  14,          1 ],
		
		'ss_summary'     => [  16,          1 ],
		'3tph'           => [  18,          1 ],
		'4tph'           => [  19,          1 ],
		'5tph'           => [  20,          1 ],
		'geo_bend'       => [  21,          1 ],
		'chirality'      => [  22,          1 ],
		'beta_br1la'     => [  23,          1 ],
		'beta_br2la'     => [  24,          1 ],

		'bb_part1nu'     => [  25,          4 ],
		'bb_part2nu'     => [  29,          4 ],
		'betash_lab'     => [  33,          1 ],
		
		'solv_acces'     => [  34,          4 ],
		
		'hb1_nh_o_p'     => [  39,          6 ],
		'hb1_nh_o_e'     => [  46,          4 ],
		
		'hb1_o_hn_p'     => [  50,          6 ],
		'hb1_o_hn_e'     => [  57,          4 ],
		
		'hb2_nh_o_p'     => [  61,          6 ],
		'hb2_nh_o_e'     => [  68,          4 ],
		
		'hb2_o_hn_p'     => [  72,          6 ],
		'hb2_o_hn_e'     => [  79,          4 ],
		
		'tco'            => [  85,          6 ],
		
		'kappa'          => [  91,          6 ],
		
		'alpha'          => [  97,          6 ],
	
		'phi'            => [ 103,          6 ],

		'psi'            => [ 109,          6 ],
		
		'x_ca'           => [ 115,          7 ],
		
		'y_ca'           => [ 122,          7 ],
		
		'z_ca'           => [ 129,          7 ] );


=head1 CONSTRUCTOR


=cut


=head2 new

 Title         : new
 Usage         : makes new object of this class
 Function      : Constructor
 Example       : $dssp_obj = Bio::DSSP:Res->new( filename or FILEHANDLE )
 Returns       : object (ref)
 Args          : filename ( must be proper DSSP output file )

=cut

sub new {
    my ( $class, @args ) = @_;
    my $self = $class->SUPER::new( @args );
    my $io = Bio::Root::IO->new( @args );
    $self->_parse( $io->_fh() );
    $io->close();
    return $self;
}

=head1 ACCESSORS


=cut

# GLOBAL FEATURES / INFO / STATS

=head2 totSurfArea

 Title         : totSurfArea
 Usage         : returns total accessible surface area in square And.
 Function      :
 Example       : $surArea = $dssp_obj->totSurfArea();
 Returns       : scalar
 Args          : none

=cut

sub totSurfArea {
    my $self = shift;
    return $self->{ 'Head' }->{ 'ProAccSurf' };
}

=head2 numResidues

 Title         : numResidues
 Usage         : returns the total number of residues in all chains or
                 just the specified chain if a chain is specified
 Function      :
 Example       : $num_res = $dssp_obj->numResidues();
 Returns       : scalar int
 Args          : none


=cut

sub numResidues {
    my $self = shift;
    my $chain = shift;
    if ( !( $chain ) ) {
	return $self->{'Head'}->{'TotNumRes'};
    }
    else {
	my ( $num_res,
	     $cont_seg );
	my $cont_seg_pnt = $self->_contSegs();
	foreach $cont_seg ( @{ $cont_seg_pnt } ) {
	    if ( $chain eq $cont_seg->[ 2 ] ) {
		# this segment is part of the chain we want
		$num_res += ( $self->_toDsspKey( $cont_seg->[ 1 ] )
			      - $self->_toDsspKey( $cont_seg->[ 0 ] )
			      + 1 ); # this works because we know the
				     # the region between the start
				     # and end of a dssp key is
				     # continuous
	    }
	}
	return $num_res;
    }
}

#  STRAIGHT FROM PDB ENTRY

=head2 pdbID

 Title         : pdbID
 Usage         : returns pdb identifier ( 1FJM, e.g.)
 Function      :
 Example       : $pdb_id = $dssp_obj->pdbID();
 Returns       : scalar string
 Args          : none


=cut

sub pdbID {
    my $self = shift;
    return $self->{'Head'}->{'PDB'};
}

=head2 pdbAuthor

 Title         : pdbAuthor
 Usage         : returns author field
 Function      :
 Example       : $auth = $dssp_obj->pdbAuthor()
 Returns       : scalar string
 Args          : none


=cut

sub pdbAuthor {
    my $self = shift;
    return $self->{'Head'}->{'AUTHOR'};
}

=head2 pdbCompound

 Title         : pdbCompound
 Usage         : returns pdbCompound given in PDB file
 Function      :
 Example       : $cmpd = $dssp_obj->pdbCompound();
 Returns       : scalar string
 Args          : none


=cut

sub pdbCompound {
    my $self = shift;
    return $self->{'Head'}->{'COMPND'};
}

=head2 pdbDate

 Title         : pdbDate
 Usage         : returns date given in PDB file
 Function      :
 Example       : $pdb_date = $dssp_obj->pdbDate();
 Returns       : scalar
 Args          : none


=cut

sub pdbDate {
    my $self = shift;
    return $self->{'Head'}->{'DATE'};
}

=head2 pdbHeader

 Title         : pdbHeader
 Usage         : returns header info from PDB file
 Function      :
 Example       : $header = $dssp_obj->pdbHeader();
 Returns       : scalar
 Args          : none


=cut

sub pdbHeader {
    my $self = shift;
    return $self->{'Head'}->{'HEADER'};
}

=head2 pdbSource

 Title         : pdbSource
 Usage         : returns pdbSource information from PDBSOURCE line
 Function      :
 Example       : $pdbSource = $dssp_obj->pdbSource();
 Returns       : scalar
 Args          : none


=cut

sub pdbSource {
    my $self = shift;
    return $self->{'Head'}->{'SOURCE'};
}


# RESIDUE SPECIFIC ACCESSORS

=head2 resAA

 Title         : resAA
 Usage         : fetches the 1 char amino acid code, given an id
 Function      :
 Example       : $aa = $dssp_obj->resAA( '20:A' ); # pdb id as arg
 Returns       : 1 character scalar string
 Args          : RESIDUE_ID


=cut

sub resAA {
    my $self = shift;
    my @args = @_;
    my $dssp_key = $self->_toDsspKey( @args );
    return $self->{ 'Res' }->[ $dssp_key ]->{ 'amino_acid' };
}

=head2 resPhi

 Title         : resPhi
 Usage         : returns phi angle of a single residue
 Function      : accessor
 Example       : $phi = $dssp_obj->resPhi( RESIDUE_ID )
 Returns       : scalar
 Args          : RESIDUE_ID


=cut

sub resPhi {
    my $self = shift;
    my @args = @_;
    my $dssp_key = $self->_toDsspKey( @args );
    return $self->{ 'Res' }->[ $dssp_key ]->{ 'phi' };
}

=head2 resPsi

 Title         : resPsi
 Usage         : returns psi angle of a single residue
 Function      : accessor
 Example       : $psi = $dssp_obj->resPsi( RESIDUE_ID )
 Returns       : scalar
 Args          : RESIDUE_ID


=cut

sub resPsi {
    my $self = shift;
    my @args = @_;
    my $dssp_key = $self->_toDsspKey( @args );
    return $self->{ 'Res' }->[ $dssp_key ]->{ 'psi' };
}

=head2 resSolvAcc

 Title         : resSolvAcc
 Usage         : returns solvent exposed area of this residue in
                 square Andstroms
 Function      :
 Example       : $solv_acc = $dssp_obj->resSolvAcc( RESIDUE_ID );
 Returns       : scalar
 Args          : RESIDUE_ID


=cut

sub resSolvAcc {
    my $self = shift;
    my @args = @_;
    my $dssp_key = $self->_toDsspKey( @args );
    return $self->{ 'Res' }->[ $dssp_key ]->{ 'solv_acces' };
}

=head2 resSurfArea

 Title         : resSurfArea
 Usage         : returns solvent exposed area of this residue in
                 square Andstroms
 Function      :
 Example       : $solv_acc = $dssp_obj->resSurfArea( RESIDUE_ID );
 Returns       : scalar
 Args          : RESIDUE_ID


=cut

sub resSurfArea {
    my $self = shift;
    my @args = @_;
    my $dssp_key = $self->_toDsspKey( @args );
    return $self->{ 'Res' }->[ $dssp_key ]->{ 'solv_acces' };
}

=head2 resSecStr

 Title         : resSecStr
 Usage         : $ss = $dssp_obj->resSecStr( RESIDUE_ID );
 Function      : returns the DSSP secondary structural designation of this residue
 Example       :
 Returns       : a character ( 'B', 'E', 'G', 'H', 'I', 'S', 'T', or ' ' )
 Args          : RESIDUE_ID
 NOTE          : The range of this method differs from that of the
    resSecStr method in the STRIDE SecStr parser.  That is because of the
    slightly different format for STRIDE and DSSP output.  The resSecStrSum
    method exists to map these different ranges onto an identical range.

=cut

sub resSecStr {
    my $self = shift;
    my @args = @_;
    my $dssp_key = $self->_toDsspKey( @args );
    my $ss_char = $self->{ 'Res' }->[ $dssp_key ]->{ 'ss_summary' };
    return $ss_char if $ss_char;
    return ' ';
}


=head2 resSecStrSum

 Title         : resSecStrSum
 Usage         : $ss = $dssp_obj->resSecStrSum( $id );
 Function      : returns what secondary structure group this residue belongs
                 to.  One of:  'H': helix ( H, G, or I )
                               'B': beta  ( B or E )
                               'T': turn  ( T or S )
                               ' ': none  ( ' ' )
                 This method is similar to resSecStr, but the information
                 it returns is less specific.
 Example       :
 Returns       : a character ( 'H', 'B', 'T', or ' ' )
 Args          : dssp residue number of pdb residue identifier


=cut

sub resSecStrSum {
    my $self = shift;
    my @args = @_;
    my $dssp_key = $self->_toDsspKey( @args );
    my $ss_char = $self->{ 'Res' }->[ $dssp_key ]->{ 'ss_summary' };
    if ( $ss_char eq 'H' || $ss_char eq 'G' || $ss_char eq 'I' ) {
	return 'H';
    }
    if ( $ss_char eq ' ' || !( $ss_char ) ) {
	return ' ';
    }
    if ( $ss_char eq 'B' || $ss_char eq 'E' ) {
	return 'B';
    }
    else {
	return 'T';
    }
}

# DSSP SPECIFIC

=head2 hBonds

 Title         : hBonds
 Usage         : returns number of 14 different types of H Bonds
 Function      :
 Example       : $hb = $dssp_obj->hBonds
 Returns       : pointer to 14 element array of ints
 Args          : none
 NOTE          : The different type of H-Bonds reported are, in order:
    TYPE O(I)-->H-N(J)
    IN PARALLEL BRIDGES
    IN ANTIPARALLEL BRIDGES
    TYPE O(I)-->H-N(I-5)
    TYPE O(I)-->H-N(I-4)
    TYPE O(I)-->H-N(I-3)
    TYPE O(I)-->H-N(I-2)
    TYPE O(I)-->H-N(I-1)
    TYPE O(I)-->H-N(I+0)
    TYPE O(I)-->H-N(I+1)
    TYPE O(I)-->H-N(I+2)
    TYPE O(I)-->H-N(I+3)
    TYPE O(I)-->H-N(I+4)
    TYPE O(I)-->H-N(I+5)

=cut

sub hBonds {
    my $self = shift;
    return $self->{ 'HBond'};
}

=head2 numSSBr

 Title         : numSSBr
 Usage         : returns info about number of SS-bridges
 Function      :
 Example       : @SS_br = $dssp_obj->numSSbr();
 Returns       : 3 element scalar int array
 Args          : none


=cut

sub numSSBr {
    my $self = shift;
    return ( $self->{'Head'}->{'TotSSBr'},
	     $self->{'Head'}->{'TotIaSSBr'},
	     $self->{'Head'}->{'TotIeSSBr'} );
}

=head2 resHB_O_HN

 Title         : resHB_O_HN
 Usage         : returns pointer to a 4 element array
                 consisting of: relative position of binding
                 partner #1, energy of that bond (kcal/mol),
                 relative positionof binding partner #2,
                 energy of that bond (kcal/mol).  If the bond
                 is not bifurcated, the second bond is reported
                 as 0, 0.0
 Function      : accessor
 Example       : $oBonds_ptr = $dssp_obj->resHB_O_HN( RESIDUE_ID )
 Returns       : pointer to 4 element array
 Args          : RESIDUE_ID


=cut

sub resHB_O_HN {
    my $self = shift;
    my @args = @_;
    my $dssp_key = $self->_toDsspKey( @args );
    return ( $self->{ 'Res' }->[ $dssp_key ]->{ 'hb1_o_hn_p' },
	     $self->{ 'Res' }->[ $dssp_key ]->{ 'hb1_o_hn_e' },
	     $self->{ 'Res' }->[ $dssp_key ]->{ 'hb2_o_hn_p' },
	     $self->{ 'Res' }->[ $dssp_key ]->{ 'hb2_o_hn_e' } );
}


=head2 resHB_NH_O

 Title         : resHB_NH_O
 Usage         : returns pointer to a 4 element array
                 consisting of: relative position of binding
                 partner #1, energy of that bond (kcal/mol),
                 relative positionof binding partner #2,
                 energy of that bond (kcal/mol).  If the bond
                 is not bifurcated, the second bond is reported
                 as 0, 0.0
 Function      : accessor
 Example       : $nhBonds_ptr = $dssp_obj->resHB_NH_O( RESIDUE_ID )
 Returns       : pointer to 4 element array
 Args          : RESIDUE_ID


=cut

sub resHB_NH_O {
    my $self = shift;
    my @args = @_;
    my $dssp_key = $self->_toDsspKey( @args );
    return ( $self->{ 'Res' }->[ $dssp_key ]->{ 'hb1_nh_o_p' },
	     $self->{ 'Res' }->[ $dssp_key ]->{ 'hb1_nh_o_e' },
	     $self->{ 'Res' }->[ $dssp_key ]->{ 'hb2_nh_o_p' },
	     $self->{ 'Res' }->[ $dssp_key ]->{ 'hb2_nh_o_e' } );
}


=head2 resTco

 Title         : resTco
 Usage         : returns tco angle around this residue
 Function      : accessor
 Example       : resTco = $dssp_obj->resTco( RESIDUE_ID )
 Returns       : scalar
 Args          : RESIDUE_ID


=cut

sub resTco {
    my $self = shift;
    my @args = @_;
    my $dssp_key = $self->_toDsspKey( @args );
    return $self->{ 'Res' }->[ $dssp_key ]->{ 'tco' };
}


=head2 resKappa

 Title         : resKappa
 Usage         : returns kappa angle around this residue
 Function      : accessor
 Example       : $kappa = $dssp_obj->resKappa( RESIDUE_ID )
 Returns       : scalar
 Args          : RESIDUE_ID ( dssp or PDB )


=cut

sub resKappa {
    my $self = shift;
    my @args = @_;
    my $dssp_key = $self->_toDsspKey( @args );
    return $self->{ 'Res' }->[ $dssp_key ]->{ 'kappa' };
}


=head2 resAlpha

 Title         : resAlpha
 Usage         : returns alpha angle around this residue
 Function      : accessor
 Example       : $alpha = $dssp_obj->resAlpha( RESIDUE_ID )
 Returns       : scalar
 Args          : RESIDUE_ID ( dssp or PDB )


=cut

sub resAlpha {
    my $self = shift;
    my @args = @_;
    my $dssp_key = $self->_toDsspKey( @args );
    return $self->{ 'Res' }->[ $dssp_key ]->{ 'alpha' };
}

# VALUE ADDED METHODS (NOT JUST PARSE/REPORT)

=head2 secBounds

 Title         : secBounds
 Usage         : gets residue ids of boundary residues in each
                 contiguous secondary structural element of specified
                 chain
 Function      : returns pointer to array of 3 element arrays.  First
                 two elements are the PDB IDs of the start and end points,
                 respectively and inclusively.  The last element is the
                 DSSP secondary structural assignment code,
                 i.e. one of : ('B', 'E', 'G', 'H', 'I', 'S', 'T', or ' ')
 Example       : $ss_elements_pts = $dssp_obj->secBounds( 'A' );
 Returns       : pointer to array of arrays
 Args          : chain id ( 'A', for example ).  No arg => no chain id


=cut

sub secBounds {
    my $self = shift;
    my $chain = shift;
    my %sec_bounds;

    $chain = '-' if ( !( $chain ) || $chain eq ' ' || $chain eq '-' );

    # if we've memoized this chain, use that
    if ( $self->{ 'SecBounds' } ) {
	# check to make sure chain is valid
	if ( !( $self->{ 'SecBounds' }->{ $chain } ) ) {
	    $self->throw( "No such chain: $chain\n" );
	}
	return $self->{ 'SecBounds' }->{ $chain };
    }

    my ( $cur_element, $i, $cur_chain, $beg, );

    #initialize
    $cur_element = $self->{ 'Res' }->[ 1 ]->{ 'ss_summary' };
    $beg = 1;

    for ( $i = 2; $i <= $self->_numResLines() - 1; $i++ ) {
	if ( $self->{ 'Res' }->[ $i ]->{ 'amino_acid' } eq '!' ) {
	    # element is terminated by a chain discontinuity
	    push( @{ $sec_bounds{ $self->_pdbChain( $beg ) } },
		  [ $self->_toPdbId( $beg ),
		    $self->_toPdbId( $i - 1 ),
		    $cur_element ] );
	    $i++;
	    $beg = $i;
	    $cur_element = $self->{ 'Res' }->[ $i ]->{ 'ss_summary' };
	}
	
	elsif ( $self->{ 'Res' }->[ $i ]->{ 'ss_summary' } ne $cur_element ) {
	    # element is terminated by beginning of a new element
	    push( @{ $sec_bounds{ $self->_pdbChain( $beg ) } },
		  [ $self->_toPdbId( $beg ),
		    $self->_toPdbId( $i - 1 ),
		    $cur_element ] );
	    $beg = $i;
	    $cur_element = $self->{ 'Res' }->[ $i ]->{ 'ss_summary' };
	}
    }
    #last residue
    if ( $self->{ 'Res' }->[ $i ]->{ 'ss_summary' } eq $cur_element ) {
	push( @{ $sec_bounds{ $self->_pdbChain( $beg ) } },
	      [ $self->_toPdbId( $beg ),
		$self->_toPdbId( $i ),
		$cur_element ] );
    }

    else {
	push( @{ $sec_bounds{ $self->_pdbChain( $beg ) } },
	      [ $self->_toPdbId( $beg ),
		$self->_toPdbId( $i - 1 ),
		$cur_element ] );
	push( @{ $sec_bounds{ $self->_pdbChain( $i ) } },
	      [ $self->_toPdbId( $i ),
		$self->_toPdbId( $i ),
		$self->{ 'Res' }->[ $i ]->{ 'ss_summary' } ] );
    }

    $self->{ 'SecBounds' } = \%sec_bounds;

    # check to make sure chain is valid
    if ( !( $self->{ 'SecBounds' }->{ $chain } ) ) {
	$self->throw( "No such chain: $chain\n" );
    }

    return $self->{ 'SecBounds' }->{ $chain };
}



=head2 chains

 Title         : chains
 Usage         : returns pointer to array of chain I.D.s (characters)
 Function      :
 Example       : $chains_pnt = $dssp_obj->chains();
 Returns       : array of characters, one of which may be ' '
 Args          : none


=cut

sub chains {
    my $self = shift;
    my $cont_segs = $self->_contSegs();
    my %chains;
    my $seg;
    foreach $seg ( @{ $cont_segs } ) {
	$chains{ $seg->[ 2 ] } = 1;
    }
    my @chains = keys( %chains );
    return \@chains;
}


=head2 residues

    Title : residues
    Usage : returns array of residue identifiers for all residues in
    the output file, or in a specific chain
    Function :
    Example : @residues_ids = $dssp_obj->residues()
    Returns : array of residue identifiers
    Args : if none => returns residue ids of all residues of all
    chains (in order); if chain id is given, returns just the residue
    ids of residues in that chain


=cut

# Can't use the standard interface for getting the amino acid,
# pdb_resnum, etc. in this method because we don't *know* the residue
# indentifiers - we are building a list of them.
sub residues {
    my $self  = shift;
    my $chain = shift;
    my @residues;
    my $num_res = $self->_numResLines();
    my $aa;
    for ( my $i = 1; $i <= $num_res; $i++ ) {
	# find what character was in the slot for tha amino acid code,
	# if it's a '!' we know this is not a *real* amino acid, it's
	# a chain discontinuity marker 
	$aa = $self->{ 'Res' }->[ $i ]->{ 'amino_acid' };
	if ( $aa ne '!' ) {
	    if ( !$chain ||
		 $chain eq $self->{ 'Res' }->[ $i ]->{ 'pdb_chain' } ) {
		push( @residues, 
		      $self->{ 'Res' }->[ $i ]->{ 'pdb_resnum' }.
		      $self->{ 'Res' }->[ $i ]->{ 'insertionco' }.
		      ":".
		      $self->{ 'Res' }->[ $i ]->{ 'pdb_chain' } );
	    }
	}
    }
    return @residues;
}


=head2 getSeq

 Title         : getSeq
 Usage         : returns a Bio::PrimarySeq object which represents a good
                 guess at the sequence of the given chain
 Function      : For most chains of most entries, the sequence returned by
                 this method will be very good.  However, it is inherently
                 unsafe to rely on DSSP to extract sequence information about
                 a PDB entry.  More reliable information can be obtained from
                 the PDB entry itself.
 Example       : $pso = $dssp_obj->getSeq( 'A' );
 Returns       : (pointer to) a PrimarySeq object
 Args          : Chain identifier.  If none given, ' ' is assumed.  If no ' '
                 chain, the first chain is used.


=cut

sub getSeq {
    my $self  = shift;
    my $chain = shift;

    my ( $pot_chain,
	 $seq,
	 $frag_num,
	 $frag,
	 $curPdbNum,
	 $lastPdbNum,
	 $gap_len,
	 $i,
	 $id,
	 );
    my @frags;

    if ( !( $chain ) ) {
	$chain = ' ';
    }

    if ( $self->{ 'Seq' }->{ $chain } ) {
	return $self->{ 'Seq' }->{ $chain };
    }

    my $contSegs_pnt = $self->_contSegs();

    # load up specified chain
    foreach $pot_chain ( @{ $contSegs_pnt } ) {
	if ( $pot_chain->[ 2 ] eq $chain ) {
	    push( @frags, $pot_chain );
	}
    }
    
    # if that didn't work, just get the first one
    if ( !( @frags ) ) {
	$chain = $contSegs_pnt->[ 0 ]->[ 2 ];
	foreach $pot_chain ( @{ $contSegs_pnt } ) {
	    if ( $pot_chain->[ 2 ] eq $chain ) {
		push( @frags, $pot_chain );
	    }
	}
    }

    # now build the sequence string
    $seq = "";
    $frag_num = 0;
    foreach $frag ( @frags ) {
	$frag_num++;
	if ( $frag_num > 1 ) {  # we need to put in some gap seq
	    $curPdbNum = $self->_pdbNum( $frag->[ 0 ] );
	    $gap_len = $curPdbNum - $lastPdbNum - 1;
	    if ( $gap_len > 0 ) {
		$seq .= 'u' x $gap_len;
	    }
	    else {
		$seq .= 'u';
	    }
	}
	for ( $i = $frag->[ 0 ]; $i <= $frag->[ 1 ]; $i++ ) {
	    $seq .= $self->_resAA( $i );
	}
	$lastPdbNum = $self->_pdbNum( $i - 1 );
    }



    $id = $self->pdbID();
    $id .= ":$chain";

    $self->{ 'Seq' }->{ $chain } =  Bio::PrimarySeq->new ( -seq => $seq,
							   -id  => $id,
							   -moltype => 'protein'
							   );
    return $self->{ 'Seq' }->{ $chain };
}

=head1 INTERNAL METHODS


=cut

=head2 _pdbChain

 Title         : _pdbChain
 Usage         : returns the pdb chain id of given residue
 Function      :
 Example       : $chain_id = $dssp_obj->pdbChain( DSSP_KEY );
 Returns       : scalar
 Args          : DSSP_KEY ( dssp or pdb )


=cut

sub _pdbChain {
    my $self = shift;
    my $dssp_key = shift;
    return $self->{ 'Res' }->[ $dssp_key ]->{ 'pdb_chain' };
}

=head2 _resAA

 Title         : _resAA
 Usage         : fetches the 1 char amino acid code, given a dssp id
 Function      :
 Example       : $aa = $dssp_obj->_resAA( dssp_id );
 Returns       : 1 character scalar string
 Args          : dssp_id


=cut

sub _resAA {
    my $self = shift;
    my $dssp_key = shift;
    return $self->{ 'Res' }->[ $dssp_key ]->{ 'amino_acid' };
}


=head2 _pdbNum

 Title        : _pdbNum
 Usage        : fetches the numeric portion of the identifier for a given
                residue as reported by the pdb entry.  Note, this DOES NOT
                uniquely specify a residue.  There may be an insertion code
                and/or chain identifier differences.
 Function     :
 Example      : $pdbNum = $self->_pdbNum( DSSP_ID );
 Returns      : a scalar
 Args         : DSSP_ID


=cut

sub _pdbNum {
    my $self = shift;
    my $dssp_key = shift;
    return $self->{ 'Res' }->[ $dssp_key ]->{ 'pdb_resnum' };
}

=head2 _pdbInsCo

 Title        : _pdbInsCo
 Usage        : fetches the Insertion Code for this residue, if it has one.
 Function     :
 Example      : $pdbNum = $self->_pdbInsCo( DSSP_ID );
 Returns      : a scalar
 Args         : DSSP_ID


=cut

sub _pdbInsCo {
    my $self = shift;
    my $dssp_key = shift;
    return $self->{ 'Res' }->[ $dssp_key ]->{ 'insertionco' };
}

=head2 _toPdbId

 Title        : _toPdbId
 Usage        : Takes a dssp key and builds the corresponding
                PDB identifier string
 Function     :
 Example      : $pdbId = $self->_toPdbId( DSSP_ID );
 Returns      : scalar
 Args         : DSSP_ID

=cut

sub _toPdbId {
    my $self = shift;
    my $dssp_key = shift;
    my $pdbId = ( $self->_pdbNum( $dssp_key ).
		  $self->_pdbInsCo( $dssp_key ) );
    my $chain = $self->_pdbChain( $dssp_key );
    $pdbId = "$pdbId:$chain" if $chain;
    return $pdbId;
}

=head2 _contSegs

 Title         : _contSegs
 Usage         : find the endpoints of continuous regions of this structure
 Function      : returns pointer to array of 3 element array.
                 Elements are the dssp keys of the start and end points of each
                 continuous element and its PDB chain id (may be blank).
                 Note that it is common to have several
                 continuous elements with the same chain id.  This occurs
                 when an internal region is disordered and no structural
                 information is available.
 Example       : $cont_seg_ptr = $dssp_obj->_contSegs();
 Returns       : pointer to array of arrays
 Args          : none


=cut

sub _contSegs {
    my $self = shift;
    if ( $self->{ 'contSegs' } ) {
	return $self->{ 'contSegs' };
    }
    else {
	# first time, so make contSegs
	my ( $cur_chain, $i, $beg );
	my @contSegs;
	#initialize
	$cur_chain = $self->_pdbChain( 1 );
	$beg = 1;
	#internal residues
	for ( $i = 2; $i <= $self->_numResLines() - 1; $i++ ) {
	    if ( $self->{ 'Res' }->[ $i ]->{ 'amino_acid' } eq '!' ) {
		push( @contSegs, [ $beg, $i - 1, $cur_chain ] );
		$beg = $i + 1;
		$cur_chain = $self->_pdbChain( $i + 1 );
	    }
	}
	# last residue must be the end of a chain
	push( @contSegs, [ $beg, $i, $cur_chain ] );

	$self->{ 'contSegs' } = \@contSegs;
	return $self->{ 'contSegs' };
    }
}

=head2 _numResLines

 Title         : _numResLines
 Usage         : returns the total number of residue lines in this
                 dssp file.
                 This number is DIFFERENT than the number of residues in
                 the pdb file because dssp has chain termination and chain
                 discontinuity 'residues'.
 Function      :
 Example       : $num_res = $dssp_obj->_numResLines();
 Returns       : scalar int
 Args          : none


=cut

sub _numResLines {
    my $self = shift;
    return ( $#{$self->{ 'Res' }} );
}

=head2 _toDsspKey

 Title         : _toDsspKey
 Usage         : returns the unique dssp integer key given a pdb residue id.
                 All accessor methods require (internally)
                 the dssp key.   This method is very useful in converting
                 pdb keys to dssp keys so the accessors can accept pdb keys
                 as argument.  PDB Residue IDs are inherently
                 problematic since they have multiple parts of
                 overlapping function and ill-defined or observed
                 convention in form.  Input can be in any of the formats
                 described in the DESCRIPTION section above.
 Function      :
 Example       : $dssp_id = $dssp_obj->_pdbKeyToDsspKey( '10B:A' )
 Returns       : scalar int
 Args          : pdb residue identifier: num[insertion code]:[chain]


=cut

sub _toDsspKey {
    # Consider adding lookup table for 'common' name (like 20:A) for
    # fast access.  Could be built during parse of input.

    my $self = shift;
 
    my ( $key_num, $chain_id, $ins_code ) = @_;

    if ( ! $chain_id) { # parse the lone argument
	( $key_num, $chain_id, $ins_code ) 
	    = $key_num =~ m/([0-9]+)
                            ([a-zA-z]?)
                            (?::([a-zA-Z]))?/xms
              ? ( $1, $2, $3 )
	      : $self->throw("Could not derive PDB key $key_num");
     }
    

    # Now find the residue which fits this description.  Linear search is
    # probably not the best way to do this, but oh well...
    for ( my $i = 1; $i <= $self->_numResLines(); $i++ ) {
	unless ( ($self->{'Res'}->[$i]->{'term_sig'} eq '*') ||
		 ($self->{'Res'}->[$i]->{'amino_acid'} eq '!') ) {
	    # chain break 'residue', doesn't match anything
	    if ( $key_num == $self->{'Res'}->[$i]->{'pdb_resnum'} ) {
		if ( $chain_id ) { # if a chain was specified
		    if ( $chain_id eq $self->{'Res'}->[$i]->{'pdb_chain'} ) {
			# and it's the right one
			if ( $ins_code ) { # if insertion code was specified
			    if ( $ins_code eq $self->{'Res'}->[$i]->{'insertionco'} ) {
				# and it's the right one
				return $i;
			    }
			}
			elsif ( $self->{'Res'}->[$i]->{'insertionco'} eq ''  ) { # no isertion code specified, but need to check that the located residue doesn't have an insertion code E.g. pdb1aye fails on this
			    return $i;
			}
		    }
		}
		else { # no chain was specified
		    return $i;
		}
	    }
	}
    }
    $self->throw( "PDB key not found." );
}

=head2 _parse

 Title         : _parse
 Usage         : parses dssp output
 Function      :
 Example       : used by the constructor
 Returns       :
 Args          : input source ( handled by Bio::Root:IO )


=cut

sub _parse {
    my $self = shift;
    my $file = shift;
    my $cur;
    my $current_chain;
    my ( @elements, @hbond );
    my ( %head, %his, );
    my $element;
    my $res_num;

    $cur = <$file>;
    unless ( $cur =~ /^==== Secondary Structure Definition/ ) {
	$self->throw( "Not dssp output" );
	return;
    }

    # REFERENCE line (always there)
    $cur = <$file>;
    ( $element ) = ( $cur =~ /^REFERENCE\s+(.+?)\s+\./ );
    $head{ 'REFERENCE' } = $element;

    $cur = <$file>;
    # Check for HEADER line (not always there)
    if ( $cur =~ /^HEADER\s/ ) {
	@elements = split( /\s+/, $cur );
	pop( @elements ); # take off that annoying period
	$head{ 'PDB' } = pop( @elements );
	$head{ 'DATE' } = pop( @elements );
	# now, everything else is "header" except for the word
	# HEADER
	shift( @elements );
	$element = shift( @elements );
	while ( @elements ) {
	    $element = $element." ".shift( @elements );
	}
	$head{ 'HEADER' } = $element;
	
	$cur = <$file>;
    }

    # Check for COMPND line (not always there)
    if ( $cur =~ /^COMPND\s/ ) {
	($element) = ( $cur =~ /^COMPND\s+(.+?)\s+\./ );
	$head{ 'COMPND' } = $element;
	
	$cur = <$file>;
    }

    # Check for SOURCE or PDBSOURCE line (not always there)
    if ( $cur =~ /^PDBSOURCE\s/ ) {
	($element) = ( $cur =~ /^PDBSOURCE\s+(.+?)\s+\./ );
	$head{ 'SOURCE' } = $element;
	
	$cur = <$file>;
    }

    elsif ( $cur =~ /^SOURCE\s/ ) {
	($element) = ( $cur =~ /^SOURCE\s+(.+?)\s+\./ );
	$head{ 'SOURCE' } = $element;
	
	$cur = <$file>;
    }

    # Check for AUTHOR line (not always there)
    if ( $cur =~ /^AUTHOR/ ) {
	($element) = ( $cur =~ /^AUTHOR\s+(.+?)\s+/ );
	$head{ 'AUTHOR' } = $element;

	$cur = <$file>;
    }

    # A B C D E TOTAL NUMBER OF RESIDUES, NUMBER ... line
    @elements = split( /\s+/, $cur );
    shift( @elements );
    $head{ 'TotNumRes' } = shift( @elements );
    $head{ 'NumChain' }  = shift( @elements );
    $head{ 'TotSSBr' }   = shift( @elements );
    $head{ 'TotIaSSBr' } = shift( @elements );
    $head{ 'TotIeSSBr' } = shift( @elements );

    $cur = <$file>;
    ( $element ) = ( $cur =~ /\s*(\d+\.\d*)\s+ACCESSIBLE SURFACE OF PROTEIN/ );
    $head{ 'ProAccSurf' } = $element;
    $self->{ 'Head' } = \%head;

    for ( my $i = 1; $i <= 14; $i++ ) {
	$cur = <$file>;
	( $element ) =
	    $cur =~ /\s*(\d+)\s+\d+\.\d+\s+TOTAL NUMBER OF HYDROGEN/;
	push( @hbond, $element );
#	$hbond{ $hBondType } = $element;
    }
    $self->{ 'HBond' } = \@hbond;

    my $histogram_finished = 0;
    while ( !($histogram_finished) && chomp( $cur = <$file> ) ) {
	if ( $cur =~ /RESIDUE AA STRUCTURE/ ) {
	    $histogram_finished = 1;
	}
    }

    while ( $cur = <$file> ) {
	if ( $cur =~ m/^\s*$/ ) {
	    next;
	}
	$res_num = substr( $cur, 0, 5 );
	$res_num =~ s/\s//g;
	$self->{ 'Res' }->[ $res_num ] = &_parseResLine( $cur );
    }
}


=head2 _parseResLine

 Title         : _parseResLine
 Usage         : parses a single residue line
 Function      :
 Example       : used internally
 Returns       :
 Args          : residue line ( string )


=cut

sub _parseResLine() {
    my $cur = shift;
    my ( $feat, $value );
    my %elements;

    foreach $feat ( keys %lookUp ) {
	$value = substr( $cur, $lookUp{ $feat }->[0],
			 $lookUp{ $feat }->[1] );
	$value =~ s/\s//g;
	$elements{$feat} = $value ;
    }

    # if no chain id, make it '-' (like STRIDE...very convenient)
    if ( !( $elements{ 'pdb_chain' } ) || $elements{ 'pdb_chain'} eq ' ' ) {
	$elements{ 'pdb_chain' } = '-';
    }
    return \%elements;
}

1;

