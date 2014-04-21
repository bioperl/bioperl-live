# $id $
#
# bioperl module for Bio::Structure::SecStr::STRIDE::Res.pm
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

Bio::Structure::SecStr::STRIDE::Res - Module for parsing/accessing stride output

=head1 SYNOPSIS

 my $stride_obj = Bio::Structure::SecStr::STRIDE::Res->new( '-file' => 'filename.stride' );

 # or

 my $stride_obj = Bio::Structure::SecStr::STRIDE::Res->new( '-fh' => \*STDOUT );

 # Get secondary structure assignment for PDB residue 20 of chain A
 $sec_str = $stride_obj->resSecStr( '20:A' );

 # same
 $sec_str = $stride_obj->resSecStr( 20, 'A' )

=head1 DESCRIPTION

STRIDE::Res is a module for objectifying STRIDE output.  STRIDE is a
program (similar to DSSP) for assigning secondary structure to
individual residues of a pdb structure file.

    ( Knowledge-Based Protein Secondary Structure Assignment,
    PROTEINS: Structure, Function, and Genetics 23:566-579 (1995) )

STRIDE is available here:
http://webclu.bio.wzw.tum.de/stride/

Methods are then available for extracting all of the infomation
present within the output or convenient subsets of it.

Although they are very similar in function, DSSP and STRIDE differ
somewhat in output format.  Thes differences are reflected in the
return value of some methods of these modules.  For example, both
the STRIDE and DSSP parsers have resSecStr() methods for returning
the secondary structure of a given residue.  However, the range of
return values for DSSP is ( H, B, E, G, I, T, and S ) whereas the
range of values for STRIDE is ( H, G, I, E, B, b, T, and C ).  See
individual methods for details.

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

=head2 MailingLists

UsUser feedback is an integral part of the evolution of this and other
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

The Rest of the documentation details each method.
Internal methods are preceded with a _.


=cut

package Bio::Structure::SecStr::STRIDE::Res;
use strict;
use Bio::Root::IO;
use Bio::PrimarySeq;

use base qw(Bio::Root::Root);

our %ASGTable = ( 'aa'         =>  0,
		  'resNum'     =>  1,
		  'ssAbbr'     =>  2,
		  'ssName'     =>  3,
		  'phi'        =>  4,
		  'psi'        =>  5,
		  'surfArea'   =>  6 );

our %AATable = ( 'ALA' => 'A', 'ARG' => 'R', 'ASN' => 'N',
		 'ASP' => 'D', 'CYS' => 'C', 'GLN' => 'Q',
		 'GLU' => 'E', 'GLY' => 'G', 'HIS' => 'H',
		 'ILE' => 'I', 'LEU' => 'L', 'LYS' => 'K',
		 'MET' => 'M', 'PHE' => 'F', 'PRO' => 'P',
		 'SER' => 'S', 'THR' => 'T', 'TRP' => 'W',
		 'TYR' => 'Y', 'VAL' => 'V' );

=head2 new

 Title         : new
 Usage         : makes new object of this class
 Function      : Constructor
 Example       : $stride_obj = Bio::Structure::SecStr::STRIDE:Res->new( '-file' =>  filename 
						     # or 
						     '-fh'   => FILEHANDLE )
 Returns       : object (ref)
 Args          : filename or filehandle( must be proper STRIDE output )

=cut

sub new {
    my ( $class, @args ) = @_;
    my $self = $class->SUPER::new( @args );
    my $io   = Bio::Root::IO->new( @args );
    $self->_parse( $io ); # not passing filehandle !
    $io->close();
    return $self;
}

# GLOBAL FEATURES / INFO / STATS

=head2 totSurfArea

 Title         : totSurfArea
 Usage         : returns sum of surface areas of all residues of all
                 chains considered.  Result is memoized.
 Function      :
 Example       : $tot_SA = $stride_obj->totSurfArea();
 Returns       : scalar
 Args          : none


=cut

sub totSurfArea {
    my $self = shift;
    my $total = 0;
    my ( $chain, $res );

    if ( $self->{ 'SurfArea' } ) {
	return $self->{ 'SurfArea' };
    }
    else {
	foreach $chain ( keys %{$self->{ 'ASG' }} ) {
	    for ( my $i = 1; $i <= $#{$self->{'ASG'}->{$chain}}; $i++ ) {
		$total += 
		    $self->{'ASG'}->{$chain}->[$i]->[$ASGTable{'surfArea'}];
	    }
	}
    }

    $self->{ 'SurfArea' } = $total;
    return $self->{ 'SurfArea' };
   
}

=head2 numResidues

 Title         : numResidues
 Usage         : returns total number of residues in all chains or
                 just the specified chain
 Function      : 
 Example       : $tot_res = $stride_obj->numResidues();
 Returns       : scalar int
 Args          : none or chain id


=cut

sub numResidues {
    my $self = shift;
    my $chain = shift;
    my $total = 0;
    my $key;
    foreach $key ( keys %{$self->{ 'ASG' }} ) {
	if ( $chain ) {
	    if ( $key eq $chain ) {
		$total += $#{$self->{ 'ASG' }{ $key }};
	    }
	}
	else {
	    $total += $#{$self->{ 'ASG' }{ $key }};
	}
    }
    return $total;
}

# STRAIGHT FROM THE PDB ENTRY

=head2 pdbID

 Title         : pdbID
 Usage         : returns pdb identifier ( 1FJM, e.g. )
 Function      : 
 Example       : $pdb_id = $stride_obj->pdbID();
 Returns       : scalar string
 Args          : none


=cut

sub pdbID {
    my $self = shift;
    return $self->{ 'PDB' };
}
=head2 pdbAuthor

 Title         : pdbAuthor
 Usage         : returns author of this PDB entry
 Function      : 
 Example       : $auth = $stride_obj->pdbAuthor()
 Returns       : scalar string
 Args          : none


=cut

sub pdbAuthor {
    my $self = shift;
    return join( ' ', @{ $self->{ 'HEAD' }->{ 'AUT' } } );
}

=head2 pdbCompound

 Title         : pdbCompound
 Usage         : returns string of what was found on the  
                 CMP lines
 Function      : 
 Example       : $cmp = $stride_obj->pdbCompound();
 Returns       : string
 Args          : none


=cut

sub pdbCompound {
    my $self = shift;
    return join( ' ', @{ $self->{ 'HEAD' }->{ 'CMP' } } );
}

=head2 pdbDate

 Title         : pdbDate
 Usage         : returns date given in PDB file
 Function      :
 Example       : $pdb_date = $stride_obj->pdbDate();
 Returns       : scalar
 Args          : none


=cut

sub pdbDate {
    my $self = shift;
    return $self->{ 'DATE' };
}

=head2 pdbHeader

 Title         : pdbHeader
 Usage         : returns string of characters found on the PDB header line
 Function      :
 Example       : $head = $stride_obj->pdbHeader();
 Returns       : scalar
 Args          : none


=cut

sub pdbHeader {
    my $self = shift;
    return $self->{ 'HEAD' }->{ 'HEADER' };
}

=head2 pdbSource

 Title         : pdbSource
 Usage         : returns string of what was found on SRC lines
 Function      : 
 Example       : $src = $stride_obj->pdbSource();
 Returns       : scalar
 Args          : none


=cut

sub pdbSource {
    my $self = shift;
    return join( ' ', @{ $self->{ 'HEAD' }->{ 'SRC' } } );
}

# RESIDUE SPECIFIC ACCESSORS

=head2 resAA

 Title         : resAA
 Usage         : returns 1 letter abbr. of the amino acid specified by
                 the arguments
 Function      : 
 Examples      : $aa = $stride_obj->resAA( RESIDUE_ID );
 Returns       : scalar character
 Args          : RESIDUE_ID


=cut

sub resAA {
    my $self = shift;
    my @args = @_;
    my ( $ord, $chain ) = $self->_toOrdChain( @args );
    return ( $AATable{$self->{'ASG'}->{$chain}->[$ord]->[$ASGTable{'aa'}]} );
}

=head2 resPhi

 Title         : resPhi
 Usage         : returns phi angle of specified residue
 Function      :
 Example       : $phi = $stride_obj->resPhi( RESIDUE_ID );
 Returns       : scaler
 Args          : RESIDUE_ID


=cut

sub resPhi {
    my $self = shift;
    my @args = @_;
    my ( $ord, $chain ) = $self->_toOrdChain( @args );
    return $self->{ 'ASG' }->{ $chain }->[ $ord ]->[ $ASGTable{ 'phi' } ];
}

=head2 resPsi

 Title         : resPsi
 Usage         : returns psi angle of specified residue
 Function      :
 Example       : $psi = $stride_obj->resPsi( RESIDUE_ID );
 Returns       : scalar
 Args          : RESIDUE_ID


=cut

sub resPsi {
    my $self = shift;
    my @args = @_;
    my ( $ord, $chain ) = $self->_toOrdChain( @args );
    return $self->{ 'ASG' }->{ $chain }->[ $ord ]->[ $ASGTable{ 'psi' } ];
}

=head2 resSolvAcc

 Title         : resSolvAcc
 Usage         : returns stride calculated surface area of specified residue
 Function      : 
 Example       : $sa = $stride_obj->resSolvAcc( RESIDUE_ID );
 Returns       : scalar
 Args          : RESIDUE_ID


=cut

sub resSolvAcc {
    my $self = shift;
    my @args = @_;
    my ( $ord, $chain ) = $self->_toOrdChain( @args );
    return $self->{ 'ASG' }->{ $chain }->[ $ord ]->[ $ASGTable{ 'surfArea' } ];
}

=head2 resSurfArea

 Title         : resSurfArea
 Usage         : returns stride calculated surface area of specified residue
 Function      : 
 Example       : $sa = $stride_obj->resSurfArea( RESIDUE_ID );
 Returns       : scalar
 Args          : RESIDUE_ID


=cut

sub resSurfArea {
    my $self = shift;
    my @args = @_;
    my ( $ord, $chain ) = $self->_toOrdChain( @args );
    return $self->{ 'ASG' }->{ $chain }->[ $ord ]->[ $ASGTable{ 'surfArea' } ];
}

=head2 resSecStr

 Title         : resSecStr 
 Usage         : gives one letter abbr. of stride determined secondary
                 structure of specified residue
 Function      : 
 Example       : $ss = $stride_obj->resSecStr( RESIDUE_ID );
 Returns       : one of: 'H' => Alpha Helix
                         'G' => 3-10 helix
                         'I' => PI-helix
                         'E' => Extended conformation
                         'B' or 'b' => Isolated bridge
                         'T' => Turn
                         'C' => Coil
                         ' ' => None
                # NOTE:  This range is slightly DIFFERENT from the
                #        DSSP method of the same name
 Args          : RESIDUE_ID


=cut

sub resSecStr {
    my $self = shift;
    my @args = @_;
    my ( $ord, $chain ) = $self->_toOrdChain( @args );
    return $self->{ 'ASG' }->{ $chain }->[ $ord ]->[ $ASGTable{ 'ssAbbr' } ];
}

=head2 resSecStrSum

 Title         : resSecStrSum
 Usage         : gives one letter summary of secondary structure of
                 specified residue.  More general than secStruc() 
 Function      :
 Example       : $ss_sum = $stride_obj->resSecStrSum( RESIDUE_ID );
 Returns       : one of: 'H' (helix), 'B' (beta), 'T' (turn), or 'C' (coil)
 Args          : residue identifier(s) ( SEE INTRO NOTE )


=cut

sub resSecStrSum {
    my $self = shift;
    my @args = @_;
    my $ss_char = $self->resSecStr( @args );

    if ( $ss_char eq 'H' || $ss_char eq 'G' || $ss_char eq 'I' ) {
	return 'H';
    }
    if ( $ss_char eq 'E' || $ss_char eq 'B' || $ss_char eq 'b' ) {
	return 'B';
    }
    if ( $ss_char eq 'T' ) {
	return 'T';
    }
    else {
	return 'C';
    }
}

# STRIDE SPECIFIC

=head2 resSecStrName

 Title         : resSecStrName
 Usage         : gives full name of the secondary structural element
                 classification of the specified residue
 Function      : 
 Example       : $ss_name = $stride_obj->resSecStrName( RESIDUE_ID );
 Returns       : scalar string
 Args          : RESIDUE_ID


=cut

sub resSecStrName {
    my $self = shift;
    my @args = @_;
    my ( $ord, $chain ) = $self->_toOrdChain( @args );
    return $self->{ 'ASG' }->{ $chain }->[ $ord ]->[ $ASGTable{ 'ssName' } ];
}

=head2 strideLocs

 Title         : strideLocs
 Usage         : returns stride determined contiguous secondary
    structural elements as specified on the LOC lines
 Function      : 
 Example       : $loc_pnt = $stride_obj->strideLocs();
 Returns       : pointer to array of 5 element arrays.
    0 => stride name of structural element
    1 => first residue pdb key (including insertion code, if app.)
    2 => first residue chain id
    3 => last residue pdb key (including insertion code, if app.)
    4 => last residue chain id
    NOTE the differences between this range and the range of SecBounds()
 Args          : none


=cut

sub strideLocs {
    my $self = shift;
    return $self->{ 'LOC' };
}

# VALUE ADDED METHODS (NOT JUST PARSE/REPORT)

=head2 secBounds

 Title         : secBounds
 Usage         : gets residue ids of boundary residues in each
                 contiguous secondary structural element of specified
                 chain 
 Function      : 
 Example       : $ss_bound_pnt = $stride_obj->secBounds( 'A' );
 Returns       : pointer to array of 3 element arrays.  First two elements
                 are the PDB IDs of the start and end points, respectively
                 and inclusively.  The last element is the STRIDE secondary
                 structural element code (same range as resSecStr).
 Args          : chain identifier ( one character ).  If none, '-' is assumed


=cut

sub secBounds {
    # Requires a chain name.  If left blank, we assume ' ' which equals '-'
    my $self  = shift;
    my $chain = shift;
    my @SecBounds;

    $chain = '-' if ( !( $chain ) || $chain eq ' ' || $chain eq '-' );

    # if we've memoized this one, use that
    if ( $self->{ 'SecBounds' }->{ $chain } ) {
	return $self->{ 'SecBounds' }->{ $chain }; 
    }

    #check to make sure chain is valid
    if ( !( $self->{ 'ASG' }->{ $chain } ) ) {
	$self->throw( "No such chain: $chain\n" );
    }

    my $cur_element = $self->{ 'ASG' }->{ $chain }->[ 1 ]->
	[ $ASGTable{ 'ssAbbr' } ];
    my $beg = 1;
    my $i;

    for ( $i = 2; $i <= $#{$self->{'ASG'}->{$chain}}; $i++ ) {
	if ( $self->{ 'ASG' }->{ $chain }->[ $i ]->[ $ASGTable{ 'ssAbbr' } ] 
	     ne $cur_element ) {
	    push( @SecBounds, [ $beg, $i -1 , $cur_element ] );
	    $beg = $i;
	    $cur_element = $self->{ 'ASG' }->{ $chain }->[ $i ]->
		[ $ASGTable{ 'ssAbbr' } ];
	}
    }
    
    if ( $self->{ 'ASG' }->{ $chain }->[ $i ]->[ $ASGTable{ 'ssAbbr' } ] 
	 eq $cur_element ) {
	push( @SecBounds, [ $beg, $i, $cur_element ] );
    }
    else {
	push( @SecBounds, [ $beg, $i - 1, $cur_element ], 
	      [ $i, $i, $self->{ 'ASG' }->{ $chain }->[ $i ]->
		[ $ASGTable{ 'ssAbbr' } ] ] );
    }
    
    $self->{ 'SecBounds' }->{ $chain } = \@SecBounds;
    return $self->{ 'SecBounds' }->{ $chain };
}

=head2 chains

 Title         : chains
 Usage         : gives array chain I.D.s (characters)
 Function      :
 Example       : @chains = $stride_obj->chains();
 Returns       : array of characters
 Args          : none


=cut

sub chains {
    my $self = shift;
    my @chains = keys ( %{ $self->{ 'ASG' } } );
    return \@chains;
}

=head2 getSeq

 Title         : getSeq
 Usage         : returns a Bio::PrimarySeq object which represents an
                 approximation at the sequence of the specified chain.
 Function      : For most chain of most entries, the sequence returned by
                 this method will be very good.  However, it it inherently 
                 unsafe to rely on STRIDE to extract sequence information about
                 a PDB entry.  More reliable information can be obtained from
                 the PDB entry itself.  If a second option is given
                 (and evaluates to true), the sequence generated will
                 have 'X' in spaces where the pdb residue numbers are
                 discontinuous.  In some cases this results in a
                 better sequence object (when the  discontinuity is
		 due to regions which were present, but could not be
		 resolved).  In other cases, it will result in a WORSE
                 sequence object (when the discontinuity is due to
		 historical sequence numbering and all sequence is
		 actually resolved).
 Example       : $pso = $dssp_obj->getSeq( 'A' );
 Returns       : (pointer to) a PrimarySeq object
 Args          : Chain identifier.  If none given, '-' is assumed.  


=cut

sub getSeq {
    my $self    = shift;
    my $chain   = shift;
    my $fill_in = shift;

    if ( !( $chain ) ) {
	$chain = '-';
    }

    if ( $self->{ 'Seq' }->{ $chain } ) {
	return $self->{ 'Seq' }->{ $chain };
    }

    my ( $seq, 
	 $num_res,
	 $last_res_num,
	 $cur_res_num, 
	 $i,
	 $step,
	 $id
	 );

    $seq = "";
    $num_res = $self->numResidues( $chain );
    $last_res_num = $self->_pdbNum( 1, $chain );
    for ( $i = 1; $i <= $num_res; $i++ ) {
	if ( $fill_in ) {
	    $cur_res_num = $self->_pdbNum( $i, $chain );
	    $step = $cur_res_num - $last_res_num;
	    if ( $step > 1 ) {
		$seq .= 'X' x ( $step - 1 );
	    }
	}
	$seq .= $self->_resAA( $i, $chain );
	$last_res_num = $cur_res_num;
    }

    $id = $self->pdbID();
    $id .= "$chain";

    $self->{ 'Seq' }->{ $chain } = Bio::PrimarySeq->new( -seq => $seq,
							 -id  => $id,
							 -moltype => 'protein'
							 );

    return $self->{ 'Seq' }->{ $chain };
}

=head1 INTERNAL METHODS

=head2 _pdbNum

 Title        : _pdbNum
 Usage        : fetches the numeric portion of the identifier for a given
                residue as reported by the pdb entry.  Note, this DOES NOT
                uniquely specify a residue.  There may be an insertion code
                and/or chain identifier differences.
 Function     : 
 Example      : $pdbNum = $self->pdbNum( 3, 'A' );
 Returns      : a scalar
 Args         : valid ordinal num / chain combination


=cut

sub _pdbNum {
    my $self  = shift;
    my $ord   = shift;
    my $chain = shift;
    if ( !( $self->{ 'ASG' }->{ $chain }->[ $ord ] ) ) {
	$self->throw( "No such ordinal $ord in chain $chain.\n" );
    }
    my $pdb_junk = $self->{ 'ASG' }->{ $chain }->[ $ord ]->[ $ASGTable{ 'resNum' } ];
    my $num_part;
    ( $num_part ) = ( $pdb_junk =~ /(-*\d+).*/ );
    return $num_part;
}

=head2 _resAA

 Title         : _resAA
 Usage         : returns 1 letter abbr. of the amino acid specified by
                 the arguments
 Function      : 
 Examples      : $aa = $stride_obj->_resAA( 3, '-' );
 Returns       : scalar character
 Args          : ( ord. num, chain )


=cut

sub _resAA {
    my $self  = shift;
    my $ord   = shift;
    my $chain = shift;
    if ( !( $self->{ 'ASG' }->{ $chain }->[ $ord ] ) ) {
	$self->throw( "No such ordinal $ord in chain $chain.\n" );
    }
    return ( $AATable{$self->{'ASG'}->{$chain}->[$ord]->[$ASGTable{'aa'}]} );
}

=head2 _pdbInsCo

 Title        : _pdbInsCo
 Usage        : fetches the Insertion code for this residue.
 Function     : 
 Example      : $pdb_ins_co = $self->_pdb_ins_co( 15, 'B' );
 Returns      : a scalar
 Args         : ordinal number and chain


=cut

sub _pdbInsCo {
    my $self  = shift;
    my $ord   = shift;
    my $chain = shift;
    if ( !( $self->{ 'ASG' }->{ $chain }->[ $ord ] ) ) {
	$self->throw( "No such ordinal $ord in chain $chain.\n" );
    }
    my $pdb_junk = $self->{ 'ASG' }->{ $chain }->[ $ord ]->[ $ASGTable{ 'resNum' } ];
    my $letter_part;
    ( $letter_part ) = ( $pdb_junk =~ /\d+(\D+)/ ); # insertion code can be any 
                                                    # non-word character(s)
    return $letter_part;
}

=head2 _toOrdChain

 Title         : _toOrdChain
 Usage         : takes any set of residue identifying parameters and
    wrestles them into a two element array:  the chain and the ordinal
    number of this residue.  This two element array can then be
    efficiently used as keys in many of the above accessor methods
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

  #  ('#:C) or ( #, 'C' )
       | |
       | -Chain ID
       ---Numeric portion of residue ID.

  If a residue is incompletely specified then the first residue that 
  fits the arguments is returned.  For example, if 19 is the argument 
  and there are three chains, A, B, and C with a residue whose number 
  is 19, then 19:A will be returned (assuming its listed first).

 Function      :
 Example       : my ( $ord, $chain ) = $self->_toOrdChain( @args );
 Returns       : two element array
 Args          : valid set of residue identifier(s) ( SEE NOTE ABOVE )


=cut

sub _toOrdChain {
    my $self = shift;
    my $arg_str;

    my ( $key_num, $chain_id, $ins_code, $key, $i );
    
    # check to see how many args are given
    if ( $#_ >= 1 ) { # multiple args
	$key_num = shift;
	if ( $#_ >= 1 ) { # still multiple args => ins. code, too
	    $ins_code = shift;
	    $chain_id = shift;
	}
	else { # just one more arg. => chain_id
	    $chain_id = shift;
	}
    }
    else { # only single arg.  Might be number or string
	$arg_str = shift;
	if ( $arg_str =~ /:/ ) {
	    # a chain is specified
	    ( $chain_id ) = ( $arg_str =~ /:(.)/);
	    $arg_str =~ s/:.//;
	}
	if ( $arg_str =~ /[A-Z]|[a-z]/ ) {
	    # an insertion code is specified
	    ( $ins_code ) = ( $arg_str =~ /([A-Z]|[a-z])/ );
	    $arg_str =~ s/[A-Z]|[a-z]//g;
	}
	#now, get the number bit-> everything still around
	$key_num = $arg_str;
    }     

    $key = "$key_num$ins_code";
    if ( !( $chain_id ) || $chain_id eq ' ' ) {
	$chain_id = '-';
    }    

    if ( !( $self->{ 'ASG' }->{ $chain_id } ) ) {
	$self->throw( "No such chain: $chain_id" );
    }
    
    for ( $i = 1; $i <= $#{$self->{ 'ASG' }->{ $chain_id }}; $i++ ) {
	if ( $self->{ 'ASG' }->{ $chain_id }->[ $i ]->[ $ASGTable{ 'resNum' } ] eq
	     $key ) {
	    return ( $i, $chain_id );
	}
    }
    
    $self->throw( "No such key: $key" );

}

=head2 _parse

 Title         : _parse
 Usage         : as name suggests, parses stride output, creating object
 Function      :
 Example       : $self->_parse( $io );
 Returns       : 
 Args          : valid Bio::Root::IO object


=cut

sub _parse {
    my $self = shift;
    my $io = shift;
    my $file = $io->_fh();

    # Parse top lines
	if ( $self->_parseTop( $io ) ) {
	$self->throw( "Not stride output" );
    }

    # Parse the HDR, CMP, SCR, and AUT lines
    $self->_parseHead( $io );

    # Parse the CHN, SEQ, STR, and LOC  lines
    $self->_parseSummary( $io ); # we're ignoring this

    # Parse the ASG lines
    $self->_parseASG( $io );
}

=head2 _parseTop

 Title         : _parseTop
 Usage         : makes sure this looks like stride output
 Function      :
 Example       : 
 Returns       :
 Args          :


=cut

sub _parseTop {
    my $self = shift;
    my $io = shift;
    my $file = $io->_fh();
    my $cur = <$file>;
    if ( $cur =~ /^REM  ---/ ) {
	return 0;
    }
    return 1;
}

=head2 _parseHead

 Title         : _parseHead
 Usage         : parses
 Function      : HDR, CMP, SRC, and AUT lines
 Example       :
 Returns       :
 Args          :


=cut

sub _parseHead {
    my $self = shift;
    my $io = shift;
    my $file = $io->_fh();
    my $cur;
    my $element;
    my ( @elements, @cmp, @src, @aut );
    my %head = {};
    my $still_head = 1;

    $cur = <$file>;
    while ( $cur =~ /^REM / ) {
	$cur = <$file>;
    }

    if ( $cur =~ /^HDR / ) {
	@elements = split( /\s+/, $cur );
	shift( @elements );
	pop( @elements );
	$self->{ 'PDB' }  = pop( @elements );
	$self->{ 'DATE' } = pop( @elements );
	# now, everything else is "header" except for the word
	# HDR
	$element = join( ' ', @elements );
	$head{ 'HEADER' } = $element;
    }

    $cur = <$file>;
    while ( $cur =~ /^CMP / ) {
	( $cur ) = ( $cur =~ /^CMP\s+(.+?)\s*\w{4}$/ );
	push( @cmp, $cur );
	$cur = <$file>;
    }

    while ( $cur =~ /^SRC / ) {
	( $cur ) = ( $cur =~ /^SRC\s+(.+?)\s*\w{4}$/ );
	push( @src, $cur );
	$cur = <$file>;
    }

    while ( $cur =~ /^AUT / ) {
	( $cur ) = ( $cur =~ /^AUT\s+(.+?)\s*\w{4}$/ );
	push( @aut, $cur );
	$cur = <$file>;
    }

    $head{ 'CMP' } = \@cmp;
    $head{ 'SRC' } = \@src;
    $head{ 'AUT' } = \@aut;
    $self->{ 'HEAD' } = \%head;
}

=head2 _parseSummary

 Title         : _parseSummary
 Usage         : parses LOC lines
 Function      :
 Example       :
 Returns       :
 Args          :


=cut

sub _parseSummary {
    my $self = shift;
    my $io = shift;
    my $file = $io->_fh();
    my $cur = <$file>;
    my $bound_set;
    my $element;
    my ( @elements, @cur );
    my @LOC_lookup = ( [ 5,  12 ],   # Element name
	# reduntdant	       [ 18, 3 ],    # First residue name
		       [ 22, 5 ],    # First residue PDB number
		       [ 28, 1 ],    # First residue Chain ID
	# redundant	       [ 35, 3 ],    # Last residue name
		       [ 40, 5 ],    # Last residue PDB number
		       [ 46, 1 ]  ); # Last residue Chain ID

    #ignore these lines
    while ( $cur =~ /^REM |^STR |^SEQ |^CHN / ) {
	$cur = <$file>;
    }

    while ( $cur =~ /^LOC / ) {
	foreach $bound_set ( @LOC_lookup ) {
	    $element = substr( $cur, $bound_set->[ 0 ], $bound_set->[ 1 ] );
	    $element =~ s/\s//g;
	    push( @cur, $element );
	}
	push( @elements, [ @cur ] );
	$cur = <$file>;
	@cur = ();
    }
    $self->{ 'LOC' } = \@elements;

}

=head2 _parseASG

 Title         : _parseASG
 Usage         : parses ASG lines
 Function      :
 Example       :
 Returns       :
 Args          :


=cut

sub _parseASG {
    my $self = shift;
    my $io = shift;
    my $file = $io->_fh();
    my $cur = <$file>;
    my $bound_set;
    my $ord_num;
    my ( $chain, $last_chain );
    my $element;
    my %ASG;
    my ( @cur, @elements );
    my @ASG_lookup = ( [ 5,  3 ],  # Residue name
		  #    [ 9,  1 ],  # Chain ID
		       [ 10, 5 ],  # PDB residue number (w/ins.code)
		  #    [ 16, 4 ],  # ordinal stride number
		       [ 24, 1 ],  # one letter sec. stru. abbr.
		       [ 26, 13],  # full sec. stru. name
		       [ 42, 7 ],  # phi angle
		       [ 52, 7 ],  # psi angle
		       [ 64, 5 ] );# residue solv. acc.

    while ( $cur =~ /^REM / ) {
	$cur = <$file>;
    }

    while ( $cur =~ /^ASG / ) {
	# get ordinal number for array key
	$ord_num = substr( $cur, 16, 4 );
	$ord_num =~ s/\s//g;

	# get the chain id
	$chain = substr( $cur, 9, 1 );
	
	if ( $last_chain && ( $chain ne $last_chain ) ) {
	    $ASG{ $last_chain } = [ @elements ];
	    @elements = ();
	}

	# now get the rest of the info on this line
	foreach $bound_set ( @ASG_lookup ) {
	    $element = substr( $cur, $bound_set->[ 0 ], 
			       $bound_set->[ 1 ] );
	    $element =~ s/\s//g;
	    push( @cur, $element );
	}
	$elements[ $ord_num ] = [ @cur ];
	$cur = <$file>;
	@cur = ();
	$last_chain = $chain;
    }

    $ASG{ $chain } = [ @elements ];

    $self->{ 'ASG' } = \%ASG;
}

1;



