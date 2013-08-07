#!perl

=head1 NAME

bp_search2gff

=head1 SYNOPSIS

Usage:

  bp_search2gff [-o outputfile] [-f reportformat] [-i inputfilename]  OR file1 file2 ..

=head1 DESCRIPTION

This script will turn a SearchIO report (BLAST, FASTP, SSEARCH, 
AXT, WABA) into GFF.

The options are:

   -i infilename      - (optional) inputfilename, will read
                        either ARGV files or from STDIN
   -o filename        - the output filename [default STDOUT]
   -f format          - search result format (blast, fasta,waba,axt)
                        (ssearch is fasta format). default is blast.
   -t/--type seqtype  - if you want to see query or hit information
                        in the GFF report
   -s/--source        - specify the source (will be algorithm name
                        otherwise like BLASTN)
   --method           - the method tag (primary_tag) of the features
                        (default is similarity)
   --scorefunc        - a string or a file that when parsed evaluates
                        to a closure which will be passed a feature
                        object and that returns the score to be printed
   --locfunc          - a string or a file that when parsed evaluates
                        to a closure which will be passed two
                        features, query and hit, and returns the
                        location (Bio::LocationI compliant) for the
                        GFF3 feature created for each HSP; the closure
                        may use the clone_loc() and create_loc()
                        functions for convenience, see their PODs
   --onehsp           - only print the first HSP feature for each hit
   -p/--parent        - the parent to which HSP features should refer
                        if not the name of the hit or query (depending
                        on --type)
   --target/--notarget - whether to always add the Target tag or not
   -h                 - this help menu
   --version          - GFF version to use (put a 3 here to use gff 3)
   --component        - generate GFF component fields (chromosome)
   -m/--match         - generate a 'match' line which is a container
                        of all the similarity HSPs
   --addid            - add ID tag in the absence of --match
   -c/--cutoff        - specify an evalue cutoff

Additionally specify the filenames you want to process on the
command-line.  If no files are specified then STDIN input is assumed.
You specify this by doing: bp_search2gff E<lt> file1 file2 file3

=head1 AUTHOR

Jason Stajich, jason-at-bioperl-dot-org

=head1 Contributors

Hilmar Lapp, hlapp-at-gmx-dot-net

=cut

use strict;
use warnings;
use Bio::Tools::GFF;
use Getopt::Long;
use Bio::SearchIO;
use Bio::Location::Simple;    # pre-declare to simplify $locfunc implementations
use Bio::Location::Atomic;    # pre-declare to simplify $locfunc implementations
use Storable qw(dclone);      # for cloning location objects
use Bio::Factory::FTLocationFactory;

my (
    $output,       # output file (if not stdout)
    $input,        # name of the input file
    $format,       # format of the input file, defauly is blast
    $type,         # 'query' or 'hit'
    $cutoff,       # cut-off value for e-value filter
    $sourcetag,    # explicit source tag (will be taken from program
                   # otherwise
    $methodtag,    # primary tag (a.k.a. method), default 'similarity'
    $gffver,       # GFF version (dialect) to write
    $scorefunc,    # closure returning the score for a passed feature
    $locfunc,      # closure returning a location object for a passed
                   # query and hit feature
    $addid,        # flag: whether to always add the ID for $match == 0
    $parent,       # the name of the parent to use; if set and $match == 0
                   # will always add the target
    $comp,         # flag: whether to print a component feature
    $addtarget,    # flag: whether to always add the Target tag, default
                   # is true
    $match,        # flag: whether to print match lines as containers
    $onehsp,       # flag: whether to consider only the first HSP for a hit
    $quiet,        # flag: run quietly
    $help          # flag: show help screen
);                 

# set defaults:
$format    = 'blast';
$type      = 'query';
$gffver    = 2;
$methodtag = "similarity";
$addtarget = 1;

GetOptions(
    'i|input:s'   => \$input,
    'component'   => \$comp,
    'm|match'     => \$match,
    'o|output:s'  => \$output,
    'f|format:s'  => \$format,
    's|source:s'  => \$sourcetag,
    'method=s'    => \$methodtag,
    'addid'       => \$addid,
    'scorefunc=s' => \$scorefunc,
    'locfunc=s'   => \$locfunc,
    'p|parent=s'  => \$parent,
    'target!'     => \$addtarget,
    'onehsp'      => \$onehsp,
    't|type:s'    => \$type,
    'c|cutoff:s'  => \$cutoff,
    'v|version:i' => \$gffver,
    'q|quiet'     => \$quiet,
    'h|help'      => sub {
        exec( 'perldoc', $0 );
        exit(0);
    },
);
$type = lc($type);
if ( $type =~ /target/ ) { $type = 'hit' }
elsif ( $type ne 'query' && $type ne 'hit' ) {
    die("seqtype must be either 'query' or 'hit'");
}

# custom or default function returning the score
$scorefunc =
  defined($scorefunc) ? parse_code($scorefunc) : sub { shift->score };

# custom or default function returning the location
$locfunc = defined($locfunc) ? parse_code($locfunc) : sub { shift->location };

# if --match is given then $addid needs to be disabled
$addid = undef if $addid && $match;

# if no input is provided STDIN will be used
my $parser = new Bio::SearchIO(
    -format  => $format,
    -verbose => $quiet ? -1 : 0,
    -file    => $input
);

my $out;
if ( defined $output ) {
    $out = new Bio::Tools::GFF(
        -gff_version => $gffver,
        -file        => ">$output"
    );
}
else {
    $out = new Bio::Tools::GFF( -gff_version => $gffver );    # STDOUT
}
my ( %seen_hit, %seen );
my $other = $type eq 'query' ? 'hit' : 'query';

while ( my $result = $parser->next_result ) {
    my $qname = $result->query_name;
    if (   $comp
        && $type eq 'query'
        && $result->query_length )
    {
        $out->write_feature(
            Bio::SeqFeature::Generic->new(
                -start       => 1,
                -end         => $result->query_length,
                -seq_id      => $qname,
                -source_tag  => 'chromosome',
                -primary_tag => 'Component',
                -tag         => {
                    'Sequence' => $qname
                }
            )
        );
    }
    while ( my $hit = $result->next_hit ) {
        next if ( defined $cutoff && $hit->significance > $cutoff );
        my $acc = $qname;
        if ( $seen{ $qname . "-" . $hit->name }++ ) {
            $acc = $qname . "-" . $seen{ $qname . '-' . $hit->name };
        }

        if (   $comp
            && $type eq 'hit'
            && $hit->length
            && !$seen_hit{ $hit->name }++ )
        {
            $out->write_feature(
                Bio::SeqFeature::Generic->new(
                    -start       => 1,
                    -end         => $hit->length,
                    -seq_id      => $hit->name,
                    -source_tag  => 'chromosome',
                    -primary_tag => 'Component',
                    -tag         => {
                        'Sequence' => $hit->name
                    }
                )
            );
        }
        my ( %min, %max, $seqid, $name, $st );
        while ( my $hsp = $hit->next_hsp ) {
            my $feature = new Bio::SeqFeature::Generic;
            my ( $proxyfor, $otherf );
            if ( $type eq 'query' ) {
                ( $proxyfor, $otherf ) = ( $hsp->query, $hsp->hit );
                $name ||= $hit->name;
            }
            else {
                ( $otherf, $proxyfor ) = ( $hsp->query, $hsp->hit );
                $name ||= $acc;
            }
            $proxyfor->score( $hit->bits ) unless ( $proxyfor->score );
            if ( ( $gffver == 3 ) && ( $match || $parent ) ) {
                $feature->add_tag_value( 'Parent', $parent || $name );
            }

            $min{$type} = $proxyfor->start
              unless defined $min{$type} && $min{$type} < $proxyfor->start;
            $max{$type} = $proxyfor->end
              unless defined $max{$type} && $max{$type} > $proxyfor->end;
            $min{$other} = $otherf->start
              unless defined $min{$other} && $min{$other} < $otherf->start;
            $max{$other} = $otherf->end
              unless defined $max{$other} && $max{$other} > $otherf->end;
            if ( $addtarget || $match ) {
                $feature->add_tag_value( 'Target', 'Sequence:' . $name );
                $feature->add_tag_value( 'Target', $otherf->start );
                $feature->add_tag_value( 'Target', $otherf->end );
            }
            if ($addid) {
                $feature->add_tag_value( 'ID', $name );
            }

            $feature->location( &$locfunc( $proxyfor, $otherf ) );

            #  strand for feature is always going to be product of
            #  query & hit strands so that target can always be just
            #  '+'
            $feature->strand( $proxyfor->strand * $otherf->strand );
            if ($sourcetag) {
                $feature->source_tag($sourcetag);
            }
            else {
                $feature->source_tag( $proxyfor->source_tag );
            }
            $feature->score( &$scorefunc($proxyfor) );
            $feature->frame( $proxyfor->frame );
            $feature->seq_id( $proxyfor->seq_id );
            $feature->primary_tag($methodtag);

            # add annotation if encoded in the query description
            my $desc = $result->query_description;
            while ( $desc =~ /\/([^=]+)=(\S+)/g ) {
                $feature->add_tag_value( $1, $2 );
            }
            $seqid ||= $proxyfor->seq_id;
            $out->write_feature($feature);
            $st ||= $sourcetag || $proxyfor->source_tag;
            last if $onehsp;
        }
        if ($match) {

            my $matchf = Bio::SeqFeature::Generic->new(
                -start       => $min{$type},
                -end         => $max{$type},
                -strand      => $hit->strand($type) * $hit->strand($other),
                -primary_tag => 'match',
                -source_tag  => $st,
                -score       => $hit->bits,
                -seq_id      => $seqid
            );
            if ( $gffver == 3 ) {
                $matchf->add_tag_value( 'ID', $name );
            }
            $matchf->add_tag_value( 'Target', "Sequence:$name" );
            $out->write_feature($matchf);
        }
    }
}

sub parse_code {
    my $src = shift;
    my $code;

    # file or subroutine?
    if ( -r $src ) {
        if ( !( ( $code = do $src ) && ( ref($code) eq "CODE" ) ) ) {
            die "error in parsing code block $src: $@" if $@;
            die "unable to read file $src: $!"         if $!;
            die "failed to run $src, or it failed to return a closure";
        }
    }
    else {
        $code = eval $src;
        die "error in parsing code block \"$src\": $@" if $@;
        die "\"$src\" fails to return a closure"
          unless ref($code) eq "CODE";
    }
    return $code;
}

=head2 clone_loc

Title   : clone_loc
Usage   : my $l = clone_loc($feature->location);
Function: Helper function to simplify the task of cloning locations
           for --locfunc closures.

          Presently simply implemented using Storable::dclone().
Example :
Returns : A L<Bio::LocationI> object of the same type and with the
          same properties as the argument, but physically different.
          All structured properties will be cloned as well.
Args    : A L<Bio::LocationI> compliant object

=cut

sub clone_loc {
    return dclone(shift);
}

=head2 create_loc

Title   : create_loc
Usage   : my $l = create_loc("10..12");
Function: Helper function to simplify the task of creating locations
          for --locfunc closures. Creates a location from a feature-
          table formatted string.

Example :
Returns : A L<Bio::LocationI> object representing the location given
          as formatted string.
Args    : A GenBank feature-table formatted string.

=cut

sub create_loc {
    return Bio::Factory::FTLocationFactory->from_string(shift);
}
