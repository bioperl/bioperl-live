#!/usr/bin/perl -w

use strict;

use Bio::OntologyIO::simpleGOparser;
use Bio::Root::IO;
use Bio::Tools::GFF;

use Getopt::Long;

our $USAGE = 'validate_gff_via_ontology.pl [-trace|t] [-type_only] [-mapping|m MYTYPE=REALTYPE] GFF-FILE ONTOLOGY-FILE';


eval { 
    require "GO/Parser.pm";
};
if ($@) {
    print <<EOM

Currently you need the go-dev perl modules to run this - go to

 http://godatabase.sf.net

Fetch the latest codebase from cvs and add go-dev/perl-api to your PERL5LIB

EOM
}

# NOTE:
# Currently requires go-dev perl modules; see www.godatabase.org/dev

my $trace;
my $help;
my $types_only;
my $type_mapping = {};
GetOptions("help|h"=>\$help,
           "trace|t"=>\$trace,
           "types_only"=>\$types_only,
           "mapping|m=s%"=>$type_mapping,
          );

if ($help) {
   exec('perldoc', $0);
   exit 0;
}

if (@ARGV != 2) {
    print $USAGE, "\n";
    die;
}

my $gfffn = shift @ARGV;
my $gffio = Bio::Tools::GFF->new(-file => $gfffn);
my @features = $gffio->features;
$gffio->close();


sub uniquify {
    my %h = ();
    grep {
        my $done = $h{$_};
        $h{$_} = 1;
        !$done;
    } @_;
}

sub trace {
    if ($trace) {
        my $fmt = shift;
        printf($fmt, @_);
    }
}

#foreach (@features) {
#    printf "%s %s\n", $_->type, $_->type_string;
#}

my $fn = shift @ARGV;

my $p = GO::Parser->new({handler=>'obj'});
$p->parse($fn);
my $graph = $p->handler->graph;

my @bad_types = ();
my @good_types = ();

my @bad_partofs = ();

foreach my $f (@features) {
    check($f);
}
@bad_types = uniquify @bad_types;
@good_types = uniquify @good_types;
@bad_partofs = uniquify @bad_partofs;
print "BAD : @bad_types\n";
print "GOOD: @good_types\n";

print "BAD PART-OFS: @bad_partofs\n" if @bad_partofs;


sub check {
    my $f = shift;
    my $type = $f->resolve_type($graph,
                                $type_mapping);
    trace "CHECKING %s of type %s [%s]\n", 
      $f->unique_id || '?', $f->type->name, $f->type->identifier || '?';
    if (!$type) {
        push(@bad_types, $f->type_string);
    }
    else {
        push(@good_types, sprintf("%s [%s]", $f->type_string, $f->type->identifier));
    }
    my @subf = $f->sub_SeqFeature;
    foreach (@subf) {
        trace "DESCENDING FEATURE GRAPH TO %s %s\n",
          $_->unique_id || '?', $_->type_string;
        check($_);
    }
    if ($types_only) {
        return;
    }
    
    # subfeatures should have types resolved by the call above;
    # let's check if the parent relationships conform to PART-OFs in SO

    my $parent_type_id = $f->type->identifier;
    if (!$parent_type_id) {
        trace "Can't check subfeatures of %s %s as type is not in ontology\n",
          $f->type_string, $f->unique_id || '?';
        return;
    }

    foreach my $subf (@subf) {
        my $ok = 0;
        my $child_type_id = $subf->type->identifier;
        trace "\n=======\nPARTOF CHECK BETWEEN %s AND %s\n",
          $subf->type->name, $f->type->name;
        if (!$child_type_id) {
            trace "Can't check subfeature %s %s as type is not in ontology\n",
              $subf->type_string, $subf->unique_id || '?';
            next;
        }
        my $child_type_all_poss = 
          $graph->get_reflexive_parent_terms_by_type($child_type_id,
                                                     'isa');
        my $parent_type_all_poss = 
          $graph->get_reflexive_parent_terms_by_type($parent_type_id,
                                                     'isa');

        # check rule:
        #  ChildTypeAllPossible part-of ParentTypeAllPossible
        foreach my $c (@$child_type_all_poss) {
            foreach my $p (@$parent_type_all_poss) {
                trace " checking if %20s is-partof %-20s...\n",
                  $c->name, $p->name;
                my $rels =
                  $graph->get_relationships_between_terms($p->acc, $c->acc);
                trace "   R:@$rels\n";
                my @po = grep {$_->type eq 'partof'} @$rels;
                if (@po) {
                    trace "   YES!\n";
                    $ok = 1;
                }
                else {
                    trace "   NO!\n";
                }
            }
        }
        if ($ok) {
            trace "**OK**\n";
        }
        else {
            trace "**NAUGHTY**\n";
            push(@bad_partofs, sprintf("[%s PARTOF %s]",
                                       $subf->type_string,
                                       $f->type_string));
        }
    }
}

__END__

=head1 NAME

  validate_gff_via_ontology.pl

=head1 SYNOPSIS

  validate_gff_via_ontology.pl -trace -m cds=coding_sequence myfeatures.gff sofa.ontology

=head1 DESCRIPTION

This script will validate a GFF (version3 is best) file against an
ontology of sequence features. Currently the ontology must be in
GO/GOBO format, support for other formats (eg DAML+OIL, OWL) are
forthcoming.

The validation is in two parts:

first of all the types in the GFF file are checked vs the terms in the
ontology. An identifier (which will be a SO accession if the SO/SOFA
ontology is used) is assigned if it exists (case insensitive matching
is used). If the type used in the GFF file does not exist in the
ontology file the script will notify you with a message like:

  BAD:  locus transcrpt geme

You can obtain the canonical SOFA ontology here:

 http://sofa.sf.net/

The second part of the validation checks to see if the 'subfeature'
relationships specified in the GFF (using Parent and ID attributes in
GFF3) are valid. For instance, making 'transcript' a subfeature of
'repeat' is obviously silly, and it is banned because there is no
'part-of' restriction between transcript and repeat in the canonical
SOFA ontology. The rules are a bit more subtle than this, as we also
have to traverse the subsumption hierarchy.

=head2 RULES

a subfeature can only be part of a superfeature
if there is a direct part-of link for the relevant feature types
OR such a link can be obtained by traversing either of the
two graphs upwards


 notes: R* is the reflexive transitive closure of R
 notes: R+ is the transitive closure of R

(reflexive includes the relationship x R x)

  WHERE ChildFeat is an object of type SeqFeature 
        ParentFeat is an object of type SeqFeature
  IF
    ChildFeat part-of ParentFeat

  THEN

  THIS MUST BE SATISFIED FOR SOME POSSIBLE BINDING OF THE VARIABLES BELOW
  (if it starts with a capital it is an unbound variable):
  
  ChildFeat has-feature-type ChildType
  ParentFeat has-feature-type ParentType

  ChildType is-a* ChildTypeAllPossible
  ParentType is-a* ParentTypeAllPossible

  # don't allow nonreflexive transitive closure on part-of
  # we *could* allow ChildTypeAllPossible part-of+ ParentTypeAllPossible
  # which would mean we could attach exons directly to genes

  ChildTypeAllPossible part-of ParentTypeAllPossible
  
  ;;;; RULES IN PROLOG:
  ;;;; (making predicates arguments so we can reason over preds)
  
  ;; closure of R is R+
  stmt(X, 'is-a+', Y):-
    stmt(X, 'is-a', Y).
  
  stmt(X, 'is-a+', Y):-
    stmt(X, 'is-a', Z),
    stmt(Z, 'is-a+', Y).
  
  ;; reflexive closure of R is R*
  stmt(X, 'is-a*', X).
  stmt(X, 'is-a*', Y):-
    stmt(X, 'is-a+', Y).
  
  ;; a subfeature can only be part of a superfeature
  ;; if there is a direct part-of link for the relevant feature types
  ;; OR such a link can be obtained by traversing either of the
  ;; two graphs upwards
  stmt(ChildFeat, 'part-of', ParentFeat):-
    stmt(ChildFeat, 'instance-of', ChildType),
    stmt(ParentFeat, 'instance-of', ParentType),
  
    stmt(ChildType, 'is-a*', ChildTypeAllPossible),
    stmt(ParentType, 'is-a*', ParentTypeAllPossible),
  
    stmt(ChildTypeAllPossible, 'is-a*', ParentTypeAllPossible).


=cut
