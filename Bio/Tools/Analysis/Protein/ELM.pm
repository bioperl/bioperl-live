#
# BioPerl module for Bio::Tools::Analysis::Protein::ELM
#
# Please direct questions and support issues to <bioperl-l@bioperl.org>
#
# Cared for by Richard Adams <richard.adams@ed.ac.uk>
#
# Copyright Richard Adams
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1     NAME

Bio::Tools::Analysis::Protein::ELM - a wrapper around the ELM server which predicts short functional motifs on amino acid sequences

=head1     SYNOPSIS

  # get a Bio::Seq object to start with, or a Bio::PrimaryI object.

  my $tool = Bio::Tools::Analysis::Protein::ELM->
      new(seq => $seqobj->primary_seq() );
  $tool->compartment(['ER', 'Golgi']);
  $tool->species(9606);
  $tool->run;
  my @fts = $tool->Result('Bio::SeqFeatureI');
  $seqobj->addSeqFeature(@fts);

=head1    DESCRIPTION

This module is a wrapper around the ELM server L<http://elm.eu.org/>
which predicts short functional motifs on amino acid sequences.

False positives can be limited by providing values for the species
and cellular compartment of the protein. To set the species attribute,
use either a L<Bio::Species> object or an NCBI taxon ID number.  To set
the cell compartment attribute (any number of compartments can be
chosen) use an array reference to a list of compartment names.

Results can be obtained either as raw text output, parsed into a
data structure, or as Bio::SeqFeature::Generic objects.

=head1 SEE ALSO

L<Bio::SimpleAnalysisI>,
L<Bio::WebAgent>

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

=head1 AUTHORS

Richard Adams, Richard.Adams@ed.ac.uk,

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

use strict;
package Bio::Tools::Analysis::Protein::ELM;
use vars qw(%cc);
use HTML::HeadParser;
use Bio::SeqFeature::Generic;
use HTTP::Request::Common qw(POST);
use IO::String;
use base qw(Bio::Tools::Analysis::SimpleAnalysisBase);

## valid cell compartments ##
%cc = (
      all            => 1,
      nucleus        => 'GO:0005634',
      extracellular  => 'GO:0005576',
      cytoplasm      => 'GO:0005737',
      peroxisome     => 'GO:0005777',
      glycosome      => 'GO:0020015',
      glyoxisome     => 'GO:0009514',
      golgi          => 'GO:0005794',
      er             => 'GO:0005783',
      lysosome       => 'GO:0005764',
      endosome       => 'GO:0005768',
      plasma_membrane=> 'GO:0005886',
		);

my $URL           = 'http://elm.eu.org/cgimodel.py';
my $ANALYSIS_NAME = 'ELM';
my $INPUT_SPEC    =
    [
     {
      'mandatory' => 'true',
      'type'      => 'Bio::PrimarySeqI',
      'name'      => 'seq',
     },
     {
      'mandatory' => 'false',
      'type'      => 'taxon_id or Bio::Species object',
      'name'      => 'species',
      'default'   => '9606',
     },
     {
      'mandatory' => 'false',
      'type'      => 'string',
      'name'      => 'compartment',
      'default'   => [1],
     },
    ];

my  $RESULT_SPEC =
    {
     ''                 => 'bulk',              # same as undef
     'Bio::SeqFeatureI' => 'ARRAY of Bio::SeqFeature::Generic',
     'parsed'           => '{motif1_name=>{locus=>[],
					   peptide=>[],
					   regexp=>[]
					  },
			    }',
    };
my $ANALYSIS_SPEC= {name        => 'ELM',
				    type        => 'Protein',
                    version     => 'n/a',
                    supplier    =>'BioComputing Unit, EMBL',
					description =>'Prediction of linear functional motifs
                                  in proteins',
					reference   => 'NAR, 31:3625-3630'};


sub _init {
    my $self = shift;
    $self->url($URL);
    $self->{'_ANALYSIS_SPEC'} = $ANALYSIS_SPEC;
    $self->{'_INPUT_SPEC'}    = $INPUT_SPEC;
    $self->{'_RESULT_SPEC'}   = $RESULT_SPEC;
    $self->{'_ANALYSIS_NAME'} = $ANALYSIS_NAME;
    return $self;
}

=head2        compartment

 name        : compartment
 usage       : $elm->compartment(['golgi', 'er']);
 purpose     : get/setter for cell compartment specifications
 arguments   : None, single compartment string or ref to array of
               compartment names.
 returns     : Array of compartment names (default if not previously set).

=cut

sub compartment {

    my ($self, $arg) = @_;
    if ($arg) {

        # convert to array ref if not one already
	if (ref ($arg) ne 'ARRAY') {
            $arg = [$arg];
	}

        ## now add params if valid
	for my $param (@$arg) {
            if (exists($cc{lc($param)})) {
                push @{$self->{'_compartment'}} , $cc{$param};
            } else {
                $self->warn("invalid argument ! Must be one of " .
                            join "\n", keys %cc );
            }
        }                       #end of for loop

    }                           #endif $arg
    return defined($self->{'_compartment'})? $self->{'_compartment'}
        : $self->input_spec()->[2]{'default'};

}

=head1      species

 name      : species
 usage     : $tool->species('9606');
 purpose   : get/setter for species selction for ELM server
 arguments : none, taxon_id or Bio::Species object
 returns   : a string of the ncbi taxon_id

=cut

sub species {
    my ($self, $arg) = @_;

    if ($arg) {
        if (ref($arg) && $arg->isa('Bio::Species')) {
            $self->{'_species'} = $arg->ncbi_taxid();
        } elsif ($arg =~ /^\d+$/) {
            $self->{'_species'} = $arg;
        } else {
            $self->warn("Argument must be a Bio::Species object or ".
                        " an integer NCBI taxon id. ");
        }
    }                           #end if $arg
    return defined($self->{'_species'})?$self->{'_species'}
        :$self->input_spec()->[1]{'default'};

}

sub  _run {
    my $self  = shift;
    $self->delay(1);
    # delay repeated calls by default by 3 sec, set delay() to change
    #$self->sleep;
    $self->status('TERMINATED_BY_ERROR');

    #### this deals with being able to submit multiple checkboxed
    #### slections

    #1st of all make param array
    my @cc_str;
    my @cmpts = @{$self->compartment()};
    for (my $i = 0; $i <= $#cmpts ; $i++) {
        splice @cc_str, @cc_str, 0,   'userCC',$cmpts[$i];
    }
    my %h = (swissprotId      => "",
             sequence         => $self->seq->seq,
             userSpecies      => $self->species,
             typedUserSpecies => '',
             fun              => "Submit");
    splice (@cc_str, @cc_str,0, ( map{$_, $h{$_}} keys %h));


    my $request = POST $self->url(),
        Content_Type => 'form-data',
            Content  => \@cc_str;
    $self->debug( $request->as_string);
    my $r1 = $self->request($request);
    if ( $r1->is_error  ) {
	$self->warn(ref($self)." Request Error:\n".$r1->as_string);
	return;
    }

    my $text = $r1->content;
    my ($url) = $text =~ /URL=\S+(fun=\S+r=\d)/s;
    #$url =~ s/amp;//g ;
    my ($resp2);
    $url = $URL . "?" .$url;
    while (1) {
	my $req2 = HTTP::Request->new(GET=>$url);
	my $r2 = $self->request ($req2);
	if ( $r2->is_error ) {
	    $self->warn(ref($self)." Request Error:\n".$r2->as_string);
	    return;
	}
	$resp2 = $r2->content();

	if ($resp2 !~ /patient/s) {
	    $self->status('COMPLETED');
	    $resp2=~ s/<[^>]+>/ /sg;
            $self->{'_result'} = $resp2;
	    return;
	} else {
	    print "." if $self->verbose > 0;
	    $self->sleep(1);
	}
    }
}

=head1      result

 name      : result
 usage     : $tool->result('Bio::SeqFeatureI');
 purpose   : parse results into sequence features or basic data format
 arguments : 1. none    (retrieves raw text without html)
             2. a value (retrieves data structure)
             3. 'Bio::SeqFeatureI' (returns array of sequence features)
                tag names are : {method => 'ELM', motif => motifname,
                                 peptide => seqeunce of match,
                                 concensus => regexp of match}.
 returns   : see arguments.

=cut

sub result {
    my ($self, $val) = @_;
    if ($val) {
        if (!exists($self->{'_parsed'}) ) {
            $self->_parse_raw();
        }
        if ($val eq 'Bio::SeqFeatureI') {
            my @fts;
            for my $motif (keys %{$self->{'_parsed'}}) {
                for (my $i = 0; $i< scalar @{$self->{'_parsed'}{$motif}{'locus'}};$i++) {
                    my ($st, $end) = split /\-/, $self->{'_parsed'}{$motif}{'locus'}[$i];
                    push @fts, Bio::SeqFeature::Generic->new
                        (
                         -start       => $st,
                         -end         => $end,
                         -primary_tag => 'Domain',
                         -source      => 'ELM',
                         -tag   => {
                                    method    => 'ELM',
                                    motif     => $motif,
                                    peptide   => $self->{'_parsed'}{$motif}{'peptide'}[$i],
                                    concensus => $self->{'_parsed'}{$motif}{'regexp'}[0],
                                   });
                }
            }
            return @fts;
        }                       #end if BioSeqFeature
        return $self->{'_parsed'};
    }                           #endif ($val)
    return $self->{'_result'};
}

## internal sub to parse raw data into internal data structure which is cached.
sub _parse_raw {
    my $self = shift;
    my $result = IO::String->new($self->{'_result'});
    my $in_results = 0;
    my $name;
    my %results;
    my $last;
    while (my $l = <$result>) {
        next unless  $in_results > 0 ||$l =~ /^\s+Elm\s+Name\s+Instances/;
        $in_results++;          #will be set whnstart of results reached.
        last if $l =~ /List of excluded/;
        next unless $in_results >1;

        my @line_parts = split /\s+/, $l;
        shift @line_parts;
        ## if result has motif name on 1 line
        if (scalar @line_parts == 1 && $line_parts[0]=~ /^\s*(\w+_\w+)/) {
            $name = $1;
            next;
        }
        ## else if is line with loci /seq matches
        elsif (@line_parts > 1) {
            my $index = 0;      ## array index
            my $read_loci = 0;  ## flag to know that loci are being read
            while ($index <= $#line_parts) {
                my $word = $line_parts[$index++];
                if ($read_loci ==0 && $word =~/_/) {
                    $name = $word;
                } elsif ($read_loci == 0 && $word =~ /^\w+$/ ) {
	            push @{$results{$name}{'peptide'}}, $word;
                } elsif ($word =~ /\d+\-\d+/) {
                    $read_loci = 1;
                    push @{$results{$name}{'locus'}}, $word;
                } else {        ## only get here if there are elements
                    last;
                }
            }                   #end of while
            push @{$results{$name}{'regexp'}}, $line_parts[$#line_parts];
        }                       #end of elsif

    }                           #end of while

    $self->{'_parsed'} = 	\%results;
}
1;
