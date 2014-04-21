#
# BioPerl module for Bio::Tools::Analysis::DNA::ESEfinder
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Heikki Lehvaslaiho <heikki-at-bioperl-dot-org>
#
# Copyright Richard Adams
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Tools::Analysis::DNA::ESEfinder - a wrapper around ESEfinder
server

=head1 SYNOPSIS

  use Bio::Tools::Analysis::DNA::ESEfinder;
  use strict;

  my $seq; # a Bio::PrimarySeqI or Bio::SeqI object

  $seq = Bio::Seq->new(
       -primary_id => 'test',
       -seq=>'atgcatgctaggtgtgtgttttgtgggttgtactagctagtgat'.
       -alphabet=>'dna');

  my $ese_finder = Bio::Tools::Analysis::DNA::ESEfinder->
      new(-seq => $seq);

  # run ESEfinder prediction on a DNA sequence
  $ese_finder->run();

  die "Could not get a result"
      unless $ese_finder->status =~ /^COMPLETED/;

  print $ese_finder->result;      # print raw prediction to STDOUT

  foreach my $feat ( $ese_finder->result('Bio::SeqFeatureI') ) {

      # do something to SeqFeature
      # e.g. print as GFF
      print $feat->gff_string, "\n";
      # or store within the sequence - if it is a Bio::SeqI
      $seq->add_SeqFeature($feat)

  }

=head1 DESCRIPTION

This class is a wrapper around the ESEfinder web server which uses
experimentally defined scoring matrices to identify possible exonic
splicing enhancers in human transcripts.

The results can be retrieved in 4 ways.

=over 4

=item 1.

C<$ese_finder-E<gt>result('')> retrieves the raw text output of the
program

=item 2.

C<$ese_finder-E<gt>result('all')> returns a Bio::Seq::Meta::Array object
with prediction scores for all residues in the sequence


=item 3.

C<$ese_finder-E<gt>result('Bio::SeqFeatureI')> returns an array of
Bio::SeqFeature objects for sequences with significant scores. Feature
tags are score, motif, SR_protein and method

=item 4.

C<$ese_finder-E<gt>result('raw')> returns an array of significant matches
with each element being a reference to [SR_protein, position, motif, 
score]

=back

See L<http://rulai.cshl.edu/tools/ESE2/>

This the second implentation of Bio::SimpleAnalysisI which hopefully
will make it easier to write wrappers on various services. This class
uses a web resource and therefore inherits from L<Bio::WebAgent>.

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
Heikki Lehvaslaiho, heikki-at-bioperl-dot-org

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...
#should have own 

package Bio::Tools::Analysis::DNA::ESEfinder;

use Data::Dumper;
use IO::String;
use Bio::SeqIO;
use HTTP::Request::Common qw (POST);
use HTML::HeadParser;
use Bio::SeqFeature::Generic;
use Bio::Seq::Meta::Array;
use Bio::WebAgent;
use strict;

#inherits directly from SimpleAnalysisBase
use base qw(Bio::Tools::Analysis::SimpleAnalysisBase);


#global vars are now file-scoped lexicals

my $URL = 'http://rulai.cshl.org/cgi-bin/tools/ESE/esefinder.cgi';
	
my $ANALYSIS_NAME = 'ESEfinder';


my $ANALYSIS_SPEC =
    {
     'name' => 'ESEfinder',
     'type' => 'DNA', #compulsory entry as is used for seq checking
     'version' => '2.0',
     'supplier' => 'Krainer lab, Cold Spring Harbor Laboratory, POBOX100, Bungtown Rd, COld Spring Harbor, NY, USA',
     'description' => 'to identify exonic splicing elements in human transcripts',
    };

my $INPUT_SPEC =
    [{
      'mandatory' => 'true',
      'type' => 'Bio::PrimarySeqI',
      'name' => 'sequence',
     }];

my $RESULT_SPEC =
    {
     '' => 'bulk',  # same as undef
     'Bio::SeqFeatureI' => 'ARRAY of Bio::SeqFeature::Generic',
     'raw' => 'Array of [ SR_protein, position, motif, score]',
     'all' => 'Bio::Seq::Meta::Array object'
    };


### unique to this module ##
sub _init {
    ## fills in fixed data for class ##
    my $self = shift;
    $self->url($URL);
    $self->{'_ANALYSIS_SPEC'} =$ANALYSIS_SPEC;
    $self->{'_INPUT_SPEC'} =$INPUT_SPEC;
    $self->{'_RESULT_SPEC'} =$RESULT_SPEC;
    $self->{'_ANALYSIS_NAME'} =$ANALYSIS_NAME;
    return $self;
}


sub _run {
    my $self  = shift;
    my $seq_fasta;
    my $stringfh = IO::String->new($seq_fasta);
    my $seqout = Bio::SeqIO->new(-fh => $stringfh,
                                -format => 'fasta');
    $seqout->write_seq($self->seq);
    $self->debug($seq_fasta);
    $self->delay(1);
    # delay repeated calls by default by 3 sec, set delay() to change
    $self->sleep;
    $self->status('TERMINATED_BY_ERROR');

    my $request = POST $self->url,
        #Content_Type => 'x-www-form-urlencoded',
            Content  => [
                         protein1 => 1,
                         protein2 => 1,
                         protein3 => 1,
                         protein4 => 1,
                         radio_sf2 => 0,
                         radio_sc35 => 0,
                         radio_srp40 => 0,
                         radio_srp55 => 0,
                         sequence =>$seq_fasta,
                        ];
    my $content = $self->request($request);
    if( $content->is_error  ) {
	$self->throw(ref($self)." Request Error:\n".$content->as_string);
    }

    my $text = $content->content; #1st reponse
    my ($tmpfile) = $text =~ /value="(tmp\/.+txt)"/;
    # now get data for all residues #
    my $rq2 = POST 'http://rulai.cshl.org/cgi-bin/tools/ESE/resultfile.txt',
        #Content_Type => 'x-www-form-urlencoded',
            Content => [
                        fname => $tmpfile,
                       ];
    my $ua2 = Bio::WebAgent->new();
    my $content2 = $ua2->request($rq2);
    if( $content2->is_error  ) {
	$self->throw(ref($self)." Request Error:\n".$content2->as_string);
    }

    my $text2 = $content2->content;
    $self->{'_result'} = $text2;		
    $self->status('COMPLETED') if $text2 ne '';

    #print Dumper $response;
}


sub result {

    #make sec feat of above threshold scores #

    my ($self,$value) = @_;

    my @sig_pdctns;
    my @fts;

    if ($value ) {
	my $result = IO::String->new($self->{'_result'});
        my $current_SR;
        my $all_st_flag = 0;
        my %all;
        while (my $line = <$result>) {
            #make array of all scores or threshold depending on $value
            last if $line =~ /^All scores/ && $value ne 'all' or $line =~ /2001,/;
            $all_st_flag++ if $line =~ /All scores/;
            next if $value eq 'all' && $all_st_flag == 0;

            #parse line
            if ($line =~ /^Protein/) {
                ($current_SR) = $line =~/:\s+(\S+)/;
                $current_SR =~ s{/}{_}; # remove unallowed charcters from hash
            }
            if ( $line =~/^\d+/ && $value ne 'all') {
                push @sig_pdctns, [$current_SR, split /\s+/, $line] ;
            } elsif ($line =~ /^\d+/) {

                push @{$all{$current_SR}}, [split /\s+/, $line];
            }
        }

        if ($value eq 'Bio::SeqFeatureI') {
            foreach (@sig_pdctns) {
                #make new ese object for each row of results
                push @fts, Bio::SeqFeature::Generic->new
                    (
                     -start => $_->[1],
                     -end => $_->[1] + length($_->[2]) -1,
                     -source => 'ESEfinder',
                     -primary => 'ESE',
                     -tag =>{
                             score =>$_->[3],
                             motif=> $_->[2],
                             SR_protein=> $_->[0],
                             method=> 'ESEfinder',
                            },
                    );
            }
            return @fts;
        }
        ## convert parsed data into a meta array format
        elsif ($value eq 'all') {
            bless ($self->seq, "Bio::Seq::Meta::Array");
            $self->seq->isa("Bio::Seq::MetaI")
                || $self->throw("$self is not a Bio::Seq::MetaI");

            for my $prot (keys %all) {
                my @meta;
                my $len =  scalar @{$all{$prot}} ;
                for (my $i = 0; $i < $len; $i++ ) {
                    $meta[$i] = $all{$prot}[$i][2];
                }

                # assign default name here so that the
                # Bio::Seq::Meta::Array can work for all classes
                # implementing it and we can avoid having to make
                # asubclass for each implementation

                $Bio::Seq::Meta::Array::DEFAULT_NAME = "ESEfinder_SRp55";
                my $meta_name = $self->analysis_spec->{'name'} . "_" . "$prot";
                $self->seq->named_meta($meta_name,\@meta );
            }
            # return  seq array object implementing meta sequence #
            return $self->seq;

        }
		#return ref to array of arrays
        return \@sig_pdctns;
    }
    return $self->{'_result'};
}


1;





