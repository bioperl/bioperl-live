# $Id: Mitoprot.pm,
#
# BioPerl module for Bio::Tools::Analysis::Protein::Mitoprot
# Copyright Richard Adams
#
# You may distribute this module under the same terms as perl itself
# POD documentation - main docs before the code

=head1 NAME

Bio::Tools::Analysis::Protein::Mitoprot - a wrapper around Mitoprot
server

=head1 SYNOPSIS

  use Bio::Tools::Analysis::Protein::Mitoprot;

  use Bio::PrimarySeq;
  my $seq = Bio::PrimarySeq->new
    (-seq=>'IKLCVHHJHJHJHJHJHJHNLAILAKAHLIELALAL',
     -primary_id=>'test'); # a Bio::PrimarySeqI object

  my $mitoprot = Bio::Tools::Analysis::Protein::Mitoprot->new
     ( -seq => $seq
     ); # sequence must be  >!5aa long and start with an M.

  # run Mitoprot prediction on a DNA sequence
  my $mitoprot->run();


  die "Could not get a result" unless $mitoprot->status =~ /^COMPLETED/;

  print $mitoprot->result;     # print raw prediction to STDOUT

  foreach my $feat ( $mitoprot->result('Bio::SeqFeatureI') ) {

      # do something to SeqFeature
      # e.g. print as GFF
      print $feat->gff_string, "\n";
      # or store within the sequence - if it is a Bio::RichSeqI
      $seq->add_SeqFeature($feat);

 }

=head1 DESCRIPTION

This class is a wrapper around the Mitoprot web server which
calculates the probability of a sequence containing a mitochondrial
targetting peptide. See http://mips.gsf.de/cgi-bin/proj/medgen/mitofilter
for more details.

The results can be obtained in 3 formats:

=over 3

=item 1

The raw text of the program output

  my $rawdata = $analysis_object->result;

=item 2

An reference to a hash of  scores :

  my $data_ref = $analysis_object->result('parsed'); print "predicted
  export prob is $data_ref->{'export_prob'}\n"; #

key values of returned hash are input_length, basic_aas, acidic_aas,
export_prob, charge, cleavage_site.

=item 3

A Bio::SeqFeature::Generic object

  my $ft = $analysis_object->result(Bio::SeqFeatureI);
  print "export prob is ", ($ft->each_tag_value('export_prob'))[0]  ,"\n";


This the second implentation of Bio::SimpleAnalysisI which hopefully
will make it easier to write wrappers on various services. This class
uses a web resource and therefore inherits from Bio::WebAgent.

=back

=head1 SEE ALSO

L<Bio::SimpleAnalysisI>,
L<Bio::Tools::Analysis::SimpleAnalysisBase>,
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


# Let the code begin...


package Bio::Tools::Analysis::Protein::Mitoprot;
use vars qw($FLOAT);
use strict;

use IO::String;
use Bio::SeqIO;
use HTTP::Request::Common qw(GET);
use Bio::SeqFeature::Generic;

use base qw(Bio::Tools::Analysis::SimpleAnalysisBase);
$FLOAT = '[+-]?\d*\.\d*';

my $URL = 'http://ihg.gsf.de/cgi-bin/paolo/mitofilter?';


my %STATUS =  map { $_ => 1 } qw(CREATED COMPLETED TERMINATED_BY_ERROR);

my $MIN_LEN = 60;               #min len for protein analysis
my $ANALYSIS_NAME = "Mitoprot";

my $ANALYSIS_SPEC =
    {
     'name'        => 'Mitoprot',
     'type'        => 'Protein',
     'version'     => '1.0a4',
     'supplier'    => 'Munich Information Center for ProteinSequences',
     'description' => 'mitochondrial sig seq prediction',
    };

my $INPUT_SPEC =
    [
     {
      'mandatory' => 'true',
      'type'      => 'Bio::PrimarySeqI',
      'name'      => 'seq',          #value must be name of method used to set value
     },
    ];

my $RESULT_SPEC =
    {
     '' => 'raw text results',  # same as undef
     'Bio::SeqFeatureI' => 'ARRAY of Bio::SeqFeature::Generic',
     'all' => 'hash of results',
    };



### unique to this module ##

=head2 result

 Usage   : $job->result (...)
 Returns : a result created by running an analysis
 Args    : various

The method returns a result of an executed job. If the job was
terminated by an error the result may contain an error message instead
of the real data.

This implementation returns differently processed data depending on
argument:

=over 3

=item undef

Returns the raw ASCII data stream but without HTML tags

=item 'Bio::SeqFeatureI'

The argument string defines the type of bioperl objects returned in an
array.  The objects are L<Bio::SeqFeature::Generic>.  Feature primary
tag is "SigSeq".  Feature tags are input_length , basic_aas,
acidic_aas, export_prob, charge, cleavage_site, method.

=item 'parsed'

hash references of parsed results { input_length =E<gt>, basic_aas=E<gt>,
acidic_aas=E<gt>, export_prob=E<gt>, charge=E<gt>, cleavage_site=E<gt>}.

=back

=cut


sub result {
    my ($self,$value) = @_;
    #make sec feat of above threshold scores #

    my @sig_pdctns;
    my @fts;

    if ($value ) {
        my $result = IO::String->new($self->{'_result'});
        my %results;
        while (my $line = <$result>) {
            #make array of all scores or threshold depending on $value
            next unless $line =~ /\d/ || $line =~ /^Cle/;
            if ($line =~ /^Net[^+\-\d]+  # Net, then anything except +,- or digit
                          ((\+|-)?\d+)/x) #then get charge with optional + or -
                              {
                               $results{'charge'} = $1;
                              } elsif ($line =~ /^Input[^\d]+(\d+)/ ) {
                                  $results{'input_length'} = $1;
                              } elsif ($line =~ /basic[^\d]+(\d+)$/ ) {
                                  $results{'basic_aas'} = $1;
                              } elsif ($line =~ /acidic[^\d]+(\d+)$/) {
                                  $results{'acidic_aas'} = $1;
                              } elsif ($line =~ /^Cleavage[^\d]+(\d+)$/) {
                                  $results{'cleavage_site'} = $1;
                              } elsif ($line =~ /^Cleavage/) {
                                  $results{'cleavage_site'} = 'not predictable';
                              } elsif ($line =~ /^of export[^\d]+((0|1)\.\d+)$/) {
                                  $results{'export_prob'} = $1;
                              }
        }

        if ($value eq 'Bio::SeqFeatureI') {
            push @fts, Bio::SeqFeature::Generic->new
                (
                 -start => 1,
                 -end => ($results{'cleavage_site'} =~
                          /^\d+$/)?$results{'cleavage_site'}:$self->seq->length,
                 -source => 'Mitoprot',
                 -primary => 'Region',
                 -tag =>{
                         export_prob   => $results{'export_prob'},
                         charge        => $results{'charge'},
                         basic_aas     => $results{'basic_aas'},
                         acid_aas      => $results{'acidic_aas'},
                         region_name   => 'Transit_peptide',
                         method        => 'MitoProt',
                         cleavage_site => $results{'cleavage_site'},
                        },
                );
            return @fts;        #return Bioseqfeature array
        }
        ## convert parsed data into a meta array format
        else  {
            return \%results;   # hash based results ref
        }
    }
    return $self->{'_result'};
}

sub _init {
    my $self = shift;
    $self->url($URL);
    $self->{'_ANALYSIS_SPEC'} =$ANALYSIS_SPEC;
    $self->{'_INPUT_SPEC'} =$INPUT_SPEC;
    $self->{'_RESULT_SPEC'} =$RESULT_SPEC;
    $self->{'_ANALYSIS_NAME'} =$ANALYSIS_SPEC->{'name'};
    return $self;
}

sub _process_arguments {
    #extra checking for sequence length
    #mitoprot specific argument testing
    my ($self, $args) = @_;
    #use base checking for existence of mandatory fields
    $self->SUPER::_process_arguments($args) ;

    #then check specifics
    $self->throw ("1st_aa must be M") if $self->seq->subseq(1,1) !~ /M/i;
    $self->throw ("sequence must be at least 15aa long") if $self->seq->length< 15;
    return;
}



sub _run {
    #request submitted by get not by post
    my $self  = shift;
    $self->delay(1);
    $self->sleep;

    $self->status('TERMINATED_BY_ERROR');
    my $url = $self->url . "seq=".lc($self->seq->seq). "&seqnam=";
    my $request = GET $url;
    my $content = $self->request($request);
    my $text = $content->content; #1st reponse

    #remove html stuff
    $text =~ s/.*<PRE>(.*)<\/PRE>.*/$1/s;
    $text =~ s/<[^>]+>//sg;

    $self->status('COMPLETED') if $text ne '' && $self->seq->length > $MIN_LEN;
    $self->{'_result'} = $text;

}

1;
