#
# BioPerl module for Bio::Tools::Analysis::Protein::NetPhos
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

Bio::Tools::Analysis::Protein::NetPhos - a wrapper around NetPhos server

=head1 SYNOPSIS

  use Bio::Tools::Analysis::Protein::NetPhos;

  my $seq; # a Bio::PrimarySeqI object
  my $threshold  = "0.90";

  my $netphos = Bio::Tools::Analysis::Protein::NetPhos->new
     ( -seq => $seq,
       -threshold => $threshold );

  # run NetPhos prediction on a sequence
  my $netphos->run();

  # alternatively you can say
  $netphos->seq($seq)->threshold($threshold)->run;

  die "Could not get a result" unless $netphos->status =~ /^COMPLETED/;

  print $netphos->result;     # print raw prediction to STDOUT

  foreach my $feat ( $netphos->result('Bio::SeqFeatureI') ) {

      # do something to SeqFeature
      # e.g. print as GFF
      print $feat->gff_string, "\n";
      # or store within the sequence - if it is a Bio::RichSeqI
      $seq->add_SeqFeature($feat)

 }

=head1 DESCRIPTION

This class is wrapper around the NetPhos 2.0 server which produces
neural network predictions for serine, threonine and tyrosine
phosphorylation sites in eukaryotic proteins.

See L<http://www.cbs.dtu.dk/services/NetPhos/>.

This the first implentation of Bio::SimpleAnalysisI which hopefully
will make it easier to write wrappers on various services. This class
uses a web resource and therefore inherits from Bio::WebAgent.

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


package Bio::Tools::Analysis::Protein::NetPhos;
use vars qw($FLOAT);
use strict;
use IO::String;
use Bio::SeqIO;
use HTTP::Request::Common qw (POST);
use Bio::SeqFeature::Generic;

use base qw(Bio::Tools::Analysis::SimpleAnalysisBase);

    $FLOAT = '[+-]?\d*\.\d*';
    my $URL = 'http://www.cbs.dtu.dk/cgi-bin/nph-webface';


    my $ANALYSIS_SPEC =
        {
         'name'        => 'NetPhos',
         'type'        => 'Protein',
         'version'     => '2.0',
         'supplier'    => 'Center for Biological Sequence Analysis,
                           Technical University of Denmark',
         'description' => 'Prediction of serine, threonine and tyrosine
                             phosphorylation sites in eukaryotic proteins',
        };

    my $INPUT_SPEC =
        [
         {
          'mandatory' => 'true',
          'type'      => 'Bio::PrimarySeqI',
          'name'      => 'seq',
         },
         {
          'mandatory' => 'false',
          'type'      => 'float',
          'name'      => 'threshold',
          'default'   => 0.8,
         }
        ];

    my $RESULT_SPEC =
        {
         ''                 => 'bulk',  # same as undef
         'Bio::SeqFeatureI' => 'ARRAY of Bio::SeqFeeature::Generic',
         'raw'              => 'Array of [ position, score, residue ]'
        };


=head2 result

 Name    : result
 Usage   : $job->result (...)
 Returns : a result created by running an analysis
 Args    : none (but an implementation may choose
           to add arguments for instructions how to process
           the raw result)

The method returns a scalar representing a result of an executed
job. If the job was terminated by an error the result may contain an
error message instead of the real data (or both, depending on the
implementation).

This implementation returns differently processed data depending on
argument:

=over 3

=item undef

Returns the raw ASCII data stream but without HTML tags

=item 'Bio::SeqFeatureI'

The argument string defined the type of bioperl objects returned in an
array.  The objects are L<Bio::SeqFeature::Generic>.

=item anything else

Array of array references of [ position, score, residue].


=back


=cut

sub result {
    my ($self,$value) = @_;

    my @predictions;
    my @fts;

    if ($value ) {

        my $result = IO::String->new($self->{'_result'});
        while (<$result>) {
            next if /^____/;
            /^\S+ +(\d+) +\w+ +(0\.\d+) +.([STY])/;
            next unless $3 and $2 > $self->threshold;
            push @predictions, [$1, $2, $3];
        }
        if ($value eq 'Bio::SeqFeatureI') {
            foreach  (@predictions) {
                push @fts, Bio::SeqFeature::Generic->new
                    (-start   => $_->[0],
                     -end     => $_->[0] ,
                     -source  => 'NetPhos',
                     -primary => 'Site',
                     -tag     => {
                               score   => $_->[1],
                               residue => $_->[2] });
            }
            return @fts;
        }
        return \@predictions;
    }

    return $self->{'_result'};
}

=head2  threshold

 Usage   : $job->threshold(...)
 Returns  : The significance threshold of a prediction
 Args     : None (retrieves value) or a value beween 0 and 1.
 Purpose  : Get/setter of the threshold to be sumitted for analysis.

=cut

sub threshold {
   my ($self,$value) = @_;
   if( defined $value) {
       if ( $value !~ /$FLOAT/ or $value < 0 or $value > 1 ) {
           $self->throw("I need a value between 0 and 1 , not  [". $value. "]")
       }
       $self->{'_threshold'} = $value;
       return $self;
   }
   return $self->{'_threshold'} || $self->input_spec->[1]{'default'} ;
}

sub _init 
		{
	my $self = shift;
	$self->url($URL);
	$self->{'_ANALYSIS_SPEC'} =$ANALYSIS_SPEC;
	$self->{'_INPUT_SPEC'} =$INPUT_SPEC;
	$self->{'_RESULT_SPEC'} =$RESULT_SPEC;
	$self->{'_ANALYSIS_NAME'} =$ANALYSIS_SPEC->{name};
	return $self;
}

sub _run {
    my $self = shift;

    # format the sequence into fasta
    my $seq_fasta;
    my $stringfh = IO::String->new($seq_fasta);
    my $seqout = Bio::SeqIO->new(-fh => $stringfh,
                                -format => 'fasta');
    $seqout->write_seq($self->seq);
    $self->debug($seq_fasta);

    # delay repeated calls by default by 3 sec, set delay() to change
    $self->sleep;

    $self->status('TERMINATED_BY_ERROR');

    my $request = POST $self->url,
            Content_Type => 'form-data',
            Content      => [configfile => '/usr/opt/www/pub/CBS/services/NetPhos-2.0/NetPhos.cf',
                             SEQPASTE   => $seq_fasta];
    my $content = $self->request($request);
    my $text    = $content->content;

    my ($result_url) = $text =~ /follow <a href="(.*?)"/;
    return 0 unless $result_url;
    $self->debug("url is $result_url\n\n");

    my $ua2      = $self->clone;
    my $content2 = $ua2->request(POST $result_url);

    my $ua3      = $self->clone;
    $result_url  =~ s/&.*//;
    $self->debug("final result url is $result_url\n");
    my $content3 = $ua3->request(POST $result_url);
    #print Dumper $content3;
    my $response = $content3->content;


    $response =~ s/.*<pre>(.*)<\/pre>.*/$1/s;
    $response =~ s/<.*?>//gs;

    $self->{'_result'} = $response;

    $self->status('COMPLETED') if $response ne '';

}

1;
