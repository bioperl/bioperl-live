# $Id$
#
# BioPerl module for Bio::Tools::Analysis::Protein::NetPhos
#
# Cared for by Heikki Lehvaslaiho <heikki@ebi.ac.uk>
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
  $netphos->seq($seq)->threshold($threshod)->run;

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

  bioperl-l@bioperl.org                       - General discussion
  http://bio.perl.org/MailList.html           - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via email
or the web:

  bioperl-bugs@bio.perl.org
  http://bugzilla.bioperl.org/

=head1 AUTHORS

Richard Adams, Richard.Adams@ed.ac.uk, 
Heikki Lehvaslaiho, heikki@ebi.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Tools::Analysis::Protein::NetPhos;
use vars qw(@ISA $URL $FLOAT %STATUS $ANALYSIS_SPEC $INPUT_SPEC $RESULT_SPEC);
use strict;

use Data::Dumper;
use Bio::WebAgent;
use IO::String;
use Bio::SeqIO;
use HTTP::Request::Common qw (POST);
use Bio::SeqFeature::Generic;
use Bio::SimpleAnalysisI;

@ISA = qw(Bio::WebAgent Bio::SimpleAnalysisI);

BEGIN {

    $URL = 'http://www.cbs.dtu.dk/cgi-bin/nph-webface';

    $FLOAT = '[+-]?\d*\.\d*';

    %STATUS =  map { $_ => 1 } qw(CREATED COMPLETED TERMINATED_BY_ERROR);

    $ANALYSIS_SPEC =
        {'name' => 'NetPhos',
         'type' => 'Protein',
         'version' => '2.0',
         'supplier' => 'Center for Biological Sequence Analysis, Technical University of Denmark',
         'description' => 'Prediction of serine, threonine and tyrosine phosphorylation sites in eukaryotic proteins',
        };

    $INPUT_SPEC =
        [
         {
          'mandatory' => 'true',
          'type' => 'Bio::PrimarySeqI',
          'name' => 'seq'
         },
         {
          'mandatory' => 'false',
          'type' => 'float',
          'name' => 'threshold'
         }
        ];

    $RESULT_SPEC =
        {
         '' => 'bulk',  # same as undef
         'Bio::SeqFeatureI' => 'ARRAY of Bio::SeqFeeature::Generic',
         'raw' => 'Array of [ position, score, residue ]'
        };

    # Setting the default wait before calling again.
    # Change by setting delay() in seconds
    sub delay_policy {return 3 ;}

}


sub new {
    my $class = shift;

    my $self = $class->SUPER::new();
    while( @_ ) {
	my $key = lc shift;
        $key =~ s/^-//;
        $self->$key(shift);
    }

    $self->url($URL);

    return $self; # success - we hope!
}

=head1 Local attribute methods

=cut

=head2 seq

 Title   : seq
 Usage   : $obj->seq($seq);
 Function: Set and read the input sequence object to be analyzed
 Returns : sequence object or undef
 Args    : sequence object (optional)

=cut

sub seq {
   my ($self,$value) = @_;
   if( defined $value) {
       $self->throw("I need a Bio::PrimarySeqI, not  [". $value. "]")
	   unless $value->isa('Bio::PrimarySeqI');
       $self->throw("I need a protein seq, not  [". $value->alphabet. "]")
	   unless $value->alphabet eq 'protein';
       $self->{'_seq'} = $value;
       return $self;
   }
   return $self->{'_seq'} ;
}

=head2 threshold

 Title   : threshold
 Usage   : $obj->threshold($threshold);
 Function: Set and read threshold value for prediction scores, defaults to 0.8
 Returns : sequence object or undef
 Args    : sequence object (optional)

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
   return $self->{'_threshold'} || 0.8 ;
}


=head1 Methods from L<Bio::SimpleAnalysisI>

=cut


=head2 analysis_name

 Usage   : $tool->analysis_name;
 Returns : a name of this analysis
 Args    : none

=cut

sub analysis_name { return 'NetPhos'; }

=head2 analysis_spec

 Usage   : $tool->analysis_spec;
 Returns : a hash reference describing this analysis
 Args    : none

The returned hash uses the following keys (not all of them always
present, perhaps others present as well): C<name>, C<type>, C<version>,
C<supplier>, C<installation>, C<description>.

=cut

sub analysis_spec { return $ANALYSIS_SPEC }


=head2 input_spec

 Usage   : $tool->input_spec;
 Returns : an array reference with hashes as elements
 Args    : none

The analysis input data are named, and can be also associated with a
default value, with allowed values and with few other attributes. The
names are important for feeding the analysis with the input data (the
inputs are given to methods C<run> and C<wait_for> as name/value
pairs).

In this implementation, the parameters can equally well be passed to
the contructor or set independetly after object creation. They will be
checked against input_spec before running the analysis.


The hash returned by the input_spec gives you (for each in/out data and
parameter) just:

=over 5

=item name

e.g. 'name' =E<gt> 'sequence'

=item type

e.g. 'type' =E<gt> 'float' (there is no predefined set of types - that
would be a bigger task..., BioMoby tries to address it, for example)

=item default value

e.g. 'default_value' =E<gt> 3

=item mandatory

true or false, e.g. 'mandatory' =E<gt> 'true'

=item allowed values

An array ref with elements representing all possible values; this is,
however, used only for discrete values.

=back

=cut

sub input_spec { $INPUT_SPEC }


=head2 result_spec

 Usage   : $tool->result_spec;
 Returns : a hash reference with result names as keys
           and result types as values
 Args    : none

=cut

sub result_spec { $RESULT_SPEC }

=head2 run

 Usage   : $tool->run ( {seq=>$my.seq', threshold=>0.9} )
 Returns : $self
 Args    : data and parameters for this execution
           (in various formats)

=cut


sub run {
    my ($self, $args) = @_;

    $self->_process_arguments ($args) if $args;

    # check input
    $self->throw("Need a sequence object as an input") unless $self->seq;
    $self->debug(Data::Dumper->Dump([$self],[$self]));

    # internal run()
    $self->_run;
    return $self;
}

sub wait_for {
    my ($self, $args) = @_;
    $self->run($args);
}

sub _process_arguments {
    my ($self, $args) = @_;

    my %spec;
    map {$spec{ $_->{'name'} } = $_ } @{$self->input_spec};

    $self->debug(Data::Dumper->Dump([\%spec, $args],[\%spec, $args]));
    foreach my $key (keys %$args) {
        my $value = $args->{$key};

        $self->throw("Unknown argument [$key]")
            unless $spec{$key};
        $self->$key($value);
    }

    foreach my $key (keys %spec) {
        $self->throw("Mandatory argument [$key] is not set")
            if $spec{$key}{'mandatory'} eq 'true' and not defined $self->$key;
    }
}


sub _run {
    my $self = shift;

    # format the sequence into fasta
    my $seq_fasta;
    my $stringfh = new IO::String($seq_fasta);
    my $seqout = new Bio::SeqIO(-fh => $stringfh,
                                -format => 'fasta');
    $seqout->write_seq($self->seq);
    $self->debug($seq_fasta);

    # delay repeated calls by default by 3 sec, set delay() to change
    $self->sleep;

    $self->status('TERMINATED_BY_ERROR');

    my $request = POST $self->url,
        Content_Type => 'form-data',
            Content => [configfile => '/usr/opt/www/pub/CBS/services/NetPhos-2.0/NetPhos.cf',
                        SEQPASTE => $seq_fasta];
    my $content = $self->request($request);
    my $text = $content->content;

    my ($result_url) = $text =~ /follow <a href="(.*?)"/;
    return 0 unless $result_url;
    $self->debug("url is $result_url\n\n");

    my $ua2 = $self->clone;
    my $content2 = $ua2->request(POST $result_url);

    my $ua3 = $self->clone;
    $result_url =~ s/&.*//;
    $self->debug("final result url is $result_url\n");
    my $content3 = $ua3->request(POST $result_url);
    #print Dumper $content3;
    my $response = $content3->content;

#    open FH, '/home/heikki/netphos/netfosresults.html'|| die "could not open resultfile";
#    while (<FH>) {
#        $response .= $_;;
#    }

    $response =~ s/.*<pre>(.*)<\/pre>.*/$1/s;
    $response =~ s/<.*?>//gs;

    $self->{'_result'} = $response;

    $self->status('COMPLETED') if $response ne '';

    #print Dumper $response;
}

=head2 status

 Usage   : $tool->status
 Returns : string describing a status of the execution
 Args    : none

It returns one of the following strings (and perhaps more if a server
implementation extended possible job states):

  CREATED              (not run yet)
  COMPLETED            (run and finished normally)
  TERMINATED_BY_ERROR  (run and finished with an error or a signal)

=cut


sub status {
   my ($self,$value) = @_;

   if( defined $value) {
       $self->throw("Not a valid status value [$value]\n".
		    "Valid values are ". join(", ", keys %STATUS ))
	   unless defined $STATUS{$value};
       $self->{'_status'} = $value;
   }
   return $self->{'_status'} || 'CREATED' ;
}


=head2 result

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

The agrument string defined the type of bioperl objects returned in an
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
            /^\w+ +(\d+) +\w+ +(0\.\d+) +.([STY])/;
            next unless $3 and $2 > $self->threshold;
            print;
            push @predictions, [$1, $2, $3];
        }
        if ($value eq 'Bio::SeqFeatureI') {
            foreach  (@predictions) {
                push @fts, Bio::SeqFeature::Generic->new
                    (-start => $_->[0],
                     -end => $_->[0] ,
                     -source => 'NetPhos',
                     -primary => 'NetPhos_p',
                     -tag => {
                              score => $_->[1],
                              residue => $_->[2] });
            }
            return @fts;
        }
        return @predictions;
    }

    return $self->{'_result'};
}



1;
