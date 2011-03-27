
=head1 NAME

Bio::FeatureIO::interpro - read features from InterPro XML

=head1 SYNOPSIS

  my $in = Bio::FeatureIO(-format=>'interpro');
  while (my $feat = $in->next_feature) {
    # do something with the Bio::SeqFeatureI object
  }

=head1 DESCRIPTION

See L<http://www.ebi.ac.uk/interpro/documentation.html>.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

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
of the bugs and their resolution. Bug reports can be submitted via
the web:

  https://redmine.open-bio.org/projects/bioperl/

=head1 AUTHOR - Allen Day

Email allenday@ucla.edu

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...

package Bio::FeatureIO::interpro;

use strict;
use base qw(Bio::FeatureIO);
use Bio::SeqFeature::Annotated;
use Bio::OntologyIO;

use Bio::Annotation::Comment;
use Bio::Annotation::DBLink;
use Bio::Annotation::OntologyTerm;
use Bio::Annotation::SimpleValue;
use Bio::Annotation::Target;

use URI::Escape;
use XML::DOM;
use XML::DOM::XPath;

sub _initialize {
  my($self,%arg) = @_;

  $self->SUPER::_initialize(%arg);
  $self->xml_parser(XML::DOM::Parser->new());
  my $buf;
  while(($buf = $self->_readline()) && $buf !~ /<protein/){
    next;
  }
  $self->_pushback($buf);
}

sub next_feature {
  my $self =shift;
  my $buf;    #line buffer
  my $ok = 0; #true if there is another <protein/> record in stream
  my $record; #holds the record to be parsed and returned.

  #try to dump buffer from last record before moving on to next record
  my $f = $self->_shift_feature_buffer();
  if($f){
    return $f;
  }

  while(my $buf = $self->_readline()){
    $ok = 1 if $buf =~ m!<protein!;
    $record .= $buf;
    last if $buf =~ m!</protein>!;
  }
  return unless $ok;

  my $dom = $self->xml_parser->parse($record);


  my ($pNode) = $dom->findnodes('/protein');

  my @iNodes = $pNode->findnodes('/protein/interpro');

  foreach my $iNode (@iNodes){
    my @cNodes = $iNode->findnodes('classification');
    my @mNodes = $iNode->findnodes('match');

    #we don't handle these
    #my @nNodes = $iNode->findnodes('contains');
    #my @fNodes = $iNode->findnodes('found_in');

    foreach my $mNode (@mNodes){
      my @lNodes = $mNode->findnodes('location');
      foreach my $lNode (@lNodes){
        my $feature = Bio::SeqFeature::Annotated->new(
                                                      -start  => $lNode->getAttribute('start'),
                                                      -end    => $lNode->getAttribute('end'),
                                                      -score  => $lNode->getAttribute('score'),
#                                                      -seq_id => $pNode->getAttribute('id'),
                                                     );
        $feature->seq_id->value($pNode->getAttribute('id'));

#warn $pNode->getAttribute('id');

        $feature->source( $lNode->getAttribute('evidence') );

        my $t = Bio::Annotation::OntologyTerm->new(-identifier => 'SO:0000417', -name => 'polypeptide_domain');
        $feature->add_Annotation('type',$t);

        my $c = Bio::Annotation::Comment->new(-tagname => 'comment', -text => $iNode->getAttribute('name'));
        $feature->add_Annotation($c);

        my $d = Bio::Annotation::DBLink->new();
        $d->database($mNode->getAttribute('dbname'));
        $d->primary_id($mNode->getAttribute('id'));
        $d->optional_id($mNode->getAttribute('name'));
        $feature->annotation->add_Annotation('dblink',$d);

        my $s = Bio::Annotation::SimpleValue->new(-tagname => 'status', -value => $lNode->getAttribute('status'));
        $feature->annotation->add_Annotation($s);

        foreach my $cNode (@cNodes){
          my $o = Bio::Annotation::OntologyTerm->new(-identifier => $cNode->getAttribute('id'));
          $feature->annotation->add_Annotation('ontology_term',$o);
        }

        $self->_push_feature_buffer($feature);
      }
    }
  }

  return $self->_shift_feature_buffer;
}

=head2 _push_feature_buffer()

 Usage   :
 Function:
 Returns : 
 Args    :


=cut

sub _push_feature_buffer {
  my ($self,$f) = @_;

  if(ref($f)){
    push @{ $self->{feature_buffer} }, $f;
  }
}

=head2 _shift_feature_buffer()

 Usage   :
 Function:
 Returns : 
 Args    :


=cut

sub _shift_feature_buffer {
  my ($self) = @_;
  return $self->{feature_buffer} ? shift @{ $self->{feature_buffer} } : undef;
}

=head2 xml_parser()

 Usage   : $obj->xml_parser($newval)
 Function: 
 Example : 
 Returns : value of xml_parser (a scalar)
 Args    : on set, new value (a scalar or undef, optional)


=cut

sub xml_parser {
  my($self,$val) = @_;
  $self->{'xml_parser'} = $val if defined($val);
  return $self->{'xml_parser'};
}

1;
