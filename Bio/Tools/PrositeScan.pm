
=head1 NAME

Bio::Tools::PrositeScan - Parser for ps_scan result

=head1 SYNOPSIS

  use Bio::Tools::PrositeScan;

  my $factory = Bio::Tools::PrositeScan->new(
      -file   => 'out.PrositeScan',
      -format => 'fasta'
  );

  while(my $match = $factory->next_prediction){
      #  $match is a Bio::SeqFeature::FeaturePair

      # Sequence ID
      my $seq_id = $match->seq_id;

      # PROSITE accession number
      my $psac = $match->hseq_id;

      # Coordinates
      my @coords = ( $match->start, $match->end );

      # Subsequence
      my $seq = $match->feature1->seq;
  }

=head1 DESCRIPTION

This is a parser of the output of the ps_scan program. It takes either a file
handle or a file name, and returns a L<Bio::SeqFeature::FeaturePair> object.

Note that the current implementation parses the entire file at once.

=head1 AUTHOR

Juguang Xiao, juguang@tll.org.sg

=head1 SEE ALSO

=over

=item * L<ps_scan software|ftp://ftp.expasy.org/databases/prosite/ps_scan>

=item * L<PROSITE User Manual|http://prosite.expasy.org/prosuser.html>

=back

=cut

# Let the code begin...

package Bio::Tools::PrositeScan;
use vars qw(@FORMATS);
use strict;
use Bio::Seq;
use Bio::SeqFeature::Generic;
use Bio::SeqFeature::FeaturePair;

use base qw(Bio::Root::Root Bio::Root::IO);
@FORMATS = qw(SCAN FASTA PSA MSA PFF MATCHLIST);

=head2 new

  Title   : new
  Usage   : Bio::Tools::PrositeScan->new(-file => 'out.PrositeScan');
            Bio::Tools::PrositeScan->new(-fh => \*FH);
  Returns : L<Bio::Tools::PrositeScan>
  Args    : -format => string representing the format type for the
                       ps_scan output, REQUIRED

The C<-format> argument must currently be set to C<fasta> since this is the
only parser implemented. This corresponds with using the ps_scan arguments
C<-o fasta>.

=cut

sub new {
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(@args);
    $self->_initialize_io(@args);
    my ($format) = $self->_rearrange([qw(FORMAT)], @args);
    $format || $self->throw("format needed");
    if(grep /^$format$/i, @FORMATS){
        $self->format($format);
    }else{
        $self->throw("Invalid format, [$format]");
    }
    return $self;
}

sub format {
    my $self = shift;
    return $self->{_format} = shift if(@_);
    return $self->{_format};
}

=head2 next_prediction

  Title   : new
  Usage   :
      while($result = $factory->next_prediction){
          ;
      }

  Returns : a Bio::SeqFeature::FeaturePair object where
            feature1 is the matched subsequence and
            feature2 is the PROSITE accession number.
            See <http://prosite.expasy.org/prosuser.html#conv_ac>.

=cut

sub next_prediction {
    my ($self) = @_;
    unless($self->_parsed){
        $self->_parse;
        $self->_parsed(1);
    }
    return shift @{$self->{_matches}};
}

sub next_result {
    return shift->next_prediction;
}

sub _parsed {
    my $self = shift;
    return $self->{_parsed} = 1 if @_ && $_[0];
    return $self->{_parsed};
}

sub _parse {
    my $self = shift;
    my $format = $self->format;
    if($self->format =~ /^fasta$/){
        $self->_parse_fasta;
    }else{
        $self->throw("the [$format] parser has not been written");
    }
}

sub _parse_fasta {
    my ($self) = @_;
    my @matches;
    my $fp;
    my $seq;
    while(defined($_ = $self->_readline)){
        chop;
        if(/^\>([^>]+)/){
            my $fasta_head = $1;
            if($fasta_head =~ /([^\/]+)\/(\d+)\-(\d+)(\s+)\:(\s+)(\S+)/){
                my $q_id = $1;
                my $q_start = $2;
                my $q_end = $3;
                my $h_id = $6;
                if(defined $fp){
                    $self->_attach_seq($seq, $fp);
                    push @matches, $fp;
                }
                $fp = Bio::SeqFeature::FeaturePair->new(
                    -feature1 => Bio::SeqFeature::Generic->new(
                        -seq_id => $q_id,
                        -start => $q_start,
                        -end => $q_end
                    ),
                    -feature2 => Bio::SeqFeature::Generic->new(
                        -seq_id => $h_id,
                        -start => 0,
                        -end => 0
                    )
                );
                $seq = '';
            }else{
                $self->throw("ERR:\t\[$_\]");
            }
        }else{ # sequence lines, ignored
            $seq .= $_;
        }
    }
    if(defined $fp){
        $self->_attach_seq($seq, $fp);
        push @matches, $fp;
    }
    push @{$self->{_matches}}, @matches;
}

sub _attach_seq {
    my ($self, $seq, $fp) = @_;
    if(defined $fp){
        my $whole_seq = 'X' x ($fp->start-1);
        $whole_seq .= $seq;
        $fp->feature1->attach_seq(
            Bio::Seq->new(-seq => $whole_seq)
        );
    }
}

1;
