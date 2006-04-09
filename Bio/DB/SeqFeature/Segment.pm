package Bio::DB::SeqFeature::Segment;

use strict;

use base 'Bio::SeqFeature::CollectionI';
use Bio::DB::GFF::Util::Rearrange;
use overload '""' => \&as_string;

###
# new()
#
# Call as Bio::DB::SeqFeature::Segment->new($seqfeature,$store)
#
# or
# Bio::DB::SeqFeature::Segment->new(-seqid=>$seqid,-start=>$start,-end=>$end,-strand=>$strand,-store=>$store)
#
sub new {
  my $class = shift;
  my ($store,$seqid,$start,$end,$strand) = @_;
  return bless {
		store => $store,
		seqid => $seqid,
		start => $start,
		end   => $end,
		strand => $strand,
	       },ref($class) || $class;
}

sub features {
  my $self = shift;
  my @args = @_ == 1 ? (-type=>shift) : @_;
  $self->{store}->features(@args,-seqid=>$self->{seqid},-start=>$self->{start},-end=>$self->{end});
}

sub get_seq_stream {
  my $self = shift;
  $self->{store}->get_seq_stream(@_,-seqid=>$self->{seqid},-start=>$self->{start},-end=>$self->{end});
}

sub get_feature_stream { shift->get_seq_stream(@_) }

sub start  { shift->{start}  }
sub end    { shift->{end}    }
sub seq_id { shift->{seqid}  }
sub strand { shift->{strand} }
sub ref    { shift->seq_id   }
sub length { my $self = shift; return $self->end-$self->start+1; }
sub factory { shift->{store} }
sub store   { shift->{store} }
sub type    { shift->primary_tag }
sub primary_tag    { 'Segment' }
sub display_name { shift->as_string }
sub abs_ref      { shift->ref}
sub abs_start    { shift->start}
sub abs_end      { shift->end}

sub as_string {
  my $self = shift;
  my $label = $self->seq_id;
  my $start = $self->start || '';
  my $end   = $self->end   || '';
  return "$label:$start..$end";
}


1;
