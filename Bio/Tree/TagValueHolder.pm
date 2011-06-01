package Bio::Tree::TagValueHolder;
use strict;

sub set_tag_value{
    my ($self,$tag,@values) = @_;
    if( ! defined $tag || ! scalar @values ) {
	$self->warn("cannot call set_tag_value with an undefined value");
    }
    $self->remove_tag ($tag);
    map { push @{$self->{'_tags'}->{$tag}}, $_ } @values;
    return scalar @{$self->{'_tags'}->{$tag}};
}

sub add_tag_value{
    my ($self,$tag,$value) = @_;
    if( ! defined $tag || ! defined $value ) {
	$self->warn("cannot call add_tag_value with an undefined value".($tag ? " ($tag)" : ''));
	$self->warn($self->stack_trace_dump,"\n");
    }
    push @{$self->{'_tags'}->{$tag}}, $value;
    return scalar @{$self->{'_tags'}->{$tag}};
}

sub remove_tag {
   my ($self,$tag) = @_;
   if( exists $self->{'_tags'}->{$tag} ) {
       $self->{'_tags'}->{$tag} = undef;
       delete $self->{'_tags'}->{$tag};
       return 1;
   }
   return 0;
}

sub remove_all_tags{
   my ($self) = @_;
   $self->{'_tags'} = {};
   return;
}

sub get_all_tags{
   my ($self) = @_;
   my @tags = sort keys %{$self->{'_tags'} || {}};
   return @tags;
}

sub get_tag_values{
   my ($self,$tag) = @_;
   return wantarray ? @{$self->{'_tags'}->{$tag} || []} :
                     (@{$self->{'_tags'}->{$tag} || []})[0];
}

sub get_tag_value {
    my ($self,$tag) = @_;
    return (@{$self->{'_tags'}->{$tag} || []})[0];
}

sub has_tag {
   my ($self,$tag) = @_;
   return exists $self->{'_tags'}->{$tag};
}

sub get_tagvalue_hash {
  my $self = shift;
  my $tags = $self->{'_tags'};
  
  my $new_tags;
  foreach my $tag (keys %$tags) {
    $new_tags->{$tag} = $tags->{$tag}[0];
  }
  return $new_tags;
}

1;
