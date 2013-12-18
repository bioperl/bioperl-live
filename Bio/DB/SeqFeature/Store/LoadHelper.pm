package Bio::DB::SeqFeature::Store::LoadHelper;

=head1 NAME

Bio::DB::SeqFeature::Store::LoadHelper -- Internal utility for Bio::DB::SeqFeature::Store

=head1 SYNOPSIS

  # For internal use only.

=head1 DESCRIPTION

For internal use only

=head1 SEE ALSO

L<bioperl>,
L<Bio::DB::SeqFeature::Store>,
L<Bio::DB::SeqFeature::Segment>,
L<Bio::DB::SeqFeature::NormalizedFeature>,
L<Bio::DB::SeqFeature::GFF2Loader>,
L<Bio::DB::SeqFeature::Store::DBI::mysql>,
L<Bio::DB::SeqFeature::Store::berkeleydb>

=head1 AUTHOR

Lincoln Stein E<lt>lstein@cshl.orgE<gt>.

Copyright (c) 2006 Cold Spring Harbor Laboratory.

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.

=cut

use strict;
use DB_File;
use File::Path 'rmtree';
use File::Temp 'tempdir';
use File::Spec;
use Fcntl qw(O_CREAT O_RDWR);

our $VERSION = '1.10';

my %DBHandles;

sub new {
    my $class   = shift;
    my $tmpdir  = shift;

    my $template = 'SeqFeatureLoadHelper_XXXXXX';

    my @tmpargs = $tmpdir ? ($template,DIR=>$tmpdir) : ($template);
    my $tmppath = tempdir(@tmpargs,CLEANUP=>1);
    my $self    = $class->create_dbs($tmppath);
    $self->{tmppath} = $tmppath;
    return bless $self,$class;
}

sub DESTROY {
    my $self = shift;
    # Destroy all filehandle references
    # before trying to delete files and folder
    %DBHandles = ();
    undef $self->{IndexIt};
    undef $self->{TopLevel};
    undef $self->{Local2Global};
    undef $self->{Parent2Child};
    rmtree $self->{tmppath};
#    File::Temp::cleanup() unless $self->{keep};
}

sub create_dbs {
    my $self = shift;
    my $tmp  = shift;
    my %self;
    # experiment with caching these handles in memory
    my $hash_options           = DB_File::HASHINFO->new();
    # Each of these hashes allow only unique keys
    for my $dbname (qw(IndexIt TopLevel Local2Global)) {
	unless ($DBHandles{$dbname}) {
	    my %h;
	    tie(%h,'DB_File',File::Spec->catfile($tmp,$dbname),
		O_CREAT|O_RDWR,0666,$hash_options);
	    $DBHandles{$dbname} = \%h;
	}
	$self{$dbname} = $DBHandles{$dbname};
	%{$self{$dbname}} = ();
    }

    # The Parent2Child hash allows duplicate keys, so we
    # create it with the R_DUP flag.
    my $btree_options           = DB_File::BTREEINFO->new();
    $btree_options->{flags}     = R_DUP;
    unless ($DBHandles{'Parent2Child'}) {
	my %h;
	tie(%h,'DB_File',File::Spec->catfile($tmp,'Parent2Child'),
	    O_CREAT|O_RDWR,0666,$btree_options);
	$DBHandles{'Parent2Child'} = \%h;
    }
    $self{Parent2Child}    = $DBHandles{'Parent2Child'};
    %{$self{Parent2Child}} = ();
    return \%self;
}

sub indexit {
    my $self = shift;
    my $id   = shift;
    $self->{IndexIt}{$id} = shift if @_;
    return $self->{IndexIt}{$id};
}

sub toplevel {
    my $self = shift;
    my $id   = shift;
    $self->{TopLevel}{$id} = shift if @_;
    return $self->{TopLevel}{$id};
}

sub each_toplevel {
    my $self = shift;
    my ($id) = each %{$self->{TopLevel}};
    $id;
}

sub local2global {
    my $self = shift;
    my $id   = shift;
    $self->{Local2Global}{$id} = shift if @_;
    return $self->{Local2Global}{$id};
}

sub add_children {
    my $self      = shift;
    my $parent_id = shift;
    # (@children) = @_;
    $self->{Parent2Child}{$parent_id} = shift while @_;
}

sub children {
    my $self = shift;
    my $parent_id = shift;

    my @children;

    my $db        = tied(%{$self->{Parent2Child}});
    my $key       = $parent_id;
    my $value     = '';
    for (my $status = $db->seq($key,$value,R_CURSOR);
	 $status    == 0 && $key eq $parent_id;
	 $status    = $db->seq($key,$value,R_NEXT)
	) {
	push @children,$value;
    }
    return wantarray ? @children: \@children;
}

# this acts like each() and returns each parent id and an array ref of children
sub each_family {
    my $self = shift;

    my $db        = tied(%{$self->{Parent2Child}});

    if ($self->{_cursordone}) {
	undef $self->{_cursordone};
	undef $self->{_parent};
	undef $self->{_child};
	return;
    }

    # do a slightly tricky cursor search
    unless (defined $self->{_parent}) {
	return unless $db->seq($self->{_parent},$self->{_child},R_FIRST) == 0;
    }

    my $parent   = $self->{_parent};
    my @children = $self->{_child};

    my $status;
    while (($status = $db->seq($self->{_parent},$self->{_child},R_NEXT)) == 0
	   && $self->{_parent} eq $parent
	) {
	push @children,$self->{_child};
    }

    $self->{_cursordone}++ if $status != 0;
    
    return ($parent,\@children);
}

sub local_ids {
    my $self = shift;
    my @ids  = keys %{$self->{Local2Global}}
                   if $self->{Local2Global};
    return \@ids;
}

sub loaded_ids {
    my $self = shift;
    my @ids  = values %{$self->{Local2Global}}
                     if $self->{Local2Global};
    return \@ids;
}

1;
