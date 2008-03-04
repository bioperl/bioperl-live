package Bio::DB::SeqFeature::Store::LoadHelper;

use strict;
use DB_File;
use File::Temp 'tempdir';
use File::Spec;
use Fcntl qw(O_CREAT O_RDWR);

sub new {
    my $class   = shift;
    my $tmpdir  = shift;

    my @tmpargs = $tmpdir ? (DIR=>$tmpdir) : ();
    my $tmppath = tempdir(@tmpargs,CLEANUP=>1);
    my $self    = $class->create_dbs($tmppath);
    return bless $self,$class;
}

sub create_dbs {
    my $self = shift;
    my $tmp  = shift;
    my %self;

    # Each of these hashes allow only unique keys
    for my $dbname qw(IndexIt TopLevel Local2Global) {
	my %h;
	tie(%h,'DB_File',File::Spec->catfile($tmp,$dbname),
	    O_CREAT|O_RDWR,0666,$DB_BTREE);
	$self{$dbname} = \%h;
    }

    # The Parent2Child hash allows duplicate keys, so we
    # create it with the R_DUP flag.
    my $btree_dups = DB_File::BTREEINFO->new();
    $btree_dups->{flags} = R_DUP;
    my %h;
    tie(%h,'DB_File',File::Spec->catfile($tmp,'Parent2Child'),
	O_CREAT|O_RDWR,0666,$btree_dups);
    $self{Parent2Child} = \%h;

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

1;
