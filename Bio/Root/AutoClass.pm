package Bio::Root::AutoClass;
use strict;

sub new {
  my($class,@args)=@_;
  $class=(ref $class)||$class;
  my @isa=$class->ISA;
  my $temp=bless({},$class);
  # Find a base class who really wants to create the object, if any
  my $new_method=$temp->_new_method($class);
  # Create the object
  # Rebless what comes from new since some BioPerl classes don't bless into subclasses
  my $self=$new_method? bless(&$new_method($class,@args),$class): bless({},$class);
  my $args=_fix_args(@args);
  my $package=__PACKAGE__;
  $self->{'_init'}{$package}++;	# Don't init AutoClass again
  $self->_init($class,$args);
  $self;
}

sub _new_method {
  my($self,$class)=@_;
  return if ($self->{'_new'}{$class});
  $self->{'_new'}{$class}++;
  my @isa=$class->ISA;
  my $new_method;
  for my $super (@isa) {
    next if $super eq __PACKAGE__;
    if (UNIVERSAL::isa($super,__PACKAGE__)) { # he's one of ours, so recurse
      $new_method=$self->_new_method($super);
#      my $new_super=$super->can('_new_method');
#      $new_method=$self->$new_super($super);
    } else {			         # an external class -- can he new?
      $new_method=$super->can('new');
    }
    last if $new_method;
  }
  $new_method;
}
     
sub _init {
  my($self,$class,$args)=@_;
  return if ($self->{'_init'}{$class});
  $self->{'_init'}{$class}++;
  my @isa=$class->ISA;
  for my $super (@isa) {
    next if $super eq __PACKAGE__;
    my $init_super=$super->can('_init');
    $self->$init_super($super,$args) if $init_super && !($self->{'_init'}{$super});
  }
  my %synonyms=SYNONYMS($class);
  my $attributes=[AUTO_ATTRIBUTES($class),OTHER_ATTRIBUTES($class),keys %synonyms];
  $self->_set_attributes($attributes,$args);
  my $init_self=$class->can('_init_self');
  $self->$init_self($class,$args) if $init_self;
}
sub set {
  my($self,@args)=@_;
  my $args=_fix_args(@args);
  while(my($key,$value)=each %$args) {
    my $func=$self->can(lc($key));
    $self->$func($value) if $func;
  }
}
sub get {
  my($self,@args)=@_;
  my @keys=_fix_keyword(@args);
  my @results;
  for my $key (@keys) {
    my $func=$self->can(lc($key));
    $func? push(@results,$self->$func()): undef;
  }
  wantarray? @results: $results[0];
}

sub class {ref $_[0];}
sub ISA {
  my($class)=@_;
  $class=$class->class if ref $class; # get class if called as object method
  no strict 'refs';
  @{$class.'::ISA'}
}
sub AUTO_ATTRIBUTES {
  my($class)=@_;
  $class=$class->class if ref $class; # get class if called as object method
  no strict 'refs';
  @{$class.'::AUTO_ATTRIBUTES'}
}
sub OTHER_ATTRIBUTES {
  my($class)=@_;
  $class=$class->class if ref $class; # get class if called as object method
  no strict 'refs';
  @{$class.'::OTHER_ATTRIBUTES'}
}
sub SYNONYMS {
  my($class)=@_;
  $class=$class->class if ref $class; # get class if called as object method
  no strict 'refs';
  %{$class.'::SYNONYMS'}
}

sub declare {
  my($class,$attributes,$synonyms,$case)=@_;
  for my $func (@$attributes) {
    my $keyword=_fix_keyword($func);
    my $sub='*'.$class.'::'.$func."=sub {\@_>1? \$_[0]->{\'$keyword\'}=\$_[1]: \$_[0]->{\'$keyword\'};}";
    eval $sub;
  }
  while(my($func,$old_func)=each %$synonyms) {
    next if $func eq $old_func;	# avoid infinite recursion if old and new are the same
    my $sub='*'.$class.'::'.$func."=sub {\$_[0]->$old_func(\@_[1..\$\#_])}";
    eval $sub;
  }
  if ($case=~/lower|lc/i) {	# create lowercase versions of each method, too
    for my $func (@$attributes) {
      my $lc_func=lc $func;
      next if $lc_func eq $func; # avoid infinite recursion if func already lowercase
      my $sub='*'.$class.'::'.$lc_func."=sub {\$_[0]->$func(\@_[1..\$\#_])}";
      eval $sub;
    }
  }
  if ($case=~/upper|uc/i) {	# create uppercase versions of each method, too
    for my $func (@$attributes) {
      my $uc_func=uc $func;
      next if $uc_func eq $func; # avoid infinite recursion if func already uppercase
      my $sub='*'.$class.'::'.$uc_func."=sub {\$_[0]->$func(\@_[1..\$\#_])}";
      eval $sub;
    }
  }
}

sub _fix_args {
  my(@args)=@_;
  @args=@{$args[0]} if @args==1 && 'ARRAY' eq ref $args[0];
  @args=%{$args[0]} if @args==1 && 'HASH' eq ref $args[0];
  @args=$args[0]->_outcast if @args==1 && UNIVERSAL::isa($args[0],__PACKAGE__);
  @args=%{$args[0]} if @args==1 && $args[0]=~/HASH/;
  die("Malformed keyword argument list (odd number of elements): @args") if @args%2;
  my(%args,%counts);
  while(@args) {
    my($key,$value)=(_fix_keyword(shift @args),shift @args);
    $args{$key}=$value if $counts{$key}==0;
    $args{$key}=[$args{$key},$value] if $counts{$key}==1;
    push(@{$args{$key}},$value) if $counts{$key}>1;
    $counts{$key}++;
  }
  \%args;
}
sub _fix_keyword {
  my @keywords=@_;		# copies input, so update-in-place doesn't munge it
  for my $keyword (@keywords) {
    next unless defined $keyword;
    $keyword=~s/^-*(.*)$/\U\1/ unless ref $keyword; # updates in place
  }
  wantarray? @keywords: $keywords[0];
}
sub _fix_keywords {_fix_keyword(@_);}

sub _set_attributes {
  my($self,$attributes,$args)=@_;
  for my $func (@$attributes) {
    my $keyword=_fix_keyword($func);
    next unless exists $args->{$keyword};
    $self->$func($args->{$keyword});
  }
}
sub _outcast {
  my($self,$target)=@_;
  %$self;
}
sub _is_positional {
  @_%2 || $_[0]!~/^-/;
}

1;
