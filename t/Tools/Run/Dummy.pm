package Dummy;
use strict;
use warnings;

use lib '.';
use lib '..';
use Dummy::Config;

use Bio::Tools::Run::WrapperBase;
use Bio::Tools::Run::WrapperBase::CommandExts;

use base qw(Bio::Tools::Run::WrapperBase Bio::Root::Root);

sub new {
    my ($class,@args) = @_;
    my $self = $class->SUPER::new(@args);
    return $self;
}

1;

