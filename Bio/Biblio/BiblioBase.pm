# $Id$
#
# BioPerl module for Bio::Biblio::BiblioBase
#
# Cared for by Martin Senger <senger@ebi.ac.uk>
# For copyright and disclaimer see below.

# POD documentation - main docs before the code

=head1 NAME

Bio::Biblio::BiblioBase - An abstract base for other biblio classes

=head1 SYNOPSIS

 # to be written

=head1 DESCRIPTION

 # to be written


=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

  bioperl-l@bioperl.org              - General discussion
  http://bioperl.org/MailList.shtml  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution. Bug reports can be submitted via
email or the web:

  bioperl-bugs@bioperl.org
  http://bioperl.org/bioperl-bugs/

=head1 AUTHORS

Heikki Lehvaslaiho (heikki@ebi.ac.uk)
Martin Senger (senger@ebi.ac.uk)

=head1 COPYRIGHT

Copyright (c) 2002 European Bioinformatics Institute. All Rights Reserved.

This module is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.

=head1 DISCLAIMER

This software is provided "as is" without warranty of any kind.

=cut


# Let the code begin...


package Bio::Biblio::BiblioBase;
use strict;
use vars qw(@ISA $AUTOLOAD);

use Bio::Root::Root;

@ISA = qw(Bio::Root::Root);

sub _abstractDeath {
  my $self = shift;
  my $package = ref $self;
  my $caller = (caller(1))[3];

  $self->throw ("Abstract method '$caller' should not been called.\n" .
		"Not your fault - author of $package should be blamed!");

}

# these methods should not be called here;
# they should be implemented by a subclass
sub _accessible { shift->_abstractDeath; }
sub _attr_type { shift->_abstractDeath; }


#
# deal with 'set_' and 'get_' methods
#
sub AUTOLOAD {
    my ($self, $newval) = @_;

    if ($AUTOLOAD =~ /.*::(\w+)/ && $self->_accessible ("_$1")) {
	my $attr_name = "_$1";
	my $attr_type = $self->_attr_type ($attr_name);
	my $ref_sub =
	    sub {
		my ($this, $new_value) = @_;
		return $this->{$attr_name} unless defined $new_value;

		# here we continue with 'set' method
		my ($newval_type) = ref ($new_value) || 'string';
		my ($expected_type) = $attr_type || 'string';
#		$this->throw ("In method $AUTOLOAD, trying to set a value of type '$newval_type' but '$expected_type' is expected.")
		$this->throw ($this->_wrong_type_msg ($newval_type, $expected_type, $AUTOLOAD))
		    unless ($newval_type eq $expected_type) or
		      UNIVERSAL::isa ($new_value, $expected_type);
                       
		$this->{$attr_name} = $new_value;
		return $new_value;
	    };

        no strict 'refs'; 
        *{$AUTOLOAD} = $ref_sub;
        use strict 'refs'; 
        return $ref_sub->($self, $newval);
    }

    $self->throw ("No such method: $AUTOLOAD");
}

# 

sub new {
    my ($caller, @args) = @_;
    my $class = ref ($caller) || $caller;

    # create and bless a new instance    
    my ($self) = $class->SUPER::new (@args);	

    # make a hashtable from @args
    my %param = @args;
    @param { map { lc $_ } keys %param } = values %param; # lowercase keys

    # set all @args into this object with 'set' values;
    # change '-key' into '_key', and making keys lowercase
    my $new_key;
    foreach my $key (keys %param) {
	($new_key = $key) =~ s/-/_/og;   # change it everywhere, why not
        my $method = lc (substr ($new_key, 1));   # omitting the first '_'
        no strict 'refs'; 
        $method->($self, $param { $key });
    }

    # done
    return $self;
}

#
# set methods test whether incoming value is of a correct type;
# here we return message explaining it
#
sub _wrong_type_msg {
    my ($self, $given_type, $expected_type, $method) = @_;
    my $msg = 'In method ';
    if (defined $method) {
	$msg .= $method;
    } else {
	$msg .= (caller(1))[3];
    }
    return ("$msg: Trying to set a value of type '$given_type' but '$expected_type' is expected.");
}

#
# probably just for debugging
# TBD: to decide...
#
sub print_me {
    my ($self) = @_;
    use Data::Dumper;
    return Data::Dumper->Dump ( [$self], ['Citation']);
}

1;
__END__
