#
# BioPerl module for Bio::Biblio::BiblioBase
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Martin Senger <senger@ebi.ac.uk>
# For copyright and disclaimer see below.

# POD documentation - main docs before the code

=head1 NAME

Bio::Biblio::BiblioBase - An abstract base for other biblio classes

=head1 SYNOPSIS

 # do not instantiate this class directly

=head1 DESCRIPTION

It is a base class where all other biblio data storage classes inherit
from. It does not reflect any real-world object, it exists only for
convenience, in order to have a place for shared code.

=head2 new()

The I<new()> class method constructs a new biblio storage object.  It
accepts list of named arguments - the same names as attribute names
prefixed with a minus sign. Available attribute names are listed in
the documentation of the individual biblio storage objects.

=head2 Accessors

All attribute names can be used as method names. When used without any
parameter the method returns current value of the attribute (or
undef), when used with a value the method sets the attribute to this
value and also returns it back. The set method also checks if the type
of the new value is correct.

=head2 Custom classes

If there is a need for new attributes, create your own class which
usually inherits from I<Bio::Biblio::Ref>. For new types of providers
and journals, let your class inherit directly from this
I<Bio::Biblio::BiblioBase> class.

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
of the bugs and their resolution. Bug reports can be submitted via the
web:

  https://redmine.open-bio.org/projects/bioperl/

=head1 AUTHOR

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
use vars qw($AUTOLOAD);


use base qw(Bio::Root::Root);

# these methods should not be called here;
# they should be implemented by a subclass
sub _accessible { shift->throw_not_implemented(); }
sub _attr_type { shift->throw_not_implemented(); }

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
    require Data::Dumper;
    return Data::Dumper->Dump ( [$self], ['Citation']);
}

1;
__END__
