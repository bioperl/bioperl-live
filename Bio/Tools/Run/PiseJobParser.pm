# $Id$
#

=head1 NAME

Bio::Tools::Run::PiseJobParser

=head1 SYNOPSIS

  #

=head1 DESCRIPTION

   Parsing of Pise XHTML output to extract results files and piping menus.

=cut

#'

package Bio::Tools::Run::PiseJobParser;

use vars qw(@ISA);
use strict;
use Bio::Root::Root;

@ISA = qw(Bio::Root::Root);

sub new {
    my ($class, $verbose) = @_;
    my $self = $class->SUPER::new();
    if ($verbose) {
	$self->{VERBOSE} = $verbose;
    } else {
	$self->{VERBOSE} = 0;
    }

    $self->{DEBUG} = 0;

    return $self;
}

sub characters {
    my ($self, $element) = @_; 
    chomp ($element->{Data});
    print STDERR $element->{Data}  if ($self->{DEBUG});
    if ($element->{Data} =~ /Results not available yet/) {
	$self->{terminated} = 0;
    }
    if ($element->{Data} =~ /Results:/) {
	$self->{output_files} = 1;
	$self->{terminated} = 1;
    }
    if ($element->{Data} =~ /this files will remain accessible for/) {
	$self->{output_files} = 0;
	$self->{result_url} = 1;
    }
    if ($element->{Data} =~ /You can save them individually/) {
	$self->{result_url} = 0;
    }
    if ($element->{Data} =~ /upon completion of the job/) {
	print STDERR "job detached\n" if $self->{DEBUG};
	$self->{result_url} = 1;
    }
    if ($element->{Data} =~ /Please wait for the results of this query before submitting/) {
	$self->{result_url} = 0;
    }
    if ($self->{check_message}) {
	$self->{error_message} .= $element->{Data};
    }
}

sub comment {
    my ($self, $element) = @_; 
    if ($element->{Data} =~ /USER ERROR/) {
	$self->{error} = 1;
    }
}

sub start_element {
    my ($self, $element) = @_;
    print STDERR "\nstart element: ",$element->{Name},"\n" if ($self->{DEBUG});

    my %attributes = %{ $element->{Attributes} };
    foreach my $attr (keys %attributes) {
	print STDERR "\t$attr $attributes{$attr}\n"  if ($self->{DEBUG});
    }
    if ($element->{Name} eq "HTML") {
	#$self->{terminated} = 1;
	$self->{terminated} = 0;
	$self->{error} = 0;
	$self->{error_message} = "";
	$self->{hrefs} = [];
    } elsif ($element->{Name} eq "A") {
	$self->{href} = $attributes{HREF};
	# so nothing could work for Pise installation where Pise is in the
	# url... :-(
	#if ($PiseJobParser::href !~ /Pise/ && $PiseJobParser::href ne "" && $PiseJobParser::output_files) {
	if ($self->{href} ne "") {
	    if ($self->{output_files}) {
		push (@{$self->{hrefs}}, $self->{href} );
		print STDERR "Bio::Tools::Run::PiseJobParser: href=",$self->{href} ,"\n" if ($self->{DEBUG});
	    } elsif ($self->{result_url}) {
		if (! $self->{bioweb_result}) {
		    $self->{bioweb_result} = $self->{href};
		}
	    }
	}
    } elsif ($element->{Name} eq "H3") {
	if ($self->{error}) {
	    $self->{check_message} = 1;
	}
    } elsif ($element->{Name} eq "FORM") {
	my $action = $attributes{action};
	if ($action =~ /connect.pl/) {
	    $self->{connected}{$self->{href}} = 1;
	    print STDERR "\t",$self->{href}, " is connected to...\n" if ($self->{DEBUG});
	} elsif ($action =~ /results.pl/) {
	    $self->{terminated} = 0;
	}
    } elsif ($element->{Name} eq "INPUT") {
	my $name=$attributes{NAME};
	my $value=$attributes{VALUE};
	if ($name eq "scratch_dir") {
	    $self->{scratch_dir} = $value;
	} else {
	    $self->{value}{$name} = $value;
	}
	if ($self->{connected}{$self->{href}}) {
	    if ($name eq "piped_file_type") {
		print STDERR "DEBUG> ",$self->{href}," = $value\n"  if ($self->{DEBUG});
		$self->{piped_file_type}{$self->{href}} = $value;
	    }
	}
    } elsif ($element->{Name} eq "OPTION") {
	my $option=1;
	my $value=$attributes{VALUE};
	my $command;
	my @with_piped_files;
	my $with_href;
	my $root_url;
	my $with_param;
	my $with_value;
	($command,@with_piped_files) = split(",",$value);
	if ($self->{connected}{$self->{href}}) {
	    push (@{$self->{pipes}{$self->{href}}}, $command);

	    #print STDERR "pipes:\n";
	    #foreach my $f (keys %{$self->{pipes}}) {
		#my @p = @{ $self->{pipes}{$f} };
		#foreach my $p (@p) {
		#    print STDERR "\tf: $f\tp: $p\n";
		#}
	    #}
	    
	    ($root_url = $self->{href}) =~ s/(.+)\/.+/$1/;
	    foreach my $with_file (@with_piped_files) {
		($with_param,$with_value) = split("=",$with_file);	
		my $with_href = "$root_url/$with_value" ;
		push @{ $self->{with_href}{$self->{href}} }, $with_href;
		if ( ! (grep {$command eq $_ } @{$self->{pipes}{$with_href}}) ) {
		    push (@{$self->{pipes}{$with_href}},$command);
		    $self->{piped_file_type}{$with_href} = $with_param;
		}
	    }
	    print STDERR "\t\t",$self->{href}," is connected to $value\n"  if ($self->{DEBUG});
	}
    } 

}

sub end_element {
    my ($self, $element) = @_;
    print STDERR "\nend element: ",$element->{Name},"\n" if ($self->{DEBUG});
    if ($element->{Name} eq "H3") {
	if ($self->{error}) {
	    $self->{check_message}=0;
	}
    }
}

sub pipes {
    my $self = shift;
    if (defined $self->{pipes}) {
	return %{$self->{pipes}};
    }
}

sub piped_file_type {
    my $self = shift;
    my $href = shift;
    return $self->{piped_file_type}{$href};
}

sub bioweb_result {
    my $self = shift;
    return $self->{bioweb_result};
}

sub scratch_dir {
    my $self = shift;
    return $self->{scratch_dir};
}

sub hrefs {
    my $self = shift;
    return @{ $self->{hrefs} };
}

sub terminated {
    my $self = shift;
    return $self->{terminated};
}

sub error {
    my $self = shift;
    return $self->{error};
}

sub error_message {
    my $self = shift;
    return $self->{error_message};
}

1;
