#
# BioPerl module for Bio::OntologyIO::simplehierarchy
#
# Please direct questions and support issues to <bioperl-l@bioperl.org>
#
# Cared for by Allen Day, allenday@ucla.edu
#
# (c) Allen Day, allenday@ucla.edu, 2003.
# (c) Department of Human Genetics, UCLA Medical School, 2003.
#
# You may distribute this module under the same terms as perl itself.
# Refer to the Perl Artistic License (see the license accompanying this
# software package, or see http://www.perl.com/language/misc/Artistic.html)
# for the terms under which you may use, modify, and redistribute this module.
#
# THIS PACKAGE IS PROVIDED "AS IS" AND WITHOUT ANY EXPRESS OR IMPLIED
# WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED WARRANTIES OF
# MERCHANTIBILITY AND FITNESS FOR A PARTICULAR PURPOSE.
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::OntologyIO::simplehierarchy - a base class parser for simple hierarchy-by-indentation
                  type formats

=head1 SYNOPSIS

  use Bio::OntologyIO;

  # do not use directly -- use via Bio::OntologyIO
  my $parser = Bio::OntologyIO->new
        ( -format        => "simplehierarchy",
          -file          => "pathology_terms.csv",
          -indent_string => ",",
          -ontology_name => "eVOC",
          -term_factory  => $fact,
        );

  my $ontology = $parser->next_ontology();

=head1 DESCRIPTION

Needs Graph.pm from CPAN.  This class is nearly identical to
OntologyIO::dagflat, see L<Bio::OntologyIO::dagflat> for details.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to the
Bioperl mailing lists  Your participation is much appreciated.

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
the bugs and their resolution.  Bug reports can be submitted via the
web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR

Allen Day

Email: allenday@ucla.edu

=head2 CONTRIBUTOR

Christian Zmasek

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package  Bio::OntologyIO::simplehierarchy;

use strict;

use Data::Dumper;
use File::Basename;
use Bio::Root::IO;
use Bio::Ontology::OBOEngine;
use Bio::Ontology::Ontology;
use Bio::Ontology::TermFactory;

use constant TRUE         => 1;
use constant FALSE        => 0;


use base qw(Bio::OntologyIO);


=head2 new

 Title   : new
 Usage   : see SYNOPSIS
 Function: Creates a new simplehierarchy parser.
 Returns : A new simplehierarchy parser object, implementing Bio::OntologyIO.
 Args    : -files         => a single ontology flat file holding the
                             term relationships, or an array ref holding
                             the file names
           -file          => if there is only a single flat file, it may
                             also be specified via the -file parameter
           -ontology_name => the name of the ontology, defaults to
                             "Gene Ontology"
           -file_is_root  => Boolean indicating whether a virtual root
                             term is to be added, the name of which will
                             be derived from the file name. Default is false.
                             Enabling this allows one to parse multiple input
                             files into the same ontology and still have
                             separately rooted.
           -engine        => the L<Bio::Ontology::OntologyEngineI> object
                             to be reused (will be created otherwise); note
                             that every L<Bio::Ontology::OntologyI> will
                             qualify as well since that one inherits from the
                             former.
           -indent_string => the string used to indent hierarchical
                             levels in the file.

                             For a file like this:

                             term0
                               subterm1A
                                 subterm2A
                               subterm1B
                               subterm1C

                             indent_string would be "  ".  Defaults to
                             one space (" ").
           -comment_char  => Allows specification of a regular
                             expression string to indicate a comment line.
                             Currently defaults to "[\|\-]".
                             Note: this is not yet implemented.

See L<Bio::OntologyIO>.

=cut

# in reality, we let OntologyIO::new do the instantiation, and override
# _initialize for all initialization work
sub _initialize {
    my ($self, @args) = @_;
    $self->SUPER::_initialize( @args );

    my ( $indent,$files,$fileisroot,$name,$eng ) =
        $self->_rearrange([qw(INDENT_STRING
                              FILES
                              FILE_IS_ROOT
                              ONTOLOGY_NAME
                              ENGINE)
                           ], @args);

    $self->_done( FALSE );
    $self->_not_first_record( FALSE );
    $self->_term( "" );
    $self->file_is_root($fileisroot) if defined($fileisroot);
    $indent = ' ' unless defined($indent); #reasonable default?
    # the indentation string may have escaped chars
    if (($indent =~ /\\/) && ($indent !~ /[\$\`]/)) {
        $indent = "\$indent = \"$indent\"";
        eval $indent;
    }
    $self->indent_string($indent);
    delete $self->{'_ontologies'};

    # ontology engine (and possibly name if it's an OntologyI)
    $eng = Bio::Ontology::OBOEngine->new() unless $eng;
    if($eng->isa("Bio::Ontology::OntologyI")) {
        $self->ontology_name($eng->name());
        $eng = $eng->engine() if $eng->can('engine');
    }
    $self->_ont_engine($eng);

    # flat files to parse
    $self->{_flat_files} = $files ? ref($files) ? $files : [$files] : [];

    # ontology name (overrides implicit one through OntologyI engine)
    $self->ontology_name($name) if $name;

} # _initialize

=head2 ontology_name

 Title   : ontology_name
 Usage   : $obj->ontology_name($newval)
 Function: Get/set the name of the ontology parsed by this module.
 Example :
 Returns : value of ontology_name (a scalar)
 Args    : on set, new value (a scalar or undef, optional)

=cut

sub ontology_name{
    my $self = shift;

    return $self->{'ontology_name'} = shift if @_;
    return $self->{'ontology_name'};
}


=head2 parse

 Title   : parse()
 Usage   : $parser->parse();
 Function: Parses the files set with "new" or with methods
           defs_file and _flat_files.

           Normally you should not need to call this method as it will
           be called automatically upon the first call to
           next_ontology().

 Returns : [Bio::Ontology::OntologyEngineI]
 Args    :

=cut

sub parse {
    my $self = shift;

    # setup the default term factory if not done by anyone yet
    $self->term_factory(Bio::Ontology::TermFactory->new(
                                             -type => "Bio::Ontology::Term"))
        unless $self->term_factory();

    # create the ontology object itself
    my $ont = Bio::Ontology::Ontology->new(-name => $self->ontology_name(),
                                           -engine => $self->_ont_engine());

    # set up the ontology of the relationship types
    foreach ($self->_part_of_relationship(),
             $self->_is_a_relationship(),
             $self->_related_to_relationship()) {
        $_->ontology($ont);
    }

    # pre-seed the IO system with the first flat file if -file wasn't provided
    if(! $self->_fh) {
        $self->_initialize_io(-file => shift(@{$self->_flat_files()}));
    }

    while($self->_fh) {
          $self->_parse_flat_file($ont);
          # advance to next flat file if more are available
          if(@{$self->_flat_files()}) {
            $self->close();
            # reset the virtual root so that the next one is generated from
            # the next file
            $self->_virtual_root(undef);
            # now re-initialize the IO object
            $self->_initialize_io(-file => shift(@{$self->_flat_files()}));
          } else {
            last; # nothing else to parse so terminate the loop
          }
    }
    $self->_add_ontology($ont);
    # not needed anywhere, only because of backward compatibility
    return $self->_ont_engine();
} # parse

=head2 next_ontology

 Title   : next_ontology
 Usage   :
 Function: Get the next available ontology from the parser. This is the
           method prescribed by Bio::OntologyIO.
 Example :
 Returns : An object implementing Bio::Ontology::OntologyI, and undef if
           there is no more ontology in the input.
 Args    :

=cut

sub next_ontology{
    my $self = shift;

    # parse if not done already
    $self->parse() unless exists($self->{'_ontologies'});
    # return next available ontology
    return shift(@{$self->{'_ontologies'}}) if exists($self->{'_ontologies'});
    return;
}

=head2 _flat_files

 Title   : _flat_files
 Usage   : $files_to_parse = $parser->_flat_files();
 Function: Get the array of ontology flat files that need to be parsed.

           Note that this array will decrease in elements over the
           parsing process. Therefore, it\'s value outside of this
           module will be limited. Also, be careful not to alter the
           array unless you know what you are doing.

 Returns : a reference to an array of zero or more strings
 Args    : none

=cut

sub _flat_files {
    my $self = shift;

    $self->{_flat_files} = [] unless exists($self->{_flat_files});
    return $self->{_flat_files};
}


# INTERNAL METHODS
# ----------------

=head2 _defs_io

 Title   : _defs_io
 Usage   : $obj->_defs_io($newval)
 Function: Get/set the Bio::Root::IO instance representing the
           definition file, if provided (see defs_file()).
 Example :
 Returns : value of _defs_io (a Bio::Root::IO object)
 Args    : on set, new value (a Bio::Root::IO object or undef, optional)

=cut

sub _defs_io{
    my $self = shift;

    return $self->{'_defs_io'} = shift if @_;
    return $self->{'_defs_io'};
}

sub _add_ontology {
    my $self = shift;
    $self->{'_ontologies'} = [] unless exists($self->{'_ontologies'});
    foreach my $ont (@_) {
        $self->throw(ref($ont)." does not implement Bio::Ontology::OntologyI")
            unless ref($ont) && $ont->isa("Bio::Ontology::OntologyI");
        push(@{$self->{'_ontologies'}}, $ont);
    }
}

# This simply delegates. See SimpleGOEngine.
sub _add_term {
    my ( $self, $term, $ont ) = @_;

    $term->ontology($ont) if $ont && (! $term->ontology);
    $self->_ont_engine()->add_term( $term );


} # _add_term

# This simply delegates. See SimpleGOEngine
sub _part_of_relationship {
    my ( $self, $term ) = @_;

    return $self->_ont_engine()->part_of_relationship();

} # _part_of_relationship



# This simply delegates. See SimpleGOEngine
sub _is_a_relationship {
    my ( $self, $term ) = @_;

    return $self->_ont_engine()->is_a_relationship();

} # _is_a_relationship


# This simply delegates. See SimpleGOEngine
sub _related_to_relationship {
    my ( $self, $term ) = @_;

    return $self->_ont_engine()->related_to_relationship();

} # _is_a_relationship


# This simply delegates. See SimpleGOEngine
sub _add_relationship {
    my ( $self, $parent, $child, $type, $ont ) = @_;

    # note the triple terminology (subject,predicate,object) corresponds to
    # (child,type,parent)
    $self->_ont_engine()->add_relationship( $child, $type, $parent, $ont );


} # _add_relationship


# This simply delegates. See SimpleGOEngine
sub _has_term {
    my ( $self, $term ) = @_;

    $term = $self->ontology_name() .'|'. $term
        unless ref($term) || !$self->ontology_name();
    return $self->_ont_engine()->has_term( $term );

} # _add_term

# This simply delegates after prefixing the namespace name if it is just a
# base identifier. See SimpleGOEngine
sub _get_terms{
    my $self = shift;
    my @args = ();

    while(@_) {
        unshift(@args, pop(@_)); # this actually does preserve the order
        $args[0] = $self->ontology_name() .'|'. $args[0]
            unless ref($args[0]) || !$self->ontology_name();
    }
    return $self->_ont_engine->get_terms(@args);
}


# This parses the relationships files
sub _parse_flat_file {
  my $self = shift;
  my $ont  = shift;

  my @stack       = ();
  my $prev_indent = -1;
  my $parent      = "";
  my $prev_term   = "";

  my $indent_string = $self->indent_string;


  while ( my $line = $self->_readline() ) {
        if ( $line =~ /^[$indent_string]*[\|\-]/ ) { #this is not yet generalized
          next;
        }

        my ($current_term) = $line =~ /^[$indent_string]*(.*)/;
        my $current_indent = $self->_count_indents( $line );
        chomp $current_term;
        # remove extraneous delimiter characters at the end of the name if any
        $current_term =~ s/[$indent_string]+$//;
        # remove double quotes surrounding the entry, if any
        $current_term =~ s/^\"(.*)\"$/$1/;
        # also, the name might contain a synonym
        my $syn = $current_term =~ s/\s+{([^}]+)}// ? $1 : undef;

         if ( ! $self->_has_term( $current_term ) ) {
           my $term = $self->_create_ont_entry($current_term);
          # add synonym(s) if any
          $term->add_synonym(split(/[;,]\s*/,$syn)) if $syn;
          # add to the machine
           $self->_add_term( $term, $ont );

          #go on to the next term if a root node.
          if($current_indent == 0) {
              # add the virtual root as parent if there is one
              if($self->_virtual_root()) {
                  $self->_add_relationship($self->_virtual_root(),
                                           $term,
                                           $self->_is_a_relationship(),
                                           $ont);
              }
              $prev_indent = $current_indent;
              $prev_term = $current_term;
              push @stack, $current_term;
              next;
          }
         }

        # note: we are ensured to see the parent first in this type of file,
        # so we never need to possibly insert the parent here

         if ( $current_indent != $prev_indent  ) {
           if ( $current_indent == $prev_indent + 1 ) {
                 push( @stack, $prev_term );
           } elsif ( $current_indent < $prev_indent ) {
                 my $n = $prev_indent -  $current_indent;
                 for ( my $i = 0; $i < $n; ++$i ) {
                   pop( @stack );
                 }
           } else {
                 $self->throw("format error: indentation level $current_indent "
                             ."is more than one higher than the previous "
                             ."level $prev_indent ('$current_term', "
                             ."file ".$self->file.")" );
           }
         }

         $parent = $stack[-1];

        if($parent ne $current_term) { #this prevents infinite recursion from a parent linking to itself
          $self->_add_relationship($self->_get_terms($parent),
                                   $self->_get_terms($current_term),
                                   $self->_is_a_relationship(),
                                   $ont);
        }

        $prev_indent = $current_indent;
        $prev_term   = $current_term;
  }
  return $ont;
} # _parse_relationships_file



# Parses the 1st term id number out of line.
sub _get_first_termid {
    my ( $self, $line ) = @_;

    if ( $line =~ /;\s*([A-Z]{1,8}:\d{7})/ ) {
        return $1;
    }
    else {
        $self->throw( "format error: no term id in line \"$line\"" );
    }

} # _get_first_termid

# Counts the indents at the beginning of a line in the relationships files
sub _count_indents {
  my ( $self, $line ) = @_;

  my $indent = $self->indent_string;

  if ( $line =~ /^($indent+)/ ) {
        return (length($1)/length($indent));
  }
  else {
        return 0;
  }
} # _count_indents


# Holds the GO engine to be parsed into
sub _ont_engine {
    my ( $self, $value ) = @_;

    if ( defined $value ) {
        $self->{ "_ont_engine" } = $value;
    }

    return $self->{ "_ont_engine" };
} # _ont_engine


# Used to create ontology terms.
# Arguments: name, id
sub _create_ont_entry {
    my ( $self, $name, $termid ) = @_;

    my $term = $self->term_factory->create_object(-name => $name,
                                                  -identifier => $termid);

    return $term;

} # _create_ont_entry


# Holds whether first record or not
sub _not_first_record {
    my ( $self, $value ) = @_;

    if ( defined $value ) {
        unless ( $value == FALSE || $value == TRUE ) {
            $self->throw( "Argument to method \"_not_first_record\" must be either ".TRUE." or ".FALSE );
        }
        $self->{ "_not_first_record" } = $value;
    }

    return $self->{ "_not_first_record" };
} # _not_first_record



# Holds whether done or not
sub _done {
    my ( $self, $value ) = @_;

    if ( defined $value ) {
        unless ( $value == FALSE || $value == TRUE ) {
            $self->throw( "Found [$value] where [" . TRUE
            ." or " . FALSE . "] expected" );
        }
        $self->{ "_done" } = $value;
    }

    return $self->{ "_done" };
} # _done


# Holds a term.
sub _term {
    my ( $self, $value ) = @_;

    if ( defined $value ) {
        $self->{ "_term" } = $value;
    }

    return $self->{ "_term" };
} # _term

=head2 indent_string

 Title   : indent_string
 Usage   : $obj->indent_string($newval)
 Function:
 Example :
 Returns : value of indent_string (a scalar)
 Args    : on set, new value (a scalar or undef, optional)

=cut

sub indent_string{
    my $self = shift;

    return $self->{'indent_string'} = shift if @_;
    return $self->{'indent_string'};
}

=head2 file_is_root

 Title   : file_is_root
 Usage   : $obj->file_is_root($newval)
 Function: Boolean indicating whether a virtual root term is to be
           added, the name of which will be derived from the file
           name.

           Enabling this allows one to parse multiple input files into the
           same ontology and still have separately rooted.

 Example :
 Returns : value of file_is_root (a scalar)
 Args    : on set, new value (a scalar or undef, optional)

=cut

sub file_is_root{
    my $self = shift;

    return $self->{'file_is_root'} = shift if @_;
    return $self->{'file_is_root'};
}

=head2 _virtual_root

 Title   : _virtual_root
 Usage   : $obj->_virtual_root($newval)
 Function:
 Example :
 Returns : value of _virtual_root (a scalar)
 Args    : on set, new value (a scalar or undef, optional)

=cut

sub _virtual_root{
    my $self = shift;

    return $self->{'_virtual_root'} = shift if @_;

    # don't return anything if not in file_is_root mode, or if we don't
    # have a file to derive the root node from
    return unless $self->file_is_root() && $self->file();

    # construct it if we haven't done this before
    if(! $self->{'_virtual_root'}) {
        my ($rt,undef,undef) = fileparse($self->file(), '\..*');
        $rt =~ s/_/ /g;
        $rt = $self->_create_ont_entry($rt);
        $self->_add_term($rt, $self->ontology_name());
        $self->{'_virtual_root'} = $rt;
    }

    return $self->{'_virtual_root'};
}

1;
