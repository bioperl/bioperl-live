# $Id$
#
# BioPerl module for Bio::SeqFeatureProducer
#
# Cared for by Jason Stajich <jason@chg.mc.duke.edu>,
# 
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::SeqFeatureProducer - add Sequence Features from an Analysis to an
sequence object in an easy fashion

=head1 SYNOPSIS

    use Bio::SeqFeatureProducer;
    my $sfp = new Bio::SeqFeatureProducer(-method=>'genscan',
					  -input =>$analysisfile,
					  );
   $sfp->add_features($seq);

   # or

   use Bio::Tools::Genscan;
   my $parser = new Bio::Tools::Genscan(-file=>$analysisfile);
   my $sfp = new Bio::SeqFeatureProducer();
   $sfp->add_features($seq,$parser);

=head1 DESCRIPTION

SeqFeatureProducer is a helper class for adding features to objects 
from a SeqAnalysisParserI compliant objects.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this
and other Bioperl modules. Send your comments and suggestions preferably
 to one of the Bioperl mailing lists.
Your participation is much appreciated.

  bioperl-l@bioperl.org                - General discussion
  http://bio.perl.org/MailList.html    - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via email
or the web:

  bioperl-bugs@bio.perl.org
  http://bio.perl.org/bioperl-bugs/

=head1 AUTHOR - Jason Stajich

Email Jason Stajich <jason@chg.mc.duke.edu>

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut

package Bio::SeqFeatureProducer;
use strict;
use Bio::Root::RootI;
use vars qw(@ISA %DRIVERVALUES);

@ISA = qw(Bio::Root::RootI); 

BEGIN {
    %DRIVERVALUES = ( 
		      'mzef' => { MODULE => "Bio::Tools::MZEF" },
		      'genscan' => { MODULE => "Bio::Tools::Genscan" }
		      );
    # would like to GFF and GAME support added here
}

sub new {
    my ($class, @args) = @_;
    my $self = bless {}, $class;
    $self->_initialize(@args);
    return $self;
}

sub _initialize {
    my ($self, @args) = @_;
    my $make = $self->SUPER::_initialize(@args);
    my ($method, $input, $params) = $self->_rearrange([qw(METHOD 
							  INPUT
							  PARAMS)],
						      @args);
    
    if( defined $method ) {
	$method = lc $method;
	if( !defined $DRIVERVALUES{$method} ) {
	    $self->throw("Asked for AnalysisParser $method which does not exist\nonly (". join(",", keys %DRIVERVALUES) . ") supported");	    
	}

	$params = [] if( !defined $params );
	my $modulename = $DRIVERVALUES{$method}->{MODULE};
	if( &_load_module($modulename) == 0 ) {
	    return undef;
	}
	my $inputmethod = '-file';
	if( ref($input) =~ /GLOB/i ) {
	    $inputmethod = '-fh';
	}
	my $parser = $modulename->new($inputmethod => $input,@$params);
	$self->analysis_parser($parser);
    }
    return $make;
}

=head2 add_features

 Title   : add_features
 Usage   : $obj->add_features($seq,$analfeatParser)
 Function: adds features to a sequence based on the analysis parser
 Example :
 Returns : void
 Args    : B<seq>    - a sequence object to receive new features
	   B<analFeatureParser> - a SeqAnalysisParserI compliant object 
             [optional]

=cut

sub add_features {
    my ($self, $seq, $parser) = @_;
    
    if( !defined $seq || !ref($seq) || ! $seq->isa('Bio::SeqI') ) {
	$self->throw("Must provide a Bio::SeqI  as input");
    }
    if( defined $parser && ref($parser) && 
	( $parser->isa('Bio::SeqAnalysisParserI') || 
	  $parser->can('next_feature') )) {
	$self->analysis_parser($parser);
    } else {
	$parser = $self->analysis_parser($parser);
    }
    
    if( !defined $parser || ! ref($parser) || 
	! $parser->isa("Bio::SeqAnalysisParserI" ) ) {
	$self->warn("Cannot get features from object " .
		    ref($parser) . "\n");
    }
    
    while( my $feature = $parser->next_feature() ) {
	$seq->add_SeqFeature($feature);
    }
}

=head2 analysis_parser

 Title   : analysis_parser
 Usage   : $obj->analysis_parser($parser); or
           my $parser = $obj->analysis_parser()
 Function: Get/Set method for analysis parser object
 Example : 
 Returns : Bio::SeqAnalysisParserI object
 Args    : B<parser> [optional] - Bio::SeqAnalysisParserI object

=cut

sub analysis_parser {
    my ($self, $value) = @_;
    if( defined $value ) {
	if( ! ref($value) || ! $value->isa('Bio::SeqAnalysisParserI') ) {
	    $self->warn("specifying a value for analysis_parser ($value) which is not a Bio::SeqAnalysisParserI");
	    return;
	}
	$self->{'_analysisparser'} = $value;
    }
    return $self->{'_analysisparser'};
}

=head2 _load_module

 Title   : _load_module
 Usage   : *INTERNAL SeqFeatureProducer stuff*
 Function: Loads up (like use) a module at run time on demand
 Example :
 Returns :
 Args    :

=cut

sub _load_module {
  my ($name) = @_;
  my ($module, $load, $m);

  $module = "_<$name.pm";
  $load = "$name.pm";
  $load =~ s/:/\//g;

  return 1 if $main::{$module};
  eval {
    require $load;
  };
  if ( $@ ) {
    print STDERR <<END;
$load: $name cannot be found
Exception $@
END
  ;
    return;
  }
  return 1;
}

1;
