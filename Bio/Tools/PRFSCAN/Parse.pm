# Perl Module for PRFSCAN results
#
#  by Evgueni.Zdobnov@ebi.ac.uk
#

=head1 NAME

Bio::Tools::PRFSCAN::Parse - Object representing PRFSCAN output results

=head1 SYNOPSIS

   # parse a PRFSCAN result file
   $res = new Bio::Tools::PRFSCAN::Parse( -file => 'output.prfscan');


   foreach $feature ( $res->each_Feature ) {
      print $feature->gff_string();
   }

=head1 DESCRIPTION

This object represents PRFSCAN output and produces one SeqFeatureSet containing all matches. 

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists.  Your participation is much appreciated.

  bioperl-l@bioperl.org             - General discussion
  http://bio.perl.org/MailList.html - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via email
or the web:

  bioperl-bugs@bio.perl.org
  http://bio.perl.org/bioperl-bugs/

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

package Bio::Tools::PRFSCAN::Parse;

use vars qw(@ISA);
use Carp;
use strict;

use Bio::Root::Object;
use Bio::SeqFeatureSet;
use FileHandle;

@ISA = ( 'Bio::SeqFeatureSet' );

sub _initialize {
  my($self,@args) = @_;

  my $make = new Bio::SeqFeatureSet;

  my($file,$fh,$parsetype) = $self->_rearrange([qw(FILE
						   FH
						   TYPE
						   )],
					       @args);

  if( !defined $file && ! defined $fh ) {
      $self->throw("No file/filehandle definition to result file");
  }

  if( defined $file ) {
    $fh = new FileHandle;
    $fh->open($file) || $self->throw("Could not open file [$file] $!");
  }

  $self->_parse_PRFSCAN($fh);

  return $make; # success - we hope!
}

=head2 _parse_PRFSCAN

 Title   : _parse_PRFSCAN
 Usage   : $res->_parse_PRFSCAN($filehandle)
 Function:
 Returns : 
 Args    :


=cut

sub _parse_PRFSCAN {
  my $self = shift;
  my $file = shift;
  my ($sequenceName);
  my $count;

  # Read a complete result set for one sequence in one go.
  {
    local ($/) = "\n>>>";
    
    while (<$file>) {
            
	# Sequence name
	unless (($sequenceName) = /^>*([a-zA-Z0-9_]+)/) {
	    $self->throw("Could not find sequence name in \n$_");
	}
	
	# store matches
        #>KRINGLE_2  13.773   1367 pos.    100 -   164 PS50070|
	while (/>(\w+).*?pos\.\s*(\d+)\s*-\s*(\d+)\s*(PS\d+)\|.*?\n/sg) {
	  # Store the results in a new SeqFeature object
	  my $feature = Bio::SeqFeature::Generic->new();
	  
	  $feature->seqname($sequenceName);
	  $feature->add_tag_value('match_id', $1);
	  $feature->primary_tag($4);
	  $feature->start($2);
	  $feature->end($3);
	  $feature->source_tag('PRFSCAN');
	  $self->add_Feature($feature);
	  $count++;
	}
    }
  }
  return $count;
}
1;				
__END__


