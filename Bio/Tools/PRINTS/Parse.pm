# Perl Module for fingerPRINTScan results
#
# Cared for by Henning Hermjakob, hhe@ebi.ac.uk
#

=head1 NAME

Bio::Tools::PRINTS::Parse - Object representing fingerPRINTScan output results

=head1 SYNOPSIS

   # parse a fingerPRINTScan result file
   $res = new Bio::Tools::PRINTS::Parse( -file => 'output.prints');


   foreach $feature ( $res->each_Feature ) {
      print $feature->gff_string();
   }

=head1 DESCRIPTION

This object represents fingerPRINTScan output and produces one SeqFeatureSet containing all significant matches. 

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

package Bio::Tools::PRINTS::Parse;

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

  $self->_parse_fingerPRINTScan($fh);

  return $make; # success - we hope!
}

=head2 _parse_fingerPRINTScan

 Title   : _parse_fingerPRINTScan
 Usage   : $res->_parse_fingerPRINTScan($filehandle)
 Function:
 Returns : 
 Args    :


=cut

sub _parse_fingerPRINTScan{
  my $self = shift;
  my $file = shift;
  my ($sequenceName, $matchName, $matchLine, $from, $to, $length);
  my $count;
  my $matchAC;

  # Read a complete result set for one sequence in one go.
  {
    local ($/) = "\n3TBF\n";
    
    while (<$file>) {
      
      # Skip if there are no significant results
      if (/NO SIGNIFICANT RESULTS ABOVE/) {
	next;
      }
      
      my %matchACs = ();

      # Sequence name
      next unless (($sequenceName) = /^Sn\;\s+(\w+)\s*/);
      #$self->throw("Could not find sequence name in \n$_");
     

      # store ACs of matches
      while (/1TBH\s+(\S+)\s+.*(PR\d+)\s*\n/gx) {
	$matchACs{$1} = $2;
      }

      # Parse the significant results
      # this should match
      # 3TBH DISINTEGRIN     1  of  2  53.72   710     4.25e-09  CAHGLCCEDCQLKPAGTACR                                    20   2    464  479 
      
      foreach $matchLine (/3TBH.*\n/g){

	if (($matchName, $length, $from) 
	    = $matchLine =~ /3TBH\s+            # match line start
                             (\w+)              # pattern name
                              \s+
                              \d+               # match number
                              \s+of\s+
                              \d+               # number of matches
                              \s+
                              \d+\.?\d*         # IdScore
                              \s+
                              \d+               # PfScore
                              \s+
                              \S+               # Pvalue
                              \s+
                              \w+               # Sequence
                              \s+
                             (\d+)              # Match length
                              \s+
                              \d+               # low ???
                              \s+
                             (\d+)              # match position
                              \s+
                              \d+               # high ???
                              \s*\n             # end of line
                            /x){
	    $to = $from + $length-1;
	    # Store the results in a new SeqFeature object
	    my $feature = Bio::SeqFeature::Generic->new();

	    $feature->seqname($sequenceName);
	    $feature->add_tag_value('match_id', $matchName);
	    # If the Prints AC has been found, use it.
	    if ($matchACs{$matchName}) {
	      $feature->primary_tag($matchACs{$matchName});
	    } else {
	      $feature->primary_tag('-');
	    }
	    $feature->start($from);
	    $feature->end($to);
	    $feature->source_tag('fingerPRINTScan');
	    $self->add_Feature($feature);
	    $count++;
	}
      } 
    }
  }
  return $count;
}
1;				# says use was ok
__END__


