
#
# BioPerl module for Bio::DB::Universal
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Ewan Birney <birney@ebi.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::DB::Universal - Artificial database that delegates to specific databases

=head1 SYNOPSIS

    $uni = Bio::DB::Universal->new();

    # by default connects to web databases. We can also
    # substitute local databases

    $embl = Bio::Index::EMBL->new( -filename => '/some/index/filename/locally/stored');
    $uni->use_database('embl',$embl);

    # treat it like a normal database. Recognises strings
    # like gb|XXXXXX and embl:YYYYYY

    $seq1 = $uni->get_Seq_by_id("embl:HSHNRNPA");
    $seq2 = $uni->get_Seq_by_acc("gb|A000012");

    # with no separator, tries to guess database. In this case the
    # _ is considered to be indicative of swissprot
    $seq3 = $uni->get_Seq_by_id('ROA1_HUMAN');

=head1 DESCRIPTION

Artificial database that delegates to specific databases, with a
"smart" (well, smartish) guessing routine for what the ids. No doubt
the smart routine can be made smarter.

The hope is that you can make this database and just throw ids at it -
for most easy cases it will sort you out. Personally, I would be
making sure I knew where each id came from and putting it into its own
database first - but this is a quick and dirty solution.

By default this connects to web orientated databases, with all the
reliability and network bandwidth costs this implies. However you can
subsistute your own local databases - they could be Bio::Index
databases (DBM file and flat file) or bioperl-db based (MySQL based)
or biocorba-based (whatever you like behind the corba interface).

Internally the tags for the databases are

   genbank - ncbi dna database
   embl    - ebi's dna database (these two share accession number space)
   swiss   - swissprot + sptrembl (EBI's protein database)

We should extend this for RefSeq and other sequence databases which 
are out there... ;)

Inspired by Lincoln Stein, written by Ewan Birney.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists.  Your participation is much appreciated.

  bioperl-l@bio.perl.org

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

=head1 AUTHOR - Ewan Birney

Email birney@ebi.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::DB::Universal;
use strict;

# Object preamble - inherits from Bio::Root::Root


use Bio::DB::GenBank;
use Bio::DB::SwissProt;
use Bio::DB::EMBL;


use base qw(Bio::DB::RandomAccessI Bio::Root::Root);
# new() can be inherited from Bio::Root::Root

sub new {
    my ($class) = @_;

    my $self = {};
    bless $self,$class;

    $self->{'db_hash'} = {};

    # default databases

    $self->use_database('embl',Bio::DB::EMBL->new);
    $self->use_database('genbank',Bio::DB::GenBank->new);
    $self->use_database('swiss',Bio::DB::GenBank->new);

    return $self;
}


=head2 get_Seq_by_id

 Title   : get_Seq_by_id
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_Seq_by_id{
   my ($self,$str) = @_;

   my ($tag,$id) = $self->guess_id($str);

   return $self->{'db_hash'}->{$tag}->get_Seq_by_id($id);
}


=head2 get_Seq_by_acc

 Title   : get_Seq_by_acc
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_Seq_by_acc {
   my ($self,$str) = @_;

   my ($tag,$id) = $self->guess_id($str);

   return $self->{'db_hash'}->{$tag}->get_Seq_by_acc($id);
}



=head2 guess_id

 Title   : guess_id
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub guess_id{
   my ($self,$str) = @_;
   
   if( $str =~ /(\S+)[:|\/;](\w+)/ ) {
       my $tag;
       my $db = $1;
       my $id = $2;
       if( $db =~ /gb/i || $db =~ /genbank/i || $db =~ /ncbi/i ) {
	   $tag = 'genbank';
       } elsif ( $db =~ /embl/i || $db =~ /emblbank/ || $db =~ /^em/i ) {
	   $tag = 'embl';
       } elsif ( $db =~ /swiss/i || $db =~ /^sw/i || $db =~ /sptr/ ) {
	   $tag = 'swiss';
       } else {
	   # throw for the moment
	   $self->throw("Could not guess database type $db from $str");
       }
       return ($tag,$id);

   } else {
       my $tag;
       # auto-guess from just the id
       if( $str =~ /_/ ) {
	   $tag = 'swiss';
       } elsif ( $str =~ /^[QPR]\w+\d$/ ) {
	   $tag = 'swiss';
       } elsif ( $str =~ /[A-Z]\d+/ ) {
	   $tag = 'genbank';
       } else {
	   # default genbank...
	   $tag = 'genbank';
       }
       return ($tag,$str);
   }

   
}


=head2 use_database

 Title   : use_database
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub use_database{
   my ($self,$name,$database) = @_;

   $self->{'db_hash'}->{$name} = $database;
}

1;
