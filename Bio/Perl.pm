#
# BioPerl module for Bio::Perl
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

Bio::Perl - Functional access to BioPerl for people who don't know objects

=head1 SYNOPSIS

  use Bio::Perl;

  # will guess file format from extension
  $seq_object = read_sequence($filename);

  # forces genbank format
  $seq_object = read_sequence($filename,'genbank');

  # reads an array of sequences
  @seq_object_array = read_all_sequences($filename,'fasta');

  # sequences are Bio::Seq objects, so the following methods work
  # for more info see Bio::Seq, or do 'perldoc Bio/Seq.pm'

  print "Sequence name is ",$seq_object->display_id,"\n";
  print "Sequence acc  is ",$seq_object->accession_number,"\n";
  print "First 5 bases is ",$seq_object->subseq(1,5),"\n";

  # get the whole sequence as a single string

  $sequence_as_a_string = $seq_object->seq();

  # writing sequences

  write_sequence(">$filename",'genbank',$seq_object);

  write_sequence(">$filename",'genbank',@seq_object_array);

  # making a new sequence from just a string

  $seq_object = new_sequence("ATTGGTTTGGGGACCCAATTTGTGTGTTATATGTA",
      "myname","AL12232");

  # getting a sequence from a database (assumes internet connection)

  $seq_object = get_sequence('swissprot',"ROA1_HUMAN");

  $seq_object = get_sequence('embl',"AI129902");

  $seq_object = get_sequence('genbank',"AI129902");

  # BLAST a sequence (assummes an internet connection)

  $blast_report = blast_sequence($seq_object);

  write_blast(">blast.out",$blast_report);


=head1 DESCRIPTION

Easy first time access to BioPerl via functions.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists.  Your participation is much appreciated.

  bioperl-l@bioperl.org

=head2 Support 

Please direct usage questions or support issues to the mailing list:

I<bioperl-l@bioperl.org>

rather than to the module maintainer directly. Many experienced and 
reponsive experts will be able look at the problem and quickly 
address it. Please include a thorough description of the problem 
with code and data examples if at all possible.

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution. Bug reports can be submitted via the web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Ewan Birney

Email birney@ebi.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

#'
# Let the code begin...


package Bio::Perl;
use vars qw(@EXPORT @EXPORT_OK $DBOKAY);
use strict;
use Carp;

use Bio::SeqIO;
use Bio::Seq;
use Bio::Root::Version '$VERSION';
BEGIN {
    eval {
        require Bio::DB::EMBL;
        require Bio::DB::GenBank;
        require Bio::DB::SwissProt;
        require Bio::DB::RefSeq;
        require Bio::DB::GenPept;
    };
    if( $@ ) {
        $DBOKAY = 0;
    } else {
        $DBOKAY = 1;
    }
}

use base qw(Exporter);

@EXPORT = qw(read_sequence read_all_sequences write_sequence
             new_sequence get_sequence translate translate_as_string
             reverse_complement revcom revcom_as_string
             reverse_complement_as_string blast_sequence write_blast);

@EXPORT_OK = @EXPORT;


=head2 read_sequence

 Title   : read_sequence
 Usage   : $seq = read_sequence('sequences.fa')
           $seq = read_sequence($filename,'genbank');

           # pipes are fine
           $seq = read_sequence("my_fetching_program $id |",'fasta');

 Function: Reads the top sequence from the file. If no format is given, it will
           try to guess the format from the filename. If a format is given, it
           forces that format. The filename can be any valid perl open() string
           - in particular, you can put in pipes

 Returns : A Bio::Seq object. A quick synopsis:
           $seq_object->display_id - name of the sequence
           $seq_object->seq        - sequence as a string

 Args    : Two strings, first the filename - any Perl open() string is ok
           Second string is the format, which is optional

For more information on Seq objects see L<Bio::Seq>.

=cut

sub read_sequence{
    my ($filename,$format) = @_;

    if( !defined $filename ) {
        confess "read_sequence($filename) - usage incorrect";
    }

    my $seqio;

    if( defined $format ) {
        $seqio = Bio::SeqIO->new( '-file' => $filename, '-format' => $format);
    } else {
        $seqio = Bio::SeqIO->new( '-file' => $filename);
    }

    my $seq = $seqio->next_seq();

    return $seq;
}


=head2 read_all_sequences

 Title   : read_all_sequences
 Usage   : @seq_object_array = read_all_sequences($filename);
           @seq_object_array = read_all_sequences($filename,'genbank');

 Function: Just as the function above, but reads all the sequences in the
           file and loads them into an array.

           For very large files, you will run out of memory. When this
           happens, you've got to use the SeqIO system directly (this is
           not so hard! Don't worry about it!).

 Returns : array of Bio::Seq objects

 Args    : two strings, first the filename (any open() string is ok)
           second the format (which is optional)

See L<Bio::SeqIO> and L<Bio::Seq> for more information

=cut

sub read_all_sequences{
    my ($filename,$format) = @_;

    if( !defined $filename ) {
        confess "read_all_sequences($filename) - usage incorrect";
    }

    my $seqio;

    if( defined $format ) {
        $seqio = Bio::SeqIO->new( '-file' => $filename, '-format' => $format);
    } else {
        $seqio = Bio::SeqIO->new( '-file' => $filename);
    }

    my @seq_array;

    while( my $seq = $seqio->next_seq() ) {
        push(@seq_array,$seq);
    }

    return @seq_array;
}


=head2 write_sequence

 Title   : write_sequence
 Usage   : write_sequence(">new_file.gb",'genbank',$seq)
           write_sequence(">new_file.gb",'genbank',@array_of_sequence_objects)

 Function: writes sequences in the specified format

 Returns : true

 Args    : filename as a string, must provide an open() output file
           format as a string
           one or more sequence objects


=cut

sub write_sequence{
    my ($filename,$format,@sequence_objects) = @_;

    if( scalar(@sequence_objects) == 0 ) {
        confess("write_sequence(filename,format,sequence_object)");
    }

    my $error = 0;
    my $seqname = "sequence1";

    # catch users who haven't passed us a filename we can open
    if( $filename !~ /^\>/ && $filename !~ /^|/ ) {
        $filename = ">".$filename;
    }

    my $seqio = Bio::SeqIO->new('-file' => $filename, '-format' => $format);

    foreach my $seq ( @sequence_objects ) {
        my $seq_obj;

        if( !ref $seq ) {
            if( length $seq > 50 ) {
                # odds are this is a sequence as a string, and someone has not figured out
                # how to make objects. Warn him/her and then make a sequence object
                # from this
                if( $error == 0 ) {
                    carp("WARNING: You have put in a long string into write_sequence.\n".
                         "I suspect this means that this is the actual sequence\n".
                         "In the future try the\n".
                         "  new_sequence method of this module to make a new sequence object.\n".
                         "Doing this for you here\n");
                    $error = 1;
                }

                $seq_obj = new_sequence($seq,$seqname);
                $seqname++;
            } else {
                confess("You have a non object [$seq] passed to write_sequence. It maybe that you".
                        "want to use new_sequence to make this string into a sequence object?");
            }
        } else {
            if( !$seq->isa("Bio::SeqI") ) {
                confess("object [$seq] is not a Bio::Seq object; can't write it out");
            }
            $seq_obj = $seq;
        }

        # finally... we get to write out the sequence!
        $seqio->write_seq($seq_obj);
    }
    1;
}

=head2 new_sequence

 Title   : new_sequence
 Usage   : $seq_obj = new_sequence("GATTACA", "kino-enzyme");

 Function: Construct a sequency object from sequence string
 Returns : A Bio::Seq object

 Args    : sequence string
           name string (optional, default "no-name-for-sequence")
           accession - accession number (optional, no default)

=cut

sub new_sequence{
    my ($seq,$name,$accession) = @_;

    if( !defined $seq ) {
        confess("new_sequence(sequence_as_string) usage");
    }

    $name ||= "no-name-for-sequence";

    my $seq_object = Bio::Seq->new( -seq => $seq, -id => $name);

    $accession && $seq_object->accession_number($accession);

    return $seq_object;
}

=head2 blast_sequence

 Title   : blast_sequence
 Usage   : $blast_result = blast_sequence($seq)
           $blast_result = blast_sequence('MFVEGGTFASEDDDSASAEDE');

 Function: If the computer has Internet accessibility, blasts
           the sequence using the NCBI BLAST server against nrdb.

           It chooses the flavour of BLAST on the basis of the sequence.

           This function uses Bio::Tools::Run::RemoteBlast, which itself
           use Bio::SearchIO - as soon as you want to know more, check out
           these modules
 Returns : Bio::Search::Result::GenericResult.pm

 Args    : Either a string of protein letters or nucleotides, or a
           Bio::Seq object

=cut

sub blast_sequence {
    my ($seq,$verbose) = @_;

    if( !defined $verbose ) {
        $verbose = 1;
    }

    if( !ref $seq ) {
        $seq = Bio::Seq->new( -seq => $seq, -id => 'blast-sequence-temp-id');
    } elsif ( !$seq->isa('Bio::PrimarySeqI') ) {
        croak("[$seq] is an object, but not a Bio::Seq object, cannot be blasted");
    }

    require Bio::Tools::Run::RemoteBlast;

    my $prog = ( $seq->alphabet eq 'protein' ) ? 'blastp' : 'blastn';
    my $e_val= '1e-10';

    my @params = ( '-prog' => $prog,
                   '-expect' => $e_val,
                   '-readmethod' => 'SearchIO' );

    my $factory = Bio::Tools::Run::RemoteBlast->new(@params);

    my $r = $factory->submit_blast($seq);
    if( $verbose ) {
        print STDERR "Submitted Blast for [".$seq->id."] ";
    }
    sleep 5;

    my $result;

    LOOP :
    while( my @rids = $factory->each_rid) {
        foreach my $rid ( @rids ) {
            my $rc = $factory->retrieve_blast($rid);
            if( !ref($rc) ) {
                if( $rc < 0 ) {
                    $factory->remove_rid($rid);
                }
                if( $verbose ) {
                    print STDERR ".";
                }
                sleep 10;
            } else {
                $result = $rc->next_result();
                $factory->remove_rid($rid);
                last LOOP;
            }
        }
    }

    if( $verbose ) {
        print STDERR "\n";
    }
    return $result;
}

=head2 write_blast

 Title   : write_blast
 Usage   : write_blast($filename,$blast_report);

 Function: Writes a BLAST result object (or more formally
           a SearchIO result object) out to a filename
           in BLAST-like format

 Returns : none

 Args    : filename as a string
           Bio::SearchIO::Results object

=cut

sub write_blast {
    my ($filename,$blast) = @_;

    if( $filename !~ /^\>/ && $filename !~ /^|/ ) {
        $filename = ">".$filename;
    }

    my $output = Bio::SearchIO->new( -output_format => 'blast', -file => $filename);

    $output->write_result($blast);

}

=head2 get_sequence

 Title   : get_sequence
 Usage   : $seq_object = get_sequence('swiss',"ROA1_HUMAN");

 Function: If the computer has Internet access this method gets
           the sequence from Internet accessible databases. Currently
           this supports Swissprot ('swiss'), EMBL ('embl'), GenBank
           ('genbank'), GenPept ('genpept'), and RefSeq ('refseq').

           Swissprot and EMBL are more robust than GenBank fetching.

           If the user is trying to retrieve a RefSeq entry from
           GenBank/EMBL, the query is silently redirected.

 Returns : A Bio::Seq object

 Args    : database type - one of swiss, embl, genbank, genpept, or
           refseq

=cut

my $genbank_db = undef;
my $genpept_db = undef;
my $embl_db = undef;
my $swiss_db = undef;
my $refseq_db = undef;

sub get_sequence{
    my ($db_type,$identifier) = @_;
    if( ! $DBOKAY ) {
        confess ("Your system does not have one of LWP, HTTP::Request::Common, IO::String\n".
                 "installed so the DB retrieval method is not available.\n".
                 "Full error message is:\n $!\n");
        return;
    }
    $db_type = lc($db_type);

    my $db;

    if( $db_type =~ /genbank/ ) {
        if( !defined $genbank_db ) {
            $genbank_db = Bio::DB::GenBank->new();
        }
        $db = $genbank_db;
    }
    if( $db_type =~ /genpept/ ) {
        if( !defined $genpept_db ) {
            $genpept_db = Bio::DB::GenPept->new();
        }
        $db = $genpept_db;
    }

    if( $db_type =~ /swiss/ ) {
        if( !defined $swiss_db ) {
            $swiss_db = Bio::DB::SwissProt->new();
        }
        $db = $swiss_db;
    }

    if( $db_type =~ /embl/ ) {
        if( !defined $embl_db ) {
            $embl_db = Bio::DB::EMBL->new();
        }
        $db = $embl_db;
    }

    if( $db_type =~ /refseq/ or ($db_type !~ /swiss/ and
                                 $identifier =~ /^\s*N\S+_/)) {
        if( !defined $refseq_db ) {
            $refseq_db = Bio::DB::RefSeq->new();
        }
        $db = $refseq_db;
    }

    my $seq;

    if( $identifier =~ /^\w+\d+$/ ) {
        $seq = $db->get_Seq_by_acc($identifier);
    } else {
        $seq = $db->get_Seq_by_id($identifier);
    }

    return $seq;
}


=head2 translate

 Title   : translate
 Usage   : $seqobj = translate($seq_or_string_scalar)

 Function: translates a DNA sequence object OR just a plain
           string of DNA to amino acids
 Returns : A Bio::Seq object

 Args    : Either a sequence object or a string of
           just DNA sequence characters

=cut

sub translate {
    my ($scalar) = shift;

    my $obj;
    if( ref $scalar ) {
        if( !$scalar->isa("Bio::PrimarySeqI") ) {
            confess("Expecting a sequence object not a $scalar");
        } else {
            $obj= $scalar;
        }
    } else {
        # check this looks vaguely like DNA
        my $n = ( $scalar =~ tr/ATGCNatgcn/ATGCNatgcn/ );
        if( $n < length($scalar) * 0.85 ) {
            confess("Sequence [$scalar] is less than 85% ATGCN, which doesn't look very DNA to me");
        }
        $obj = Bio::PrimarySeq->new(-id => 'internalbioperlseq',-seq => $scalar);
    }
    return $obj->translate();
}


=head2 translate_as_string

 Title   : translate_as_string
 Usage   : $seqstring = translate_as_string($seq_or_string_scalar)

 Function: translates a DNA sequence object OR just a plain
           string of DNA to amino acids
 Returns : A string of just amino acids

 Args    : Either a sequence object or a string of
           just DNA sequence characters

=cut

sub translate_as_string {
    my ($scalar) = shift;
    my $obj = Bio::Perl::translate($scalar);
    return $obj->seq;
}


=head2 reverse_complement

 Title   : reverse_complement
 Usage   : $seqobj = reverse_complement($seq_or_string_scalar)

 Function: reverse complements a string or sequence argument
           producing a Bio::Seq - if you want a string, you
           can use reverse_complement_as_string
 Returns : A Bio::Seq object

 Args    : Either a sequence object or a string of
           just DNA sequence characters

=cut

sub reverse_complement {
    my ($scalar) = shift;

    my $obj;

    if( ref $scalar ) {
        if( !$scalar->isa("Bio::PrimarySeqI") ) {
            confess("Expecting a sequence object not a $scalar");
        } else {
            $obj= $scalar;
        }

    } else {

        # check this looks vaguely like DNA
        my $n = ( $scalar =~ tr/ATGCNatgcn/ATGCNatgcn/ );

        if( $n < length($scalar) * 0.85 ) {
            confess("Sequence [$scalar] is less than 85% ATGCN, which doesn't look very DNA to me");
        }

        $obj = Bio::PrimarySeq->new(-id => 'internalbioperlseq',-seq => $scalar);
    }

    return $obj->revcom();
}

=head2 revcom

 Title   : revcom
 Usage   : $seqobj = revcom($seq_or_string_scalar)

 Function: reverse complements a string or sequence argument
           producing a Bio::Seq - if you want a string, you
           can use reverse_complement_as_string

           This is an alias for reverse_complement
 Returns : A Bio::Seq object

 Args    : Either a sequence object or a string of
           just DNA sequence characters

=cut

sub revcom {
    return &Bio::Perl::reverse_complement(@_);
}


=head2 reverse_complement_as_string

 Title   : reverse_complement_as_string
 Usage   : $string = reverse_complement_as_string($seq_or_string_scalar)

 Function: reverse complements a string or sequence argument
           producing a string
 Returns : A string of DNA letters

 Args    : Either a sequence object or a string of
           just DNA sequence characters

=cut

sub reverse_complement_as_string {
    my ($scalar) = shift;
    my $obj = &Bio::Perl::reverse_complement($scalar);
    return $obj->seq;
}


=head2 revcom_as_string

 Title   : revcom_as_string
 Usage   : $string = revcom_as_string($seq_or_string_scalar)

 Function: reverse complements a string or sequence argument
           producing a string
 Returns : A string of DNA letters

 Args    : Either a sequence object or a string of
           just DNA sequence characters

=cut

sub revcom_as_string {
    my ($scalar) = shift;
    my $obj = &Bio::Perl::reverse_complement($scalar);
    return $obj->seq;
}


1;
