# BioPerl module for Bio::Tools::Profile
#
# Cared for by Balamurugan Kumarasamy
#
# You may distribute this module under the same terms as perl itself
# POD documentation - main docs before the code

=head1 NAME

Bio::Tools::Profile - parse Profile output

=head1 SYNOPSIS

 use Bio::Tools::Profile;
 my $profile_parser = new Bio::Tools::Profile(-fh =>$filehandle );
 while( my $profile_feat = $profile_parser->next_result ) {
       push @profile_feat, $profile_feat;
}

=head1 DESCRIPTION

 Parser for Profile output

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
 http://bugzilla.bioperl.org/
=head1 AUTHOR - Balamurugan Kumarasamy

 Email: fugui@worf.fugu-sg.org

=head1 APPENDIX

 The rest of the documentation details each of the object methods.
 Internal methods are usually preceded with a _


=cut


package Bio::Tools::Profile;
use vars qw(@ISA);
use strict;

use Bio::Root::Root;
use Bio::SeqFeature::FeaturePair;
use Bio::Root::IO;
use Bio::SeqFeature::Generic;

@ISA = qw(Bio::Root::Root Bio::Root::IO );



=head2 new

 Title   : new
 Usage   : my $obj = new Bio::Tools::Profile();
 Function: Builds a new Bio::Tools::Profile object
 Returns : Bio::Tools::Profile
 Args    : -filename
           -fh ($filehandle)

=cut

sub new {
    my($class,@args) = @_;

    my $self = $class->SUPER::new(@args);
    $self->_initialize_io(@args);

    return $self;
}

=head2 next_result

 Title   : next_result
 Usage   : my $feat = $profile_parser->next_result
 Function: Get the next result set from parser data
 Returns : L<Bio::SeqFeature::FeaturePair>
 Args    : none


=cut

sub next_result {
        my ($self) = @_;

        my %printsac;
        my $line;
        my @features;
        while ($_=$self->_readline()) { 
                $line = $_;
                chomp $line;
                my ($nscore,$rawscore,$from,$to,$hfrom,$hto,$ac) = $line =~ /(\S+)\s+(\d+)\s*pos.\s+(\d*)\s*-\s+(\d*)\s*\[\s+(\d*),\s+(\S*)\]\s*(\w+)/;
                #for example in this output line 
                #38.435   2559 pos.  19958 - 20212 [    1,    -1] PS50011|PROTEIN_KINASE_DOM Protein kinase domain profile.
                #$nscore = 38.435
                #$rawscore = 2559
                #$from = 19958
                #$end = 20212
                #$hfrom = 1
                #$hto =-1
                #$ac = PS50011
                my $feat = "$ac,$from,$to,$hfrom,$hto,$nscore";
                my $new_feat= $self->create_feature($feat);
                return $new_feat

        } 
}


=head2 create_feature

 Title   : create_feature
 Usage   : my $feat= $profile_parser->create_feature($feature)
 Function: creates a Bio::SeqFeature::FeaturePair object
 Returns : L<Bio::SeqFeature::FeaturePair>
 Args    :


=cut

sub create_feature {
        my ($self, $feat) = @_;

        my @f = split (/,/,$feat);


        my $hto = $f[4];

        if ($f[4] =~ /-1/) {
 
           $hto = $f[2] - $f[1] + 1;
        
        }


        my $feat1 = new Bio::SeqFeature::Generic ( -start => $f[1],
                                                   -end => $f[2],
                                                   -score => $f[5],
                                                   -source=>'pfscan',
                                                   -primary=>$f[0]);  
        
        my $feat2 = new Bio::SeqFeature::Generic (-start => $f[3],
                                                  -end => $hto,
                                                  );  

        my $feature = new Bio::SeqFeature::FeaturePair(-feature1 => $feat1,
                                                    -feature2 => $feat2);

        return $feature;
        
}
1;
