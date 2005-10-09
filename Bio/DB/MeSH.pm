# $Id$
#
# BioPerl module for Bio::DB::MeSH
#
# Cared for by Heikki Lehvaslaiho, heikki@ebi.ac.uk
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::DB::MeSH - Term retrieval from a Web MeSH database

=head1 SYNOPSIS

 my $mesh = new Bio::DB::MeSH();
 my $term=$mesh->get_exact_term('Butter');
 print $term->description;


=head1 DESCRIPTION

This class retrieves a term from the Medical Subject Headings database
by the National Library of Medicine of USA.  See
http://www.nlm.nih.gov/mesh/meshhome.html.

This class implements Bio::SimpleAnalysisI and wraps its methods under
get_exact_term().

The web access uses my favorite, WWW::Mechanize, but in its absense
falls back to bioperl module Bio::WebAgent which is a subclass of
LWP::UserAgent. If not even that is not installed, it uses
Bio::Root::HTTPget.


=head1 SEE ALSO

L<Bio::Phenotype::MeSH::Term>

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to the
Bioperl mailing lists Your participation is much appreciated.

  bioperl-l@bioperl.org                       - General discussion
  http://bio.perl.org/MailList.html           - About the mailing lists

=head2 Reporting Bugs

report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via the
web:

  http://bugzilla.bioperl.org/

=head1 AUTHOR

Heikki Lehvaslaiho, heikki@ebi.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::DB::MeSH;
use vars qw( @ISA );
use strict;

use Bio::Tools::Analysis::SimpleAnalysisBase;
use Bio::Phenotype::MeSH::Term;
use Bio::Phenotype::MeSH::Twig;

@ISA = qw(Bio::Tools::Analysis::SimpleAnalysisBase );


my $URL = 'http://www.nlm.nih.gov/mesh/MBrowser.html';
my $ANALYSIS_SPEC= {name => 'MeSH term retrival',
                    type => 'Entry retrieval'};
my $INPUT_SPEC = [
                  {mandatory=>'true',
                   type => 'scalar',
                   'name'=> 'value',
                  },
                 ];

my  $RESULT_SPEC =
    {
     '' => 'Bio::Phenotype::MeSH::Term',
     raw => 'raw output',
    };

sub _init {
    my $self = shift;
    $self->url($URL);
    $self->{'_ANALYSIS_SPEC'} =$ANALYSIS_SPEC;
    $self->{'_INPUT_SPEC'} =$INPUT_SPEC;
    $self->{'_RESULT_SPEC'} =$RESULT_SPEC;
    $self->{'_ANALYSIS_NAME'} = $ANALYSIS_SPEC->{'name'};
    $self->_webmodule;
    return $self;
}

sub _webmodule {
    my ($self) = shift;
    $self->{'_webmodule'} = '';
    eval {
        require WWW::Mechanize;
    };
    unless ($@) {
        $self->{'_webmodule'} = 'WWW::Mechanize';
        return;
    }
    eval {
        require LWP::UserAgent;
    };
    unless ($@) {
        $self->{'_webmodule'} = 'Bio::WebAgent';
        return;
    }
    require Bio::Root::HTTPget;
    $self->{'_webmodule'} = 'Bio::Root::HTTPget';
    1;
}

=head2 get_exact_term

  Title   : get_exact_term
  Usage   : $s = $db->get_exact_term($value);
  Function: Retrive a single MeSH term using a unique ID or exact name.
  Example :
  Returns : a Bio::Phenotype::MeSH::Term object
  Args    : scalar, UID or name of a MeSH term

The returned term object contains information about the immediate
vincinity of the term in the terminology hierarchy. See
L<Bio::Phenotype::MeSH::Twig>.

=cut


sub get_exact_term {
    my ($self, $value) = @_;
    $self->{'_term'} = undef;
    $self->run($value) if $value;
    $self->throw("Could not connect to the server")
        unless $self->status eq 'COMPLETED';
    return $self->result;
}


sub run {
    my ($self, $value) = @_;

    # check input
    $self->throw("Need a MeSH name or ID  as an input [$value]") if ref $value;

    # internal run()
    $self->_run($value);
}


sub  _run {
    my ($self, $value)  = @_;
    $self->throw('Need a value [$value]')
        unless $value;;
    # delay repeated calls by default by 3 sec, set delay() to change
#    $self->sleep;

    $self->status('TERMINATED_BY_ERROR');

    if ($self->{'_webmodule'} eq 'WWW::Mechanize') {
        print "using WWW::Mechanize...\n" if $self->verbose > 0;
        my $agent = WWW::Mechanize->new();
        $agent->get($self->url);
        $agent->status == 200
            or print STDERR "Could not connect to the server\n" and return;

        $agent->form_name('MB');

        $agent->field("term", $value);
        if ($value =~ /\w\d{6}/) {
            $agent->field("field", 'uid');
        } else {
            $agent->field("field", 'entry');
        }
        $agent->click("exact");

        $self->{'_content'} = $agent->content();
        $self->status('COMPLETED');
        return;
    }
    elsif ($self->{'_webmodule'} eq 'Bio::WebAgent') {
        print "using LWP::UserAgent...\n" if $self->verbose > 0;
        my $response;
        if ($value =~ /\w\d{6}/) {
            $self->{'_content'} =
                $response = eval {
                    $self->get("http://www.nlm.nih.gov/cgi/mesh/2003/MB_cgi?field=uid&term=$value")
                };
            print STDERR "Could not connect to the server\n" and return
                if $@;
        } else {
            $self->{'_content'} =
                eval {
                    $response = $self->get("http://www.nlm.nih.gov/cgi/mesh/2003/MB_cgi?field=entry&term=$value")
                };
            print STDERR "Could not connect to the server\n" and return
                if $@;
        }
        if ($response->is_success) {
            $self->{'_content'} =  $response->content;
            $self->status('COMPLETED');
        }
        return;
    } else {
        print "using Bio::Root::HTTPget...\n" if $self->verbose > 0;
        my $agent = new Bio::Root::HTTPget;
        if ($value =~ /\w\d{6}/) {
            $self->{'_content'} =
                eval {
                    $agent->get("http://www.nlm.nih.gov/cgi/mesh/2003/MB_cgi?field=uid&term=$value")
                };
            print STDERR "Could not connect to the server\n" and return
                if $@;
        } else {
            $self->{'_content'} =
                eval {
                    $agent->get("http://www.nlm.nih.gov/cgi/mesh/2003/MB_cgi?field=entry&term=$value")
                };
            print STDERR "Could not connect to the server\n" and return
                if $@;
        }
        $self->status('COMPLETED');
    }
}

sub result {
    my ($self,$value) = @_;

    $self->throw("Could not retrive results") unless $self->status('COMPLETED');

    # no processing
    return $self->{'_content'} if $value && $value eq 'raw';


    # create a MeSH::Term object
    $_ = $self->{'_content'};
    print substr ($_, 0, 100), "\n" if $self->verbose > 0;
    my ($id) = m|Unique ID</TH><TD>(.*?)</TD>|i;
    my ($name) = m|MeSH Heading</TH><TD>([^<]+)|i;
    my ($desc) = m|Scope Note</TH><TD>(.*?)</TD>|i;
    $desc =~ s/<.*?>//sg;

    my $term = Bio::Phenotype::MeSH::Term->new(-id => $id,
                                               -name => $name,
                                               -description => $desc
                                              );
    my ($trees) = $self->{'_content'} =~ /MeSH Tree Structures(.*)/s;

    while (m|Entry Term</TH><TD>([^<]+)|ig) {
        $term->add_synonym($1);
        print "Synonym: |$1|\n" if $self->verbose > 0;
    }

    foreach (split /<HR>/i, $trees ) {
        next unless /$name/;
        s/<TD.*?>/ /sgi;
        s/<.*?>//sg;
        s/&nbsp;/ /sg;
        #print "|$_|";
        my ($treeno) = /$name \[([^]]+)]/;
        my ($parent_treeno) = $treeno =~ /(.*)\.\d{3}/;
        my ($parent) =  /\n +(\w.+) \[$parent_treeno\]/;

        my $twig = Bio::Phenotype::MeSH::Twig->new(-parent => $parent);
        $term->add_twig($twig);

        print "Parent: |$parent|\n" if $self->verbose > 0;
        while (/\n +(\w.+) \[$treeno\./g ) {
            $twig->add_child($1);
            print "Child: |$1|\n" if $self->verbose > 0;
        }

        while (/\n +(\w.+) \[$parent_treeno\./g ) {
            next if $name eq $1;
            $twig->add_sister($1);
            print "Sister: |$1|\n" if $self->verbose > 0;
        }
    }
    return $term;
}

1;
