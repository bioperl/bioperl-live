#
# BioPerl module for Bio::DB::MeSH
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Heikki Lehvaslaiho, heikki-at-bioperl-dot-org
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::DB::MeSH - Term retrieval from a Web MeSH database

=head1 SYNOPSIS

 my $mesh = Bio::DB::MeSH->new();
 my $term = $mesh->get_exact_term('Butter');
 print $term->description;

=head1 DESCRIPTION

This class retrieves a term from the Medical Subject Headings database
by the National Library of Medicine of USA. 
See L<http://www.nlm.nih.gov/mesh/meshhome.html>.

This class implements L<Bio::SimpleAnalysisI> and wraps its methods under
L<get_exact_term>.

By default, web access uses L<WWW::Mechanize>, but in its absence
falls back to bioperl module L<Bio::WebAgent> which is a subclass of
L<LWP::UserAgent>. If not even that is not installed, it uses
L<Bio::Root::HTTPget>.

=head1 SEE ALSO

L<Bio::Phenotype::MeSH::Term>

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to the
Bioperl mailing lists Your participation is much appreciated.

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

report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via the
web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR

Heikki Lehvaslaiho, heikki-at-bioperl-dot-org

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::DB::MeSH;
use strict;

use Bio::Phenotype::MeSH::Term;
use Bio::Phenotype::MeSH::Twig;

use base qw(Bio::Tools::Analysis::SimpleAnalysisBase);


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


sub _cgi_url {
  my($self, $field, $term) = @_;
  # we don't bother to URI::Escape $field and $term as this is an untainted private sub
  return 'http://www.nlm.nih.gov/cgi/mesh/2003/MB_cgi?field='.$field.'&term='.$term;
}


sub  _run {
    my ($self, $value)  = @_;
    $self->throw('Need a value [$value]')
        unless $value;;
    # delay repeated calls by default by 3 sec, set delay() to change
#    $self->sleep;

    $self->status('TERMINATED_BY_ERROR');

    if ($self->{'_webmodule'} eq 'WWW::Mechanize') {
        $self->debug("using WWW::Mechanize...\n");
        my $agent = WWW::Mechanize->new();
        $agent->get($self->url);
        $agent->status == 200
            or $self->warn("Could not connect to the server\n") and return;

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
        $self->debug("using LWP::UserAgent...\n");
        my $response;
        if ($value =~ /\w\d{6}/) {
            $self->{'_content'} =
                $response = eval {
                    $self->get( $self->_cgi_url('uid', $value) )
                };
            $self->warn("Could not connect to the server\n") and return
                if $@;
        } else {
            $self->{'_content'} =
                eval {
                    $response = $self->get( $self->_cgi_url('entry', $value) )
                };
            $self->warn("Could not connect to the server\n") and return
                if $@;
        }
        if ($response->is_success) {
            $self->{'_content'} =  $response->content;
            $self->status('COMPLETED');
        }
        return;
    } else {
        $self->debug("using Bio::Root::HTTPget...\n");
        my $agent = Bio::Root::HTTPget->new();
        if ($value =~ /\w\d{6}/) {
            $self->{'_content'} =
                eval {
                    $agent->get( $self->_cgi_url('uid', $value) )
                };
            $self->warn("Could not connect to the server\n") and return
                if $@;
        } else {
            $self->{'_content'} =
                eval {
                    $agent->get( $self->_cgi_url('entry', $value) )
                };
            $self->debug("Could not connect to the server\n") and return
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
    $self->debug( substr($_, 0, 100) . "\n");
    my ($id)   = m|Unique \s+ ID </TH>
                   <TD (?: \s+ [^>]+ )? >
                   (.*?)                    # id
                   </TD> |ix;
    my ($name) = m|MeSH \s+ Heading </TH>
                   <TD (?: \s+ [^>]+ )? >
                   (.*?)                    # name
                   </TD> |ix;
    my ($desc) = m|Scope \s+ Note </TH>
                   <TD (?: \s+ [^>]+ )? >
                   (.*?)                    # desc
                   </TD>|isx;
    $self->throw("No description returned: $_") unless defined $desc;
    $desc =~ s/<.*?>//sg;
    $desc =~ s/(?:\r)?\n/ /g;
    $desc =~ s/\( +/\(/g;
    $desc =~ s/ {2,}/ /g;

    my $term = Bio::Phenotype::MeSH::Term->new(-id => $id,
                                               -name => $name,
                                               -description => $desc
                                              );
    my ($trees) = $self->{'_content'} =~ /MeSH Tree Structures(.*)/s;

    while (m|Entry Term</TH><TD(?: [^>]+)?>(.*?)</TD>|ig) {
        $term->add_synonym($1);
        $self->debug("Synonym: |$1|\n");
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

        $self->debug("Parent: |$parent|\n");
        while (/\n +(\w.+) \[$treeno\./g ) {
            $twig->add_child($1);
            $self->debug("Child: |$1|\n");
        }

        while (/\n +(\w.+) \[$parent_treeno\./g ) {
            next if $name eq $1;
            $twig->add_sister($1);
            $self->debug("Sister: |$1|\n");
        }
    }
    return $term;
}

1;
