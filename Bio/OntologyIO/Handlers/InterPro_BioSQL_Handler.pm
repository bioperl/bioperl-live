#
#
#
#

=head1 NAME

Bio::OntologyIO::Handlers::InterPro_BioSQL_Handler

=head1 SYNOPSIS

    see load_interpro.pl in bioperl-db/scripts/biosql/

=head1 DESCRIPTION

This module take the advantage of SAX, a stream-based XML parser technology, to
keep the used memory as small as possible. InterProHandler in the same 
directory use another scheme to keep all interpro record in the engine and keep
the ontology data in memory, almost impossible to operate it! 

=head1 AUTHOR

Juguang Xiao, juguang@tll.org.sg

=cut

package Bio::OntologyIO::Handlers::InterPro_BioSQL_Handler;
use strict;
use vars qw(@ISA);
use Bio::Root::Root;
use Bio::Ontology::Ontology;
use Bio::Ontology::Term;
use Bio::Ontology::RelationshipType;
use Bio::Ontology::Relationship;
@ISA = qw(Bio::Root::Root);

my $is_a_rel;
my $count=0;
sub new {
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(@args);
    my ($db, $version) = $self->_rearrange(
        [qw(DB  VERSION)], @args);
    defined $db or $self->throw('db must be assigned');
    $self->_db($db);
    my $ontology = Bio::Ontology::Ontology->new(
        -name => 'InterPro',
        -definition => "InterPro, $version"
    );
    $self->_ontology($ontology);
    $self->{_states} = {};
    $self->{_contents} = {};
    $is_a_rel = Bio::Ontology::RelationshipType->get_instance('IS_A');
    $is_a_rel->ontology($ontology);
    return $self;
}

sub start_document {
    my $self = shift;
    my $ont = $self->_ontology;
    my $db = $self->_db;
    $ont->add_term($self->create_term(-identifier=>'Family', -name=>'Family'));
    $ont->add_term($self->create_term(-identifier=>'Domain', -name=>'Domain'));
    $ont->add_term($self->create_term(-identifier=>'Repeat', -name=>'Repeat'));
    $ont->add_term($self->create_term(-identifier=>'PTM', 
        -name=>'post-translational modification'));
    $ont->add_term($self->create_term(
        -identifier=>'Active_site', -name=>'Active_site'));
    $ont->add_term($self->create_term(
        -identifier=>'Binding_site', -name=>'Binding_site'));

    
}

sub start_element {
    my ($self, $element) = @_;
    my $tag=$element->{Name};
    my %args=%{$element->{Attributes}};
    my $ont = $self->_ontology;

    if($tag eq 'interpro'){
        my $id = $args{id};
        my $term = $self->create_term(-identifier=>$id);
        $term->ontology($ont);
        $term->add_synonym($args{short_name});
        $term->definition();
        
        my $object_term =
            ($ont->engine->get_term_by_identifier($args{type}))[0];
        
        my $rel = Bio::Ontology::Relationship->new(
            -subject_term => $term,
            -predicate_term => $is_a_rel,
            -object_term => $object_term,
            -ontology => $ont
        );
        $self->_relationship($rel);
    }
    $self->_states->{$tag}=1;
}


sub end_element {
    my ($self, $element)=@_;
    my $tag=$element->{Name};
    if($tag eq 'interpro'){
        my $rel = $self->_relationship;
        my $db = $self->_db;
        my $prel = $db->create_persistent($rel);
        eval {$prel->create};
        if($@){
            my $interpro = $rel->subject_term;
            print STDERR "$@\n";
            print STDERR "\n\n\n";
            print STDERR $interpro->identifier, "\t", $interpro->name, "\t";
            print STDERR $interpro->description, "\n\n\n";
            print '-'x80, "\n";
            
        }
#        print $count++, "\n";
    }elsif($tag eq 'name'){
        my $rel = $self->_relationship;
        $rel->subject_term->name($self->_name);
    }elsif($tag eq 'abstract'){
        my $rel = $self->_relationship;
        my $abstract = $self->_contents->{abstract};
        $abstract =~ s/\n/ /g;
        $rel->subject_term->definition($abstract);
        $self->_contents->{abstract} = '';
    }
    $self->_states->{$tag}=0;
}


sub characters{
    my ($self, $characters) = @_;
    my $text = $characters->{Data};
    $self->_name($text) if $self->_states->{name};
    $self->_contents->{abstract} .= $text if $self->_states->{abstract};
}

sub create_term {
    my ($self, @args) = @_;
    my $term = Bio::Ontology::Term->new(@args);
    return $term;
}

sub _db {
    my $self=shift;
    if(@_){
        my $db = shift;
        if(exists $self->{_db}){ $self->throw('db can be set only once'); }
        $self->{_db}=$db;
    }
    return $self->{_db};
}

sub _ontology {
    my $self = shift;
    $self->{_ontology}=shift if @_;
    return $self->{_ontology};
}

sub _relationship {
    my $self =shift;
    $self->{_relationship}=shift if @_;
    return $self->{_relationship};
}

sub _name {
    my $self=shift;
    $self->{_name}=shift if @_;
    return $self->{_name};
}

sub _contents {
    return shift->{_contents};
}
sub _states {
    return shift->{_states};
}
1;
