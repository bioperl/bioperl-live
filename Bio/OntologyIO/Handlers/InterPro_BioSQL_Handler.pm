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

=head2 APPENDIX

The rest of the documentation details each of the object methods.
Interal methods are usually preceded with a _

=cut

package Bio::OntologyIO::Handlers::InterPro_BioSQL_Handler;
use strict;
use vars qw(@ISA);
use Bio::OntologyIO::Handlers::BaseSAXHandler;
use Bio::Ontology::Ontology;
use Bio::Ontology::Term;
use Bio::Ontology::RelationshipType;
use Bio::Ontology::Relationship;
use Bio::Annotation::DBLink;

@ISA = qw(Bio::OntologyIO::Handlers::BaseSAXHandler);

my $is_a_rel;
my $count=0;

sub _initialize {
    my($self,@args)=@_;
    $self->SUPER::_initialize(@args);
    my ($db, $version) = $self->_rearrange(
        [qw(DB  VERSION)], @args);
    defined $db or $self->throw('db must be assigned');
    $self->_db($db);
    my $ontology = Bio::Ontology::Ontology->new(
        -name => 'InterPro',
        -definition => "InterPro, $version"
    );
    $self->_ontology($ontology);
    $is_a_rel = Bio::Ontology::RelationshipType->get_instance('IS_A');
    $is_a_rel->ontology($ontology);
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
    my $self=shift;
    my $tag=$_[0]->{Name};
    my %args=%{$_[0]->{Attributes}};
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
    }elsif($tag eq 'example'){
        my $example = Bio::Annotation::DBLink->new;
        $self->_current_hash->{example} = $example;
    }elsif($tag eq 'db_xref'){
        my $top = $self->_top_tag;
        if($top eq 'example'){
            my $example = $self->_current_hash->{example};
            $example->database($args{db});
            $example->primary_id($args{dbkey});
            print "EXAmPLE:\t", $example->database, '|', $example->primary_id, "\n";
        }else{
            print STDERR "Possible??\n";
        }
    }
    $self->_visited_count_inc($tag);
    $self->_push_tag($tag);
}


sub end_element {
    my $self=shift;
    my $tag=shift->{Name};
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
        $rel->subject_term->name($self->_chars_hash->{name});
        $self->_chars_hash->{name}='';
    }elsif($tag eq 'abstract'){
        my $rel = $self->_relationship;
        my $abstract = $self->_chars_hash->{abstract};
        $abstract =~ s/\n/ /g;
        $rel->subject_term->definition($abstract);
        $self->_chars_hash->{abstract} = '';
    }elsif($tag eq 'example'){
        my $example = $self->_current_hash->{example};
        my $comment = $self->_chars_hash->{example};
        $comment =~ s/^(\s+)//; $comment =~ s/(\s+)$//;
        $example->comment($comment);
        $self->_relationship->subject_term->add_dblink($example);
        $self->_chars_hash->{example}='';
    }
    $self->_pop_tag;
    $self->_visited_count_dec($tag);
}


sub characters {
    my $self=shift;
    my $text = shift->{Data};

    my $top_tag =$self->_top_tag;
    $self->_chars_hash->{$top_tag} .= $text;
    
#    $self->_chars_hash->{abstract} .= $text if $self->_visited_count('abstract');
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

1;
