#
# BioPerl module for InterPro_BioSQL_Handler
#
# Please direct questions and support issues to <bioperl-l@bioperl.org>
#
# Cared for by Juguang Xiao, juguang@tll.org.sg
#
# Copyright Juguang Xiao
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::OntologyIO::Handlers::InterPro_BioSQL_Handler - parse an InterPro XML file and persist the resulting terms to a Biosql database

=head1 SYNOPSIS

   # see load_interpro.pl in bioperl-db/scripts/biosql/

=head1 DESCRIPTION

This module is for parsing an InterPro XML file and persist the
resulting terms to a Biosql database as soon as the term is complete
as signaled by the appropriate xml tag. This parser takes advantage of
SAX, a stream-based XML parser technology, to keep the used memory as
small as possible. The alternative parser for InterPro, module
InterProHandler, builds up the entire ontology in memory, which given
the size of the latest InterPro releases requires a huge amount of
memory.

This module takes the following non-standard arguments upon
instantiation.

   -db           the adaptor factory as returned by a call to
                 Bio::DB::BioDB->new()
   -version      the InterPro version (not available as property!)
   -term_factory the object factory to use for creating terms

Note that there are two alternatives for how to persist the terms and
relationships to the database. The default is using the adaptor
factory passed as -db or set as a property to create persistent
objects and store them in the database. The alternative is to specify
a term persistence and a relationship persistence handler; if one or
both have been set, the respective handler will be called with each
term and relationship that is to be stored. See properties
persist_term_handler and persist_relationship_handler.

=head1 AUTHOR

Juguang Xiao, juguang@tll.org.sg

=head1 Contributors

Hilmar Lapp, hlapp at gmx.net

=head2 APPENDIX

The rest of the documentation details each of the object methods.
Interal methods are usually preceded with a _

=cut

package Bio::OntologyIO::Handlers::InterPro_BioSQL_Handler;
use strict;
use Bio::Ontology::Ontology;
use Bio::Ontology::Term;
use Bio::Ontology::TermFactory;
use Bio::Ontology::RelationshipType;
use Bio::Ontology::Relationship;
use Bio::Annotation::DBLink;
use Bio::Annotation::Reference;

use base qw(Bio::OntologyIO::Handlers::BaseSAXHandler);

my $is_a_rel;
my $count=0;

sub _initialize {
    my($self,@args)=@_;
    $self->SUPER::_initialize(@args);
    my ($db, $version, $fact) = $self->_rearrange(
        [qw(DB  VERSION TERM_FACTORY)], @args);
    $self->db($db) if $db; # this is now a property and may be set later
    if (!$fact) {
        $fact = Bio::Ontology::TermFactory->new(-type=>"Bio::Ontology::Term");
    }
    $self->term_factory($fact);
    my $ontology = Bio::Ontology::Ontology->new(-name => 'InterPro');
    if (defined($version)) {
        $version = "InterPro version $version";
        $ontology->definition($version);
    }
    $self->_ontology($ontology);
    $is_a_rel = Bio::Ontology::RelationshipType->get_instance('IS_A');
    $is_a_rel->ontology($ontology);
}

=head2 term_factory

 Title   : term_factory
 Usage   : $obj->term_factory($newval)
 Function: Get/set the ontology term factory to use.

           As a user of this module it is not necessary to call this
           method as there will be default. In order to change the
           default, the easiest way is to instantiate
           L<Bio::Ontology::TermFactory> with the proper -type
           argument. Most if not all parsers will actually use this
           very implementation, so even easier than the aforementioned
           way is to simply call
           $ontio->term_factory->type("Bio::Ontology::MyTerm").

 Example :
 Returns : value of term_factory (a Bio::Factory::ObjectFactoryI object)
 Args    : on set, new value (a Bio::Factory::ObjectFactoryI object, optional)


=cut

sub term_factory{
    my $self = shift;

    return $self->{'term_factory'} = shift if @_;
    return $self->{'term_factory'};
}

=head2 db

 Title   : db
 Usage   : $obj->db($newval)
 Function: Sets or retrieves the database adaptor factory.

           The adaptor factory is a Bio::DB::DBAdaptorI compliant
           object and will be used to obtain the persistence adaptors
           necessary to serialize terms and relationships to the
           database.

           Usually, you will obtain such an object from a call to
           Bio::DB::BioDB. You *must* set this property before
           starting the parse.

           Note that this property is immutable once set, except that
           you may set it to undef. Therefore, be careful not to set
           to undef before setting the desired real value.

 Example :
 Returns : value of db (a Bio::DB::DBAdaptorI compliant object)
 Args    : on set, new value (a Bio::DB::DBAdaptorI compliant object
           or undef, optional)


=cut

sub db {
    my $self=shift;
    if(@_){
        my $db = shift;
        if ($db && exists($self->{_db}) && ($self->{_db} != $db)) {
            $self->throw('db may not be modified once set');
        }
        $self->{_db}=$db;
    }
    return $self->{_db};
}

=head2 persist_term_handler

 Title   : persist_term_handler
 Usage   : $obj->persist_term_handler($handler,@args)
 Function: Sets or retrieves the persistence handler for terms along
           with the constant set of arguments to be passed to the
           handler.

           If set, the first argument will be treated as a closure and
           be called for each term to persist to the database. The
           term will be passed as a named parameter (-term), followed
           by the other arguments passed to this setter. Note that
           this allows one to pass an arbitrary configuration to the
           handler.

           If not set, terms will be persisted along with their
           relationships using the respective persistence adaptor
           returned by the adaptor factory (see property db).

 Example :
 Returns : an array reference with the values passed on set, or an empty
           array if never set
 Args    : On set, an array of values. The first value is the handler
           as a closure; all other values will be passed to the handler
           as constant argument.


=cut

sub persist_term_handler{
    my $self = shift;

    return $self->{'persist_term_handler'} = [@_] if @_;
    return $self->{'persist_term_handler'} || [];
}

=head2 persist_relationship_handler

 Title   : persist_relationship_handler
 Usage   : $obj->persist_relationship_handler($handler,@args)
 Function: Sets or retrieves the persistence handler for relationships
           along with the constant set of arguments to be passed to
           the handler.

           If set, the first argument will be treated as a closure and
           be called for each relationship to persist to the database. The
           relationship will be passed as a named parameter (-rel), followed
           by the other arguments passed to this setter. Note that
           this allows one to pass an arbitrary configuration to the
           handler.

           If not set, relationships will be persisted along with their
           relationships using the respective persistence adaptor
           returned by the adaptor factory (see property db).

 Example :
 Returns : an array reference with the values passed on set, or an empty
           array if never set
 Args    : On set, an array of values. The first value is the handler
           as a closure; all other values will be passed to the handler
           as constant argument.


=cut

sub persist_relationship_handler{
    my $self = shift;

    return $self->{'persist_relationship_handler'} = [@_] if @_;
    return $self->{'persist_relationship_handler'} || [];
}

=head2 _persist_term

 Title   : _persist_term
 Usage   :
 Function: Persists a term to the database, using either a previously
           set persistence handler, or the adaptor factory directly.
 Example :
 Returns :
 Args    : the ontology term to persist


=cut

sub _persist_term {
    my $self = shift;
    my $term = shift;

    my ($handler,@args) = @{$self->persist_term_handler};
    if ($handler) {
        &$handler('-term' => $term, @args);
    } else {
        # no handler; we'll do this ourselves straight and simple
        my $db = $self->db();
        my $pterm = $db->create_persistent($term);
        eval {
            $pterm->create();
            $pterm->commit();
        };
        if ($@) {
            $pterm->rollback();
            $self->warn("failed to store term '".$term->name."': ".$@);
        }
    }
}

=head2 _persist_relationship

 Title   : _persist_relationship
 Usage   :
 Function: Persists a relationship to the database, using either a
           previously set persistence handler, or the adaptor factory
           directly.

 Example :
 Returns :
 Args    : the term relationship to persist


=cut

sub _persist_relationship {
    my $self = shift;
    my $rel = shift;

    my ($handler,@args) = @{$self->persist_relationship_handler};
    if ($handler) {
        &$handler('-rel' => $rel, @args);
    } else {
        # no handler; we'll do this ourselves straight and simple
        my $db = $self->db();
        my $prel = $db->create_persistent($rel);
        eval {
            $prel->create();
            $prel->commit();
        };
        if ($@) {
            $prel->rollback();
            $self->warn("failed to store relationship of subject '"
                        .$rel->subject_term->name."' to object '"
                        .$rel->object_term->name.": ".$@);
        }
    }
}

=head2 _persist_ontology

 Title   : _persist_ontology
 Usage   :
 Function: Perists the ontology itself to the database, by either
           inserting or updating it.

           Note that this will only create or update the ontology as
           an entity, not any of its terms, relationships, or
           relationship types.

 Example :
 Returns : the ontology as a peristent object with primary key
 Args    : the ontology to persist as a Bio::Ontology::OntologyI
           compliant object


=cut

sub _persist_ontology{
    my $self = shift;
    my $ont = shift;
    my $db = $self->db();

    # do a lookup first; chances are we have this already in the database
    my $adp = $db->get_object_adaptor($ont);
    # to avoid clobbering this ontology's properties with possibly older ones
    # from the database we'll need an object factory
    my $ontfact =
        Bio::Factory::ObjectFactory->new(-type=>"Bio::Ontology::Ontology");
    # do the lookup:
    my $found = $adp->find_by_unique_key($ont, '-obj_factory' => $ontfact);
    # make a persistent object of the ontology
    $ont = $db->create_persistent($ont);
    # transfer primary key if found in the lookup
    $ont->primary_key($found->primary_key) if $found;
    # insert or update
    my $result;
    eval {
        $result = $ont->store();
    };
    if ($@ || !$result) {
        $adp->rollback();
        $self->throw("failed to update ontology '"
                     .$ont->name."' in database".($@ ? ": $@" : ""));
    }

    # done - we don't commit here
    return ref($result) ? $result : $ont;
}

sub start_document {
    my $self = shift;
    my $ont = $self->_ontology;
    my @iprtypes = (
                    $self->create_term(-identifier=>'IPR:Family',
                                       -name=>'Family',
                                       -ontology => $ont),
                    $self->create_term(-identifier=>'IPR:Domain',
                                       -name=>'Domain',
                                       -ontology => $ont),
                    $self->create_term(-identifier=>'IPR:Repeat',
                                       -name=>'Repeat',
                                       -ontology => $ont),
                    $self->create_term(-identifier=>'IPR:PTM',
                                       -name=>'post-translational modification',
                                       -ontology => $ont),
                    $self->create_term(-identifier=>'IPR:Active_site',
                                       -name=>'Active_site',
                                       -ontology => $ont),
                    $self->create_term(-identifier=>'IPR:Binding_site',
                                       -name=>'Binding_site',
                                       -ontology => $ont),
                    );
    foreach my $iprtype (@iprtypes) {
        $self->_persist_term($iprtype);
        $ont->add_term($iprtype);
    }
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
        #$term->definition();

        my ($object_term) =
            ($ont->engine->get_term_by_identifier("IPR:".$args{type}));

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
            #print "EXAmPLE:\t", $example->database, '|', $example->primary_id, "\n";
        }elsif($top eq 'child'){
            ;
        }elsif($top eq 'member_list'){
            my $dblink=Bio::Annotation::DBLink->new(
                -dbname => $args{id},
                -primary_id => $args{dbkey},
                -comment => $args{name}
            );
        }elsif($top eq 'external_doc_list'){
            ;
        }elsif($top eq 'publication'){
            if($args{db} eq 'MEDLINE'){
                $self->_current_hash->{medline} =$args{dbkey};
            } elsif($args{db} eq 'PUBMED'){
                $self->_current_hash->{pubmed} =$args{dbkey};
            }else{
                $self->warn("'".$args{dbkey}."' is not a MEDLINE publication, "
                            ."don't know how to handle");
            }
        }elsif($top eq 'structure_db_links'){
            ;
        }elsif($top eq 'abstract'){
            ;
        }
        #else{
        #    $self->warn("unrecognized element '$top' in element '$tag', ignoring");
        #}
    }elsif($tag eq 'publication'){
        my $publication = Bio::Annotation::Reference->new();
        $self->_current_hash->{publication} = $publication;
    }elsif($tag eq 'author_list'){
        ;
    }elsif($tag eq 'journal'){
        ;
    }elsif($tag eq 'location'){
        ;
    }elsif($tag eq 'year'){
        ;
    } elsif (($tag eq 'dbinfo') && ($self->_top_tag eq 'release')) {
        my $entrydate = $args{file_date} || '';
        $entrydate =~ s/ \d{2}:\d{2}:\d{2}//;
        my $def = $ont->definition() || '';
        $def .= "\n" if length($def) > 0;
        $def .= $args{dbname}." version ".$args{version}.", "
            .$args{entry_count}." entries, ".$entrydate;
        $ont->definition($def);
    }
    #else{
    #    $self->warn("unrecognized element '$tag', ignoring");
    #}

    $self->_visited_count_inc($tag);
    $self->_push_tag($tag);
}


sub end_element {
    my $self=shift;
    my $tag=shift->{Name};
    my $chars_in=$self->_chars_hash->{$tag};
    if($tag eq 'interpro'){
        my $rel = $self->_relationship;
        # store subject term first in order to give the handler a chance to
        # apply whatever custom behaviour
        # (note that the object term is the InterPro type and has been stored
        # at the start of the whole document)
        $self->_persist_term($rel->subject_term);
        # the store the relationship to the InterPro type
        $self->_persist_relationship($rel);
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
        $self->_relationship->subject_term->add_dbxref(-dbxrefs => [$example]);
        $self->_chars_hash->{example}='';
    }elsif($tag eq 'publication'){
        my $publication = $self->_create_publication;
        $self->_relationship->subject_term->add_reference($publication);
    }elsif($tag eq 'author_list'){
        $self->_current_hash->{author} =$chars_in;
    }elsif($tag eq 'title'){
        $self->_current_hash->{title}=$chars_in;
    } elsif ($tag eq 'release') {
        my $ont = $self->_persist_ontology($self->_ontology);
        $self->_ontology($ont) if $ont;
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
    return shift->term_factory->create_object(@_);
}

sub _ontology {
    my $self = shift;
    return $self->{_ontology}=shift if @_;
    return $self->{_ontology};
}

sub _relationship {
    my $self =shift;
    $self->{_relationship}=shift if @_;
    return $self->{_relationship};
}
sub _create_publication {
    my $self=shift;
    my $publ = $self->_current_hash->{publication};
    my $journal = $self->_current_hash->{journal} || '<no journal>';
    my $year = $self->_current_hash->{year} || '<no year>';
    my $page_location = $self->_current_hash->{page_location} || '<no pages>';
    my $volumn = $self->_current_hash->{volumn} || '<no volume>';
    my $medline =
        $self->_current_hash->{medline} || $self->_current_hash->{pubmed};

    $publ->authors($self->_current_hash->{author});
    $publ->location("$journal, $year, V $volumn, $page_location");
    $publ->title($self->_current_hash->{title});
    $publ->medline($medline);
    if ($self->_current_hash->{pubmed}
        && ($self->_current_hash->{pubmed} != $medline)) {
        $publ->pubmed($self->_current_hash->{pubmed});
    }

# Clear the above in current hash
    $self->_current_hash->{publication} = undef;
    $self->_current_hash->{author}      = undef;
    $self->_current_hash->{journal}     = undef;
    $self->_current_hash->{year}        = undef;
    $self->_current_hash->{page_location}=undef;
    $self->_current_hash->{volumn}      = undef;
    $self->_current_hash->{title}       = undef;
    $self->_current_hash->{medline}     = undef;
    $self->_current_hash->{pubmed}      = undef;
}
1;
