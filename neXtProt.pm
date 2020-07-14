=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2020] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=head1 CONTACT

 Ensembl <http://www.ensembl.org/info/about/contact/index.html>
    
=cut

=head1 NAME

 neXtProt

=head1 SYNOPSIS

 mv neXtProt.pm ~/.vep/Plugins
 ./vep -i variations.vcf --plugin neXtProt

=head1 DESCRIPTION

=cut

package neXtProt;

use strict;
use warnings;
use JSON::XS;

use Data::Dumper;

use Bio::EnsEMBL::Variation::Utils::BaseVepPlugin;

use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepPlugin);

my $default_output = {
  neXtProt_MatureProtein => 1,
  neXtProt_NucleotidePhosphateBindingRegion => 1,
  neXtProt_Variant => 1,
  neXtProt_Domain => 1,
  neXtProt_MiscellaneousRegion => 1,
  neXtProt_InteractingRegion => 1
};

sub new {
  my $class = shift;
  
  my $self = $class->SUPER::new(@_);

  my $param_hash = $self->params_to_hash();

  if(defined($param_hash->{isoform})) {
    my $disease = $param_hash->{isoform};
    $default_output->{'neXtProt_url'} = 1;
    $self->{isoform} = $disease;
  }

  return $self;
}

sub feature_types {
  return ['Transcript'];
}

sub get_header_info {
  my $self = shift;

  my %header;

  if($self->{isoform}) {
    $header{'neXtProt_url'} = 'neXtProt URL';
  }

  # 'neXtProt_AntibodyMapping' => 'provides information as to whether an antibody mapping is “unique”, “pseudo-unique” or “not unique” in a manner analogous to that of peptide mappings',
  # 'neXtProt_TopologicalDomain' => 'Location of non-membrane regions of membrane-spanning proteins',
  # 'neXtProt_SequenceConflict' => 'Sequence discrepancies of unknown origin',
  # 'neXtProt_TransmembraneRegion' => 'Extent of a membrane-spanning region',
  # 'neXtProt_CompositionallyBiasedRegion' => 'Region of compositional bias in the protein',
  # 'neXtProt_ModifiedResidue' => 'Modified residues',
  # 'neXtProt_Repeat' => 'Positions of repeated sequence motifs or repeated domains',
  # 'neXtProt_Mutagenesis' => 'Site which has been experimentally altered by mutagenesis',
  # 'neXtProt_ModifiedResidue' => 'Modified residues',
  # $header{'neXtProt_PdbMapping'} = 'Protein 3D structure';
  # DisulfideBond
  # GlycosylationSite
  # neXtProt_DisulfideBond -> 152,177,Or C-152 with C-183
  # neXtProt_GlycosylationSite -> 167,167,O-linked (GalNAc...) threonine
  # neXtProt_CompositionallyBiasedRegion -> 671,677,Gly/Ser-rich

  $header{'neXtProt_MatureProtein'} = 'Extent of an active peptide or a polypetide chain in the mature protein';
  $header{'neXtProt_NucleotidePhosphateBindingRegion'} = 'Nucleotide phosphate binding region';
  $header{'neXtProt_Variant'} = 'Natural variant of the protein';
  $header{'neXtProt_Domain'} = 'Position and type of each modular protein domain';
  $header{'neXtProt_MiscellaneousRegion'} = 'Region of interest in the sequence';
  $header{'neXtProt_InteractingRegion'} = 'Region interacting with another macromolecule';

  return \%header;
}

sub run {
  my ($self, $tva) = @_;

  return {} unless $tva->transcript->translation;
  my $tv = $tva->transcript_variation;

  my $peptide_start = defined($tv->translation_start) ? $tv->translation_start : undef;
  my $peptide_end = defined($tv->translation_end) ? $tv->translation_end : undef;
  my $transcript_id = $tva->transcript->translation->stable_id;

  return {} unless defined($transcript_id) && defined($peptide_start) && defined($peptide_end); 

  # print "\n$transcript_id, $peptide_start, $peptide_end\n";

  my $query = $self->get_sparql_query($peptide_start,$transcript_id);

  # run command
  my $query_output = `curl -X POST -H "Accept:application/sparql-results+json" --data-urlencode "query=$query" https://sparql.nextprot.org/ 2> /dev/null`;

  my $output = decode_json ($query_output);

  # print "AFTER: ", Dumper($output);

  my %result_hash;
  my %result_hash_final;

  # Output format: 'iso','spos','epos','annot_type','callret-4'
  # 'iso' -> isoform URL to neXtProt page; 'spos' -> start position; 'epos' -> end position; 'annot_type' -> annotation type (e.g. PdbMapping, Variant, etc.);
  # 'callret-4' -> data
  my $output_list = $output->{results}->{bindings};
  return {} if (@$output_list == 0);

  foreach my $results (@$output_list) {
    my $isoform_url = $results->{iso}->{value};
    my $start_pos = $results->{spos}->{value};
    my $end_pos = $results->{epos}->{value};
    my $annot_type = $results->{annot_type}->{value};
    $annot_type =~ s/.*#//;
    my $data = $results->{'callret-4'}->{value};
    # PdbMapping values contain ';'
    if($data =~ /; /) {
      $data =~ s/; / /g;
    }
    $data =~ s/\.$//;

    # There is only one URL
    if($self->{isoform} && !$result_hash{'neXtProt_url'}) {
      my @isoform_value = ($isoform_url);
      $result_hash{'neXtProt_url'} = \@isoform_value;
    }

    # Some annot_type have more than one value
    # Need to check if it's not duplicate
    if($result_hash{'neXtProt_'.$annot_type}) {
      # $result_hash{'neXtProt_'.$annot_type} .= "|".$start_pos.','.$end_pos.','.$data;
      my $annot_type_data = $start_pos.','.$end_pos.','.$data;
      push @{$result_hash{'neXtProt_'.$annot_type}}, $annot_type_data unless grep{$_ eq $annot_type_data} @{$result_hash{'neXtProt_'.$annot_type}};
    }
    else {
      # $result_hash{'neXtProt_'.$annot_type} = $start_pos.','.$end_pos.','.$data;
      my @list_of_data;
      push @list_of_data, $start_pos.','.$end_pos.','.$data;
      $result_hash{'neXtProt_'.$annot_type} = \@list_of_data;
    }
  }

  foreach my $key (keys %result_hash) {
    my $data_to_return = $result_hash{$key};
    my $join_data = join('|', @$data_to_return);
    $result_hash_final{$key} = $join_data;
    # print $key, " -> ", $join_data, "\n";
  }

  # foreach my $key (keys %$default_output) {
  #   if($result_hash{$key}) {
  #     my $data_to_return = $result_hash{$key};
  #     my $join_data = join('|', @$data_to_return);
  #     $result_hash_final{$key} = $join_data;
  #   }
  #   else {
  #     $result_hash_final{$key} = '-';
  #   }
  # }

  return \%result_hash_final;
}

sub get_sparql_query {
  my ($self, $peptide_start, $transcript_id) = @_;

  my $query = "PREFIX : <http://nextprot.org/rdf#>
               PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
               PREFIX up: <http://purl.uniprot.org/core/>
               PREFIX isoform: <http://nextprot.org/rdf/isoform/>
               select distinct ?iso ?spos ?epos ?annot_type str(?txt)
               where {
                 values ?poi {$peptide_start}
                 values ?ensp {'$transcript_id'} # for e101 mapping is for canoncila - somtimes our cnonincal is different than the uniprot canonical
                 bind (IRI(CONCAT('http://rdf.ebi.ac.uk/resource/ensembl.protein/',?ensp)) as ?ENSP_IRI)
                 SERVICE <http://sparql.uniprot.org/sparql> {
                   SELECT * WHERE {
                     ?enst up:translatedTo ?ENSP_IRI .
                     ?enst rdfs:seeAlso  ?upiso .
                   }
                 }
               BIND(IRI(replace(str(?upiso),'http://purl.uniprot.org/isoforms/','http://nextprot.org/rdf/isoform/NX_')) AS ?iso) .
                 ?entry :isoform ?iso .
                 ?iso :positionalAnnotation ?statement .
                 ?statement rdfs:comment ?txt .
                 ?statement a ?annot_type .
                 ?statement :start ?spos; :end ?epos .
                 filter((?spos <= ?poi) && (?epos >= ?poi))
               } order by ?spos";

  return $query;
}

1;

