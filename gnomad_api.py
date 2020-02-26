import sys
import requests

# !!! IMPORTANT !!!
# GnomAD API is constantly evolving
# If this code does not work, then gnomAD API has changed and it has to be updated.
# Check gnomAD API on the website to figure out what has to be changed:
# https://gnomad.broadinstitute.org/api

GNOMAD_API_URL = 'https://gnomad.broadinstitute.org/api?query='

VARIANT_QUERY_VIEW = """query GnomadVariant($variantId: String!, $datasetId: DatasetsSupportingFetchVariantDetails!) {
  variant(variantId: $variantId, dataset: $datasetId) {
    alt
    chrom
    pos
    ref
    variantId
    xpos
    ... on GnomadVariantDetails {
      age_distribution {
        het {
          bin_edges
          bin_freq
          n_smaller
          n_larger
        }
        hom {
          bin_edges
          bin_freq
          n_smaller
          n_larger
        }
      }
      colocatedVariants
      multiNucleotideVariants {
        combined_variant_id
        changes_amino_acids
        n_individuals
        other_constituent_snvs
      }
      exome {
        ac
        an
        ac_hemi
        ac_hom
        faf95 {
          popmax
          popmax_population
        }
        filters
        populations {
          id
          ac
          an
          ac_hemi
          ac_hom
          subpopulations {
            id
            ac
            an
            ac_hom
          }
        }
        qualityMetrics {
          alleleBalance {
            alt {
              bin_edges
              bin_freq
              n_smaller
              n_larger
            }
          }
          genotypeDepth {
            all {
              bin_edges
              bin_freq
              n_smaller
              n_larger
            }
            alt {
              bin_edges
              bin_freq
              n_smaller
              n_larger
            }
          }
          genotypeQuality {
            all {
              bin_edges
              bin_freq
              n_smaller
              n_larger
            }
            alt {
              bin_edges
              bin_freq
              n_smaller
              n_larger
            }
          }
          siteQualityMetrics {
            BaseQRankSum
            ClippingRankSum
            DP
            FS
            InbreedingCoeff
            MQ
            MQRankSum
            pab_max
            QD
            ReadPosRankSum
            RF
            SiteQuality
            SOR
            VQSLOD
          }
        }
        reads {
          bamPath
          category
          indexPath
          readGroup
        }
      }
      genome {
        ac
        an
        ac_hemi
        ac_hom
        faf95 {
          popmax
          popmax_population
        }
        filters
        populations {
          id
          ac
          an
          ac_hemi
          ac_hom
          subpopulations {
            id
            ac
            an
            ac_hom
          }
        }
        qualityMetrics {
          alleleBalance {
            alt {
              bin_edges
              bin_freq
              n_smaller
              n_larger
            }
          }
          genotypeDepth {
            all {
              bin_edges
              bin_freq
              n_smaller
              n_larger
            }
            alt {
              bin_edges
              bin_freq
              n_smaller
              n_larger
            }
          }
          genotypeQuality {
            all {
              bin_edges
              bin_freq
              n_smaller
              n_larger
            }
            alt {
              bin_edges
              bin_freq
              n_smaller
              n_larger
            }
          }
          siteQualityMetrics {
            BaseQRankSum
            ClippingRankSum
            DP
            FS
            InbreedingCoeff
            MQ
            MQRankSum
            pab_max
            QD
            ReadPosRankSum
            RF
            SiteQuality
            SOR
            VQSLOD
          }
        }
        reads {
          bamPath
          category
          indexPath
          readGroup
        }
      }
      flags
      rsid
      sortedTranscriptConsequences {
        canonical
        gene_id
        gene_symbol
        hgvs
        hgvsc
        hgvsp
        lof
        lof_flags
        lof_filter
        lof_info
        major_consequence
        polyphen_prediction
        sift_prediction
        transcript_id
      }
    }
  }
}&operationName=GnomadVariant&variables=
"""

VARIANT_QUERY_PARAMS = '"datasetId": "{dataset_id}", "variantId": "{variant_id}"'

GENE_ID_VARIANTS_QUERY_VIEW = """
query GnomadGene($geneId: String!, $datasetId: DatasetsSupportingFetchVariantsByGene!) {
  gene(gene_id: $geneId) {
    chrom
    gene_id
    canonical_transcript
    gene_name
    gene_name_upper
    other_names

    clinvar_variants {
      alleleId
      clinicalSignificance
      goldStars
      majorConsequence
      pos
      variantId
    }
    variants(dataset: $datasetId) {
      consequence
      isCanon: consequence_in_canonical_transcript
      flags
      hgvs
      hgvsc
      hgvsp
      pos
      rsid
      variant_id: variantId
      xpos
      exome {
        ac
        ac_hemi
        ac_hom
        an
        af
        filters
        populations {
          id
          ac
          an
          ac_hemi
          ac_hom
        }
      }
      genome {
        ac
        ac_hemi
        ac_hom
        an
        af
        filters
        populations {
          id
          ac
          an
          ac_hemi
          ac_hom
        }
      }
    }
  }
}&operationName=GnomadGene&variables=
"""

GENE_ID_VARIANTS_QUERY_PARAMS = '"datasetId": "{dataset_id}", "geneId": "{gene_id}"'

GENE_NAME_VARIANTS_QUERY_VIEW = """
query GnomadGene($geneName: String!, $datasetId: DatasetsSupportingFetchVariantsByGene!) {
  gene(gene_name: $geneName) {
    chrom
    gene_id
    canonical_transcript
    gene_name
    gene_name_upper
    other_names

    clinvar_variants {
      alleleId
      clinicalSignificance
      goldStars
      majorConsequence
      pos
      variantId
    }
    variants(dataset: $datasetId) {
      consequence
      isCanon: consequence_in_canonical_transcript
      flags
      hgvs
      hgvsc
      hgvsp
      pos
      rsid
      variant_id: variantId
      xpos
      exome {
        ac
        ac_hemi
        ac_hom
        an
        af
        filters
        populations {
          id
          ac
          an
          ac_hemi
          ac_hom
        }
      }
      genome {
        ac
        ac_hemi
        ac_hom
        an
        af
        filters
        populations {
          id
          ac
          an
          ac_hemi
          ac_hom
        }
      }
    }
  }
}&operationName=GnomadGene&variables=
"""

GENE_NAME_VARIANTS_QUERY_PARAMS = '"datasetId": "{dataset_id}", "geneName": "{gene_name}"'


TRANSCRIPT_EXONS_QUERY_VIEW = """
query GnomadTranscript($transcriptId: String!)
{
    transcript(transcript_id: $transcriptId) {
      gene_id
      gene_name
      transcript_id
      start
      stop
      xstart
      xstop
      chrom
      strand
      exons{
          feature_type
          start
          stop
          strand
    }
  }
}&operationName=GnomadTranscript&variables=
"""

TRANSCRIPT_EXONS_QUERY_PARAMS = '"transcriptId": "{transcript_id}"'

VARIANT_GNOMAD_3_QUERY_VIEW = """
query GnomadVariant($variantId: String!, $datasetId: DatasetId!, $referenceGenome: ReferenceGenomeId!) {
  variant(variantId: $variantId, dataset: $datasetId) {
    variantId
    reference_genome
    chrom
    pos
    ref
    alt
    colocatedVariants
    multiNucleotideVariants {
      combined_variant_id
      changes_amino_acids
      n_individuals
      other_constituent_snvs
    }
    exome {
      ac
      an
      ac_hemi
      ac_hom
      faf95 {
        popmax
        popmax_population
      }
      filters
      populations {
        id
        ac
        an
        ac_hemi
        ac_hom
      }
      age_distribution {
        het {
          bin_edges
          bin_freq
          n_smaller
          n_larger
        }
        hom {
          bin_edges
          bin_freq
          n_smaller
          n_larger
        }
      }
      qualityMetrics {
        alleleBalance {
          alt {
            bin_edges
            bin_freq
            n_smaller
            n_larger
          }
        }
        genotypeDepth {
          all {
            bin_edges
            bin_freq
            n_smaller
            n_larger
          }
          alt {
            bin_edges
            bin_freq
            n_smaller
            n_larger
          }
        }
        genotypeQuality {
          all {
            bin_edges
            bin_freq
            n_smaller
            n_larger
          }
          alt {
            bin_edges
            bin_freq
            n_smaller
            n_larger
          }
        }
        siteQualityMetrics {
          BaseQRankSum
          ClippingRankSum
          DP
          FS
          InbreedingCoeff
          MQ
          MQRankSum
          pab_max
          QD
          ReadPosRankSum
          RF
          SiteQuality
          SOR
          VQSLOD
        }
      }
    }
    genome {
      ac
      an
      ac_hemi
      ac_hom
      faf95 {
        popmax
        popmax_population
      }
      filters
      populations {
        id
        ac
        an
        ac_hemi
        ac_hom
      }
      age_distribution {
        het {
          bin_edges
          bin_freq
          n_smaller
          n_larger
        }
        hom {
          bin_edges
          bin_freq
          n_smaller
          n_larger
        }
      }
      qualityMetrics {
        alleleBalance {
          alt {
            bin_edges
            bin_freq
            n_smaller
            n_larger
          }
        }
        genotypeDepth {
          all {
            bin_edges
            bin_freq
            n_smaller
            n_larger
          }
          alt {
            bin_edges
            bin_freq
            n_smaller
            n_larger
          }
        }
        genotypeQuality {
          all {
            bin_edges
            bin_freq
            n_smaller
            n_larger
          }
          alt {
            bin_edges
            bin_freq
            n_smaller
            n_larger
          }
        }
        siteQualityMetrics {
          BaseQRankSum
          ClippingRankSum
          DP
          FS
          InbreedingCoeff
          MQ
          MQRankSum
          pab_max
          QD
          ReadPosRankSum
          RF
          SiteQuality
          SOR
          VQSLOD
        }
      }
    }
    flags
    rsid
    sortedTranscriptConsequences {
      canonical
      gene_id
      gene_symbol
      hgvs
      hgvsc
      hgvsp
      lof
      lof_flags
      lof_filter
      major_consequence
      polyphen_prediction
      sift_prediction
      transcript_id
    }
  }
  clinvar_variant(variant_id: $variantId, reference_genome: $referenceGenome) {
    clinvar_variation_id
  }
}&operationName=GnomadVariant&variables=
"""

VARIANT_GNOMAD_3_QUERY_PARAMS = '"variantId": "{variant_id}", "datasetId": "{dataset_id}", "referenceGenome": "{build}"'


def get_variant_details(variant_id, dataset_id='gnomad_r2_1', timeout=None):
  '''
  GET variant details from gnomAD (json, see VARIANT_QUERY_VIEW for structure)

  Parameters
  ----------
  variant_id: str 
    Variant identifier in format CHROM-POS-REF-ALT (e.g. 1-55505463-C-T)
  dataset_id: str (default 'gnomad_r2_1')
    GnomAD dataset id, has to be one of the following:
    gnomad_r2_1
    gnomad_r2_1_controls
    gnomad_r2_1_non_neuro
    gnomad_r2_1_non_cancer
    gnomad_r2_1_non_topmed
  '''
  variant_query = VARIANT_QUERY_PARAMS.format(dataset_id=dataset_id, variant_id=variant_id)
  response = requests.post(GNOMAD_API_URL + VARIANT_QUERY_VIEW + '{' + variant_query + '}', timeout=timeout)
  return response.json()['data']['variant']


def get_gene_variants_by_gene_id(gene_id, dataset_id='gnomad_r2_1', timeout=None):
  '''
  GET gene variants from gnomAD (json, see GENE_NAME_VARIANTS_QUERY_VIEW for structure)

  Parameters
  ----------
  gene_id: str 
    Ensembl gene id from build 37 (e.g. "ENSG00000169174")
  dataset_id: str (default 'gnomad_r2_1')
    GnomAD dataset id, has to be one of the following:
    gnomad_r2_1
    gnomad_r2_1_controls
    gnomad_r2_1_non_neuro
    gnomad_r2_1_non_cancer
    gnomad_r2_1_non_topmed
  '''
  gene_variants_query = GENE_ID_VARIANTS_QUERY_PARAMS.format(dataset_id=dataset_id, gene_id=gene_id)
  response = requests.post(GNOMAD_API_URL + GENE_ID_VARIANTS_QUERY_VIEW + '{' + gene_variants_query + '}', timeout=timeout)
  return response.json()['data']['gene']


def get_gene_variants_by_gene_name(gene_name, dataset_id='gnomad_r2_1', timeout=None):
  '''
  GET gene variants from gnomAD (json, see GENE_VARIANTS_QUERY_VIEW for structure)

  Parameters
  ----------
  gene_name: str 
    Gene name, preferably from Ensembl build 37 (e.g. "ELN")
  dataset_id: str (default 'gnomad_r2_1')
    GnomAD dataset id, has to be one of the following:
    gnomad_r2_1
    gnomad_r2_1_controls
    gnomad_r2_1_non_neuro
    gnomad_r2_1_non_cancer
    gnomad_r2_1_non_topmed
  '''
  gene_variants_query = GENE_NAME_VARIANTS_QUERY_PARAMS.format(dataset_id=dataset_id, gene_name=gene_name)
  response = requests.post(GNOMAD_API_URL + GENE_NAME_VARIANTS_QUERY_VIEW + '{' + gene_variants_query + '}', timeout=timeout)

  return response.json()['data']['gene']


def get_transcript_exons(transcript_id, timeout=None):
  '''
  GET transcript exons data (json, see TRANSCRIPT_EXONS_QUERY_VIEW for structure)

  Parameters
  ----------
  transcript_id: str 
    Ensembl transcript id from build 37 (e.g. "ENST00000302118")
  '''
  transcript_exons_query = TRANSCRIPT_EXONS_QUERY_PARAMS.format(transcript_id=transcript_id)
  response = requests.post(GNOMAD_API_URL + TRANSCRIPT_EXONS_QUERY_VIEW + '{' + transcript_exons_query + '}', timeout=timeout)

  return response.json()['data']['transcript']


def get_variant_details_gnomad_3(variant_id, dataset_id='gnomad_r3', build='GRCh38', timeout=None):
  '''
  GET variant details from gnomAD (json, see VARIANT_QUERY_VIEW for structure)

  Parameters
  ----------
  variant_id: str 
    Variant identifier in format CHROM-POS-REF-ALT (e.g. 1-55505463-C-T)
  dataset_id: str (default 'gnomad_r2_1')
    GnomAD dataset id, has to be one of the following:
    gnomad_r3
  '''
  variant_query = VARIANT_GNOMAD_3_QUERY_PARAMS.format(dataset_id=dataset_id, variant_id=variant_id, build=build)
  response = requests.post(GNOMAD_API_URL + VARIANT_GNOMAD_3_QUERY_VIEW + '{' + variant_query + '}', timeout=timeout)
  return response.json()['data']['variant']


def main():
  # EXAMPLES
  #print get_variant_details('1-55505463-C-T')
  #print get_gene_variants_by_gene_id('ENSG00000169174')
  #print get_gene_variants_by_gene_name('ELN')

  #print get_variant_details_gnomad_3("1-55063514-G-A", timeout=None)
  pass

if __name__ == "__main__":
  sys.exit(main())