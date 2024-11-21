include { FILE_GUNZIP as FASTA_GUNZIP   } from '../../subworkflows/local/file_gunzip'
include { ORTHOFINDER                   } from '../../modules/nf-core/orthofinder/main'

workflow FASTA_ORTHOFINDER {
    take:
    ch_pep_fasta                // [ meta, fasta ]
    ch_external_pep_fasta       // [ meta, fasta ]

    main:
    ch_versions                 = Channel.empty()

    // SUBWORKFLOW: FILE_GUNZIP as FASTA_GUNZIP
    FASTA_GUNZIP ( ch_external_pep_fasta )

    ch_fasta_unzipped           = FASTA_GUNZIP.out.gunzip
    ch_versions                 = ch_versions.mix(FASTA_GUNZIP.out.versions)

    // MODULE: ORTHOFINDER
    ch_orthofinder_peps         = ch_fasta_unzipped
                                | map { meta, fasta -> fasta }
                                | mix(
                                    ch_pep_fasta.map { meta, fasta -> fasta }
                                )
                                | collect
                                | filter { it.size() > 1 }

    ORTHOFINDER(
        ch_orthofinder_peps.map { fastas -> [ [ id: 'genepal' ], fastas ] },
        [ [], [] ]
    )

    ch_orthofinder_statistics   = ORTHOFINDER.out.statistics
    ch_orthofinder_hogs         = ORTHOFINDER.out.hogs
    ch_versions                 = ch_versions.mix(ORTHOFINDER.out.versions)

    emit:
    orthofinder_statistics      = ch_orthofinder_statistics // [ val(meta), Comparative_Genomics_Statistics ]
    orthofinder_hogs            = ch_orthofinder_hogs       // [ val(meta), Phylogenetic_Hierarchical_Orthogroups ]
    versions                    = ch_versions               // [ versions.yml ]
}
