include { STAR_ALIGN    } from '../../modules/nf-core/star/align'
include { SAMTOOLS_CAT  } from '../../modules/nf-core/samtools/cat'

workflow ALIGN_RNASEQ {
    take:
    reads_target                // channel: [ meta, assembly_id ]
    trim_reads                  // channel: [ meta, [ fq ] ]
    rna_bam_by_assembly         // channel: [ meta, [ bam ] ]
    assembly_index              // channel: [ meta2, star_index ]

    main:
    ch_versions                 = Channel.empty()

    // MODULE: STAR_ALIGN
    ch_star_inputs              = reads_target
                                | combine(trim_reads, by:0)
                                | map { meta, assembly, fastq ->
                                    [
                                        assembly,
                                        [
                                            id: "${meta.id}.on.${assembly}",
                                            single_end: meta.single_end,
                                            target_assembly: assembly
                                        ],
                                        fastq
                                    ]
                                }
                                | combine(
                                    assembly_index.map { meta, index -> [ meta.id, index ] },
                                    by:0
                                )
                                | map { assembly, meta, fastq, index -> [ meta, fastq, index ] }

    def star_ignore_sjdbgtf     = true
    def seq_platform            = false
    def seq_center              = false

    STAR_ALIGN(
        ch_star_inputs.map { meta, fastq, index -> [ meta, fastq ] },
        ch_star_inputs.map { meta, fastq, index -> [ [ id: meta.target_assembly ], index ] },
        ch_star_inputs.map { meta, fastq, index -> [ [ id: meta.target_assembly ], [] ] },
        star_ignore_sjdbgtf,
        seq_platform,
        seq_center
    )

    ch_star_bam                 = STAR_ALIGN.out.bam_sorted_aligned
    ch_versions                 = ch_versions.mix(STAR_ALIGN.out.versions.first())

    // MODULE: SAMTOOLS_CAT
    ch_star_bam_branch          = ch_star_bam
                                | map { meta, bam ->
                                    [
                                        [ id: meta.target_assembly ],
                                        bam instanceof List ? bam.find { it =~ /Aligned/ } : bam
                                    ]
                                }
                                | mix ( rna_bam_by_assembly )
                                | groupTuple
                                | map { meta, bams -> [ meta, bams.flatten() ] }
                                | branch { meta, bamList ->
                                    bams: bamList.size() > 1
                                    bam: bamList.size() <= 1
                                }

    SAMTOOLS_CAT(
        ch_star_bam_branch.bams.map { meta, bams -> [ meta, bams.toSorted() ] }
    )

    ch_samtools_bam             = SAMTOOLS_CAT.out.bam
                                | map { meta, bam -> [ meta, [ bam ] ] }
                                | mix(
                                    ch_star_bam_branch.bam
                                )

    ch_versions                 = ch_versions.mix(SAMTOOLS_CAT.out.versions.first())

    emit:
    bam                         = ch_samtools_bam           // channel: [ [ id: target_assembly, single_end ], [ bam ] ]
    star_log_final              = STAR_ALIGN.out.log_final  // channel: [ [ id: sample.on.assembly, single_end ], log ]
    versions                    = ch_versions               // channel: [ versions.yml ]
}
