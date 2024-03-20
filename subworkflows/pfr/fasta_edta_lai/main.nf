include { CUSTOM_SHORTENFASTAIDS    } from '../../../modules/pfr/custom/shortenfastaids/main'
include { EDTA_EDTA                 } from '../../../modules/pfr/edta/edta/main'
include { LTRRETRIEVER_LAI          } from '../../../modules/pfr/ltrretriever/lai/main'
include { CUSTOM_RESTOREGFFIDS      } from '../../../modules/pfr/custom/restoregffids/main'

workflow FASTA_EDTA_LAI {

    take:
    ch_fasta                        // channel: [ val(meta), fasta ]
    ch_monoploid_seqs               // channel: [ val(meta), txt ]; Optional: Set to [] if not needed
    skip_lai                        // val; true|false

    main:
    ch_versions                     = Channel.empty()

    // MOUDLE: CUSTOM_SHORTENFASTAIDS
    CUSTOM_SHORTENFASTAIDS ( ch_fasta )

    ch_short_ids_fasta              = ch_fasta
                                    | join(CUSTOM_SHORTENFASTAIDS.out.short_ids_fasta, by:0, remainder:true)
                                    | map { meta, fasta, short_ids_fasta ->
                                        if ( fasta ) { [ meta, short_ids_fasta ?: fasta ] }
                                    }

    ch_short_ids_tsv                = CUSTOM_SHORTENFASTAIDS.out.short_ids_tsv
    ch_short_monoploid_seqs         = ch_short_ids_tsv
                                    | join(
                                        ch_monoploid_seqs ?: Channel.empty()
                                    )
                                    | map { meta, short_ids_tsv, monoploid_seqs ->
                                        map_monoploid_seqs_to_new_ids(meta, short_ids_tsv, monoploid_seqs)
                                    }
                                    | collectFile(newLine:true)
                                    | map { seqs ->
                                        def id = seqs.name.split('.mapped.monoploid.seqs.txt')[0]

                                        [ [ id: id ], seqs ]
                                    }
    ch_versions                     = ch_versions.mix(CUSTOM_SHORTENFASTAIDS.out.versions.first())

    // MODULE: EDTA_EDTA
    EDTA_EDTA(
        ch_short_ids_fasta,
        [],
        [],
        [],
        []
    )

    ch_te_lib_fasta                 = EDTA_EDTA.out.te_lib_fasta
    ch_pass_list                    = EDTA_EDTA.out.pass_list
    ch_out_file                     = EDTA_EDTA.out.out_file
    ch_te_anno_gff3                 = EDTA_EDTA.out.te_anno_gff3
    ch_versions                     = ch_versions.mix(EDTA_EDTA.out.versions.first())

    // MODULE: LAI
    ch_lai_inputs                   = skip_lai
                                    ? Channel.empty()
                                    : ch_short_ids_fasta
                                    | join(ch_pass_list)
                                    | join(ch_out_file)
                                    | map { meta, fasta, pass, out ->
                                        [ meta.id, meta, fasta, pass, out ]
                                    }
                                    | join(
                                        ch_short_monoploid_seqs
                                        | map { meta, mono -> [ meta.id, mono ] },
                                        by:0,
                                        remainder: true
                                    )
                                    | map { id, meta, fasta, pass, out, mono ->
                                        [ meta, fasta, pass, out, mono ?: [] ]
                                    }
    LTRRETRIEVER_LAI(
        ch_lai_inputs.map { meta, fasta, pass, out, mono -> [ meta, fasta ] },
        ch_lai_inputs.map { meta, fasta, pass, out, mono -> pass },
        ch_lai_inputs.map { meta, fasta, pass, out, mono -> out },
        ch_lai_inputs.map { meta, fasta, pass, out, mono -> mono }
    )

    ch_lai_log                      = LTRRETRIEVER_LAI.out.log
    ch_lai_out                      = LTRRETRIEVER_LAI.out.lai_out
    ch_versions                     = ch_versions.mix(LTRRETRIEVER_LAI.out.versions.first())

    // MODULE: CUSTOM_RESTOREGFFIDS
    ch_restorable_gff_tsv           = ch_te_anno_gff3.join(ch_short_ids_tsv)

    CUSTOM_RESTOREGFFIDS (
        ch_restorable_gff_tsv.map { meta, gff, tsv -> [ meta, gff ] },
        ch_restorable_gff_tsv.map { meta, gff, tsv -> tsv }
    )

    ch_restored_gff                 = ch_te_anno_gff3
                                    | join(CUSTOM_RESTOREGFFIDS.out.restored_ids_gff3, by:0, remainder:true)
                                    | map { meta, gff, restored_gff -> [ meta, restored_gff ?: gff ] }
    ch_versions                     = ch_versions.mix(CUSTOM_RESTOREGFFIDS.out.versions.first())

    emit:
    te_lib_fasta                    = ch_te_lib_fasta   // channel: [ val(meta), fasta ]
    te_anno_gff3                    = ch_restored_gff   // channel: [ val(meta), gff ]
    lai_log                         = ch_lai_log        // channel: [ val(meta), log ]
    lai_out                         = ch_lai_out        // channel: [ val(meta), out ]
    versions                        = ch_versions       // channel: [ versions.yml ]
}

def map_monoploid_seqs_to_new_ids(meta, short_ids_tsv, monoploid_seqs) {

    def short_ids_head = short_ids_tsv.text.split('\n')[0]

    if (short_ids_head == "IDs have acceptable length and character. No change required.") {
        return [ "${meta.id}.mapped.monoploid.seqs.txt" ] + monoploid_seqs.text.split('\n')
    }

    def orig_to_new_ids = [:]
    short_ids_tsv.text.eachLine { line ->
        def (original_id, renamed_id) = line.split('\t')
        orig_to_new_ids[original_id] = renamed_id
    }

    def mapped_ids = []
    monoploid_seqs.text.eachLine { original_id ->
        if (!orig_to_new_ids[original_id]) {
            error "Faild to find $original_id in ${monoploid_seqs}" +
            "The monoploid_seqs file is malformed!"
        }

        mapped_ids.add(orig_to_new_ids[original_id])
    }

    return [ "${meta.id}.mapped.monoploid.seqs.txt" ] + mapped_ids
}
