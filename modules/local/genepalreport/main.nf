process GENEPALREPORT {
    tag 'genepal'

    label 'process_single'
    cache false

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/52/5250888c07f78131a0ad6944289ce0910b0c255291c86b56e0b7b7b994fceacb/data':
        'community.wave.seqera.io/library/r-base_r-dplyr_r-gt_r-gtextras_pruned:0a84d5aec45783b9' }"

    input:
    path(splicing_marked_gff3, stageAs: 'results/etc/splicing_marked/*')
    path(pipeline_info, stageAs: 'results/pipeline_info/*')

    output:
    path("*.html")
    path("genepal_data")

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    ln -s \\
        \$(which genepal_report.Rmd)

    ln -s \\
        \$(which genepal_report.R)

    Rscript \\
        -e "rmarkdown::render('genepal_report.Rmd', output_file='genepal_report.html')"
    """

    stub:
    """
    touch genepal_report.html

    mkdir genepal_data
    touch genepal_data/summary_stats.json
    """
}
