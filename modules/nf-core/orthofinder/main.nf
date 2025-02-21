process ORTHOFINDER {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/orthofinder:2.5.5--hdfd78af_2':
        'biocontainers/orthofinder:2.5.5--hdfd78af_2' }"

    input:
    tuple val(meta), path(fastas, stageAs: 'input/')
    tuple val(meta2), path(prior_run)

    output:
    tuple val(meta), path("$prefix")                                , emit: orthofinder
    tuple val(meta), path("$prefix/WorkingDirectory")               , emit: working
    tuple val(meta), path("Comparative_Genomics_Statistics")        , emit: statistics      , optional: true
    tuple val(meta), path("Gene_Duplication_Events")                , emit: duplications    , optional: true
    tuple val(meta), path("Orthogroups")                            , emit: orthogroups     , optional: true
    tuple val(meta), path("Phylogenetic_Hierarchical_Orthogroups")  , emit: hogs            , optional: true
    tuple val(meta), path("Species_Tree")                           , emit: species_tree    , optional: true
    path "versions.yml"                                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def include_command = prior_run   ? "-b $prior_run" : ''

    """
    mkdir temp_pickle

    orthofinder \\
        -t $task.cpus \\
        -a $task.cpus \\
        -p temp_pickle \\
        -f input \\
        -n $prefix \\
        $include_command \\
        $args

    if [ -e input/OrthoFinder/Results_$prefix ]; then
        mv input/OrthoFinder/Results_$prefix $prefix
    fi

    if [ -e ${prior_run}/OrthoFinder/Results_$prefix ]; then
        mv ${prior_run}/OrthoFinder/Results_$prefix $prefix
    fi

    cp -r \\
        $prefix/Comparative_Genomics_Statistics \\
        Comparative_Genomics_Statistics \\
        || echo "Comparative_Genomics_Statistics was not produced"
    
    cp -r \\
        $prefix/Gene_Duplication_Events \\
        Gene_Duplication_Events \\
        || echo "Gene_Duplication_Events was not produced"
    
    cp -r \\
        $prefix/Orthogroups \\
        Orthogroups \\
        || echo "Orthogroups was not produced"
    
    cp -r \\
        $prefix/Phylogenetic_Hierarchical_Orthogroups \\
        Phylogenetic_Hierarchical_Orthogroups \\
        || echo "Phylogenetic_Hierarchical_Orthogroups was not produced"
    
    cp -r \\
        $prefix/Species_Tree \\
        Species_Tree \\
        || echo "Species_Tree was not produced"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        orthofinder: \$(orthofinder -h | sed -n 's/.*version \\(.*\\) Copy.*/\\1/p')
    END_VERSIONS
    """

    stub:
    def args   = task.ext.args   ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def include_command = prior_run   ? "-b $prior_run" : ''

    """
    mkdir -p    $prefix/Comparative_Genomics_Statistics
    mkdir       $prefix/Gene_Duplication_Events
    mkdir       $prefix/Gene_Trees
    mkdir       $prefix/Orthogroup_Sequences
    mkdir       $prefix/Orthogroups
    mkdir       $prefix/Orthologues
    mkdir       $prefix/Phylogenetic_Hierarchical_Orthogroups
    mkdir       $prefix/Phylogenetically_Misplaced_Genes
    mkdir       $prefix/Putative_Xenologs
    mkdir       $prefix/Resolved_Gene_Trees
    mkdir       $prefix/Single_Copy_Orthologue_Sequences
    mkdir       $prefix/Species_Tree
    mkdir       $prefix/WorkingDirectory

    touch       $prefix/Log.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        orthofinder: \$(orthofinder -h | sed -n 's/.*version \\(.*\\) Copy.*/\\1/p')
    END_VERSIONS
    """
}
