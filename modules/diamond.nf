/*
Build the DIAMOND database once per genome. The database is then re-used for all
pairwise comparisons against that genome instead of being rebuilt for every single
pair (N instead of N*(N-1) database builds).
*/
process diamond_makedb {
    label 'diamond'

    input:
      tuple val(name), path(fasta)

    output:
      tuple val(name), path("${fasta}.dmnd"), emit: db

    script:
    """
    diamond makedb --in ${fasta} -d ${fasta}.dmnd
    """

    stub:
    """
    touch ${fasta}.dmnd
    """
}

/*
Run diamond and format output for downstream POCP calculations.
https://www.nature.com/articles/s41592-021-01101-x
*/
process diamond {
    label 'diamond'
    publishDir "${params.output}/diamond", mode: 'copy', pattern: "*.diamond", enabled: params.keep_alignments

    input:
      tuple val(pair_id), val(name), path(fasta), val(name2), path(db)

    output:
      tuple val(pair_id), path(fasta), path("${name}-query-${name2}-db.diamond.hits"), emit: hits
      tuple val(pair_id), path("${name}-query-${name2}-db.diamond"), emit: alignments

    script:
    """
    diamond blastp --ultra-sensitive -p ${task.cpus} -q ${fasta} -d ${db} -e ${params.evalue} --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend qlen sstart send evalue bitscore slen | awk '{if(\$3>${params.seqidentity*100} && \$4>(\$9*${params.alnlength})){print \$0}}' > ${name}-query-${name2}-db.diamond
    awk '{print \$1}' ${name}-query-${name2}-db.diamond | sort -u | wc -l | tr -d '[:space:]' > ${name}-query-${name2}-db.diamond.hits
    echo "\t${name}:\tFound \$(cat ${name}-query-${name2}-db.diamond.hits) matches with an E value of less than ${params.evalue}, a sequence identity of more than ${params.seqidentity*100}%, and an alignable region of the query protein sequence of more than ${params.alnlength*100}%."
    """

    stub:
    """
    touch ${name}-query-${name2}-db.diamond
    echo "1" > ${name}-query-${name2}-db.diamond.hits
    """
}
