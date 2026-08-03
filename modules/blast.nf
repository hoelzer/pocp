/*
Build the BLAST database once per genome. The database is then re-used for all
pairwise comparisons against that genome instead of being rebuilt for every single
pair (N instead of N*(N-1) database builds).
*/
process blast_makedb {
    label 'blast'

    input:
      tuple val(name), path(fasta)

    output:
      tuple val(name), path("${fasta}.p*"), emit: db

    script:
    """
    makeblastdb -in ${fasta} -dbtype prot #-parse_seqids
    """
}

/*
Run blastp and format output for downstream POCP calculations.
*/
process blast {
    label 'blast'
    publishDir "${params.output}/blast", mode: 'copy', pattern: "*.blast", enabled: params.keep_alignments

    input:
      tuple val(pair_id), val(name), path(fasta), val(name2), path(fasta2), path(index)

    output:
      tuple val(pair_id), path(fasta), path("${name}-query-${name2}-db.blast.hits"), emit: hits
      tuple val(pair_id), path("${name}-query-${name2}-db.blast"), emit: alignments

    script:
    """
    blastp -task blastp -num_threads ${task.cpus} -query ${fasta} -db ${fasta2} -evalue ${params.evalue} -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend qlen sstart send evalue bitscore slen" | awk '{if(\$3>${params.seqidentity*100} && \$4>(\$9*${params.alnlength})){print \$0}}' > ${name}-query-${name2}-db.blast
    awk '{print \$1}' ${name}-query-${name2}-db.blast | sort -u | wc -l | tr -d '[:space:]' > ${name}-query-${name2}-db.blast.hits
    echo "\t${name}:\tFound \$(cat ${name}-query-${name2}-db.blast.hits) matches with an E value of less than ${params.evalue}, a sequence identity of more than ${params.seqidentity*100}%, and an alignable region of the query protein sequence of more than ${params.alnlength*100}%."
    """
}

/* Comments:
I removed the -parse_seqids parameter from the makeblastdb command because of an error with fasta IDs that are longer than 50 chars. strange.
*/
