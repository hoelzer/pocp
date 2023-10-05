/*
Run blastp and format output for downstream POCP calculations.
*/
process blast {
    label 'blast'
    publishDir "${params.output}/blast", mode: 'copy', pattern: "${name}-query-${name2}-db.blast*"

    input:
      tuple val(name), path(fasta), val(name2), path(fasta2) 
    
    output:
	    tuple env(genome_ids_sorted), path(fasta), file("${name}-query-${name2}-db.blast"), emit: blast 
	    tuple env(genome_ids_sorted), path(fasta), file("${name}-query-${name2}-db.blast.hits"), emit: hits
    
    script:
    """
    diamond makedb --in ${fasta2} -d ${fasta2}.dmnd 
    diamond blastp --ultra-sensitive -p ${task.cpus} -q ${fasta} -d ${fasta2}.dmnd -e ${params.evalue} --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend qlen sstart send evalue bitscore slen | awk '{if(\$3>${params.seqidentity*100} && \$4>(\$9*${params.alnlength})){print \$0}}' > ${name}-query-${name2}-db.blast
    awk '{print \$1}' ${name}-query-${name2}-db.blast | sort | uniq | wc -l | awk '{print \$1}' > ${name}-query-${name2}-db.blast.hits
    echo "\t${name}:\tFound \$(cat ${name}-query-${name2}-db.blast.hits) matches with an E value of less than ${params.evalue}, a sequence identity of more than ${params.seqidentity*100}%, and an alignable region of the query protein sequence of more than ${params.alnlength*100}%."

    genome_ids_sorted='${name} ${name2}'
    genome_ids_sorted=\$(echo \$genome_ids_sorted | xargs -n1 | sort | xargs | sed 's/ /-vs-/g')

    #ruby pocp.rb ${fasta} ${fasta2} ${params.evalue} ${params.seqidentity} ${params.alnlength} ./ ${task.cpus}
    """
}

/* Comments:
I removed the -parse_seqids parameter from the makeblastdb command because of an error with fasta IDs that are longer than 50 chars. strange.
*/