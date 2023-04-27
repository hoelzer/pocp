/*
Run blastp and format output for downstream POCP calculations.
*/
process blast {
    label 'blast'
    publishDir "${params.output}/blast", mode: 'copy', pattern: "${fasta.baseName}.blast"

    input:
      tuple val(name), path(fasta), val(name2), path(fasta2) 
    
    output:
	    tuple val(name), file("${name}.blast") 
    
    script:
    """
    ruby pocp.rb ${fasta} ${fasta2} ${params.evalue} ${params.seqidentity} ${params.alnlength} ./ ${task.cpus}
    """
}

/* Comments:
I removed the -parse_seqids parameter from the makeblastdb command because of an error with fasta IDs that are longer than 50 chars. strange.

    %makeblastdb -in ${fasta} -dbtype prot %-parse_seqids
    %blastp -task blastp -num_threads ${tasks.cpus} -query ${fasta2} -db ${fasta} -evalue ${params.evalue} -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend qlen sstart send evalue bitscore slen" | awk '{if(\$3>${params.seqidentity*100} && \$4>(\$9*${params.alnlength})){print \$0}}' > ${name}.blast`
    %c=`awk '{print $1}' #{out_dir}/#{species1_bn}.blast | sort | uniq | wc -l`.split(' ')[0].to_i
    %puts "\t#{species1_bn}:\tFound #{c} matches with an E value of less than #{EVALUE}, a sequence identity of more than #{SEQ_IDENTITY*100}%, and an alignable region of the query protein sequence of more than #{ALN_LENGTH*100}%."
*/