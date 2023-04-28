/*
Calculate POCP value:     
  pocp = (((c1.to_f + c2)/(t1 + t2)) * 100).round(4)
  puts "\n\nC1:\t#{c1}\nC2:\t#{c2}\n\nT1:\t#{t1}\nT2:\t#{t2}\n\nPOCP = [(C1+C2)/(T1+T2)]*100% = #{pocp}"
*/
process pocp {
    label 'blast'
    publishDir "${params.output}/pocp", mode: 'copy', pattern: "${genome_names}.txt"

    input:
      tuple val(genome_names), path(proteins), path(blast_hits) 
    
    output:
	    tuple val(genome_names), path("${genome_names}.txt"), env(POCP), emit: pocp
    
    script:
    """
    t1=\$(grep ">" ${proteins[0]} | wc -l)
    t2=\$(grep ">" ${proteins[1]} | wc -l)
    c1=\$(cat ${blast_hits[0]})
    c2=\$(cat ${blast_hits[1]})
    POCP=\$(python -c "print((((\$c1 + \$c2) / (\$t1 + \$t2)) * 100))")
    echo \$POCP > ${genome_names}.txt
    """
}
