/*Comment section: */

process prokka {
  label 'prokka'
  publishDir "${params.output}/prokka", mode: 'copy', pattern: "*/*.faa"
  publishDir "${params.output}/prokka", mode: 'copy', pattern: "*/*.gff"

  input: 
    tuple val(name), path(fasta)

  output:
    tuple val(name), path("${name}/${name}.faa"), emit: proteins
    tuple val(name), path("${name}/${name}.gff"), emit: annotations

  script:
    """
    prokka --gcode ${params.gcode} --cpus ${task.cpus} --outdir ${name} --prefix ${name} ${fasta}
    """

  stub:
    """
    mkdir -p ${name}
    printf '>${name}_00001 hypothetical protein\\nMSKV\\n' > ${name}/${name}.faa
    touch ${name}/${name}.gff
    """
}