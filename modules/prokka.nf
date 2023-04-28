/*Comment section: */

process prokka {
  label 'prokka'
  publishDir "${params.output}/prokka", mode: 'copy', pattern: "${name}/${name}.faa" 

  input: 
    tuple val(name), path(fasta)

  output:
    tuple val(name), path("${name}/${name}.faa"), emit: proteins

  script:
    """
    prokka --gcode ${params.gcode} --cpus ${task.cpus} --outdir ${name} --prefix ${name} ${fasta}
    """
}