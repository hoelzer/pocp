#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Author: hoelzer.martin@gmail.com

// terminal prints
println " "
println "\u001B[32mProfile: $workflow.profile\033[0m"
println " "
println "\033[2mCurrent User: $workflow.userName"
println "Nextflow-version: $nextflow.version"
println "Starting time: $nextflow.timestamp"
println "Workdir location:"
println "  $workflow.workDir\u001B[0m"
println " "
if (params.help) { exit 0, helpMSG() }
if (params.genomes == '') {exit 1, "input missing, use [--genomes]"}

// genome fasta input & --list support
if (params.genomes && params.list) { genome_input_ch = Channel
  .fromPath( params.genomes, checkIfExists: true )
  .splitCsv()
  .map { row -> [row[0], file("${row[1]}", checkIfExists: true)] }
  //.view() 
  }
  else if (params.genomes) { genome_input_ch = Channel
    .fromPath( params.genomes, checkIfExists: true)
    .map { file -> tuple(file.baseName, file) }
}

// load modules
include {prokka} from './modules/prokka'
include {blast} from './modules/blast'

// main workflow
workflow {
    proteins_ch = prokka(genome_input_ch).proteins
    proteins_ch.view()
    allvsall_ch = proteins_ch.combine(proteins_ch)
    blast(allvsall_ch)
}

// --help
def helpMSG() {
    c_green = "\033[0;32m";
    c_reset = "\033[0m";
    c_yellow = "\033[0;33m";
    c_blue = "\033[0;34m";
    c_dim = "\033[2m";
    log.info """
    ____________________________________________________________________________________________

    P.O.C.P - calculate percentage of conserved proteins

    A prokaryotic genus can be defined as a group of species with all pairwise POCP values higher than 50%.    
    
    ${c_yellow}Usage example:${c_reset}
    nextflow run hoelzer/pocp -r 0.0.1 --genomes '*.fasta' 

    ${c_yellow}Input:${c_reset}
    ${c_green} --genomes ${c_reset}           '*.fasta'         -> one strain per file
    ${c_dim}  ..change above input to csv:${c_reset} ${c_green}--list ${c_reset}   

    ${c_yellow}Options:${c_reset}
    --gcode             genetic code for Prokka annotation [default: $params.gcode]
    --evalue            Evalue for blastp [default: $params.evalue]
    --seqidentity       Sequence identity for blastp [default: $params.seqidentity]
    --alnlength         Alignment length for blastp [default: $params.alnlength]
    --cores             max cores per process for local use [default: $params.cores]
    --max_cores         max cores (in total) for local use [default: $params.max_cores]
    --memory            max memory for local use [default: $params.memory]
    --output            name of the result folder [default: $params.output]

    ${c_dim}Nextflow options:
    -with-report rep.html    cpu / ram usage (may cause errors)
    -with-dag chart.html     generates a flowchart for the process tree
    -with-timeline time.html timeline (may cause errors)

    ${c_yellow}Caching:${c_reset}
    --condaCacheDir         Location for storing the conda environments [default: $params.condaCacheDir]
    --singularityCacheDir   Location for storing the Singularity images [default: $params.condaCacheDir]
    -w                      Working directory for all intermediate results [default: work] 

    ${c_yellow}Execution/Engine profiles:${c_reset}
    The pipeline supports profiles to run via different ${c_green}Executers${c_reset} and ${c_blue}Engines${c_reset} e.g.: -profile ${c_green}local${c_reset},${c_blue}conda${c_reset}
    
    ${c_green}Executer${c_reset} (choose one):
      local
      slurm
    
    ${c_blue}Engines${c_reset} (choose one):
      conda
    
    Per default: -profile local,conda is executed.
    
    ${c_reset}
    """.stripIndent()
}