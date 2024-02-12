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
// error codes
if (params.profile) { exit 1, "--profile is WRONG use -profile" }
if ( !workflow.revision ) { 
  println "\033[0;33mWARNING: It is recommended to use a stable relese version via -r." 
  println "Use 'nextflow info hoelzer/pocp' to check for available release versions.\033[0m\n"
}
// help
if (params.help) { exit 0, helpMSG() }
// input error codes
if (params.genomes == '' && params.proteins == '') {exit 1, "input missing, use either [--genomes] or [--proteins]"}
if (params.genomes != '' && params.proteins != '') {exit 1, "provide one input, use either [--genomes] or [--proteins]"}
if (params.genome && params.protein) {exit 1, "provide only one input, use either [--genome] or [--protein]"}
// print if default params are used:
if (params.evalue == '1e-5' && params.seqidentity == 0.4 && params.alnlength == 0.5 && !params.blastp ) {
  if ( !workflow.revision ) { 
  println "\u001B[32mPOCP-nf was executed \033[0;31mwithout a stable release version\033[0m\u001B[32m and default paramters according to the original publication by Qin et al. (2014)."
  } else {
  println "\u001B[32mPOCP-nf was executed in version '${workflow.revision}' with default paramters according to the original publication by Qin et al. (2014)."
  }
  println ""
  println "e-value:\t\t${params.evalue}"
  println "Sequence identity:\t${params.seqidentity}"
  println "Alignment length:\t${params.alnlength}"
  println "\033[0m"
  // print if NO default params are used
} else {
  if ( !workflow.revision ) { 
  println "\033[0;31mPOCP-nf was executed without a stable release version and non-default paramters in comparison to the original publication by Qin et al. (2014)."
  println ""
  println "e-value used:\t\t${params.evalue}\t(original definition: 1e-5)"
  println "Sequence identity used:\t${params.seqidentity}\t(original definition: 0.4)"
  println "Alignment length used:\t${params.alnlength}\t(original definition: 0.5)"
  } else {
  println "\033[0;31mPOCP-nf was executed in version '${workflow.revision}' with non-default paramters in comparison to the original publication by Qin et al. (2014)."
  println ""
  println "e-value used:\t\t${params.evalue}\t(original definition: 1e-5)"
  println "Sequence identity used:\t${params.seqidentity}\t(original definition: 0.4)"
  println "Alignment length used:\t${params.alnlength}\t(original definition: 0.5)"
  }
  println ""
  println "This will change your POCP results."
  println "If you really want to use adjusted parameters, you must report them together with the used version of POCP-nf to ensure reproducibility!"
  println ""
  println "\033[0m"
}

// genomes fasta input & --list support
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

// proteins fasta input & --list support
if (params.proteins && params.list) { proteins_input_ch = Channel
  .fromPath( params.proteins, checkIfExists: true )
  .splitCsv()
  .map { row -> [row[0], file("${row[1]}", checkIfExists: true)] }
  //.view() 
  }
  else if (params.proteins) { proteins_input_ch = Channel
    .fromPath( params.proteins, checkIfExists: true)
    .map { file -> tuple(file.baseName, file) }
}

// optional input for one-vs-all comparisons instead of all-vs-all
// genome fasta input
if (params.genome) { genome_single_input_ch = Channel
    .fromPath( params.genome, checkIfExists: true)
    .map { file -> tuple(file.baseName, file) }
}
// protein fasta input
if (params.protein) { protein_single_input_ch = Channel
    .fromPath( params.protein, checkIfExists: true)
    .map { file -> tuple(file.baseName, file) }
}

// load modules
include {prokka as prokka; prokka as prokka_single} from './modules/prokka'
include {blast} from './modules/blast'
include {diamond} from './modules/diamond'
include {pocp; pocp_matrix} from './modules/pocp'

// main workflow
workflow {
    if (params.genomes) {
      proteins_ch = prokka(genome_input_ch)
    } else {
      proteins_ch = proteins_input_ch
    }

    // switch to one-vs-all
    if (params.genome) { protein_ch = prokka_single(genome_single_input_ch) }
    if (params.protein) { protein_ch = protein_single_input_ch }

    if (params.genome || params.protein) {
      // one-vs-all
      one_vs_all_ch_part1 = proteins_ch.combine(protein_ch).branch { id1, faa1, id2, faa2 ->
          controls: id1 != id2
              return tuple( id1, faa1, id2, faa2 )
      }
      one_vs_all_ch_part2 = proteins_ch.combine(protein_ch).branch { id1, faa1, id2, faa2 ->
          controls: id1 != id2
              return tuple( id2, faa2, id1, faa1 )
      }
      comparisons_ch = one_vs_all_ch_part1.concat(one_vs_all_ch_part2)
    } else {
      // all-vs-all
      all_vs_all_ch = proteins_ch.combine(proteins_ch).branch { id1, faa1, id2, faa2 ->
          controls: id1 != id2
              return tuple( id1, faa1, id2, faa2 )
      }
      comparisons_ch = all_vs_all_ch
    }

    // use either BLASTP or DIAMOND (default)
    if (params.blastp) {
      hits_ch = blast(comparisons_ch).hits.groupTuple()
    } else {
      hits_ch = diamond(comparisons_ch).hits.groupTuple()
    }

    pocp_matrix(
      pocp(hits_ch).map {comparison, pocp_file, pocp_value -> [pocp_file]}.collect()
    )
}

// --help
def helpMSG() {
    c_green = "\033[0;32m";
    c_reset = "\033[0m";
    c_yellow = "\033[0;33m";
    c_blue = "\033[0;34m";
    c_red = "\033[0;31m";
    c_dim = "\033[2m";
    log.info """
    ____________________________________________________________________________________________

    P.O.C.P - calculate percentage of conserved proteins.

    A prokaryotic genus can be defined as a group of species with all pairwise POCP values higher than 50%.    
    
    ${c_yellow}Usage example:${c_reset}
    nextflow run hoelzer/pocp -r 2.3.0 --genomes '*.fasta' 
    or
    nextflow run hoelzer/pocp -r 2.3.0 --proteins '*.faa' 

    Use the following commands to check for latest pipeline versions:
    
    nextflow pull hoelzer/pocp
    nextflow info hoelzer/pocp

    ${c_yellow}Input${c_reset}
    ${c_yellow}All-vs-all comparisons (default):${c_reset}
    ${c_green} --genomes ${c_reset}           '*.fasta'         -> one genome per file
    or
    ${c_green} --proteins ${c_reset}           '*.faa'          -> one protein multi-FASTA per file
    ${c_dim}  ..change above input to csv:${c_reset} ${c_green}--list ${c_reset}   

    ${c_yellow}Perform one-vs-all comparison against the additionally defined genome or protein FASTA (optional):${c_reset}
     --genome            genome.fasta         -> one genome FASTA
    or
     --protein           proteins.faa         -> one protein multi-FASTA

    ${c_yellow}General Options:${c_reset}
    --gcode             Genetic code for Prokka annotation [default: $params.gcode]
    --cores             Max cores per process for local use [default: $params.cores]
    --max_cores         Max cores (in total) for local use [default: $params.max_cores]
    --memory            Max memory for local use [default: $params.memory]
    --output            Name of the result folder [default: $params.output]

    ${c_yellow}Special Options${c_reset} ${c_red}(Danger Zone!)${c_yellow}:${c_reset}
    ATTENTION: changing these parameters will lead to different POCP values. 
    If you have good reasons to do that, you must report the changed parameters together with the used pipeline version.
    
    --evalue            Evalue for DIAMOND protein search [default: $params.evalue]
    --seqidentity       Sequence identity for DIAMOD alignments [default: $params.seqidentity]
    --alnlength         Alignment length for DIAMOND hits [default: $params.alnlength]
    --blastp            Use BLASTP instead of DIAMOND for protein alignment (slower but as in the original 2014 publication) [default: $params.blastp]

    ${c_dim}Nextflow options:
    -with-report rep.html    cpu / ram usage (may cause errors)
    -with-dag chart.html     generates a flowchart for the process tree
    -with-timeline time.html timeline (may cause errors)
    -resume                  resume a previous calculation w/o recalculating everything (needs the same run command and work dir!)

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
      mamba
      docker
      singularity
    
    Per default: -profile local,conda is executed.
    
    ${c_reset}
    """.stripIndent()
}