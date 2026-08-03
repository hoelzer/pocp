#!/usr/bin/env nextflow

// Author: hoelzer.martin@gmail.com

// load modules
include { prokka; prokka as prokka_single } from './modules/prokka'
include { blast_makedb; blast } from './modules/blast'
include { diamond_makedb; diamond } from './modules/diamond'
include { pocp; pocp_matrix } from './modules/pocp'
include { plot } from './modules/plot'

// main workflow
workflow {
    printIntro()

    if( params.help ) {
        log.info(helpMSG())
        System.exit(0)
    }

    validateParams()
    printParams()

    // one-vs-all needs at least one input in --genomes/--proteins, all-vs-all at least two
    def min_inputs = (params.genome || params.protein) ? 1 : 2

    // annotate genomes if needed, otherwise use the provided protein FASTA files
    def proteins_ch = params.genomes
        ? prokka(checkIds(readInput(params.genomes), min_inputs)).proteins
        : checkIds(readInput(params.proteins), min_inputs)

    // optional single target to switch from all-vs-all to one-vs-all
    def target_ch = channel.empty()
    if( params.genome )
        target_ch = prokka_single(readSingle(params.genome)).proteins
    else if( params.protein )
        target_ch = readSingle(params.protein)

    // pairwise comparisons, both directions are needed for the POCP formula
    def comparisons_ch = null
    if( params.genome || params.protein ) {
        // one-vs-all
        def one_way_ch = proteins_ch
            .combine(target_ch)
            .filter { id1, _faa1, id2, _faa2 -> id1 != id2 }
        comparisons_ch = one_way_ch
            .mix( one_way_ch.map { id1, faa1, id2, faa2 -> tuple(id2, faa2, id1, faa1) } )
    }
    else {
        // all-vs-all, combine() already yields both directions of each pair
        comparisons_ch = proteins_ch
            .combine(proteins_ch)
            .filter { id1, _faa1, id2, _faa2 -> id1 != id2 }
    }

    // key each comparison by the sorted pair of IDs, so that both directions group together
    def pairs_ch = comparisons_ch
        .map { id1, faa1, id2, faa2 -> tuple(id2, [id1, id2].sort().join('-vs-'), id1, faa1, faa2) }

    // build the search database once per genome instead of once per comparison
    def db_input_ch = proteins_ch.mix(target_ch).unique { entry -> entry[0] }

    // use either BLASTP or DIAMOND (default)
    def hits_ch = null
    if( params.blastp ) {
        hits_ch = blast(
            pairs_ch
                .combine( blast_makedb(db_input_ch).db, by: 0 )
                .map { id2, pair_id, id1, faa1, faa2, index -> tuple(pair_id, id1, faa1, id2, faa2, index) }
        ).hits
    }
    else {
        hits_ch = diamond(
            pairs_ch
                .combine( diamond_makedb(db_input_ch).db, by: 0 )
                .map { id2, pair_id, id1, faa1, _faa2, db -> tuple(pair_id, id1, faa1, id2, db) }
        ).hits
    }

    // exactly two hit counts per pair, one per direction
    def pocp_ch = pocp( hits_ch.groupTuple(size: 2) )

    pocp_matrix( pocp_ch.pocp.map { _pair_id, pocp_file -> pocp_file }.collect() )

    plot( pocp_matrix.out )
}

// read a FASTA glob or, with --list, a CSV of "id,path" rows
def readInput(String input) {
    return params.list
        ? channel.fromPath(input, checkIfExists: true)
            .splitCsv()
            .map { row -> tuple(row[0], file(row[1], checkIfExists: true)) }
        : channel.fromPath(input, checkIfExists: true)
            .map { fasta -> tuple(fasta.baseName, fasta) }
}

// read the single FASTA used for one-vs-all comparisons
def readSingle(String input) {
    return channel.fromPath(input, checkIfExists: true)
        .map { fasta -> tuple(fasta.baseName, fasta) }
}

/*
Sample IDs are derived from file base names and are used to build the "A-vs-B"
keys that group both alignment directions and that name the matrix rows and
columns. Catch ambiguous IDs up front instead of silently dropping comparisons.
*/
def checkIds(ch, int min_inputs) {
    return ch.toList().flatMap { entries ->
        def ids = entries.collect { entry -> entry[0] }

        if( ids.size() < min_inputs )
            error "Found ${ids.size()} input file(s), but at least ${min_inputs} are needed for a POCP calculation."

        def duplicates = ids.countBy { id -> id }.findAll { _id, count -> count > 1 }.keySet()
        if( duplicates )
            error "Input file base names must be unique, but these IDs occur more than once: ${duplicates.join(', ')}"

        def invalid = ids.findAll { id -> id.contains('-vs-') || id ==~ /.*\s.*/ }
        if( invalid )
            error "Input file base names must not contain whitespace or '-vs-': ${invalid.join(', ')}"

        return entries
    }
}

def validateParams() {
    if( params.profile )
        error "--profile is WRONG, use -profile"
    if( !params.genomes && !params.proteins )
        error "input missing, use either [--genomes] or [--proteins]"
    if( params.genomes && params.proteins )
        error "provide one input, use either [--genomes] or [--proteins]"
    if( params.genome && params.protein )
        error "provide only one input, use either [--genome] or [--protein]"
}

// terminal prints
def printIntro() {
    println " "
    println "\033[32mProfile: $workflow.profile\033[0m"
    println " "
    println "\033[2mCurrent User: $workflow.userName"
    println "Nextflow-version: $nextflow.version"
    println "Starting time: $nextflow.timestamp"
    println "Workdir location:"
    println "  $workflow.workDir\033[0m"
    println " "

    if( !workflow.revision ) {
        println "\033[0;33mWARNING: It is recommended to use a stable release version via -r."
        println "Use 'nextflow info hoelzer/pocp' to check for available release versions.\033[0m\n"
    }
}

def printParams() {
    def tool = params.blastp ? "BLASTP" : "DIAMOND"
    def version = workflow.revision ? "in version '${workflow.revision}'" : "\033[0;31mwithout a stable release version\033[0m\033[32m"

    // print if default params are used
    if( params.evalue == '1e-5' && params.seqidentity == 0.4 && params.alnlength == 0.5 ) {
        println "\033[32mPOCP-nf was executed ${version} and default paramters according to the original publication by Qin et al. (2014). ${tool} was used for protein alignments."
        println ""
        println "e-value:\t\t${params.evalue}"
        println "Sequence identity:\t${params.seqidentity}"
        println "Alignment length:\t${params.alnlength}"
        println "\033[0m"
    }
    // print if NO default params are used
    else {
        def revision = workflow.revision ? "in version '${workflow.revision}'" : "without a stable release version"
        println "\033[0;31mPOCP-nf was executed ${revision} and non-default paramters in comparison to the original publication by Qin et al. (2014). ${tool} was used for protein alignments."
        println ""
        println "e-value used:\t\t${params.evalue}\t(original definition: 1e-5)"
        println "Sequence identity used:\t${params.seqidentity}\t(original definition: 0.4)"
        println "Alignment length used:\t${params.alnlength}\t(original definition: 0.5)"
        println ""
        println "This might change your POCP results slightly."
        println "If you really want to use adjusted parameters, you must report them together with the used version of POCP-nf to ensure reproducibility!"
        println ""
        println "\033[0m"
    }
}

// --help
def helpMSG() {
    def c_green = "\033[0;32m"
    def c_reset = "\033[0m"
    def c_yellow = "\033[0;33m"
    def c_blue = "\033[0;34m"
    def c_red = "\033[0;31m"
    def c_dim = "\033[2m"
    return """
    ____________________________________________________________________________________________

    P.O.C.P - calculate percentage of conserved proteins.

    A prokaryotic genus can be defined as a group of species with all pairwise POCP values higher than 50%.

    ${c_yellow}Usage example:${c_reset}
    nextflow run hoelzer/pocp -r 2.3.4 --genomes '*.fasta'
    or
    nextflow run hoelzer/pocp -r 2.3.4 --proteins '*.faa'

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
    --width             Width in inches for the heatmap POCP plot [default: $params.width]
    --heigth            Height in inches for the heatmap POCP plot [default: $params.height]
    --keep_alignments   Also publish the raw pairwise alignment tables, one per comparison [default: $params.keep_alignments]

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
    --singularityCacheDir   Location for storing the Singularity images [default: $params.singularityCacheDir]
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
