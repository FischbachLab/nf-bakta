#!/usr/bin/env nextflow
nextflow.enable.dsl=1
// If the user uses the --help flag, print the help text below
params.help = false

// Function which prints help message text
def helpMessage() {
  log.info"""
    Genome Annotation

    Required Arguments:
      --seedfile              Location where 1 or more fasta files are stored.
      --project             Folder to place analysis outputs (default: )
      --prefix             Folder to place analysis outputs (default: )

    Options
      --outdir              Base directory for output files (default: )
    """.stripIndent()
}

log.info"""Starting""".stripIndent()

// Show help message if the user specifies the --help flag at runtime
if (params.help) {
  // Invoke the function above which prints the help message
  helpMessage()
  // Exit out and do not run anything else
  exit 0
}

// // Show help message if the user specifies a fasta file but not makedb or db
if (params.seedfile == "") {
  // Invoke the function above which prints the help message
  helpMessage()
  // Exit out and do not run anything else
  exit 0
}

/*
 * Defines the pipeline inputs parameters (giving a default value for each for them)
 * Each of the following parameters can be specified as command line options
 */

// TODO: @sunitj #1 Use seedfiles
// def fnaGlob = "${params.fastas}/*.${params.ext}"
// log.info"""Searching for file at this location: $fnaGlob""".stripIndent()
genomes_ch = Channel.from(
        file(
            params.seedfile
        ).splitCsv(
            header: true,
            sep: ","
        )
    ).filter {
        r -> (r.genome_id != null)
    }.ifEmpty {
        exit 1, "Cannot find values in the 'genome_id' column: ${params.seedfile}"
    }.filter {
        r -> (r.genome_path != null)
    }.ifEmpty {
        exit 1, "Cannot find values in the 'genome_path' column: ${params.seedfile}"
    }.map {
        r -> [r["genome_id"], [r["genome_path"]]]
    }

def outputBase = "${params.outdir}/${params.project}/${params.prefix}"

// Bakta
process BAKTA {
  tag "${id}"

  container params.docker_container_bakta

  publishDir "$outputBase/${id}"

  input:
    tuple val(id), path(assembly) from genomes_ch

  output:
    path "${id}*.gbff", optional: true
    path "${id}*.faa", optional: true
    path "${id}.sha256", optional: true

  script:
  """
  run_bakta.sh ${id} $assembly ${task.cpus} ${params.bakta_db}
  """
}
