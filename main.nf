#!/usr/bin/env nextflow
nextflow.enable.dsl=2

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

def outputBase = "${params.outdir}/${params.project}/${params.prefix}"

// Bakta
process BAKTA {
  tag "${id}"
  
  debug true

  container params.docker_container_bakta
  containerOptions "-v ${params.bakta_db}:/db:ro"

  publishDir "$outputBase/${id}", mode: 'copy', pattern: "${id}.*"

  input:
    tuple val(id), path(assembly)

  output:
    tuple val("${id}"), path("${id}.faa"), emit: proteins_tuple
    path "${id}*"

  script:
  """
  run_bakta.sh ${id} ${assembly} ${task.cpus} /db/db
  echo "BAKTA finished"
  """
}

// PFam
process PFAMS {
  tag "${id}"

  container params.docker_container_hmmer
  containerOptions "-v ${params.pfam_db}:/db:ro"

  publishDir "$outputBase/${id}/", mode: 'copy', pattern: "${id}.hmmsearch_tbl.csv"
  publishDir "$outputBase/${id}/intermediate_pfam_outputs", mode: 'copy', pattern: "${id}.*tblout.tsv"
  publishDir "$outputBase/${id}/intermediate_pfam_outputs", mode: 'copy', pattern: "${id}.*_hits.csv"

  input:
    tuple val(id), path(proteins)

  output:
    tuple val("domtblout"), val("${id}"), path("${id}.domtbl.tsv")
    tuple val("tblout"), val("${id}"), path("${id}.tbl.tsv")

  script:
  """
  hmmsearch.py -i ${proteins} -o ${id} -d /db/Pfam-A.hmm --threads ${task.cpus} --max_evalue ${params.max_evalue} --min_domain_coverage ${params.min_domain_coverage} --min_overlap ${params.min_overlap}

  echo "PFAM finished"
  """
}

workflow {
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

  genomes_ch | BAKTA
  BAKTA.out.proteins_tuple | PFAMS
}