#!/usr/bin/env nextflow



def helpMessage() {
  c_reset = params.monochrome ? '' : "\033[0m";
  c_dim = params.monochrome ? '' : "\033[2m";

  log.info"\n"+"""
  Usage:

  The typical command for running the pipeline is as follows:

  nextflow run Circrotate --in data.fasta --out out_folder

  Mandatory arguments:
    --in                      Path to input fasta file(s).

  Other options:
    --out                     Output folder name. Defaults to "pipe_out"
    --motif                   A JASPAR sites file with replication origins.
    --prokkaOpts              Extra prokka options. Must be wrapped in quotes.
    --seqBatch                Number of sequences to process per thread.
  """.stripIndent()
}

params.in = null
params.motif = "${workflow.projectDir}/testData/rep_orig_circo.sites"
params.out = "pipe_out"
params.help = null
params.prokkaOpts = ""
params.seqBatch = 1
params.monochrome = false

def summary = [:]
if (workflow.revision)
  summary['Pipeline Release'] = workflow.revision
summary['Run Name']           = workflow.runName
if (workflow.containerEngine)
  summary['Container']        = "$workflow.containerEngine - $workflow.container"
summary['FASTA input']        = params.in
if (params.prokkaOpts)
  summary['Prokka options:']  = params.prokkaOpts
summary['Seqs per Thread:']   = params.seqBatch
summary['Output dir']         = params.out
summary['Launch dir']         = workflow.launchDir
summary['Working dir']        = workflow.workDir
summary['Script dir']         = workflow.projectDir
summary['User']               = workflow.userName
log.info summary.collect { k,v -> "${k.padRight(18)}: $v" }.join("\n")

// Show help message
if (params.help) {
  helpMessage()
  exit 0
}


if (params.in == null) {
  log.info"You must specify an input FASTA file with \"--in <file>\"."
  exit 1
}


// Split the sequences into their own processes.
sequences_ch = Channel.fromPath(params.in)
          .splitFasta(by: params.seqBatch, file: true)
          .dump()
motif_ch = Channel.fromPath(params.motif)

// appends the sequence ID to the filename. Only takes the ID up until the first non-alphanumeric character that isn't '.', '_', or '-'.
// This behaviour is to prevent illegal filenames.
// The rest of the filename is maintained to ensure the output won't have any overlapping names.
process giveFileNameFastaID {
  input:
  file inp from sequences_ch// <--- --in

  output:
  file "*_${inp}" into renamedSequences_ch // ---> performProkka

  """
  fastaID=\$(awk -F "[^a-zA-Z0-9\\._-]" '/^>/ {print \$2; exit}' ${inp})
  cp $inp \${fastaID}_${inp}
  """
}

// performs Prokka on the fasta files
process performProkka {
  // Publish to folders by extension, not by sequence.
  publishDir "${params.out}/prokka-annotations",
    mode: 'copy',
    saveAs: {filename ->
      fileOut = file(filename)
      "${fileOut.getExtension()}/${fileOut.getName()}"
    }

  input:
  file inp from renamedSequences_ch // <--- giveFileNameFastaID

  output:
  file "**/*.gff" into annotatedSeqs_ch // ---> rotateSeq
  file("**/*")

  """
  prokka --outdir prokka-out --force --prefix ${inp.baseName} --cpus ${task.cpus} --gcode 1 --kingdom virus ${inp} ${params.prokkaOpts}
  """
}

// rotates the sequences to the correct nucleotide motif
process rotateSeqs {
  publishDir "${params.out}/sequences/", 
    mode: 'copy',
    saveAs: {filename ->
      fileOut = file(filename)
      "${fileOut.getExtension()}/${fileOut.getName()}"
    }

  input:
  file inp from annotatedSeqs_ch // <--- performProkka

  output:
  file "*.fasta" into combineFastas_ch // ---> combineControlSeqs
  file "*.txt" into combineFails_ch // ---> combineFails
  file "*.gff" into combineAnnotations_ch // ---> combineAnnotations

  script:
  pref = "${inp.baseName}_rotated"

  """
  rotate_seq.py --input ${inp} --output ${pref} --motif ${params.motif}
  """
}

// concatenates the sequences into a single file
process combineFastas {
  publishDir "${params.out}/sequences", mode: 'copy'

  input:
  file "*.fasta" from combineFastas_ch.collect() // <--- rotateSeqs

  output:
  file "all_sequences.fasta" into combinedFastas_ch // ---> combineFastaGff

  // strip all the blank lines
  // if none are found, suppress the error exit code by piping through cat
  // It doesn't like doing it in one stream from multiple files, so it concatenates to a single file first.
  // This way does take more disk space and time but also for some reason it would hang 
  // at 0% memory usage and just fill the work/ folder with massive (>100GB) files otherwise
  """
  cat *.fasta > tmp.f
  grep -v "^\$" tmp.f | cat - > all_sequences.fasta
  """
}

// concatenates the failed IDs into a single file
process combineFails {
  publishDir "${params.out}/sequences", mode: 'copy'

  input:
  file "*.txt" from combineFails_ch.collect() // <--- rotateSeqs

  output:
  file "all_fails.txt"

  // strip all the blank lines
  """
  cat *.txt > tmp.t
  grep -v "^\$" tmp.t | cat - > all_fails.txt
  """
}

// concatenates the gff annotations into a single file
process combineAnnotations {
  publishDir "${params.out}/sequences", mode: 'copy'

  input:
  file "*.gff" from combineAnnotations_ch.collect() // <--- rotateSeqs

  output:
  file "all_annotations.gff" into combinedAnnotations_ch // ---> combineFastaGff

  // strip the headers and place a single gff header in the final file.
  """
  cat *.gff > tmp.g
  (echo "##gff-version 3" ; grep -v "^##gff-" tmp.g ) > all_annotations.gff
  """
}

// concatenates the control sequences into a single file
process combineFastaGff {
  publishDir "${params.out}/sequences", mode: 'copy'

  input:
  file "all_sequences.fasta" from combinedFastas_ch // <--- combineFastas
  file "all_annotations.gff" from combinedAnnotations_ch // <--- combineAnnotations

  output:
  file "annotatedFastas.gff"

  """
  (cat all_annotations.gff ; echo "##FASTA" ; cat all_sequences.fasta) > annotatedFastas.gff
  """
}

workflow.onComplete {
  log.info "${workflow.runName} complete"
  log.info "Pipeline completed at: $workflow.complete"
  log.info "Execution status: ${ workflow.success ? "OK" : "failed" }" //"
}

