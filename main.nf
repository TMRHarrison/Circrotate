#!/usr/bin/env nextflow

def helpMessage() {
  c_reset = params.monochrome ? '' : "\033[0m";
  c_dim = params.monochrome ? '' : "\033[2m";

  log.info"\n"+"""
  Circ${c_dim}ovirus${c_reset}
  Rotate
  v1.1.0

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
    --genBank                 Switch for using genbank file(s) with --in.

  Debug/misc:
    --help                    Show this message.
  """.stripIndent()
}

params.in = null
params.out = "pipe_out"
params.motif = "${workflow.projectDir}/testData/rep_orig_circo.sites"
params.prokkaOpts = ""
params.seqBatch = 1
params.genBank = false

params.help = null
params.monochrome = false

// Show help message
if (params.help) {
  helpMessage()
  exit 0
}

def summary = [:]
if (workflow.revision)
  summary['Pipeline Release'] = workflow.revision
summary['Run Name']           = workflow.runName
if (workflow.containerEngine)
  summary['Container']        = "$workflow.containerEngine - $workflow.container"
summary['FASTA input']        = params.in
summary['GenBank file?']      = params.genBank
if (params.prokkaOpts)
  summary['Prokka options:']  = params.prokkaOpts
summary['Seqs per Thread:']   = params.seqBatch
summary['Output dir']         = params.out
summary['Launch dir']         = workflow.launchDir
summary['Working dir']        = workflow.workDir
summary['Script dir']         = workflow.projectDir
summary['User']               = workflow.userName
log.info summary.collect { k,v -> "${k.padRight(18)}: $v" }.join("\n")

if (params.in == null) {
  log.info"You must specify an input FASTA file with \"--in <file>\"."
  exit 1
}

motif_vch = Channel.value(file(params.motif))

genBank_vch = params.genBank ? Channel.value(file(params.in)) : Channel.value(false)

if (params.genBank) {
  log.info"Genbank file in use"
   genBankIn_ch = Channel.fromPath(params.in)
   fastaSequences_ch = Channel.empty()
}
else {
  fastaSequences_ch = Channel.fromPath(params.in)
  genBankIn_ch = Channel.empty()
}

// if there is a genbank file, make it into fasta files and pass it to giveFileNameFastaID
// Otherwise, the input should be fasta and sent straight to giveFileNameFastaID
process makeFasta {
  input:
  file inp from genBankIn_ch // <--- --in file when --genbank specified 

  output:
  file "${inp.baseName}.fasta" into gbFastaSequences_ch // ---> giveFileNameFastaID

  when:
  params.genBank

  script:
  """
  gb2fasta.py ${inp} ${inp.baseName}.fasta
  """
}

// appends the sequence ID to the filename. Only takes the ID up until the first non-alphanumeric character that isn't '.', '_', or '-'.
// This behaviour is to prevent illegal filenames.
// The rest of the filename is maintained to ensure the output won't have any overlapping names.
process giveFileNameFastaID {
  input:
  file inp from fastaSequences_ch.mix(gbFastaSequences_ch).collect().splitFasta(by: params.seqBatch, file: true) // <--- --in for fasta files, makeFasta for genbank files
  // mixing the channels allows one or the other to be optional

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
  file motifs from motif_vch // <--- --motif
  file gbinp from genBank_vch // <--- --in iff params.genbank is set

  output:
  file "*.fasta" into combineFastas_ch // ---> combineControlSeqs
  file "*.txt" optional true into combineFails_ch // ---> combineFails
  file "*.gff" into combineAnnotations_ch // ---> combineAnnotations
  file "*.csv" into combineCsv_ch // ---> combineAnnotations

  script:
  pref = "${inp.baseName}_rotated"
  gb = gbinp ? "--genbank ${gbinp}" : ""

  """
  rotate_seq.py --input ${inp} --output ${pref} --motif ${motifs} ${gb}
  """
}

// concatenates the sequences into a single file
process combineFastas {
  publishDir "${params.out}/sequences", mode: 'copy'

  input:
  file "f*.fasta" from combineFastas_ch.collect() // <--- rotateSeqs

  output:
  file "all_sequences.fasta" into combinedFastas_ch // ---> combineFastaGff

  /* strip all the blank lines
  if none are found, suppress the error exit code by piping through cat
  It doesn't like doing it in one stream from multiple files, so it concatenates to a single file first.
  This way does take more disk space and time but also for some reason it would hang 
  at 0% memory usage and just fill the work/ folder with massive (>100GB) files otherwise
  
  Hey, I found more weirdness: .collect() will give files sequential names to globs (1.fa, 2.fa, 3.fa, etc.)
  But only if there are multiple files in the queue. If there is one file, the glob gets removed completely.
  So if this is "*.fasta", the single collected file is ".fasta" which is a hidden file. Which cat can't find.
  So I have to call it f*.fasta so it goes to [f1.fasta, f2.fasta, ...] or "f.fasta"

  I guess it makes sense to do it like that, but it would have been nice if it said this was anywhere in the docs :^)
  */
  """
  cat *.fasta > tmp.f
  grep -v "^\$" tmp.f | cat - > all_sequences.fasta
  """
}

// concatenates the failed IDs into a single file
process combineFails {
  publishDir "${params.out}/sequences", mode: 'copy'

  input:
  file "t*.txt" from combineFails_ch.collect() // <--- rotateSeqs

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
  file "g*.gff" from combineAnnotations_ch.collect() // <--- rotateSeqs

  output:
  file "all_annotations.gff" into combinedAnnotations_ch // ---> combineFastaGff

  // strip the headers and place a single gff header in the final file.
  """
  cat *.gff > tmp.g
  HEADER=\$(head -n1 tmp.g)
  (echo \${HEADER} ; grep -v "^\${HEADER}" tmp.g ) > all_annotations.gff
  """
}

// concatenates the csv results into a single file
process combineResults {
  publishDir "${params.out}/sequences", mode: 'copy'

  input:
  file "c*.csv" from combineCsv_ch.collect() // <--- rotateSeqs

  output:
  file "all_results.csv"

  // strip the headers and place a single csv header in the final file.
  """
  cat *.csv > tmp.c
  HEADER=\$(head -n1 tmp.c)
  (echo \${HEADER} ; grep -v "^\${HEADER}" tmp.c ) > all_results.csv
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

