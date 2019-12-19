#!/usr/bin/env nextflow

def helpMessage() {
  c_reset = params.monochrome ? '' : "\033[0m";
  c_dim = params.monochrome ? '' : "\033[2m";

  log.info"\n"+"""
  Circ${c_dim}ovirus${c_reset}
  Rotate
  v1.1.1

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

motif_file = file(params.motif)

genBank_file = file(params.in)

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
  path inp from genBankIn_ch // <--- --in file when --genbank specified 

  output:
  path "${inp.baseName}.fasta" into gbFastaSequences_ch // ---> giveFileNameFastaID

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
  path inp from fastaSequences_ch.mix(gbFastaSequences_ch).splitFasta(by: params.seqBatch, file: true) // <--- --in for fasta files, makeFasta for genbank files
  // mixing the channels allows one or the other to be optional

  output:
  path "*_${inp}" into renamedSequences_ch // ---> performProkka

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
  path inp from renamedSequences_ch // <--- giveFileNameFastaID

  output:
  tuple path("**/*.gff"), path(sed_pat_fn) into annotatedSeqs_ch // ---> rotateSeq
  path "**/*"

  script:
  sed_pat_fn = "${inp.baseName}-sedpat.spf" // sed pattern file

  """
  prokka --outdir prokka-out --force --prefix ${inp.baseName} --cpus ${task.cpus} --gcode 1 --kingdom virus ${inp} ${params.prokkaOpts}

  awk -F" " -e '/^>/ {st = index(\$0," ") ; print "s>^\\\\(\\\\"\$1".*\\\\)>\\\\1 "substr(\$0,st+1)">"}' ${inp} > ${sed_pat_fn}
  """
  // The awk command constructs a sed pattern, each of the \\\\ escapes to one \ in the final sed (escape once for nextflow, again for bash
  // "s>^\(\"\$1"\)>\1 "substr(\$0,st+1)">"
  // s>^(\>seqID)>\1 seqDesc>
  // Basically replace the fasta header with just the ID with the entire description
  // It uses > as the delimiter because it usually doesn't appear in fasta headers except the start character.
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
  tuple path(inp), path(sedPat) from annotatedSeqs_ch // <--- performProkka
  file motifs from motif_file // <--- --motif
  file gbinp from genBank_file // <--- --in iff params.genbank is set

  output:
  path "*.fasta" optional true into combineFastas_ch // ---> combineControlSeqs
  path "*.txt" optional true into combineFails_ch // ---> combineFails
  path "*.gff" optional true into combineAnnotations_ch // ---> combineAnnotations
  path "*.csv" into combineCsv_ch // ---> combineAnnotations

  script:
  pref = "${inp.baseName}_rotated"
  gb = gbinp.name != 'NO_FILE' ? "--genbank ${gbinp}" : ""

  // Rotate the sequence, then if it puts a fasta out, put the description back on the fasta
  // piping to echo just suppresses the error code for sed if rotate_seq.py doesn't make a fasta file
  """
  rotate_seq.py --input ${inp} --output ${pref} --motif ${motifs} ${gb}
  sed -i -f ${sedPat} *.fasta | echo
  """
}

// concatenates the sequences into a single file
process combineFastas {
  publishDir "${params.out}/sequences", mode: 'copy'

  input:
  path "f*.fasta" from combineFastas_ch.collect() // <--- rotateSeqs

  output:
  path "all_sequences.fasta" into combinedFastas_ch // ---> combineFastaGff

  /*
  Hey, I found more weirdness: .collect() will give files sequential names to globs (1.fa, 2.fa, 3.fa, etc.)
  But only if there are multiple files in the queue. If there is one file, the glob gets removed completely.
  So if this is "*.fasta", the single collected file is ".fasta" which is a hidden file. Which cat can't find.
  So I have to call it f*.fasta so it goes to [f1.fasta, f2.fasta, ...] or "f.fasta"

  I guess it makes sense to do it like that, but it would have been nice if it said this was anywhere in the docs :^)
  */
  """
  cat *.fasta > all_sequences.fasta
  """
}

// concatenates the failed IDs into a single file
process combineFails {
  publishDir "${params.out}/sequences", mode: 'copy'

  input:
  path "t*.txt" from combineFails_ch.collect() // <--- rotateSeqs

  output:
  path "all_fails.txt"

  """
  cat *.txt > all_fails.txt
  """
}

// concatenates the gff annotations into a single file
process combineAnnotations {
  publishDir "${params.out}/sequences", mode: 'copy'

  input:
  path "g*.gff" from combineAnnotations_ch.collect() // <--- rotateSeqs

  output:
  path "all_annotations.gff" into combinedAnnotations_ch // ---> combineFastaGff

  /* 
  strip the headers and place a single gff header in the final file.
  It doesn't like doing it in one stream from multiple files, so it concatenates to a single file first.
  This way does take more disk space and time but also for some reason it would hang 
  at 0% memory usage and just fill the work/ folder with massive (>100GB) files otherwise
  */
  """
  strip_headers.sh *.gff > all_annotations.gff
  """
}

// concatenates the csv results into a single file
process combineResults {
  publishDir "${params.out}/sequences", mode: 'copy'

  input:
  path "c*.csv" from combineCsv_ch.collect() // <--- rotateSeqs

  output:
  path "all_results.csv"

  // strip the headers and place a single csv header in the final file.
  """
  strip_headers.sh *.csv > all_results.csv
  """
}

// concatenates the control sequences into a single file
process combineFastaGff {
  publishDir "${params.out}/sequences", mode: 'copy'

  input:
  path "all_sequences.fasta" from combinedFastas_ch // <--- combineFastas
  path "all_annotations.gff" from combinedAnnotations_ch // <--- combineAnnotations

  output:
  path "annotatedFastas.gff"

  """
  (cat all_annotations.gff ; echo "##FASTA" ; cat all_sequences.fasta) > annotatedFastas.gff
  """
}

workflow.onComplete {
  log.info "${workflow.runName} complete"
  log.info "Pipeline completed at: $workflow.complete"
  log.info "Execution status: ${ workflow.success ? "OK" : "failed" }" //"
}

