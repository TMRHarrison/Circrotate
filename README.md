# Circrotate
Circovirus Rotator

This software takes circovirus genomes, finds replication and capsid proteins through Prokka, and rotates the genome so the replication origin is the first feature in the linearized record. Support for pre-annotated genbank sequences is forthcoming.

## Dependencies:
To run properly, this pipeline needs:
 * [Nextflow](https://www.nextflow.io/)
 * [Python 3](https://www.python.org/downloads/)
 * [Prokka](https://github.com/tseemann/prokka)
 * [BioPython](https://biopython.org/)
 * [bcbiogff](https://github.com/chapmanb/bcbb/tree/master/gff)

Recommended versions of these dependencies are in the ```environment.yml``` file, for use with [Conda](https://docs.conda.io/en/latest/) virtual environments.

## Options:
| Option     | Use                                                  | Default                                                         |
|:-----------|:-----------------------------------------------------|:----------------------------------------------------------------|
| --in       | Specify the input file(s) in fasta format            | You must specify an input                                       |
| --motif    | Specify the JASPAR sites file to use                 | circovirus replication origins in testData/rep_orig_circo.sites |
| --out      | Specify the ouput folder name                        | pipe_out                                                        |
| --seqBatch | Specify the number of sequences processed per thread | 1 sequence per thread                                           |

## Output files:
* (named output)/
	* prokka-annotations/
		* One folder for each annotation format, each containing annotation files for each processed sequence.
	* sequences/
		* fasta/ organized by batch
			* the fasta files of the rotated sequences.
		* gff/ organized by batch
			* the gff annotation files from Prokka.
		* txt/ organized by batch
			* the txt files containing the lists of names of sequences that could not be rotated.
		* "all_sequences.fasta" contains all the rotated sequences in one fasta file.
		* "all_fails.txt" contains all the names of sequences that could not be rotated.
		* "all_annotations.gff" contains all the prokka annotations of rotated sequences.
		* "annotatedFastas.gff" contains all the annotations in gff format, plus all the sequences in a ##FASTA section.

## Test command
This pipeline comes with some genomes from NCBI as test material, as well as some motifs to search for.

### Test commands:
Normal usage:
```
nextflow run Circrotate --in Circrotate/testData/testFasta.fasta
```
This command can also be run by using ```-profile test```.

## Notes:
The program tries to use the start of the rep or cap gene as an anchor for the expected replication origin, in case there is a siimlar sequence elsehwere in the genome. However, this is not perfect, andthe following errors can occur:

- No annotations:
If Prokka doesn't find either rep or cap genes, the program discards the sequence and doesn't attempt to rotate it.
- No motif:
If the rep origin motif can't be located anywhere, the sequence is discarded.
- Incorrect sequence:
Since the rep origin location is based on proximity to one gene alone, it is possible that a similar sequence can exist in the gene closer to the start than the true replication origin, meaning the genome would be rotated to the sequence in the gene. It is also possible that the only seuqence recognised is in the incorrect place.

All sequences have to be re-annotated to get a consistent naming convention for the gene products. Using (inconsistent) manual annotations will probably take some regex tricks that the python script currently isn't set up to handle. However, it bears investigating, as preserving the genbnk annotations automatically makes the sequences easier to handle. Alternatively, I guess I could just load them in, apply an offset, and export them again...

## Versioning:

versioning: X.Y.Z

X: Major release version
  - Anything that adds an analysis step.
  - Things like adding a BLAST search, adding options for genome annotation programs, etc.
  - Changes the core functionality of the program.

Y: Minor release version
  - Changing the interface or adding smaller features, things that don't affect the core functionality of the program.
  - Giving a command line option to add Prokka options, changing output folder strucutre, etc.
  - Optionally different, but produces the same type of output and accepts the same inputs.

Z: Bugfix version
  - For when I make mistakes and have to push another version to change something that otherwise breaks the program.
  - Also used for back-end improvements that don't actually affect the usage or output.

## License:

This software is licenseed under the MIT license. See license.txt for details.