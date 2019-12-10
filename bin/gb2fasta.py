#!/usr/bin/env python

"""
gb2fasta.py input.gb output.fasta

converts genbank to fasta
"""

# command line arguments, GFF parser
import argparse
from Bio import SeqIO

def get_params():
    """Returns the command line arguments."""

    parser = argparse.ArgumentParser(description="""
        Uses biopython to convert genbank files to fasta.
        """.strip())
    parser.add_argument("input", help="The genbank file to be converted to fasta")
    parser.add_argument("output", help="the fasta output file location")

    return parser.parse_args()

def main():
    """Main CLI entry point for gb2fasta.py"""
    args = get_params()

    SeqIO.convert(args.input, "genbank", args.output, "fasta")

if __name__ == '__main__':
    main()
