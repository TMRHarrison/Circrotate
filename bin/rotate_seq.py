#!/usr/bin/env python

"""
rotate_seq.py --input <infile>.gff > <outfile>.fasta

This script takes an input of a gff file (usually generated by Proka) with a
##FASTA section containing the sequence data, and outputs a FASTA file of all
the control regions found. The script is capable of finding control regions in
reverse complemented and rotated sequences, unless the feature is split across
the start/end of the sequence.

The general process is as follows:
    + Annotations for the rep protein are located.
    + If necessary, adjustments are made for reverse complemented sequences.
    + The FASTA information gets rotated to start at the first occurence of the
        replication origin sequence.
"""

# command line arguments, GFF parser, BioPython tools, csv writer
import argparse
from BCBio import GFF
from Bio import motifs, SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
import csv

# Type hinting
from typing import TYPE_CHECKING, Dict, Tuple, Iterable, List, Optional
if TYPE_CHECKING:
    from Bio.SeqFeature import SeqFeature, FeatureLocation
    from Bio.SeqRecord import SeqRecord

def get_params():
    """Returns the command line arguments."""

    parser = argparse.ArgumentParser(description="""
        This script takes an input of a gff file (usually generated by Prokka) with a
        ##FASTA section containing the sequence data, rotates the sequences to their
        replication origins, and outputs the sequences, their annotations, the names
        of any sequences that failed to rotate, and a summary of the actions taken.
        The script is capable of finding replication origins in reverse complemented
        and rotated sequences.
        """.strip()
    )
    parser.add_argument("--input", help="The file to be worked on. GFF format with a ##FASTA section.")
    parser.add_argument("--motif", help="The JASPAR sites file of possible replication origins.")
    parser.add_argument("--output", help="The prefix for the output files")
    parser.add_argument("--force", action="store_true", help="Overwrites the output if it already exists.")
    parser.add_argument("--genbank", help="optional genbank annotation file. Overrides the prokka annotations when exporting.")

    return parser.parse_args()

def find_anchor(bound_start_anchor, start_anchor_annots):
    """
    Finds the first available anchor gene and returns its position

    arguments:
        bound_start_anchor: a tuple containing tuples of: (feature name, "start"/"end")
        start_anchor_annots: a list containing all the possible start anchors (even ones that are low-priority)
    return:
        the annotation of the anchor
        the anchor's specification
    """
    for s_anch in bound_start_anchor:
        for annot in start_anchor_annots:
            if annot.qualifiers["product"][0] == s_anch[0]:
                anch_pos = annot

                return annot, s_anch # Bonus: this breaks both loops
    return None, None # otherwise, return the default (which is None)

def get_anchor_pos(anchor_spec, anchor_annot):
    """
    gets the specification for an anchor annotation and the actual anchor annotation and returns the position of the target

    arguments:
        anchor_spec: 
        anchor_annot: 
    return:
        location of the anchor (int)
    """

    # all annotations are start < end, so ones on the reverse strand start at the "end" of their annotated range.
    # Use the start bound if it's anchored to the start of the annotation and on the + strand, or
    # if it's the end of an annotation on the - strand
    if ((anchor_spec[2] == 1 and anchor_spec[1] == "start")
            or (anchor_spec[2] == -1 and anchor_spec[1] == "end")):
        return anchor_annot.location.start
    else:
        return anchor_annot.location.end

def get_positions(rec, motif):
    """
    Takes a list of motifs, a sequence record, and finds all occurences of the motifs on the sequence, then returns the list
    
    arguments:
        rec: a biopython sequence record object
        motif: a biopython motif object
    return:
        list of starting positions of found motifs (int)
    """
    posses = []

    # overhang: grab a couple nt from the start we we find motifs broken up by the linearization
    # If there are len(motif) nt on one side that allows us to locate to motif, then it wasn't broken up.
    # if we take len(motif)-1 from one side, then we also need at least 1 from the other side, so it must have been split.
    overhang = len(motif)-1
    temp_rec = rec+rec[:overhang]
    for pos, seq in motif.instances.search(temp_rec.seq):
        posses.append(pos)    

    return posses

def find_offset(rec, positions, anchor_pos):
    """
    Finds the offset needed to rotate the sequence to its origin of replication
    
    arguments:
        rec: a biopython sequence record
        positions: a list of possible positions to rotate the sequence to
        anchor_pos: a "target" location to get as close as possible to
        rcomp: boolean for reverse complement
        motif: BioPython motif object
    return:
        the offset needed to rotate to replication origin
    """

    if len(positions) != 0:
        shift_dist = len(rec)
        shift_offset = 0
        
        for n in positions:
            shift_dist_new = min(shift_dist, circular_distance(n, anchor_pos, len(rec)))
            if shift_dist_new < shift_dist:
                shift_offset = n
                shift_dist = shift_dist_new
            
        return shift_offset

    return None

def rotate_seq(rec, shift_offset, rcomp):
    """
    Rotates a sequence record by an offset, and reverse complements it if necessary.

    arguments:
        rec: a biopython sequence record
        positions: a list of possible positions to rotate the sequence to
        anchor_pos: a "target" location to get as close as possible to
    return:
        the rotated sequence record object, or None if shift offset is None
    """

    if not shift_offset is None:
        out_rec = rec[shift_offset:] + rec[:shift_offset]
        if rcomp:
            out_rec = out_rec.reverse_complement(id=True)

        return out_rec

    return None

def circular_distance(a: int, b: int, C: int) -> int:
    """
    Finds the shortest distance between two points along the perimeter of a circle.

    arguments:
        a: a point on a circle's circumference.
        b: another point on the cicrle.
        C: the total circumference of the circle.
    return:
        The shortest distance along the circumference of the circle between the two points

    >>> circular_distance(2,5,10)
    3
    >>> circular_distance(12,3,15)
    6
    >>> # It even works with numbers >C or <0
    ... circular_distance(-20, 37, 10)
    3
    """
    arc = abs(a - b) % C # the distance between these in one direction -- not necessarily the shortest distance
    return min(C - arc, arc) # the arc and the complement of the arc, one of which will be shorter than the other.

def get_names(rec_list):
    """
    Gets the sequnece IDs from a sequence record list, returns a comma separated string of IDs
    """
    return ", ".join([r["rec"].id for r in rec_list])

def main():
    """Main CLI entry point for rotate_seq.py"""
    args = get_params()

    rec_succ = []
    # contains a dict with:
    #   ["rec"]: the record
    #   ["shift"]: the distance shifted
    #   ["rev_comp"]: T/F reverse complement necessary
    #   ["genbank"]: true if genbank record was successfully located.
    #   ["posses"]: positions on the original record that match the motif(s)

    rec_fail = []
    # contains a dict of:
    #   ["rec"]: the record
    #   ["step"]: failure at either "annotation" or at "motif"

    # read the input files
    with open(args.input, "r") as gff_file, \
            open(args.motif, "r") as motif_file:
        
        # if the genbank file is given, load it in as an index, otherwise, give an empty dict
        gb_ind = SeqIO.index(args.genbank, "genbank") if args.genbank else {}

        # highest to lowest priority, e.x. cap is preferred over rep
        # [0] product name
        # [1] start or end of gene
        # [2] 
        anchors = [
            (
                "Capsid protein",
                "start",
                -1
            ),
            (
                "Replication-associated protein",
                "start",
                1
            )
        ]

        motif = motifs.read(motif_file, "sites")

        for rec in GFF.parse(gff_file):
            cur_rec = rec
            cur_motif = motif
            cur_anchor_pos = 0
            rcomp = False
            
            # Make a dictionary called prod_features containing SeqFeature objects with products
            prod_features = [f for f in cur_rec.features if 'product' in f.qualifiers]

            anchor_annot, anchor_spec = find_anchor(anchors, prod_features)

            # if no anchors were found, stop looking for the rep origin
            if anchor_annot is None:
                rec_fail.append(
                    {
                        "rec": rec,
                        "step": "annotation"
                    }
                )
                continue

            # If the anchor is on the wrong strand, flip the motifs. This is less involved than flipping the record.
            if anchor_annot.location.strand != anchor_spec[2]:
                cur_motif = cur_motif.reverse_complement()
                rcomp = True
            
            # get the position of the anchor
            cur_anchor_pos = get_anchor_pos(anchor_spec, anchor_annot)
            
            # if there is a sequence with the corect id, pull it from the genbank list, otherwise, just keep using the prokka annotations.
            cur_rec = gb_ind.get(cur_rec.id, cur_rec)

            # find rep origin(s) from list
            posses = get_positions(cur_rec, cur_motif)

            # if the sequence needs to be reversed, get the "end" of the motif, which will be the start when it gets flipped
            if rcomp:
                posses = [(i + len(cur_motif)) % len(cur_rec) for i in posses]

            # rotate to location of origin closest to the anchor point
            shift_offset = find_offset(cur_rec, posses, cur_anchor_pos)
            out_rec = rotate_seq(cur_rec, shift_offset, rcomp)

            if not (out_rec is None):
                # just pop in a couple extra things:

                # sequence is circular
                out_rec.features.append(
                    SeqFeature(
                        FeatureLocation(0,len(out_rec)),
                        type="region",
                        strand=1,
                        qualifiers={"Is_circular": 1}
                    )
                )

                # if there's a source listed for a sequence, put it in as an annotation
                if "source" in out_rec.annotations:
                    out_rec.features.append(
                        SeqFeature(
                            FeatureLocation(0,len(out_rec)),
                            type="source",
                            strand=1,
                            qualifiers={"Name": out_rec.annotations["source"]}
                        )
                    )

                # put the successful record into the table to get printed out later
                rec_succ.append(
                    {
                        "rec": out_rec,
                        "shift": shift_offset,
                        "rev_comp": rcomp,
                        "genbank": cur_rec.id in gb_ind,
                        "posses": posses
                    }
                )

            else:
                rec_fail.append(
                    {
                        "rec": rec,
                        "step": "motif"
                    }
                )

    # successes/failures
    print("Successful rotation: "+get_names(rec_succ))
    print("Unsuccessful: "+get_names(rec_fail))

    # make a new file, either forcing overwrite of the old file or not, depending on the setting.

    # get the list of records from the successful records
    with open(args.output+"-annotations.gff", "w" if args.force else "x") as gffout_file:
        GFF.write([r["rec"] for r in rec_succ], gffout_file)

    # same as above but print the fasta sequences
    with open(args.output+"-seqs.fasta", "w" if args.force else "x") as seqout_file:
        SeqIO.write([r["rec"] for r in rec_succ], seqout_file, "fasta")

    # Write the fail file iff 1 or more things failed
    if rec_fail:
        # Print out all the sequence IDs of the sequences that failed to rotate
        with open(args.output+"-fails.txt", "w" if args.force else "x") as fail_file:
            fail_file.write("\n".join([r["rec"].id for r in rec_fail])+"\n")

    # the results csv with failures specified, shift distance, etc.
    with open(args.output+"-results.csv", "w" if args.force else "x") as csv_file:
        out_writer = csv.writer(csv_file)

        # header row
        out_writer.writerow(
            [
                "Name",
                "Success",
                "Shift",
                "Reversed",
                "Genbank retrieved",
                "Step failed",
                "Possible rep origins"
            ]
        )

        # rows for successes
        for r in rec_succ:
            out_writer.writerow(
                [
                    r["rec"].id,
                    "yes",
                    r["shift"],
                    r["rev_comp"],
                    r["genbank"],
                    None,
                    ", ".join(str(i) for i in r["posses"])
                ]
            )

        # rows for failures
        for r in rec_fail:
            out_writer.writerow(
                [
                    r["rec"].id,
                    "no",
                    "None",
                    False,
                    False,
                    r["step"],
                    None
                ]
            )

if __name__ == '__main__':
    main()
