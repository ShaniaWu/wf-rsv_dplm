#!/usr/bin/env python
"""Functions to help us parse files for flu."""
import csv
import json
import os

import pysam

from .util import get_named_logger, wf_parser  # noqa: ABS101


def process_typing(typing_json):
    """Process abricate typing string."""
    ## get rsv_type from typing_json
    typing = json.load(open(typing_json))
    rsv_type = typing['type']
    return(rsv_type)


def make_consensus(consensus, rsv_type, output):
    fasta = pysam.FastaFile(consensus)

    if rsv_type in ("RSV_TYPE_A", "RSV_TYPE_B"):
        ## isolate consensus for rsv type
        if rsv_type == "RSV_TYPE_A":
            sequence = fasta.fetch("typeA_EPI_ISL_412866")
            header = ">typeA"
        elif rsv_type == "RSV_TYPE_B":
            sequence = fasta.fetch("typeB_EPI_ISL_1653999")
            header = ">typeB"

        # f = open('typed.consensus.fasta', "x")
        with open(output, "w") as f:
            f.write(f"{header}\n{sequence}")
        
    elif rsv_type == "mixedAB":
        with open(output, "w") as f:
            f.write("mixedAB")
        print("Warning: Mixed typing information.")
    
    elif rsv_type == "None": # rsv_type = "None"
        with open(output, "w") as f:
            f.write("NA")
        print("Warning: No typing information.")



def main(args):
    """Run the entry point."""
    logger = get_named_logger("nextclade_helper")
    
    typing = process_typing(args.typing)
    make_consensus(args.consensus, typing, args.output)
    logger.info("nextclade helping done.")
    
    
def argparser():
    """Argument parser for entrypoint."""
    parser = wf_parser("nextclade_helper")
    parser.add_argument(
        "--typing",
        help="Typing json from abricate.")
    parser.add_argument(
        "--consensus",
        help="Consensus FASTA from sample.")
    parser.add_argument(
        "--output", default=None,
        help="Typed consensus fasta.")
    return parser

