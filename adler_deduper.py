#!/usr/bin/env python

# library import
import argparse
import sys
from sam_funcs import *


# define arguments
def parse_args(args):
    parser = argparse.ArgumentParser(description="Reference-based PCR duplicate removal from mapped SAM file")
    parser.add_argument(
        "-f",
        "--input_file", 
        help="Input SAM file - sorted, mapped reads")
    parser.add_argument(
        "-p", 
        "--paired-end", 
        action="store_true", 
        help="Designates file is paired end. If flag is not set, the program will assume single end reads.")
    parser.add_argument(
        "-u", 
        "--umi", 
        help="Designates UMI file. If set, UMI file should follow flag. If flag is not set, the program will assume randomer UMI")
    parser.add_argument(
        "-q",
        "--keep_highest_qscore",
        action="store_true",
        help="When flag is set, if the program finds a PCR duplicate, it will keep the duplicate with the highest qscore."
    )
    return parser.parse_args(args)

args = parse_args(sys.argv[1:])

# define globals
paired_end = False
input_filename = args.input_file
umi_filename = args.umi

# argparser
if args.input_file is None:
    exit("ATTENTION: Exiting script!\nPlease specify a sorted SAM file utilizing the input_file flag (-f, --input_file).\n")
if args.paired_end:
    paired_end = True
    exit("ATTENTION: Exiting script!\nApologies! Paired-end functionality not yet included. Please remove paired_end flag (-p, --paired_end).\n")

# keep highest qscore
if args.keep_highest_qscore:
    keep_highest_qscore=args.keep_highest_qscore
else:
    keep_highest_qscore=False

# instantiate Dedupe class
dedupe_me = Dedupe(
    input_filename=input_filename,
    umi_filename=umi_filename,
    retention_filename="C1_SE_uniqAlign.retain.sam", # TODO: remove declaration here - shift to argparse
    # duplicate_filename="C1_SE_uniqAlign.retain.sam",
    keep_highest_qscore=keep_highest_qscore
)
dedupe_me.dedupe()