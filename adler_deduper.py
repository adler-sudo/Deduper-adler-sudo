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
        "-o",
        "--output_retain_file",
        help="Output retain SAM file containing uniquely mapped reads."
    )
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
retain_filename = args.output_retain_file

# argparser
if args.input_file is None:
    exit("ATTENTION: Exiting script!\nPlease specify a sorted SAM file utilizing the input_file flag (-f, --input_file).")
if args.output_retain_file is None:
    exit("ATTENTION: Exiting script!\nPlease specify a filename to which the uniquely mapped reads will be written.")
if args.paired_end:
    paired_end = True
    print('Paired end SAM file indicated.')

# keep highest qscore
if args.keep_highest_qscore:
    keep_highest_qscore=args.keep_highest_qscore
else:
    keep_highest_qscore=False

print("Entering dedupe function.")

# instantiate Dedupe class
dedupe_me = Dedupe(
    input_filename=input_filename,
    umi_filename=umi_filename,
    retention_filename=retain_filename,
    keep_highest_qscore=keep_highest_qscore,
    paired_end=paired_end
)
dedupe_me.dedupe()