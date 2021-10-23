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
    return parser.parse_args(args)

args = parse_args(sys.argv[1:])

# define globals
paired_end = False
input_file = args.input_file

# argparser
if args.input_file is None:
    exit("ATTENTION: Exiting script!\nPlease specify a sorted SAM file utilizing the input_file flag (-f, --input_file).\n")
if args.paired_end:
    paired_end = True
    exit("ATTENTION: Exiting script!\nApologies! Paired-end functionality not yet included. Please remove paired_end flag (-p, --paired_end).\n")
if args.umi:
    umis = args.umi
else:
    exit("ATTENTION: Exiting script!\nApologies! Randomer UMI functionality not yet included. Please utilize umi flag (-u, --umi) and specify UMI file.\n")


