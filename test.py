#!/usr/bin/env python

from sam_funcs import *
import sys
import argparse
import os
import filecmp


# define globals
test_dir = "test/"
test_input_dir = test_dir + "inputs/"
test_output_actual_dir = test_dir + "outputs/actual/"
test_output_expected_dir = test_dir + "outputs/expected/"

input_filename = "./test.sam"
umi_filename = "test/inputs/test.umi.txt"
retention_filename = test_output_actual_dir + "/helloworld.retain.sam"
duplicate_filename = test_output_actual_dir + "helloworld.discard.sam"
read_length = 150

test_SAM_read = "NS500451:154:HWKTMBGXX:1:11101:21621:1145:AGGTTGCT	0	2	93022350	36	2S89M	*	0	0	TTCCACTGTTGCTTCATAACTGCAGTCCTAACATAAATGTCTGACATGTAGGATGATCTTAAGCAACCCCT	6AEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE<EEAAAEE"

# parse args for testing different test scenarios
# TODO: if we move this to sam_funcs then we can just import it
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

# define additional arguments
if args.input_file:
    base_filename = os.path.basename(args.input_file).split(".")[0]
    input_filename=test_input_dir + "{}.sam".format(base_filename)
    retention_filename=test_output_actual_dir + "{}.retain.sam".format(base_filename)
    duplicate_filename=test_output_actual_dir + "{}.discard.sam".format(base_filename)

    retention_comp_filename=test_output_expected_dir+"{}.retain.sam".format(base_filename)
    duplicate_comp_filename=test_output_expected_dir+"{}.discard.sam".format(base_filename)

# define test classes
test_samfile_class = SamFile(
    input_filename=input_filename,
    umi_filename=umi_filename,
    retention_filename=retention_filename,
    duplicate_filename=duplicate_filename,
    read_length=read_length
)

test_samread_class = SamRead(
    test_SAM_read
)

# assert that attributes are being assigned correctly
assert test_samfile_class.input_filename == input_filename, "Issue with test_samfile_class input_filename attribute"
assert test_samfile_class.umi_filename == umi_filename, "Issue with test_samfile_class umi_filename attribute"


def test_read_umis():
    eval_dict = {
        'AACGCCAT': [], 
        'AAGGTACG': [], 
        'AATTCCGG': [], 
        'ACACAGAG': [], 
        'ACACTCAG': [], 
        'ACACTGTG': [], 
        'ACAGGACA': []}
    test_samfile_class.read_umis()
    assert test_samfile_class.eval_dict == eval_dict, "Error in read_umis function"
test_read_umis()

def test_dedupe():
    test_samfile_class.dedupe()
    if args.input_file:
        assert filecmp.cmp(retention_filename,retention_comp_filename), "Error in {}".format(retention_filename)
        assert filecmp.cmp(duplicate_filename,duplicate_comp_filename), "Error in {}".format(duplicate_filename)
test_dedupe()

def test_parse_columns():
    assert test_samread_class.qname == test_SAM_read.split("\t")[0], "Issue with SamRead method parse_columns"
test_parse_columns()

def test_correct_start_position():
    test_samread_class.correct_start_position(10,"2S","+")
    assert test_samread_class.pos == 8, "SamRead correct_start_position method making error on forward strand."

    test_samread_class.correct_start_position(10,"2M","+")
    assert test_samread_class.pos == 10, "SamRead correct_start_position method incorrectly adjusting start position."

    test_samread_class.correct_start_position(10,"2S","-")
    assert test_samread_class.pos == 8, "SamRead correct_start_position method making error on reverse strand."
test_correct_start_position()

def test_generate_postrand():
    postrand = test_samread_class.generate_postrand(10,"+") == (10,"+"), "SamRead generate_postrand method making error in tuple generation."
test_generate_postrand()

print("All tests passed!")