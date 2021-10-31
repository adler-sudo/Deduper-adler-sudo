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
read_length = 150


# parse args for testing different test scenarios
# TODO: if we move this to sam_funcs then we can just import it
def parse_args(args):
    parser = argparse.ArgumentParser(description="Reference-based PCR duplicate removal from mapped SAM file")
    parser.add_argument(
        "-f",
        "--input_file", 
        help="Input SAM file - sorted, mapped reads"
    )
    parser.add_argument(
        "-p", 
        "--paired-end", 
        action="store_true", 
        help="Designates file is paired end. If flag is not set, the program will assume single end reads."
    )
    parser.add_argument(
        "-u", 
        "--umi", 
        help="Designates UMI file. If set, UMI file should follow flag. If flag is not set, the program will assume randomer UMI"
    )
    parser.add_argument(
        "-q",
        "--keep_highest_qscore",
        action="store_true",
        help="When flag is set, if the program finds a PCR duplicate, it will keep the duplicate with the highest qscore."
    )
    return parser.parse_args(args)

args = parse_args(sys.argv[1:])

# define additional arguments
if args.input_file:
    base_filename = os.path.basename(args.input_file).split(".")[0]
    input_filename=test_input_dir + "{}.sam".format(base_filename)
    retention_filename=test_output_actual_dir + "{}.retain.sam".format(base_filename)
    # duplicate_filename=test_output_actual_dir + "{}.discard.sam".format(base_filename)

    retention_comp_filename=test_output_expected_dir+"{}.retain.sam".format(base_filename)
    # duplicate_comp_filename=test_output_expected_dir+"{}.discard.sam".format(base_filename)



if args.keep_highest_qscore:
    keep_highest_qscore=args.keep_highest_qscore
else:
    keep_highest_qscore=False
print(keep_highest_qscore)
# define test classes
test_samfile_class = Dedupe(
    input_filename=input_filename,
    umi_filename=umi_filename,
    retention_filename=retention_filename,
    read_length=read_length,
    keep_highest_qscore=keep_highest_qscore
)


# assert that attributes are being assigned correctly
assert test_samfile_class.input_filename == input_filename, "Issue with test_samfile_class input_filename attribute"
assert test_samfile_class.umi_filename == umi_filename, "Issue with test_samfile_class umi_filename attribute"

test_samread_class = SamRead()

def test_evaluate_pair_existence():
    paired_end_dict = {
        'samplernext':{
            'key': (10, '+', 'sampleumi', 1),
            'partner_key': (20, '-', 'sampleumi', 1)
        }
    }
    eval_dict = {
        (10,'+','sampleumi',1):{
            'qscore':39,
            'line': 'sampleline\tsampleline'
        }
    }
    assert test_samfile_class.evaluate_pair_existence(
        qname='samplernext',
        rname='samplercurrent',
        paired_end_dict=paired_end_dict,
        eval_dict=eval_dict
    ) == None, "Error in evaluating pair existence."
test_evaluate_pair_existence()

def test_write_to_paired_end_dict():
    paired_end_dict = {}
    paired_end_dict_check = {
        'samplernext':{
            'key': (10, '+', 'sampleumi', 1),
            'partner_key': (20, '-', 'sampleumi', 1)
        }
    }
    assert test_samfile_class.write_to_paired_end_dict(
        paired_end_dict=paired_end_dict,
        rnext='samplernext',
        key=(10,'+','sampleumi',1),
        partner_key=(20,'-','sampleumi',1)
    ) == paired_end_dict_check, "Paired end dict not being written correctly."
test_write_to_paired_end_dict()

def test_determine_read_num():
    assert test_samread_class.determine_read_num(64) == 'R1', "SamRead method determine_read_num making error in bit 64 of bitwise flag."
    assert test_samread_class.determine_read_num(128) == 'R2', "SamRead method determine_read_num making error in bit 128 of bitwise flag."
test_determine_read_num()

def test_parse_cigar_string():
    assert test_samread_class.parse_cigar_string(cigar='2S') == ([2],['S']), "Parse cigar string function not separating ints and strings correctly."
test_parse_cigar_string()

def test_read_umis():
    umis = [
        'AACGCCAT',
        'AAGGTACG',
        'AATTCCGG', 
        'ACACAGAG', 
        'ACACTCAG', 
        'ACACTGTG', 
        'ACAGGACA']
    test_samfile_class.read_umis()
    assert test_samfile_class.umis == umis, "Error in read_umis function"
test_read_umis()

def test_dedupe():
    test_samfile_class.dedupe()
    if args.input_file:
        assert filecmp.cmp(retention_filename,retention_comp_filename), "Error in {}".format(retention_filename)
test_dedupe()

def test_reset_eval_dict():
    test_eval_dict = {
    }
    test_eval_dict2 = {
        (1,"+"): {'umi':'ATGC','quality':40},
        (2,"+"): {'umi':'ATGC','quality':40},
    }
    assert test_samfile_class.reset_eval_dict(test_eval_dict) == test_eval_dict, "SamFile reset_eval_dict method not properly resetting dictionary."
    assert test_samfile_class.reset_eval_dict(test_eval_dict2) == test_eval_dict, "SamFile reset_eval_dict method not properly resetting dictionary."
test_reset_eval_dict()

def test_parse_columns():
    sample_sam_line = "NS500451:154:HWKTMBGXX:1:11101:21621:1145:AGGTTGCT	0	2	93022350	36	2S89M	*	0	0	TTCCACTGTTGCTTCATAACTGCAGTCCTAACATAAATGTCTGACATGTAGGATGATCTTAAGCAACCCCT	6AEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE<EEAAAEE"
    test_samread_class.set_read(sample_sam_line)
    assert test_samread_class.qname == sample_sam_line.split("\t")[0], "Issue with SamRead method parse_columns"
test_parse_columns()

def test_correct_start_position():
    assert test_samread_class.correct_start_position(10,"2S","+") == 8, "SamRead correct_start_position method making error on forward strand."
    assert test_samread_class.correct_start_position(10,"2M","+") == 10, "SamRead correct_start_position method incorrectly adjusting start position."
    assert test_samread_class.correct_start_position(10,"2S","-") == 10, "SamRead correct_start_position method making error on reverse strand."
    assert test_samread_class.correct_start_position(10,"2S2M2D2M2S","-") == 18, "SamRead correct_start_position method making error on the reverse strand correction."
    assert test_samread_class.correct_start_position(10,"10M","-") == 20, "SamRead correct_start-position method making error on the reverse strand correction."
test_correct_start_position()

def test_generate_postrand():
    assert test_samread_class.generate_postrand(10,"+","ATAT",1) == (10,"+","ATAT",1), "SamRead generate_postrand method generating tuple incorrectly."
test_generate_postrand()

def test_determine_strandedness():
    assert test_samread_class.determine_strandedness(16) == "-", "SamRead determine_strandedness marking strand as forward when it should be reverse."
    assert test_samread_class.determine_strandedness(15) == "+", "SamRead determine_strandedness marking strand as reverse when it should be forward."
test_determine_strandedness()

def test_evaluate_existence():
    test_eval_dict = {
        (10,"+","A"): {
            "umi": "ACACTGTG",
            "qscore": 39
        }
    }
    
    existence = test_samfile_class.evaluate_existence(
        postrand=(10,"+","A"),
        umi="ACACTGTG",
        qscore=39,
        eval_dict=test_eval_dict,
        randomer_umi=True
    )

    assert existence == True, "Error in Dedupe method evaluate_existence. Not recognizing entry already exists in dictionary."

    existence = test_samfile_class.evaluate_existence(
        postrand=(1,"+","A"),
        umi="ACACTGTG",
        qscore=39,
        eval_dict=test_eval_dict,
        randomer_umi=True
    )

    assert existence == False, "Error in Dedupe method evaluate_existence. Making error in position part of eval_dict key."

    existence = test_samfile_class.evaluate_existence(
        postrand=(10,"-","A"),
        umi="ACACTGTG",
        qscore=39,
        eval_dict=test_eval_dict,
        randomer_umi=True
    )

    assert existence == False, "Error in Dedupe method evaluate_existence. Making error in strand evaluation."

    existence = test_samfile_class.evaluate_existence(
        postrand=(10,"-","A"),
        umi="ACACTGTG",
        qscore=39,
        eval_dict=test_eval_dict,
        randomer_umi=False
    )

    assert existence == True, "Error in Dedupe method evaluate_existence. Making error in randomer_umi."

    existence = test_samfile_class.evaluate_existence(
        postrand=(10,"-","A"),
        umi="ACACTGT",
        qscore=39,
        eval_dict=test_eval_dict,
        randomer_umi=False
    )

    assert existence == True, "Error in Dedupe method evaluate_existence. Making error in randomer_umi. May not be recognizing umi is not an umi designated by the umi file."

    existence = test_samfile_class.evaluate_existence(
        postrand=(10,"+","A"),
        umi="ACACTGTG",
        qscore=40,
        eval_dict=test_eval_dict,
        randomer_umi=True
    )

    assert existence == False, "Error in Dedupe method evaluate_existence. Making error in quality evaluation."
    
    existence = test_samfile_class.evaluate_existence(
        postrand=(10,"+","A"),
        umi="ACACTGTG",
        qscore=40,
        eval_dict=test_eval_dict,
        randomer_umi=True,
        keep_highest_qscore=True
    )

    assert existence == True, "Error in Dedupe method evaluate_existence. Incorporating quality check when it should not be."
test_evaluate_existence()

def test_determine_if_mapped():
    assert test_samread_class.determine_if_mapped(4) == False, "SamRead method determine_if_mapped incorrectly assigning bool."
test_determine_if_mapped()



print("All tests passed!")


