#!/usr/bin/env python

# library import
import argparse
import sys
import re
from Bioinfo import convert_phred, qual_score


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
input_filename:str = args.input_file
umi_filename:str = args.umi
retain_filename:str = args.output_retain_file

paired_end:bool = False
keep_highest_qscore:bool = False
randomer_umi:bool = False

# argparser
if args.input_file is None:
    exit("ATTENTION: Exiting script!\nPlease specify a sorted SAM file utilizing the input_file flag (-f, --input_file).")
if args.output_retain_file is None:
    exit("ATTENTION: Exiting script!\nPlease specify a filename to which the uniquely mapped reads will be written.")
if args.paired_end:
    paired_end = True
    print('Paired end SAM file indicated.')
if args.keep_highest_qscore:
    keep_highest_qscore = True
    print('Keeping duplicates with highest qscore.')
if args.umi is None:
    randomer_umi = True
    print('Randomer umi indicated.')
if args.paired_end and args.keep_highest_qscore:
    exit("Paired end highest qscore functionality not yet implemented.")

paired_end_dict:dict = {}
eval_dict:dict = {}
summary_dict:dict = {}

num_reads_retained:int = 0
incorrect_umi:int = 0
total_reads:int = 0

# define functions
def read_umis(
    umi_filename:str=None):
    """
    Reads UMI file and stores in umis instance attribute.

    Parameters:
    -----------
    umi_filename : str
        Filename of UMI. Should exist as .txt file with each UMI on newline.

    Returns:
    --------
    self.randomer_umi : bool
        Returns True when UMI file is NOT defined, and therefore UMIs will be treated as randomers.
        Returns False when UMI file IS defined.

    """
    if umi_filename is not None:
        open_umis = open(umi_filename,'r')
        umis = open_umis.read().split('\n')
        open_umis.close()
    return umis

if args.umi is not None:
    umis = read_umis(umi_filename=umi_filename)


def dedupe(
    input_filename:str=input_filename,
    retention_filename:str=retain_filename,
    summary_dict:dict=summary_dict,
    eval_dict:dict=eval_dict,
    paired_end_dict:dict=paired_end_dict,
    paired_end:bool=paired_end,
    total_reads:int=total_reads,
    num_reads_retained:int=num_reads_retained,
    incorrect_umi:int=incorrect_umi):
    """
    Deduplicates the input SAM file.
    """
    open_input_sam = open(input_filename,'r')
    open_retain_sam = open(retention_filename,'w')
    header_length = 0
    current_chr = '1'

    while True:
        line = open_input_sam.readline()
        if line == '': # break if end of file
            num_reads_retained += len(eval_dict)
            # write last chromosome to summary dict
            summary_dict[current_chr] = len(eval_dict)
                # dump last chromosome to sam
            if paired_end:
                for record in eval_dict:
                    raw_line1 = eval_dict[record]['raw_line1']
                    raw_line2 = eval_dict[record]['raw_line2']
                    open_retain_sam.write(raw_line1)
                    if raw_line2 is not None:
                        open_retain_sam.write(raw_line2)
            elif not paired_end:
                for record in eval_dict:
                    raw_line1 = eval_dict[record]['raw_line1']
                    open_retain_sam.write(raw_line1)
            break
        if line[0] == '@': # store header lines
            open_retain_sam.write(line)
            header_length += 1 
            continue
        total_reads += 1 # made it through all evalutation steps

        # generate all of our values for the current line
        split_line = line.split('\t')
        qname, umi, flag, rname, pos, cigar, rnext, qscore = parse_columns(split_line)
        strand = determine_strandedness(flag)
        cigar_nums, cigar_letters = parse_cigar_string(cigar=cigar)
        corrected_pos = correct_start_position(
            pos=pos,
            cigar=cigar,
            strand=strand,
            cigar_nums=cigar_nums,
            cigar_letters=cigar_letters
        )
        postrand = generate_postrand(
            corrected_pos=corrected_pos,
            strand=strand,
            umi=umi
        )

        # new chromosome endcountered
        if rname != current_chr:
            summary_dict[current_chr] = len(eval_dict)
            current_chr = rname
            
            if paired_end:
                # write unpaired reads to sam
                for unpaired_read in paired_end_dict:
                    open_retain_sam.write(paired_end_dict[unpaired_read]['raw_line1'])
                # write eval dict to sam
                for record in eval_dict:
                    raw_line1 = eval_dict[record]['raw_line1']
                    raw_line2 = eval_dict[record]['raw_line2']
                    open_retain_sam.write(raw_line1)
                    if raw_line2 is not None:
                        open_retain_sam.write(raw_line2)
                num_reads_retained += len(eval_dict) + len(paired_end_dict)
            elif not paired_end:
                for record in eval_dict:
                    raw_line1 = eval_dict[record]['raw_line1']
                    open_retain_sam.write(raw_line1)
                num_reads_retained += len(eval_dict)
            # reset the dictionaries
            eval_dict = {}
            paired_end_dict = {}

        # processing each read
        if paired_end:
            # combine current read with partner if they match
            if rnext in paired_end_dict:
                # then combine and write to eval_dict
                pair_postrand = paired_end_dict[rnext]['postrand']
                # place positive read first in combo
                if postrand[1] == '+':
                    combo_postrand = postrand + pair_postrand
                elif pair_postrand[1] == '+':
                    combo_postrand = pair_postrand + postrand
                # combine and write to eval_dict if doesn't already exist
                if combo_postrand not in eval_dict:
                        qscore2=paired_end_dict[rnext]['qscore']
                        raw_line2=paired_end_dict[rnext]['raw_line']
                        # write to eval_dict
                        eval_dict[postrand] = {
                            'qscore1':qscore,
                            'qscore2':qscore2,
                            'raw_line1':line,
                            'raw_line2':raw_line2
                        }
                # remove pair from paired end dict now that it has been transferred to eval_dict
                paired_end_dict.pop(qname,None)
                paired_end_dict.pop(rnext,None)
            # write instance to the paired_end_dict if its partner doesn't exist
            else:
                if randomer_umi:
                    paired_end_dict[qname] = {
                        'postrand':postrand,
                        'qscore1':qscore,
                        'raw_line1':line
                    }
                elif not randomer_umi:
                    if umi in umis:
                        paired_end_dict[qname] = {
                        'postrand':postrand,
                        'qscore1':qscore,
                        'raw_line1':line
                        }
                    else:
                        incorrect_umi += 1
        elif not paired_end:
            if randomer_umi:
                if postrand not in eval_dict:
                    eval_dict[postrand] = {
                        'qscore1':qscore,
                        'raw_line1':line
                    }
                elif keep_highest_qscore and eval_dict[postrand]['qscore1'] < qscore:
                    eval_dict[postrand] = {
                        'qscore1':qscore,
                        'raw_line1':line
                    }
            elif not randomer_umi:
                if umi in umis:
                    if postrand not in eval_dict:
                        # TODO: this should really be its own function
                        eval_dict[postrand] = {
                            'qscore1':qscore,
                            'raw_line1':line
                        }
                    elif keep_highest_qscore and eval_dict[postrand]['qscore1'] < qscore:
                        eval_dict[postrand] = {
                            'qscore1':qscore,
                            'raw_line1':line
                        }
                        print("higher quality read encountered")
                else:
                    incorrect_umi += 1 # TODO: add unittest for this

    open_input_sam.close()
    open_retain_sam.close()
    
    # finish up summary dict 
    # TODO: turn this into a function
    summary_dict['header_length'] = header_length
    summary_dict['total_reads'] = total_reads
    summary_dict['incorrect_umi'] = incorrect_umi
    summary_dict['reads_retained'] = num_reads_retained
    for entry in summary_dict:
        print(entry,summary_dict[entry],sep=": ")


def parse_columns(
    split_line:list=None):
    qname = split_line[0]
    umi = qname.split(':')[-1]
    flag = int(split_line[1])
    rname = split_line[2]
    pos = int(split_line[3])
    cigar = split_line[5]
    rnext = split_line[6]
    qscore = qual_score(split_line[10])
    return qname, umi, flag, rname, pos, cigar, rnext, qscore

def correct_start_position(
    pos:int=None, 
    cigar:str=None,
    strand:str=None,
    cigar_nums:list=None,
    cigar_letters:list=None):
    """
    Corrects start positiong using cigar string.

    Parameters:
    -----------
    current_position : int
        The position from the SAM file for the current read.
    cigar_string : str
        The cigar string from the SAM file for the current read.
    strand : str
        The strand position from the SAM file for the current read.

    Returns:
    --------
    self.corrected_pos : int
        The corrected position based on the input parameters
        
    """

    # initiate corrected postiion
    corrected_pos = pos

    # unittesting functionality
    if pos is not None:
        corrected_pos = pos
    if cigar is not None:
        cigar = cigar
    if strand is not None:
        strand = strand

    # evaluate if soft clipping present
    cigar_nums, cigar_letters = parse_cigar_string(
        cigar=cigar
    )
    minus_strand_map = {
        'D': +1,
        'I': 0,
        'M': +1,
        'S': +1,
        'N': +1, # TODO: check to make sure this is right
    }

    # amake position correction based on cigar string
    if strand == "+":
        # just have to subtract on the forward
        if cigar_letters[0] == 'S':
            correction = cigar_nums[0]
            corrected_pos = corrected_pos - correction
    elif strand == "-":
        # TODO: need to actually take the whole cigar string into account on the reverse
        if cigar_letters[0] == 'S':
            cigar_letters.pop(0) # ignore the first soft clip if reverse cause this is referring to the 3' end
            cigar_nums.pop(0)
        if len(cigar_letters) > 0:
            cigar_letters = [*map(minus_strand_map.get,cigar_letters)]
            cigar_nums = [map_val * distance for map_val, distance in zip(cigar_letters,cigar_nums)]
            correction = sum(list(map(int,cigar_nums))) # convert to ints before summing
            corrected_pos = corrected_pos + correction # could put minus one here, but as long as it's consistent it doesn't matter
    return corrected_pos

def parse_cigar_string(
    cigar:str=None):
    """
    Parse cigar string into corresponding numbers and letters.

    Parameters:
    -----------
    cigar : str
        Cigar string from current SAM file read.

    Returns:
    --------
    cigar_nums : list
        The parsed out cigar numbers
    cigar_letters : list
        The parsed out cigar letters
    """
    # parse into components
    cigar_nums = re.split("[A-Z]",cigar)
    cigar_nums.remove("")
    cigar_nums = [int(num) for num in cigar_nums]
    cigar_letters = re.split("[0-9]+",cigar)
    cigar_letters.remove("")
    return cigar_nums, cigar_letters

def generate_postrand(
    corrected_pos:int=None,
    strand:str=None,
    umi:str=None):
    """
    generates position, read combo for current read. To be used
    in evaluating whether duplicate or not.
    """
    postrand = (corrected_pos, strand, umi)
    return postrand

def determine_strandedness(
    flag:int=None):
    """
    determines strandedness by looking at 16 of bitwise flag
    """
    if ((flag & 16) == 16):
        strand = "-" # reverse
    else:
        strand = "+" # forward
    return strand

dedupe(
    input_filename=input_filename,
    retention_filename=retain_filename,
    summary_dict=summary_dict,
    eval_dict=eval_dict,
    paired_end_dict=paired_end_dict,
    paired_end=paired_end
)

