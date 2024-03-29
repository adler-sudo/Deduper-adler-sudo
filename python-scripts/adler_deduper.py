#!/usr/bin/env python

# library import
import argparse
import sys
import re
from Bioinfo import convert_phred, qual_score
import csv
import logging


# define arguments
def parse_args():
    parser = argparse.ArgumentParser(description="Reference-based PCR duplicate removal from mapped SAM file")
    parser.add_argument(
        "infile",
        nargs = '?',
        type = argparse.FileType('r'),
        help = "Input SAM file - sorted, mapped reads",
        default = sys.stdin
    )
    parser.add_argument(
        "outfile",
        nargs = '?',
        type = argparse.FileType('w'),
        help="Output retain SAM file containing uniquely mapped reads.",
        default = sys.stdout
    )
    parser.add_argument(
        "-p",
        "--paired",
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
    parser.add_argument(
        "-s",
        "--stats_file",
        help="Path and file name of duplicate stats file"
    )
    parser.add_argument(
        "-b",
        "--barcode_position",
        help = "Position of barcode in qname. 1 = barcode before colon. 2 = barcode after colon.",
        type = int
    )
    args = parser.parse_args()
    return args

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

def dedupe(
    input_filename:str,
    retention_filename:str,
    summary_dict:dict,
    eval_dict:dict,
    paired_end_dict:dict,
    duplicate_dict:dict,
    paired_end:bool,
    total_reads:int,
    num_reads_retained:int,
    incorrect_umi:int,
    randomer_umi:bool = True,
    barcode_position:int = 1,
    keep_highest_qscore:bool = False):
    """
    Deduplicates the input SAM file.

    Parameters:
    -----------
    input_filename : str
        The sorted input SAM file to be deduped.
    retention_filename : str
        The SAM file to write the uniquely mapped reads to.
    summary_dict : dict
        A dictionary to hold summary items. Chromosome, num_headers, lines retained, incorrect umis, and total reads.
    eval_dict : dict
        The dictionary used to hold reads that will be written to retnetion file at the end of each chromosome.

        Structure:
        ----------
        SINGLE-END:
                {
                    (self.correct_pos, self.strand, self.umi): {
                        'qscore1': float
                        'raw_line1': raw_line
                    }
                }

        PAIRED-END
        {
            (+ strand position, +, + strand umi, - strand position, -, - strand umi):{
                'qscore1': qscore + strand
                'qscore2': qscore - strand
                'raw_line1': raw_line + strand
                'raw_line2': raw_line - strand
            }
        }

    paired_end_dict : dict
        Only utilized in paired end deduplication. Holds on to reads that have not yet been paired.
        We add a '1' to the end of the qname for use in downstream differentiation of the pair in the case that they have the same qname.

        Structure:
        ----------
        {
            qname + '1': {
                'postrand':postrand,
                'qscore1':qscore,
                'raw_line1':raw_line
            }

    duplicate_dict : dict
        Tracks number of duplicates encountered for each umi

        Structure:
        ----------
        {
            umi: {
                'duplicate_reads': int,
                'total_reads': int
            }
        }
    paired_end : bool
        When marked True, indicates that input file is from paired end sequencing.
    total_reads : int
        Total reads in the file.
    num_reads_retained : int
        Number of reads retained after deduplication is complete.
    incorrect_umi : int
        Number of incorrect umis found and discarded.
    """
    # open_input_sam = open(input_filename,'r')
    # open_retain_sam = open(retention_filename,'w')
    header_length = 0
    current_chr = '1'

    while True:
        line = input_filename.readline()
        if line == '': # break if end of file
            num_reads_retained += len(eval_dict)
            # write last chromosome to summary dict
            summary_dict[current_chr] = len(eval_dict)
            # dump last chromosome to sam
            if paired_end:
                for record in eval_dict:
                    raw_line1 = eval_dict[record]['raw_line1']
                    raw_line2 = eval_dict[record]['raw_line2']
                    retention_filename.write(raw_line1)
                    if raw_line2 is not None:
                        retention_filename.write(raw_line2)
            elif not paired_end:
                for record in eval_dict:
                    raw_line1 = eval_dict[record]['raw_line1']
                    retention_filename.write(raw_line1)
            break
        # process header lines
        if line[0] == '@':
            retention_filename.write(line)
            header_length += 1
            continue
        total_reads += 1 # made it through all evalutation steps

        # generate all of our values for the current read
        split_line = line.split('\t')
        qname, umi, flag, rname, pos, cigar, rnext, qscore = parse_columns(split_line, barcode_position)
        strand = determine_strandedness(flag)
        cigar_nums, cigar_letters = parse_cigar_string(cigar=cigar)
        corrected_pos = correct_start_position(
            pos=pos,
            cigar=cigar,
            strand=strand,
            cigar_nums=cigar_nums,
            cigar_letters=cigar_letters
        )
        # generate the unique key for this read
        postrand = generate_postrand(
            corrected_pos=corrected_pos,
            strand=strand,
            umi=umi
        )

        # add umi if not in duplicate dict
        if umi not in duplicate_dict:
            duplicate_dict = barcode_to_duplicate_dict(
                duplicate_dict = duplicate_dict,
                umi = umi)

        # new chromosome encountered
        if rname != current_chr:
            summary_dict[current_chr] = len(eval_dict)
            current_chr = rname

            if paired_end:
                # # write unpaired reads to sam
                # for unpaired_read in paired_end_dict:
                #     open_retain_sam.write(paired_end_dict[unpaired_read]['raw_line1'])
                # write eval dict to sam
                for record in eval_dict:
                    raw_line1 = eval_dict[record]['raw_line1']
                    raw_line2 = eval_dict[record]['raw_line2']
                    retention_filename.write(raw_line1)
                    if raw_line2 is not None:
                        retention_filename.write(raw_line2)
                num_reads_retained += len(eval_dict) + len(paired_end_dict)
            elif not paired_end:
                for record in eval_dict:
                    raw_line1 = eval_dict[record]['raw_line1']
                    retention_filename.write(raw_line1)
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
                else:
                    combo_postrand = pair_postrand + postrand
                # write pair to eval_dict if combo key doesn't already exist
                if combo_postrand not in eval_dict:
                    qscore2=paired_end_dict[rnext]['qscore1']
                    raw_line2=paired_end_dict[rnext]['raw_line1']
                    # write to eval_dict
                    eval_dict[combo_postrand] = {
                        'qscore1':qscore,
                        'qscore2':qscore2,
                        'raw_line1':line,
                        'raw_line2':raw_line2
                    }
                    duplicate_dict = increment_duplicate_dict(
                            duplicate_dict = duplicate_dict,
                            umi = umi,
                            duplicate = False,
                            paired_end = paired_end
                        )
                else:
                    duplicate_dict = increment_duplicate_dict(
                            duplicate_dict = duplicate_dict,
                            umi = umi,
                            duplicate = True,
                            paired_end = paired_end
                        )
                # remove pair from paired end dict now that it has been transferred to eval_dict
                paired_end_dict.pop(qname,None)
                paired_end_dict.pop(rnext,None)
            # write instance to the paired_end_dict if its partner doesn't exist
            else:
                qname = qname + '1'
                if randomer_umi:
                    paired_end_dict[qname] = {
                        'postrand':postrand,
                        'qscore1':qscore,
                        'raw_line1':line,
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
                    duplicate_dict = increment_duplicate_dict(
                            duplicate_dict = duplicate_dict,
                            umi = umi,
                            duplicate = False,
                            paired_end = paired_end
                        )
                # evaluate whether qscore is higher than current stored qscore
                elif keep_highest_qscore and eval_dict[postrand]['qscore1'] < qscore:
                    eval_dict[postrand] = {
                        'qscore1':qscore,
                        'raw_line1':line
                    }
                    duplicate_dict = increment_duplicate_dict(
                            duplicate_dict = duplicate_dict,
                            umi = umi,
                            duplicate = True,
                            paired_end = paired_end
                        )
                else:
                    duplicate_dict = increment_duplicate_dict(
                            duplicate_dict = duplicate_dict,
                            umi = umi,
                            duplicate = True,
                            paired_end = paired_end
                        )
            elif not randomer_umi:
                if umi in umis:
                    if postrand not in eval_dict:
                        eval_dict[postrand] = {
                            'qscore1':qscore,
                            'raw_line1':line
                        }
                        duplicate_dict = increment_duplicate_dict(
                            duplicate_dict = duplicate_dict,
                            umi = umi,
                            duplicate = False,
                            paired_end = paired_end
                        )
                    # evaluate whether qscore is higher than current stored qscore
                    elif keep_highest_qscore and eval_dict[postrand]['qscore1'] < qscore:
                        eval_dict[postrand] = {
                            'qscore1':qscore,
                            'raw_line1':line
                        }
                        duplicate_dict = increment_duplicate_dict(
                            duplicate_dict = duplicate_dict,
                            umi = umi,
                            duplicate = True,
                            paired_end = paired_end
                        )
                    else:
                        duplicate_dict = increment_duplicate_dict(
                            duplicate_dict = duplicate_dict,
                            umi = umi,
                            duplicate = True,
                            paired_end = paired_end
                        )
                else:
                    incorrect_umi += 1

    # open_input_sam.close()
    # open_retain_sam.close()

    # finish up and print summary dict
    summary_dict['header_length'] = header_length
    summary_dict['total_reads'] = total_reads
    summary_dict['incorrect_umi'] = incorrect_umi
    summary_dict['reads_retained'] = num_reads_retained
    for entry in summary_dict:
        # print(entry,summary_dict[entry],sep=": ")
        statement = "{}: {}".format(entry, summary_dict[entry])
        logging.info(statement)

def parse_columns(
    split_line:list=None,
    barcode_position:int=1):
    """"
    Splits the current SAM file line into its component pieces.

    Parameters:
    -----------
    split_line : list
        The split line represented as a list. Each column of the current read is represented as a component of the list.

    Returns:
    --------
    qname : str
        The qname of the current SAM file read.
    umi : str
        The umi of the current SAM file read.
    flag : int
        The bitwise flag of the current SAM file read.
    rname : str
        The chromosome of the current SAM file read.
    pos : int
        The leftmost mapped position for the current SAM file read.
    cigar : str
        The CIGAR string for the current SAM file read.
    rnext : str
        The name of the pair in the case paired end reads.
    qscore : float
        The average qscore fo the current SAM file read.
    """
    qname = split_line[0]
    umi = qname.split(':')[barcode_position - 1] # -1 accounts for python zero indexing
    flag = int(split_line[1])
    rname = split_line[2]
    pos = int(split_line[3])
    cigar = split_line[5]
    # rnext sometimes equal to qname
    rnext = split_line[6]
    if rnext == '=':
        rnext = qname + '1' # add a one for downstream differentiation
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

    # split the cigar string into nums and letters
    cigar_nums, cigar_letters = parse_cigar_string(
        cigar=cigar
    )

    # TODO: add hard clipping
    minus_strand_map = {
        'D': +1,
        'I': 0,
        'M': +1,
        'S': +1,
        'N': +1,
        'H': +1
    }

    # amake position correction based on cigar string
    if strand == "+":
        # subtract the first soft clip on the forward
        if cigar_letters[0] == 'S' or cigar_letters[0] == 'H': # TODO: need to evaluate the hard clipping
            correction = cigar_nums[0]
            corrected_pos = corrected_pos - correction
    elif strand == "-":
        if cigar_letters[0] == 'S':
            cigar_letters.pop(0) # ignore the first soft clip if reverse cause this is referring to the 3' end
            cigar_nums.pop(0)
        if len(cigar_letters) > 0:
            try:
                # determine the correction quanitity from the cigar string
                cigar_letters = [*map(minus_strand_map.get,cigar_letters)]
                cigar_nums = [map_val * distance for map_val, distance in zip(cigar_letters,cigar_nums)]
                correction = sum(list(map(int,cigar_nums))) # convert to ints before summing
                corrected_pos = corrected_pos + correction # could put minus one here, but as long as it's consistent it doesn't matter
            except TypeError:
                logging.info('TypeError in CIGAR string. Please investigate.')
                pass
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
    if "" in cigar_nums:
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

    Parameters:
    -----------
    corrected_pos : int
        The cigar corrected start position of the current read.
    strand : str
        + for forward strand, - for the reverse strand
    umi : str
        The umi for the current read.

    Returns:
    --------
    postrand : tuple
        A tuple to be used as the reads unique key.

        Structure:
        ----------
        (corrected_pos, strand, umi)
    """
    postrand = (corrected_pos, strand, umi)
    return postrand

def determine_strandedness(
    flag:int=None):
    """
    Determines strandedness by looking at 16 of bitwise flag.

    Parameters:
    -----------
    flag : in
        The bitwise flag for the current read

    Returns:
    --------
    strand : str
        + if forward strand, - if reverse strand
    """
    # use the bit 16 to determine strandednes
    if ((flag & 16) == 16):
        strand = "-" # reverse
    else:
        strand = "+" # forward
    return strand

def barcode_to_duplicate_dict(
    duplicate_dict:dict,
    umi:str) -> dict:
    """
    Adds new umi to duplicate dict
    """
    duplicate_dict[umi] = {'duplicate_reads': 0, 'total_reads': 0}
    return duplicate_dict

def increment_duplicate_dict(
    duplicate_dict:dict,
    umi:str,
    duplicate:bool,
    paired_end:bool) -> dict:
    """
    Updates the duplicate counter dictionary.

    Note: anytime the eval_dict is updated, the duplicate_dict must be updated too

    {'umi': {'duplicate_reads': 0, 'total_reads': 0}}
    """
    if paired_end:
        increment = 2
    elif not paired_end:
        increment = 1

    duplicate_dict[umi]['total_reads'] += increment
    if duplicate:
        duplicate_dict[umi]['duplicate_reads'] += increment
    return duplicate_dict

# def output_duplicate_dict(
#     duplicate_dict:dict,
#     duplicate_dict_file:str) -> None:
#     """
#     Writes duplicate dict to file.
#     """
#     with open(duplicate_dict_file, 'w') as f:
#         f.write("barcode\tduplicates\ttotal\tperc_duplicates\n")
#         for barcode in duplicate_dict:
#             duplicate_reads = duplicate_dict[barcode]['duplicate_reads']
#             total_reads = duplicate_dict[barcode]['total_reads']
#             try:
#                 perc_dupe = round(duplicate_reads / total_reads * 100,2)
#             except ZeroDivisionError:
#                 perc_dupe = 0.00
#             f.write("%s\t%s\t%s\t%s\n" % (barcode, duplicate_reads, total_reads, perc_dupe))
#     return None

def main():

    args = parse_args()

    # initiate files
    umi_filename:str = args.umi

    # initatie optional input arguments
    paired_end:bool = False
    keep_highest_qscore:bool = False
    randomer_umi:bool = False

    # argparser
    if args.paired:
        paired_end = True
        logging.info('Paired end SAM file indicated.')
    if args.keep_highest_qscore:
        keep_highest_qscore = True
        logging.info('Keeping duplicates with highest qscore.')
    if args.umi is None:
        randomer_umi = True
        logging.info('Randomer umi indicated.')
    if args.paired and args.keep_highest_qscore:
        sys.exit("ATTENTION: Exiting script!\nPaired end highest qscore functionality not yet implemented.")
    # if args.stats_file is not None:
    #     duplicate_dict_file = args.stats_file
    #     print("Duplicate stats file:", duplicate_dict_file)

    # initiate all dictionaries
    paired_end_dict:dict = {} # paired_end_dict = {qname + '1': {'postrand': (pos, strand, umi), 'qscore1': qscore, 'raw_line1': rawline}}
    eval_dict:dict = {} # eval_dict = {(pos,strand,umi): {'qscore1': qscore R1, 'qscore2': qscore R2, 'raw_line1': rawline R1, 'raw_line2': rawline R2}}
    # eval_dict key changes for PE - combo of R1 and R2 postrand: (pos1, strand1, umi1, pos2, strand2, umi2)

    summary_dict:dict = {} # summary_dict counts number of occurrences of each chromosome, headers, lines retained, incorrect umis, and total reads
    # summary dict is printed to the console at completion of the run

    duplicate_dict:dict = {} # duplicate_dict keeps track of duplicatesw by barcode
    # structure {umi: {'duplicate_reads': int, 'total_reads': int}}

    # initiate counters for summary
    num_reads_retained:int = 0
    incorrect_umi:int = 0
    total_reads:int = 0

    randomer_umi = True

    # read in umi file if umis specified
    if args.umi is not None:
        umis = read_umis(umi_filename=umi_filename)

    # run dedupe
    dedupe(
        input_filename=args.infile,
        retention_filename=args.outfile,
        summary_dict=summary_dict,
        eval_dict=eval_dict,
        paired_end_dict=paired_end_dict,
        duplicate_dict=duplicate_dict,
        paired_end=paired_end,
        total_reads = total_reads,
        num_reads_retained = num_reads_retained,
        incorrect_umi = incorrect_umi,
        barcode_position = args.barcode_position
    )

    # output_duplicate_dict(
    #     duplicate_dict = duplicate_dict,
    #     duplicate_dict_file = duplicate_dict_file
    # )

if __name__ == '__main__':
    logging.basicConfig(
        format = '%(asctime)s: %(levelname)s: %(message)s',
        level = logging.INFO
    )
    main()

