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

def write_to_paired_end_dict(
    paired_end_dict:dict=None,
    qname:str=None,
    postrand:tuple=None,
    qscore:float=None,
    raw_line:str=None):
            
    """
    Write to the paired end dictionary.

    Parameters:
    -----------
    paired_end_dict : dict
        The paired-end dictionary.

        Structure:
        ----------
        {
            qname: {
                'postrand':postrand,
                'qscore1':qscore,
                'raw_line1':raw_line 
            }

    qname : str
        The name of the current read.
    qscore : float
        The quality score of the current read.
    raw_line : str
        The raw sam file line of the current read.

    Returns:
    --------
    paired_end_dict : dict
        The paired-end dictionary. See structure above.

    """
    # unittest functionality
    if paired_end_dict is not None:
        paired_end_dict = paired_end_dict
    # write to paired end
    paired_end_dict[qname] = {
        'postrand':postrand,
        'qscore1':qscore,
        'raw_line1':raw_line
    }
    return paired_end_dict

def reset_paired_end_dict(
    paired_end_dict:dict=None):
    """
    Empties the dictionary.

    Parameters:
    -----------
    paired_end_dict : dict
        The dict to be emptied.

    Returns:
    --------
    self.paired_end_dict : dict
        The empty dictionary.

        The structure of the dictionary is:
        {}
    """
    # unittest functionality
    if paired_end_dict is not None:
        paired_end_dict = paired_end_dict
    paired_end_dict = {}
    return paired_end_dict

def write_to_eval_dict(
    postrand:tuple=None,
    qscore1:float=None,
    qscore2:float=None,
    raw_line1:str=None,
    raw_line2:str=None,
    eval_dict:dict=None):

    """
    Writes the current read to the eval_dict. Could be making a new entry or overwriting an existing entry.

    Parameters:
    -----------
    postrand : tuple
        A tuple containing cigar-corrected-position, stand, and umi.
    qscore : float
        The quality score of the current read.
    raw_line : str
        The raw, tab-delimited line.
    eval_dict : dict
        Evaluation dict containing record of reads that have been written to the retention file.

        The structure of the dictionary depends on whether the reads were produced from paired or single end sequencing:
            
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

    Returns:
    --------
    self.eval_dict : dict
        Returns the eval_dict. Really only for the purpose of unittesting.S
    """
    # unittest functionality
    if eval_dict is not None:
        eval_dict = eval_dict

    # determine actions based on single or paired end
    if paired_end:
        eval_dict[postrand] = {
            'qscore1':qscore1,
            'qscore2':qscore2,
            'raw_line1':raw_line1,
            'raw_line2':raw_line2
        }
    elif not paired_end:
        # make eval_dict addition
        eval_dict[postrand] = {'qscore1':qscore1,'raw_line1':raw_line1}
    return eval_dict

def reset_eval_dict(
    eval_dict:dict=None):
    """
    Empties the dictionary.

    Parameters:
    -----------
    eval_dict : dict
        The dict to be emptied.

        The structure of the dictionary is:
            {
                (self.correct_pos, self.strand, self.umi): {
                    'qscore1': float
                    'raw_line1': raw_line
                }
            }

    Returns:
    --------
    self.eval_dict : dict
        The empty dictionary.

        The structure of the dictionary is:
        {}
    """
    # unittest functionality
    if eval_dict is not None:
        eval_dict = eval_dict
    eval_dict = {}
    return eval_dict

def evaluate_pair_existence(
    paired_end_dict:dict=None,
    umi:str=None,
    qname:str=None,
    rnext:str=None,
    postrand:str=None,
    qscore:float=None,
    raw_line:str=None,
    eval_dict:dict=None,
    randomer_umi:str=None):
    """
    Evaluates if pair has already been found. Performs proper steps in reading and writing to the paired_end_dict.
    
    Parameters:
    -----------
    paired_end_dict : dict
        The paired-end dictionary that holds the connection between the two pairs

        Structure:
        ----------
        See 'write_to_paired_end_dict' for structure.
    qname : str
        The qname of the current read being evaluated, which will be the KEY in the paired-end dict because the rname of the pair was used to generate the key.
    rnext : str
        The qname of the PAIR.=
    postrand : tuple
        The postrand of the current read.
    qscore : float
        The qscore of the current read.
    raw_line : str
        The raw_line of the current read.
    eval_dict : dict
        The evaluation dictionary. NOTE: may not need this, as the eval_dict functions will be carried out in the 'evaluate_existence' function.

    Returns:
    --------
    None
    """
    # unittest functionality
    if paired_end_dict is not None:
        paired_end_dict = paired_end_dict
    if eval_dict is not None:
        eval_dict = eval_dict
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
            write_to_eval_dict(
                postrand=combo_postrand,
                qscore1=qscore,
                qscore2=paired_end_dict[rnext]['qscore'],
                raw_line1=raw_line,
                raw_line2=paired_end_dict[rnext]['raw_line'],
                eval_dict=eval_dict
            )
        # remove pair from paired end dict now that it has been transferred to eval_dict
        paired_end_dict.pop(qname,None)
        paired_end_dict.pop(rnext,None)
    # write instance to the paired_end_dict if its partner doesn't exist
    elif rnext not in paired_end_dict:
        write_to_paired_end_dict(
            paired_end_dict=paired_end_dict,
            qname=qname,
            postrand=postrand,
            qscore=qscore,
            raw_line=raw_line
        )
    return paired_end_dict, eval_dict

def evaluate_existence(
    postrand:tuple=None,
    umi:str=None,
    qscore:float=None,
    raw_line:str=None,
    qname:str=None,
    rname:str=None,
    eval_dict:dict=None,
    randomer_umi:bool=None,
    keep_highest_qscore:bool=None,
    incorrect_umi:int=None):
    """
    Evaluates if read already exists in eval_dict. Writes read to dictionary and output file if it does not already exist.

    Parameters:
    -----------
    postrand : tuple
        A tuple containing cigar-corrected-position, stand, and umi.

        Example:
        --------

        (100, "+", "ATGCATGC", 1)

    umi : str
        Contains the UMI sequence for the current read.
    qscore : float
        The quality score of the current read.
    raw_line : str
        The raw, tab-delimited line.
    eval_dict : dict
        Evaluation dict containing record of reads that have been written to the retention file.

        The structure of the dictionary is:
            {
                (self.correct_pos, self.strand, self.umi): {
                    'qscore1': float,
                    'raw_line1': raw_line
                }
            }

    randomer_umi : bool
        Set to True when no UMI file is specified.
        Set to False when UMI file is specified.
    keep_highest_score : bool
        When set to True, will keep the highest qscore in the case of PCR duplicates.
        When set to False, will keep the first read encountered in the case of PCR duplicates.

    Returns:
    --------
    """
    if eval_dict is not None:
        eval_dict = eval_dict
    if randomer_umi is not None:
        randomer_umi = randomer_umi
    if keep_highest_qscore is not None:
        keep_highest_qscore = keep_highest_qscore
    
    # TODO: think i can make this simpler but not seeing it now
    if randomer_umi:
        if postrand not in eval_dict:
            # TODO: this should really be its own function
            eval_dict = write_to_eval_dict(
                    postrand=postrand,
                    qscore1=qscore,
                    raw_line1=raw_line)
        elif keep_highest_qscore and eval_dict[postrand]['qscore1'] < qscore:
            eval_dict = write_to_eval_dict(
                    postrand=postrand,
                    qscore1=qscore,
                    raw_line1=raw_line)
    elif not randomer_umi:
        if umi in umis:
            if postrand not in eval_dict:
                # TODO: this should really be its own function
                eval_dict = write_to_eval_dict(
                    postrand=postrand,
                    qscore1=qscore,
                    raw_line1=raw_line,
                    eval_dict=eval_dict)
            elif keep_highest_qscore and eval_dict[postrand]['qscore1'] < qscore:
                eval_dict = write_to_eval_dict(
                    postrand=postrand,
                    qscore1=qscore,
                    raw_line1=raw_line,
                    eval_dict=eval_dict)
                print("higher quality read encountered")
        elif umi not in umis:
            incorrect_umi += 1 # TODO: add unittest for this
    return eval_dict, incorrect_umi

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
            dump_dict_to_sam(
                open_retain_sam=open_retain_sam,
                eval_dict=eval_dict
            )
            break
        if line[0] == '@': # store header lines
            write_to_retain(
                open_retain_sam=open_retain_sam,
                line=line
            )
            header_length += 1 
            continue
        total_reads += 1 # made it through all evalutation steps

        # TODO: all of the SamRead class functions will exist here
        split_line = set_read(
            raw_line=line
        )
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

        # reset dictionary if we've reached a new chromosome
        if rname != current_chr:
            print("Initiating chromosome {} processing.".format(rname))
            summary_dict[current_chr] = len(eval_dict)
            current_chr = rname
            # write unpaired reads to eval dict
            if paired_end:
                for unpaired_read in paired_end_dict:
                    eval_dict[unpaired_read]['rawline1'] = paired_end_dict[unpaired_read]['raw_line1']
            dump_dict_to_sam(
                open_retain_sam=open_retain_sam,
                paired_end=paired_end,
                eval_dict=eval_dict
            )
            # TODO: need to fix how this counts for paired end
            num_reads_retained += len(eval_dict)
            # reset the dictionaries
            eval_dict = {}
            paired_end_dict = {}

        # if paired-end follow paired-end existence
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
            elif rnext not in paired_end_dict:
                    paired_end_dict[qname] = {
                        'postrand':postrand,
                        'qscore1':qscore,
                        'raw_line1':line
                    }
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
                elif umi not in umis:
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

def dump_dict_to_sam(
    open_retain_sam:None,
    paired_end:bool=None,
    eval_dict:dict=None):
    """
    write the eval_dict to the retain sam file
    """
    # determine action based on single or paired end
    if paired_end:
        for record in eval_dict:
            raw_line1 = eval_dict[record]['raw_line1']
            raw_line2 = eval_dict[record]['raw_line2']
            write_to_retain(raw_line1)
            if raw_line2 is not None:
                write_to_retain(
                    open_retain_sam=open_retain_same,
                    line=raw_line2
                )
    elif not paired_end:
        for record in eval_dict:
            line = eval_dict[record]['raw_line1']
            write_to_retain(
                open_retain_sam=open_retain_sam,
                line=line
            )

def write_to_retain(
    open_retain_sam=None,
    line:str=None):
    """
    Writes a line to the retain file.

    Parameters:
    -----------
    open_retain_sam : wrapper object for open sam retain file
        Opened in dedupe function.
    line : str
        The raw line to be written.

    Returns:
    --------
    None
    """
    open_retain_sam.write(line)
    return None

def set_read(
    raw_line:str=None):
    """
    sets the sam read
    """
    if raw_line is not None:
        line = raw_line.split('\t')
        return line

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

    # unittest functionality
    if corrected_pos is not None:
        corrected_pos = corrected_pos
    if strand is not None:
        strand = strand
    if umi is not None:
        umi = umi

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

