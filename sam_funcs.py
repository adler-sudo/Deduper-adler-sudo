from Bioinfo import convert_phred, qual_score
import re


# define clases
class Dedupe():
    """
    A class to remove PCR duplicates from a SAM file.

    Parameters
    ----------
    input_filename : str
        Defines the input filename. Must be sorted SAM format.

    umi_filename : str
        Defines the input file containing all UMIs. Txt file with each line representing a new umi.

    retention_filename : str
        Defines the output file with PCR duplicates removed. Should be SAM file. 

    read_length : int
        Length of each read.

    keep_highest_qscore : bool, default False
        When marked as True, will keep highest quality read in the case of PCR duplicate. Default activity is to keep first read encountered.
    """
    def __init__(
        self, 
        input_filename:str=None, 
        umi_filename:str=None, 
        retention_filename:str=None,
        read_length:int=150,
        keep_highest_qscore:bool=False):

        self.input_filename = input_filename
        self.umi_filename = umi_filename
        self.retention_filename = retention_filename
        self.read_length = int(read_length)
        
        self.total_reads = 0
        self.num_reads_retained = 0
        self.num_reads_discarded = 0
        self.current_chr = 1
        self.eval_dict = {}

        self.keep_highest_qscore = keep_highest_qscore

        self.umis_read = False
        self.randomer_umi = True
        self.umis = None

        self.deduped = False
        
        self.read_umis()

    def read_umis(
        self,
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
        if self.umi_filename is not None:
            open_umis = open(self.umi_filename,'r')
            self.umis = open_umis.read().split('\n')
        # randomer umis
        self.randomer_umi = True
        return self.randomer_umi

    def write_to_retain(
        self,
        line:str=None):
        """
        Writes a line to the retain file.

        Parameters:
        -----------
        
        retain_filename : str
            Defines the path to the deduped SAM file.

        line : str
            The raw line to be written.

        Returns:
        --------

        None

        """


        self.open_retain_sam.write(line)
        return None

    def evaluate_existence(
        self,
        postrand:tuple=None,
        umi:str=None,
        qscore:float=None,
        raw_line:str=None,
        eval_dict:dict=None,
        randomer_umi:bool=None,
        keep_highest_qscore:bool=None):
        """
        Evaluates if read already exists in eval_dict. Writes read to dictionary and output file if it does not already exist.

        Parameters:
        -----------

        postrand : tuple
            A tuple containing cigar-corrected-position, stand, and umi.

            Example:
            --------

            (100, "+", "ATGCATGC")

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
                        'qscore': float
                        'line': raw_line
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
        True when read already exists in eval_dict OR if UMI is not a valid UMI in the case that the UMI file is specified.
        False when read does not exist in eval_dict OR if current read is a PCR duplicate BUT has a higher qscore than the current resident.

        """
        if eval_dict is not None:
            self.eval_dict = eval_dict
        if randomer_umi is not None:
            self.randomer_umi = randomer_umi
        if keep_highest_qscore is not None:
            self.keep_highest_qscore = keep_highest_qscore
        
        # TODO: think i can make this simpler but not seeing it now
        if self.randomer_umi:
            if postrand not in self.eval_dict:
                # TODO: this should really be its own function
                self.write_to_eval_dict(
                        postrand=postrand,
                        qscore=qscore,
                        raw_line=raw_line)
            elif self.keep_highest_qscore and self.eval_dict[postrand]['qscore'] < qscore:
                self.write_to_eval_dict(
                        postrand=postrand,
                        qscore=qscore,
                        raw_line=raw_line)
            else:
                return True
        else:
            if umi in self.umis:
                if postrand not in self.eval_dict:
                    # TODO: this should really be its own function
                    self.write_to_eval_dict(
                        postrand=postrand,
                        qscore=qscore,
                        raw_line=raw_line)
                elif self.keep_highest_qscore and eval_dict[postrand]['qscore'] < qscore:
                    self.write_to_eval_dict(
                        postrand=postrand,
                        qscore=qscore,
                        raw_line=raw_line)
                    print("higher quality read encountered")
                else:
                    return True
            else:
                return True
        return False

    def dump_dict_to_sam(
        self,
        eval_dict=None):
        """
        write the eval_dict to the retain sam file
        """

        if eval_dict is not None:
            self.eval_dict = eval_dict

        for record in self.eval_dict:
            line = self.eval_dict[record]['line']
            self.write_to_retain(line)

    def write_to_eval_dict(
        self,
        postrand:tuple=None,
        qscore:float=None,
        raw_line:str=None,
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

            The structure of the dictionary is:
                {
                    (self.correct_pos, self.strand, self.umi): {
                        'qscore': float
                        'line': raw_line
                    }
                }

        Returns:
        --------

        self.eval_dict : dict
            Returns the eval_dict. Really only for the purpose of unittesting.S


        """
        # unittest functionality
        if eval_dict is not None:
            self.eval_dict = eval_dict
        # make eval_dict addition
        self.eval_dict[postrand] = {'qscore':qscore,'line':raw_line}
        return self.eval_dict

    def dedupe(self):
        """
        Deduplicataes the input SAM file.
        """
        self.open_input_sam = open(self.input_filename,'r')
        self.open_retain_sam = open(self.retention_filename,'w')
        # self.open_discard_sam = open(self.duplicate_filename,'w')

        while True:
            line = self.open_input_sam.readline()
            if line == '': # break if end of file
                self.dump_dict_to_sam() # dump last chromosome to sam
                self.num_reads_retained += len(self.eval_dict)
                break
            if line[0] == '@': # store header lines
                self.write_to_retain(
                    line=line
                )
                # self.write_to_discard(line)
                continue
            self.total_reads += 1 # made it through all evalutation steps

            read = SamRead(
                raw_line=line
            )
            read.correct_start_position()
            read.generate_postrand()

            # reset dictionary if we've reached a new chromosome
            if int(read.rname) != int(self.current_chr):
                self.current_chr = int(read.rname)
                self.dump_dict_to_sam(
                    eval_dict=self.eval_dict
                )
                self.reset_eval_dict(
                    eval_dict=self.eval_dict
                )

            # if erroneous umi, write to discard 
            # TODO: challenge is to make correction 
            self.evaluate_existence(
                postrand=read.postrand, 
                umi=read.umi, 
                qscore=read.qscore,
                raw_line=read.raw_line,
                eval_dict=self.eval_dict, 
                randomer_umi=self.randomer_umi
            )

        self.open_input_sam.close()
        self.open_retain_sam.close()
        # self.open_discard_sam.close()

        self.deduped = True

    def evaluate_qscore(
        self):
        return None

    def reset_eval_dict(
        self,
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
                        'qscore': float
                        'line': raw_line
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
            self.eval_dict = eval_dict

        # self.eval_dict = {k:[] for k in self.eval_dict}
        # return self.eval_dict
        self.eval_dict = {}
        return self.eval_dict

class SamRead:
    def __init__(
        self,
        raw_line=None):

        self.raw_line = None
        self.line = None
        self.qname = None
        self.umi = None
        self.flag = None
        self.rname = None
        self.pos = None
        self.corrected_pos = None
        self.cigar = None
        self.qscore = None
        
        # unitest fuctionality
        if raw_line is not None:
            self.raw_line = raw_line
            self.line = self.raw_line.split('\t')
            self.parse_columns()
            self.determine_strandedness()

    def set_read(
        self,
        raw_line=None):
        """
        sets the sam read
        """
        if raw_line is not None:
            self.line = raw_line.split('\t')

            self.parse_columns()
            self.determine_strandedness()
            return self.line

    def parse_columns(self):
        self.qname = self.line[0]
        self.umi = self.qname.split(':')[-1]
        self.flag = int(self.line[1])
        self.rname = self.line[2]
        self.pos = int(self.line[3])
        self.cigar = self.line[5]
        self.rnext = self.line[6]
        self.qscore = qual_score(self.line[10])

    def correct_start_position(
        self,
        current_position:int=None, 
        cigar_string:str=None,
        strand:str=None):
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
        self.corrected_pos = self.pos

        # unittesting functionality
        if current_position is not None:
            self.corrected_pos = current_position
        if cigar_string is not None:
            self.cigar = cigar_string
        if strand is not None:
            self.strand = strand

        # evaluate if soft clipping present
        cigar_nums, cigar_letters = self.parse_cigar_string(
            cigar=self.cigar
        )
        minus_strand_map = {
            'D': +1,
            'I': -1,
            'M': +1,
            'S': +1,
            'N': +1,
        }

        # amake position correction based on cigar string
        if self.strand == "+":
            # just have to subtract on the forward
            if cigar_letters[0] == 'S':
                correction = cigar_nums[0]
                self.corrected_pos = self.corrected_pos - correction
        elif self.strand == "-":
            # TODO: need to actually take the whole cigar string into account on the reverse
            if cigar_letters[0] == 'S':
                cigar_letters.pop(0) # ignore the first soft clip if reverse cause this is referring to the 3' end
                cigar_nums.pop(0)
            if len(cigar_letters) > 0:
                cigar_letters = [*map(minus_strand_map.get,cigar_letters)]
                cigar_nums = [map_val * distance for map_val, distance in zip(cigar_letters,cigar_nums)]
                correction = sum(list(map(int,cigar_nums))) # convert to ints before summing
                self.corrected_pos = self.corrected_pos + correction # could put minus one here, but as long as it's consistent it doesn't matter
        return self.corrected_pos

    def generate_postrand(
        self,
        corrected_pos=None,
        strand=None,
        umi=None):
        """
        generates position, read combo for current read. To be used
        in evaluating whether duplicate or not.
        """

        # unittest functionality
        if corrected_pos is not None:
            self.corrected_pos = corrected_pos
        if strand is not None:
            self.strand = strand
        if umi is not None:
            self.umi = umi

        self.postrand = (self.corrected_pos, self.strand, self.umi)
        return self.postrand

    def determine_strandedness(
        self,
        flag=None):
        """
        determines strandedness by looking at 16 of bitwise flag
        """
        # unittest functionality
        if flag is not None:
            self.flag = flag
        if ((self.flag & 16) == 16):
            self.strand = "-" # reverse
        else:
            self.strand = "+" # forward
        return self.strand

    def determine_if_mapped(
        self,
        flag:int=None):
        """
        Determines if read is mapped using bit 4 of the bitwise flag

        Parameters:
        -----------

        flag : int
            Bitwise flag of current read.

        Returns:
        --------

        mapped : bool
            True if read is mapped.
            False is read is unmapped.

        """
        # unittest functionality
        if flag is not None:
            self.flag = flag
        if ((self.flag) & 4):
            self.mapped = False
        else:
            self.mapped = True
        return self.mapped

    def parse_cigar_string(
        self,
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
        # unittest functionality
        if cigar is not None:
            self.cigar = cigar
        # parse into components
        cigar_nums = re.split("[A-Z]",self.cigar)
        cigar_nums.remove("")
        cigar_nums = [int(num) for num in cigar_nums]
        cigar_letters = re.split("[0-9]+",self.cigar)
        cigar_letters.remove("")
        return cigar_nums, cigar_letters


