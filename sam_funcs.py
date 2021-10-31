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
        keep_highest_qscore:bool=False,
        paired_end:bool=None):

        self.input_filename = input_filename
        self.umi_filename = umi_filename
        self.retention_filename = retention_filename
        self.read_length = int(read_length)
        
        self.total_reads = 0
        self.num_reads_retained = 0
        self.num_reads_discarded = 0
        self.current_chr = '1'
        self.eval_dict = {}

        self.keep_highest_qscore = keep_highest_qscore

        self.umis_read = False
        self.randomer_umi = False
        self.umis = None
        self.incorrect_umi = 0

        self.deduped = False

        self.paired_end = False
        self.paired_end_dict = {}

        self.summary_dict = {}
        
        if not self.randomer_umi:
            self.read_umis(
                umi_filename=self.umi_filename
            )

    def write_to_paired_end_dict(
        self,
        paired_end_dict:dict=None,
        rnext:dict=None,
        key:tuple=None,
        partner_key:tuple=None):
        """
        Write to the paired end dictionary.

        Parameters:
        -----------

        paired_end_dict : dict
            The paired-end dictionary.

            Structure:
            ----------

            {
                rnext: {
                    'key': eval_dict key (self.correct_pos, self.strand, self.umi, instance),
                    'partner_key': eval_dict key belonging to mate (self.correct_pos, self.strand, self.umi, instance)
                }

        rnext : str
            The name of the mate read.

        key : tuple
            The key of the first encountered read.

        partner_key : tuple
            The key belonging to the second encountered read.

        Returns:
        --------

        paired_end_dict : dict
            The paired-end dictionary. See structure above.

        """
        # unittest functionality
        if paired_end_dict is not None:
            self.paired_end_dict = paired_end_dict
        # write to paired end
        paired_end_dict[rnext] = {
            'key': key,
            'partner_key': partner_key
        }
        return self.paired_end_dict

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
            open_umis.close()
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

    # TODO: write this function
    def evaluate_pair_existence(
        self,
        paired_end_dict:dict=None,
        qname:str=None,
        rname:str=None,
        postrand:str=None,
        qscore:float=None,
        raw_line:str=None,
        eval_dict:dict=None):
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

        rname : str
            The qname of the PAIR.

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
            self.paired_end_dict = paired_end_dict
        if eval_dict is not None:
            self.eval_dict = eval_dict
        # run evaluation
        if qname in self.paired_end_dict: # if qname in dict, pair has already been found
            # add second encountered read key to instance of first read
            self.paired_end_dict[qname]['partner_key'] = postrand
            # grab partner key
            partner_key = self.paired_end_dict[rname]['key']
        else:
            partner_key = ''

        write_to_paired_end_dict(
            paired_end_dict=self.paired_end_dict,
            rnext=rname,
            key=postrand,
            partner_key=partner_key
        )

        return None

    def evaluate_existence(
        self,
        postrand:tuple=None,
        umi:str=None,
        qscore:float=None,
        raw_line:str=None,
        qname:str=None,
        rname:str=None,
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
                    (self.correct_pos, self.strand, self.umi, instance): {
                        'qscore': float,
                        'line': raw_line,
                        'qname': qname,
                        'rname': rname

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
                self.incorrect_umi += 1 # TODO: add unittest for this
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
                    (self.correct_pos, self.strand, self.umi, instance): {
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
        Deduplicates the input SAM file.
        """
        self.open_input_sam = open(self.input_filename,'r')
        self.open_retain_sam = open(self.retention_filename,'w')
        # self.open_discard_sam = open(self.duplicate_filename,'w')
        self.header_length = 0

        while True:
            line = self.open_input_sam.readline()
            if line == '': # break if end of file
                self.num_reads_retained += len(self.eval_dict)
                # write last chromosome to summary dict
                self.summary_dict[self.current_chr] = len(self.eval_dict)
                 # dump last chromosome to sam
                self.dump_dict_to_sam(
                    eval_dict=self.eval_dict
                )
                break
            if line[0] == '@': # store header lines
                self.write_to_retain(
                    line=line
                )
                # self.write_to_discard(line)
                self.header_length += 1 
                continue
            self.total_reads += 1 # made it through all evalutation steps

            read = SamRead(
                raw_line=line
            )
            read.correct_start_position()
            read.generate_postrand()

            # reset dictionary if we've reached a new chromosome
            if read.rname != self.current_chr:
                self.summary_dict[self.current_chr] = len(self.eval_dict)
                self.current_chr = read.rname
                self.dump_dict_to_sam(
                    eval_dict=self.eval_dict
                )
                self.num_reads_retained += len(self.eval_dict)
                self.reset_eval_dict(
                    eval_dict=self.eval_dict
                )

            # check if pair has already been found
            # TODO: i'm thinking it would be best if a pair has already been found, to make a change to the first and second reads eval_dict value
            # TODO: this is great to perform some kind of operation here and then pass that info into evaluate_existence
            if self.paired_end:
                print('paired-end')

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
        
        # finish up summary dict 
        # TODO: turn this into a function
        self.summary_dict['header_length'] = self.header_length
        self.summary_dict['total_reads'] = self.total_reads
        self.summary_dict['incorrect_umi'] = self.incorrect_umi
        self.summary_dict['reads_retained'] = self.num_reads_retained
        for entry in self.summary_dict:
            print(entry,self.summary_dict[entry],sep=": ")

        self.deduped = True

    # TODO: write this function
    def evaluate_qscore(
        self):
        """
        Evaluates if qscore is higher than current placeholder and performs actions accordingly.
        """
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

    # TODO: write this function
    def write_out_summary_dict(
        self):
        """
        Finishes writing to and prints out the summary dict.
        """
        # TODO: this needs to include header length, incorrect umis, reads per chromosome, reads discarded, total reads 
        return None

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
        self.read_num = None
        
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
            'N': +1, # TODO: check to make sure this is right
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
        corrected_pos:int=None,
        strand:str=None,
        umi:str=None,
        instance:int=None):
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

        self.postrand = (self.corrected_pos, self.strand, self.umi, 1)
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

    def determine_read_num(
        self,
        flag:int=None):
        """
        Determines if first or second read in paired-end sam file. Bit 64 is R1 and bit 128 is R2.

        Parameters:
        -----------

        flag : int
            Bitwise flag.

        Returns:
        --------

        read_num : str
            'R1' or 'R2' depending on if it's the first or second read in the sequence.

        """
        # unittest functionality
        if flag is not None:
            self.flag = flag
        if ((self.flag & 64) == 64):
            self.read_num = 'R1'
        elif ((self.flag & 128) == 128):
            self.read_num = 'R2'
        return self.read_num
