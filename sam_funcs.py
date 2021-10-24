from Bioinfo import convert_phred, qual_score


# define clases
class Dedupe():
    def __init__(
        self, 
        input_filename=None, 
        umi_filename=None, 
        retention_filename=None,
        read_length=150,
        keep_highest_qscore=False):

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
        self):
        """
        reads and stores umis
        """
        # umi file specified
        # if self.umi_filename is not None:
        #     open_umis = open(self.umi_filename,'r')
        #     umis = open_umis.read().split('\n')
        #     # add umis to the eval dictionary
        #     for umi in umis:
        #         if umi != '': # in case there is trailing newline
        #             self.eval_dict[umi] = {}
        #     open_umis.close()
        #     self.umis_read = True
        if self.umi_filename is not None:
            open_umis = open(self.umi_filename,'r')
            self.umis = open_umis.read().split('\n')
        # randomer umis
        self.randomer_umi = True

        return self.randomer_umi

    def write_to_retain(
        self, 
        line):
        """
        writes line to retain
        """
        self.open_retain_sam.write(line)

    def evaluate_existence(
        self,
        postrand=None,
        umi=None,
        qscore=None,
        raw_line=None,
        eval_dict=None,
        randomer_umi=None,
        keep_highest_qscore=None):
        """
        evaluates if read has already been retained

        postrand, umi, and qscore will be passed from samread class

        returns True if unique
        returns False if duplicate
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
        postrand=None,
        qscore=None,
        raw_line=None,
        eval_dict=None):
        
        if eval_dict is not None:
            self.eval_dict = eval_dict

        self.eval_dict[postrand] = {'qscore':qscore,'line':raw_line}

        return self.eval_dict

    def dedupe(self):
        """
        deduplicats the input sam file
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
        eval_dict=None):
        """
        resets the values of the eval dict to empty lists
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
        self.qscore = qual_score(self.line[10])

    def correct_start_position(
        self,
        current_position=None, 
        cigar_string=None,
        strand=None):
        """
        corrects start position using cigar string
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
        temp_cigar = [c for c in self.cigar if not c.isdigit()]

        # make correction if there is a soft clip
        if temp_cigar[0] == 'S':
            correction = int(self.cigar.split('S')[0])
        
            # add if reverse, subtract if forward
            if self.strand == "+":
                self.corrected_pos = self.corrected_pos - correction
            else:
                self.corrected_pos = self.corrected_pos + correction
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
