# define clases
class SamFile():
    def __init__(
        self, 
        input_filename:str='samp.sort.sam', 
        umi_filename:str='STL96.txt', 
        retention_filename:str='helloworld.retain.sam',
        duplicate_filename:str='helloworld.discard.sam',
        read_length:int=150
        ):
        self.input_filename = input_filename
        self.umi_filename = umi_filename
        self.retention_filename = retention_filename
        self.duplicate_filename = duplicate_filename
        self.read_length = int(read_length)
        self.total_reads = 0
        self.num_reads_retained = 0
        self.num_reads_discarded = 0
        self.eval_dict = {}
        self.umis_read = False
        self.deduped = False
        self.current_chr = 1

        self.read_umis()

    def read_umis(self):
        """
        reads and stores umis
        """
        open_umis = open(self.umi_filename,'r')
        umis = open_umis.read().split('\n')
        
        # add umis to the eval dictionary
        for umi in umis:
            if umi != '': # in case there is trailing newline
                self.eval_dict[umi] = []
        
        open_umis.close()

        self.umis_read = True

    def write_to_retain(self, line):
        """
        writes line to retain
        """
        self.open_retain_sam.write(line)

    def write_to_discard(self, line):
        """
        writes line to discard
        """
        self.open_discard_sam.write(line)

    # TODO: still need to complete this
    def write_to_eval_dict(self):
        return None

    def dedupe(self):
        """
        deduplicats the input sam file
        """
        self.open_input_sam = open(self.input_filename,'r')
        self.open_retain_sam = open(self.retention_filename,'w')
        self.open_discard_sam = open(self.duplicate_filename,'w')

        while True:
            line = self.open_input_sam.readline()
            
            if line == '': # break if end of file
                break

            if line[0] == '@': # store header lines
                self.write_to_retain(line)
                self.write_to_discard(line)
                continue
            
            self.total_reads += 1 # made it through all evalutation steps

            read = SamRead(line)
            read.correct_start_position()
            read.generate_postrand()

            # reset dictionary if we've reached a new chromosome
            if int(read.rname) != int(self.current_chr):
                self.current_chr = int(read.rname)
                self.reset_eval_dict()

            # if erroneous umi, write to discard 
            # TODO: challenge is to make correction 
            if read.umi not in self.eval_dict:
                self.write_to_discard(line)
                self.num_reads_discarded += 1
            # if new read, add to dictionary and write to retain file
            elif read.postrand not in self.eval_dict[read.umi]:
                self.eval_dict[read.umi].append(read.postrand)
                self.write_to_retain(line)
                self.num_reads_retained += 1
            else:
                self.write_to_discard(line)
                self.num_reads_discarded += 1

        self.open_input_sam.close()
        self.open_retain_sam.close()
        self.open_discard_sam.close()

        self.deduped = True

    def reset_eval_dict(self):
        """
        resets the values of the eval dict to empty lists
        """
        self.eval_dict = {k:[] for k in self.eval_dict}

class SamRead:
    def __init__(self,line:list):
        self.line = line.split('\t')
        self.qname = None
        self.umi = None
        self.rname = None
        self.pos = None
        self.cigar = None
        self.flag = None
        self.strand = None
        self.postrand = None

        self.parse_columns()
        self.determine_strandedness()

    def parse_columns(self):
        self.qname = self.line[0]
        self.umi = self.qname.split(':')[-1]
        self.flag = int(self.line[1])
        self.rname = self.line[2]
        self.pos = int(self.line[3])
        self.cigar = self.line[5]

    def correct_start_position(
        self,
        current_position=None, 
        cigar_string=None,
        strand=None):
        """
        corrects start position using cigar string
        """

        # unittesting functionality
        if current_position is not None:
            self.pos = current_position
        
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
            if self.strand == 1:
                self.pos = self.pos + correction
            else:
                self.pos = self.pos - correction

    def generate_postrand(
        self,
        position=None,
        strand=None):
        """
        generates position, read combo for current read. To be used
        in evaluating whether duplicate or not.
        """

        # unittest functionality
        if position is not None:
            self.pos = position

        if strand is not None:
            self.strand = strand

        self.postrand = (self.pos, self.strand)

        return(self.postrand)

    def determine_strandedness(self):
        """
        determines strandedness by looking at 16 of bitwise flag
        """
        
        if ((self.flag & 16) == 16):
            self.strand = "-" # reverse
        else:
            self.strand = "+" # forward
