#!/usr/bin/env python

# library import
import argparse


# define arguments
parser = argparse.ArgumentParser(description="Reference-based PCR duplicate removal from mapped SAM file")
parser.add_argument("-f", "--file", help="Input SAM file - sorted, mapped reads")
parser.add_argument("-p", "--paired-end", action="store_true", help="Designates file is paired end. If flag is not set, the program will assume single end reads.")
parser.add_argument("-u", "--umi", help="Designates UMI file. If set, UMI file should follow flag. If flag is not set, the program will assume randomer UMI")
args = parser.parse_args()

# define globals
paired_end = False
input_file = args.file

# argparser
if args.paired_end:
    paired_end = True
    print(input_file)
if args.umi:
    umis = args.umi
    print(umis)

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

    def dedupe(self):
        """
        deduplicats the input sam file
        """
        open_input_sam = open(self.input_filename,'r')
        open_retain_sam = open(self.retention_filename,'w')
        open_discard_sam = open(self.duplicate_filename,'w')

        while True:
            line = open_input_sam.readline()
            
            if line == '': # break if end of file
                break

            if line[0] == '@': # continue header lines
                continue
            
            self.total_reads += 1 # made it through all evalutation steps

            read = SamRead(line)

            # reset dictionary if we've reached a new chromosome
            if int(read.rname) != int(self.current_chr):
                self.current_chr = int(read.rname)
                self.reset_eval_dict()

            # if erroneous umi, write to discard 
            # TODO: challenge is to make correction 
            if read.umi not in self.eval_dict:
                open_discard_sam.write(line)
                self.num_reads_discarded += 1
            # if new read, add to dictionary and write to retain file
            elif read.postrand not in self.eval_dict[read.umi]:
                self.eval_dict[read.umi].append(read.postrand)
                open_retain_sam.write(line)
                self.num_reads_retained += 1
            else:
                open_discard_sam.write(line)
                self.num_reads_discarded += 1

        open_input_sam.close()
        open_retain_sam.close()
        open_discard_sam.close()

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
        self.correct_start_position()
        self.generate_postrand()

    def parse_columns(self):
        self.qname = self.line[0]
        self.umi = self.qname.split(':')[-1]
        self.flag = int(self.line[1])
        self.rname = self.line[2]
        self.pos = int(self.line[3])
        self.cigar = self.line[5]

    def correct_start_position(self):
        """
        corrects start position using cigar string
        """
        
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

    def generate_postrand(self):
        """
        generates position, read combo for current read. To be used
        in evaluating whether duplicate or not.
        """

        self.postrand = (self.pos, self.strand)

    def determine_strandedness(self):
        """
        determines strandedness by looking at 16 of bitwise flag
        """
        
        if ((self.flag & 16) == 16):
            self.strand = 1 # reverse
        else:
            self.strand = 0 # forward
