# Deduper pseudocode

## Objective

Generate pseudocode to remove PCR duplicates from a SAM file. 

## Problem

When performing an RNA-seq experiment, it is important to remove duplicates for the purpose of removing PCR bias and generating results that can be accurately compared across experimental treatments. We will produce a script that retains a single copy of each read.

## Generate examples

[Path to tests](test/)

## Algorithm

1. Open input SAM file
2. Samtool sort
3. Evaluate each read and retain single copy
    * Bitwise flag - FLAG col 2
        * 4: mapped - don't care about unmapped
        * 16: strand
    * start position - POS col 4
    * read length - this will be the same for each read
    * CIGAR string
        * soft clipping - S
            * use to correct start position
        * alternative splicing - M
    * UMI
    * chromosome - RNAME col 3
    * UMI - QNAME col 1
4. Output to SAM file


## High level functions

In order to get some OOP practice, I'm going to approach this problem from an OOP mindset. Therefore, I will create classes for each object at each step of the algorithm. 

I still need to evaluate how and when I want to sort the SAM file. This will be dependent upon dictionary constraints. I would like to utilize a dictionary in tracking which reads have been written, but may be able to utilize my samtools to generate checkpoints for a dictionary reset. For example, if I sort the SAM file by chromosome position, then there's no need to include chromosome in my evaluation dictionary. Instead I can just reset my dictionary once I reach the next chromosome in the sorted SAM file. Similarly, I could do this for UMI. At this point we would just be using the bitwise flag and read position in determining if it has already been written or not. 

* Class1: SamFile
    * input:
        * input_filename
        * umi_filename
        * retention_filename
        * duplicate_filename
        * read_length
    * attributes:
        * input_filename
        * umi_filename
        * retention_filename
        * duplicate_filename
        * read_length
        * total_reads
        * num_mapped_reads
        * num_unmapped_reads
        * num_reads_maintained
            * how many reads did we retain
        * num_reads_discarded
            * how many reads have been thrown out
        * eval_dict
            * keeps track of which reads already exist in SAM so that we can discard the duplicates
        * umis_read
            * boolean to confirm umis have been read
        * deduped
            * boolean to tell us whether dedupe has been run
    * functions:
        * samtool_sort(input_filename,sorted_filename)
            * note: may sort prior to input but just putting here for now or could use bash commands here if we want to just keep it all within python
            * description: sort the sam file
            * test examples:
                * input:
                    * input_filename="test/inputs/test1.sam"
                    * sorted_filename="test/inputs/test1.sorted.sam"

                    * defined in test example output
                * output:
                    * 
            * return statement:
                * return None
        * dedupe(input_filename,retention_filename,duplicate_filename)
            * description: 
                * opens file
                * reads line by line
                    * utilizing Class2 properly parse
                * closes file
            * test examples:
                * inputs:
                    * input_filename="test/inputs/test1.sam"
                    * retention_filename="test/outputs/test1.retain.sam"
                    * duplicate_filename="test/outputs/test1.discard.sam"
                * output: 
                    * "Deduped test/inputs/test1.sam"
                    * "Umis pulled from STL96.txt"
                    * "Retained reads output to test/outputs/actual/test1.retain.sam"
                    * "Discarded reads output to test/outputs/actual/test1.discard.sam"
                    * "2 reads processed"
                    * "2 mapped reads"
                    * "0 unmapped reads"
                    * "1 read maintained"
                    * "1 read discarded"
            * return statement:
                * return print_statment
        * read_umis(umi_filename)
            * description: read in the umi file and write to eval_dict
            * test examples:
                * input: 
                    * umi_filename="STL96.txt"
                * output: 
                    * eval_dict={"AACGCCAT":"","AAGGTACG":"",...,"TTCGTTCG":""}
            * return statement:
                * return eval_dict

* Class2: SamRead
    * input: 
        * next line - already split list by '\t'
    * attributes:
        * qname - will need this for UMI and ID
        * umi - will obtain from QNAME
        * rname - chromosome
        * pos - start position
        * cigar - contains cigar string
        * flag - bitwise flag
        * strand - plus or minus
    * functions:
        * parse_columns(line)
            * description: pull umi, rname, pos, cigar, flag, and strand from strand and assign to corresponding attributes
            * test examples:
                * input: 
                    * line=["NS500451:154:HWKTMBGXX:1:11101:21621:1145:AAGGTACG", "0", "2", "93022350", "36", "71M", "*", "0", "0", "TTCCACTGTTGCTTCATAACTGCAGTCCTAACATAAATGTCTGACATGTAGGATGATCTTAAGCAACCCCT","6AEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE<EEAAAEE"]
                * output: 
                    * None - just writes to attributes as listed above
            * return statement:
                * return None
        * extract_umi(qname)
            * description: pulls umi out of the qname
            * test examples:
                * input:
                    * qname="NS500451:154:HWKTMBGXX:1:11101:21621:1145:AAGGTACG"
                * output:
                    * None - just writes umi to attribute
            * return statement:
                * return None
        * correct_position(pos,cigar)
            * description: corrects the start position using pos and cigar
            * test examples:
                * input:
                    * pos=10
                    * cigar="2S10M"
                * output:
                    * None - just corrects self.pos to 8
            * return statement:
                * return None
        * evaluate_existence(eval_dict,identity)
            * description: evaluates whether read exists using evaluation dictionary from SamFile class, make additions to many of the attributes contained in the SamFile class
            * test examples:
                * input:
                    * eval_dict = {}
                    * identity = {self.umi:{self.chr,self.pos,self.strand}}
                * output: 
                    * False
            * return statement:
                * return False
        * if evaluate_existence returns False
            * write_to_retain(line,retention_filename)
                * descriptoin: writes unique information to retention output sam
                * test examples:
                    * input:
                        * line=["NS500451:154:HWKTMBGXX:1:11101:21621:1145:AAGGTACG", "0", "2", "93022350", "36", "71M", "*", "0", "0", "TTCCACTGTTGCTTCATAACTGCAGTCCTAACATAAATGTCTGACATGTAGGATGATCTTAAGCAACCCCT","6AEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE<EEAAAEE"]
                        * retention_filename="test/outputs/test1.retain.sam"
                    * output:
                        * None
                * return statement:
                    * return None
            * write_to_dict(eval_dict,identity)
                * description: writes the unique info to the master dict 
                * test examples:
                    * input:
                        * eval_dict = {}
                        * identity = {self.umi:{self.chr,self.pos,self.strand}}
                    * output:
                        * None
                * return statement:
                    * return None
        * if evaluate_existence returns True
            * write_to_discard(line)
                * description: writes nonunique information to duplicate outputsam
                * test examples:
                    * input:
                        * line=["NS500451:154:HWKTMBGXX:1:11101:21621:1145:AAGGTACG", "0", "2", "93022350", "36", "71M", "*", "0", "0", "TTCCACTGTTGCTTCATAACTGCAGTCCTAACATAAATGTCTGACATGTAGGATGATCTTAAGCAACCCCT","6AEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE<EEAAAEE"]
                    * output:
                        * None
                * return statmenet
                    * return None
    * note: the object for this class will need to be overwritten for each read so that we don't end up with everything in memory
    




