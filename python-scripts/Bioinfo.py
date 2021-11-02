#!/usr/bin/env python

# import modules
from math import sqrt
from typing import Any
import os
import filecmp


# define constants
DNAbases = set('ATGCNatcgn')
RNAbases = set('AUGCNaucgn')

def convert_phred(letter: str) -> int:
    """Converts a single character into a phred score"""
    
    score = ord(letter) - 33
    
    return score

# convert_phred unit tests
if __name__ == "__main__":
    assert convert_phred('I') == 40, "validate ascii conversion"
    print("convert_phred functioning correctly")


def qual_score(phred_score: str, RNAflag: bool=False) -> float:
    """calculates average phred score of a read"""

    # size to be used for average
    size = len(phred_score)
    
    # establish empty total
    total = 0
    
    # score each letter and add to total
    for letter in phred_score:
        score = convert_phred(letter)
        total += score
    
    # calc average
    avg_score = total / size
    
    return avg_score

# qual_score tests
if __name__ == "__main__":
    assert qual_score('I') == 40, "validate ascii to score is correct"
    assert qual_score('GGII') == 39, "validate qual_score summing and dividing correctly"
    print("passed all qual score tests")



def validate_base_seq(seq: str,RNAflag: bool=False) -> bool:
    '''Returns True if string is composed of only As, Ts (or Us if RNAflag),
    Gs, Cs. False otherwise. Case insensitive.'''

    # determine if additional characters exist in the input string
    check = set(seq) <= (RNAbases if RNAflag else DNAbases)
    
    return check

# validate_base_seq unit tests
if __name__ == "__main__":
    assert validate_base_seq("AATAGAT") == True, "Validate base seq does not contain non DNA values"
    assert validate_base_seq("AAUAGAU", True) == True, "validate base seq contians only RNA values"
    assert validate_base_seq("Hi there!") == False, "validate base seq fails the DNA/RNA test"
    assert validate_base_seq("Hi there!", True) == False, "validate base seq fails on non-DNA/RNA test"
    print("Passed DNA and RNA tests")
    


def init_list(lst: list, value: Any, length: int=1) -> list:
    '''
    This function takes a value and length as input and will return a list of
    length "length" with each element of the list containing the value "value".
    
    "value" can be of Any type, including another list
    '''
    
    lst = [value for _ in range(length)]
    
    return lst


def populate_list(file: str) -> (list,float):
    """Given a fastq file, calculate total phred score at each position and
    return a tuple that contains a list of sums at each base and a count of the number 
    of total lines in the file"""
    
    # initiate list
    lst = init_list([])
    
    # initiate score sum
    linecount = 0
    
    # open file
    with open(file) as f:
        
        # loop through each line
        while True:
            
            # read next line
            newline = f.readline()
            
            # break if empty (end of file)
            if newline == '':
                break
                
            # grab quality score line of record
            if linecount % 4 == 3:
                
                # loop through each letter in quality score
                for i in range(len(newline) - 1):

                    # convert phred quality score
                    score = Bioinfo.convert_phred(newline[i])
                    
                    # assign to position w/i list
                    lst[i] += score

            linecount+=1

    return lst, linecount


# define some functions
def median(l: list) -> float:
    """calculate the median value of a list"""
    
    # determine length of list
    lengthoflist = len(l)
    
    # sort list
    l = sorted(l)
    
    # calculate median
    if lengthoflist % 2 == 1:
        medindex = lengthoflist // 2
        
        med = float(l[medindex])

    else:
        medindexhigh = int(lengthoflist / 2)
        medindexlow = medindexhigh - 1 # because 0-indexed
        
        med = float((l[medindexlow] + l[medindexhigh]) / 2)
        
    return med


def mean_calc(l: list) -> float:
    """calculate the mean of a list"""
    
    # determine length of list
    lengthoflist = len(l)
    
    # sum all values of list
    total = sum(l)
    
    # determine average
    avg = float(total / lengthoflist)
    
    return avg
    

def var_calc(l: list) -> float:
    """calculate the variance of a list"""
    
    # determine length of list
    lengthoflist = len(l)
    
    # determine average of list
    avgscore = mean_calc(l)
    
    # generate diff from mean at each pos
    square_mean_diff = [(pos - avgscore)**2 for pos in l]
    
    # average of square mean diff
    var = sum(square_mean_diff) / lengthoflist
    
    return var


def stdev_calc(l: list) -> float:
    """calculate the standard deviation of a list"""
    
    # determine length of list
    var = var_calc(l)
    
    # calc stdev
    stdev = sqrt(var)
    
    return stdev


# needs to be generalized
def list_of_lists(file: str, all_qscores: list) -> (list,int):
    """Generate a list of lists of quality scores, given a file and empty list of lists.
    
    The index of the outer list will represent the 0-based position of each read. The inner
    list will hold the qscore for each observation at the position specified by the index
    of the outer list"""
    
#     # create individual instance for each list in all_qscores
#     all_qscores = [x[:] for x in all_qscores]
    
    # initiate score sum
    linecount = 0
    
    # avoid 'copied' lists
    all_qscores = [x[:] for x in all_qscores]
    
    # open file
    with open(file) as f:
        
        
        # loop through each line
        while True:
            
            # read next line
            newline = f.readline().rstrip()
            
            # break if empty (end of file)
            if newline == '':
                break
            
            # update statement
            if linecount % 100000 == 0:
                print(linecount)
            
            # grab quality score line of record
            if linecount % 4 == 3:
                
                # loop through each letter in the newline and append phred score to correct list in all_qscores
                for i, letter in enumerate(newline):

                    phred_score = Bioinfo.convert_phred(letter)
                                       
                    all_qscores[i].append(phred_score)

            # increment linecounter
            linecount+=1
    
    return all_qscores, linecount

# define our fasta line sequence concatenator
def oneline_fasta(input_file: str, output_file: str):
    """
    removes newline characters from fasta file and writes to new fasta file

    Parameters:
    -----------
    input_file: name of input fasta file

    output_file: name of output fasta file

    Outputs:
    --------
    writes new fasta file with one line sequences

    """

    # overwrite output file if exists
    with open(output_file,'w'):
        pass

    # read each line of the input file
    with open(input_file) as f:

        # write first line
        firstline = f.readline().rstrip()
        with open(output_file,'a') as nf:
            nf.write(firstline+'\n')
        
        # write to output file
        with open(output_file,'a') as nf:

            while True:

                # read and strip next line
                line = f.readline().rstrip()

                # break at end of file
                if line == '':
                    break

                # identify and write headers         
                if line.startswith('>'):
                    nf.write('\n'+line+'\n')
                
                # concat nucleotides
                else:
                    nf.write(line)


# oneline_fasta tests
if __name__ == "__main__":
    oneline_fasta(input_file='unittests/oneline_fasta_test.fasta',output_file='unittests/oneline_fasta_test.output.fasta')
    assert os.path.exists('unittests/oneline_fasta_test.output.fasta'), "validate oneline_fasta writing files"
    assert filecmp.cmp('unittests/oneline_fasta_test.output.fasta','unittests/oneline_fasta_test.gold.fasta') == True, "validate oneline_fasta concatenating sequences correctly"
    os.remove('unittests/oneline_fasta_test.output.fasta') # removes the test output if fasta concatenated correctly
    print("oneline_fasta functioning as expected")


# define gc content calculator function
def gc_content(sequence: str) -> float:
    """
    Parameters:
    -----------
    sequence: sequence for which to calculate gc content

    Returns:
    --------
    gc_per: float that represents percentage gc content

    """
    # validate nucleic acid sequence
    validate_base_seq(sequence)

    # calculate 
    c: int = sequence.count('C') + sequence.count('c')
    g: int = sequence.count('G') + sequence.count('g')

    # total length
    seq_len: int = len(sequence)

    # gc calc
    gc_per: float = (c + g) / seq_len * 100

    return gc_per

# gc_content unit tests
if __name__ == "__main__":
    assert gc_content('cgcg') == 100.0, "validaate summing and dividing of gc_content"
    assert gc_content('CGcg') == 100.0, "validaate summing and dividing of gc_content"
    assert gc_content('ATATCGCG') == 50.0, "validaate summing and dividing of gc_content"
    assert gc_content('ATATcgCG') == 50.0, "validaate summing and dividing of gc_content"
    print("gc_content functioning correctly")


# define protein length function
def protein_length(input_protein_file: str, input_biomart_file: str, output_fasta_file: str) -> int:
    """
    Parameters:
    -----------
        input_protein_file: a fasta file of CONCATENATED (removed newlines from sequences) peptide sequences - ftp from ensembl
            -> example file
                Homo_sapiens.GRCh38.pep.all.concat.fa
        
        input_biomart_file: biomart file of gene_id, protein_id, and gene_name

    Outputs:
    --------
        -> prints length of dictionary
        -> generates an output fasta file with new header and longest associated protein for the given gene_id
        -> example header:
            >protein_id gene_id gene_name

        -> example record:
            >ENSP00000488240 ENSG00000282253 TRBD1
            GTGG



    ATTENTION: newline characters WITHIN protein sequences must already have been removed from this input file
    """

    # initiate dictionary for tracking gene lengths
    gene_dict: dict = {}

    # open file and extract contents
    with open(input_protein_file) as f:

        while True:
            # grab header line
            header: str = f.readline().rstrip()

            # grab amino acid sequence
            aa_seq: str = f.readline().rstrip()

            # break at end of file
            if header == '':
                print(aa_seq)
                break

            # split the header
            split_header: list = header.split()

            # find protein and gene
            protein_id: str = split_header[0][1:].split('.')[0] # removing > and .          
            gene_id: str = split_header[3].split(':')[1].split('.')[0] # removing label and .

            # calculate length of protein
            prot_length: int = len(aa_seq)

            # add protein if not in dictionary
            if not gene_id in gene_dict:
                
                gene_dict[gene_id] = {
                    'prot_length': prot_length,
                    'amino_acid_seq': aa_seq,
                    'protein_id': protein_id,
                    'gene_name': ''
                 }

            else:
                
                # extract current longest length
                max_length: int = gene_dict[gene_id]['prot_length']

                # replace current max values if the length is greater
                if prot_length > max_length:

                    # replace max length, associated aa_seq, and associated protein_id
                    gene_dict[gene_id]['prot_length'] = prot_length
                    gene_dict[gene_id]['amino_acid_seq'] = aa_seq
                    gene_dict[gene_id]['protein_id'] = protein_id



    # print length of gene dictionary to console
    print('length of gene_dict:',len(gene_dict))

    # grab biomart lines to include in dictionary
    with open(input_biomart_file) as bf:
        
        # toss header
        bf.readline()

        # loop through each line
        while True:

            # initialize variables
            gene_name: str = ''

            # read in next line
            line: str = bf.readline()
            
            # break statement
            if line == '':
                break
            
            # split into component pieces
            line = line.rstrip().split('\t')
            
            # extract info
            len_line = len(line)
            
            if len_line != 3:
                continue
            
            gene_id = line[0]
            gene_name = line[1]
            protein_id = line[2]

            # add gene name and protein id to dictionary
            if gene_id in gene_dict:

                # confirm protein id matches
                if gene_dict[gene_id]['protein_id'] == protein_id:

                    # if gene name has not yet been found, add gene name to dict
                    if gene_dict[gene_id]['gene_name'] == '':
                        gene_dict[gene_id]['gene_name'] = gene_name
                    # TODO: i found that some don't have names. is this correct?
            
    # write to new longest protein per gene fasta file
    with open(output_fasta_file,'w') as off:
        for k, sub_d in gene_dict.items():

            # create and write header line
            # TODO: may need to add tabs instead of spaces
            newheader = '>' + sub_d['protein_id'] + ' ' + k + ' ' + sub_d['gene_name'] + '\n'
            off.write(newheader)

            # write amino acid sequence
            aa_seq = sub_d['amino_acid_seq']
            off.write(aa_seq)
            off.write('\n') # newline


# define contig length adjust function
def adjust_contig_length(kmer_contig_length: int, kmer_length: int) -> int:
    """
    calculates physical contig length given a kmer contig length and the kmer length
    """
    adj_contig_length = int(kmer_contig_length) - 1 + kmer_length

    return adj_contig_length


# calculate nucleotide wise coverage of contig
def nucleotide_wise_coverage(contig_length: int, kmer_length: int, kmer_coverage: float) -> int:
    """
    Parameters:
    -----------
        contig_length: physical length of the contig
        kmer_length: length of the kmer
        kmer_coverage: reported kmer coverage 
            - like from contigs.fa velvet

    Returns:
    --------
        nucleotide-wise coverage


    Ck = C * (L - k + 1) / L



    """

    nt_coverage = kmer_coverage * contig_length / (contig_length - kmer_length + 1)
    
    return nt_coverage


# define kmer_normalizer
def kmer_normalizer(input_file: str, kmer_coverage_limit: int, kmer_length: int, output_file: str):
    """
    normalize kmer coverage of input fastq file and write to new fastq file

    Parameters:
    -----------
    input_file: input fastq file to be normalized

    kmer_coverage_limit: desired coverage limit to normalize to (normalization is inclusive)

    kmer_length: length of kmer

    output_file: output normalized fastq file
    """
    
    # initialize empty dictionary
    kmer_dict: dict = {}

    # open output file 
    with open(input_file) as f:
        
        # initialize record counter
        recordcount: int = 0

        # each loop handles a SINGLE RECORD
        while True:
            
            # store next record
            record: list = [f.readline().rstrip() for _ in range(4)]
            seq: str = record[1]

            # break if end of file
            if seq == '':
                break

            # kmerize current read
            for i in range(len(seq)):
                kmer = seq[i:i+kmer_length]
                if kmer in kmer_dict:
                    kmer_dict[kmer] += 1
                else:
                    kmer_dict[kmer] = 1
            
            # initialize empty list
            kmer_coverage: list = []

            # append kmer coverage to list
            for i in range(len(seq)):
                kmer = seq[i:i+kmer_length]
                v = kmer_dict[kmer]
                kmer_coverage.append(v)

            # calc median score for current read
            read_median: float = median(kmer_coverage)

            # determine if below threshold
            if read_median <= kmer_coverage_limit:
                
                # if so, write to file
                with open(output_file,'a') as o: # TODO: need to remove this from loop
                    for line in record:
                        o.write(str(line))
                        o.write("\n")
        
            # update statement
            if recordcount % 100000 == 0:
                print('records read:',recordcount)

            # increment record count
            recordcount += 1


# use matplotlib to plot a emperical kmer spectra plot
def kmer_spectra_plot(kmer_freq_dict: dict,kmer_size: int,xmin: int,xmax: int,yscale: str,graph_title: str,output_img_file: str):
    """
    creates emperical kmer spectra plot and saves to output_img_file
    
    Parameters:
    -----------
        kmer_freq_dict:
            - dictionary where keys are frequency of occurence and values are number of kmers,
        
        kmer_size:
            - kmer size used to create the kmer_freq_dict

        xmin:
            - minimum x value on plot
        
        xmax:
            - maximum x value on plot

        yscale:
            - yscale of plot (i.e. "log")

        output_file:
            - file to write emperical kmer spectra plot to

    """

    # generate and save plot
    plt.figure(figsize=(15,4))
    plt.bar(kmer_freq_dict.keys(),kmer_freq_dict.values())
    plt.title("emperical kmer spectra\nk={}, read length={}".format(kmer_size,read_length))
    plt.xlabel("frequency of occurrence")
    plt.ylabel("number of kmers (log scale)")
    plt.yscale(yscale)
    plt.xlim([xmin,xmax])
    plt.title(graph_title)
    plt.savefig(output_img_file)

