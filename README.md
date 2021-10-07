# Deduper

## Part 1
Use this repo template to create your own Deduper repo - you should do all your work in your own repository. Please name it `Deduper-<github-user-name>`.

Write up a strategy for writing a Reference Based PCR Duplicate Removal tool. That is, given a sam file of uniquely mapped reads, remove all PCR duplicates (retain only a single copy of each read). Develop a strategy that avoids loading everything into memory. You should not write any code for this portion of the assignment. Be sure to:
- Define the problem
- Write examples:
    - Include a properly formated input sam file
    - Include a properly formated expected output sam file
- Develop your algorithm using pseudocode
- Determine high level functions
    - Description
    - Function headers
    - Test examples (for individual functions)
    - Return statement
    
For this portion of the assignment, you should design your algorithm for single-end data, with 96 UMIs. UMI information will be in the QNAME, like so: ```NS500451:154:HWKTMBGXX:1:11101:15364:1139:GAACAGGT```. Discard any UMIs with errors (or error correct, if you're feeling ambitious).

## Part 2
An important part of writing code is reviewing code - both your own and other's. In this portion of the assignment, you will be assigned 3 students' algorithms to review. Be sure to evaluate the following points:
- Does the proposed algorithm make sense to you? Can you follow the logic?
- Does the algorithm do everything it's supposed to do? (see part 1)
- Are proposed functions reasonable? Are they "standalone" pieces of code?

You can find your assigned reviewees on Canvas. You can find your fellow students' repositories at 
```
github.com/<user>/Deduper-<github-user-name>
```
Be sure to leave comments on their repositories by creating issues or by commenting on the pull request.

## Part 3
Write your deduper function!

Given a SAM file of uniquely mapped reads, remove all PCR duplicates (retain only a single copy of each read). Remember:
- Samtools sort
- Account for all possible CIGAR strings (including adjusting for soft clipping, etc.)
- Strand
- Single-end reads
- Known UMIs
- Considerations:
    - Millions of reads – avoid loading everything into memory!
    - Be sure to utilize functions appropriately
    - Appropriately comment code and include doc strings
- **CHALLENGE**: Include options for
    - Single-end vs paired-end
    - Known UMIs vs randomers
    - Choice of duplicate written to file
    
You MUST:
- Write Python 3.9 compatible code
- Include the following argparse options
    - ```-f```, ```--file```: required arg, absolute file path
    - ```-p```, ```--paired```: optional arg, designates file is paired end (not single-end)
    - ```-u```, ```--umi```: optional arg, designates file containing the list of UMIs (unset if randomers instead of UMIs)
    - ```-h```, ```--help```: optional arg, prints a USEFUL help message (see argparse docs)
        - If your script is not capable of dealing with a particular option (ex: no paired-end functionality), your script should print an error message and quit
- Output the first read encountered if duplicates are found
    - You may include an additional argument to designate output of a different read (highest quality or random or ???)
- Output a properly formatted SAM file with “_deduped” appended to the filename
- Name your python script ```<your_last_name>_deduper.py```

