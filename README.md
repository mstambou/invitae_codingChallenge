# Github Repository for Invitae Coding Challenge

## Problem Specification

Given a Transcript to Chromosome alignment:

```
                0    5    10   15   20     25   30   35   40   45   50    
Chromosome:CHR1 ACTGTCATGTACGTTTAGCTAGCC--TAGCTAGGGACCTAGATAATTTAGCTAG
TR1                GTCATGTA-------CTAGCCGGTA-----------AGATAAT
                   |    |           |    |               |   |
                   0    5           10   15              20  24
```

We can compactly express this alignment in the same way that we compactly represent a read alignment in the ​ SAM/BAM format​ : using a position and CIGAR string. In this case, the (0-based) position is CHR1:3, and the CIGAR string is ​ 8M7D6M2I2M11D7M​ . For this exercise, you may assume that the transcript is always mapped from genomic 5’ to 3’


The objective is then to translate a (0-based) transcript coordinate to a (0 based) genome coordinate. For example the fifth base in TR1 (i.e. TR1:4) maps to genome coordinate CHR1:7. Similarly, TR1:13 maps to CHR1:23 and TR1:14 maps to an insertion immediately before CHR1:24. 

The software should be implemented in a language of your choice (with some preference towards python), and should conform to community-preferred style guidelines (e.g., python should be “pythonic” and conform to PEP-8). Code will be evaluated both on correctness and overall quality. Your solution should take the following inputs:
1. A four column (tab-separated) file containing the transcripts. The first column is the transcript name, and the remaining three columns indicate it’s genomic mapping: chromosome name, 0-based starting position on the chromosome, and CIGAR string indicating the mapping.
2. A two column (tab-separated) file indicating a set of queries. The first column is a transcript name, and the second column is a 0-based transcript coordinate.Your solution should handle errors appropriately. 

The output is a four column tab separated file with one row for each of the input queries. The first two columns are exactly the two columns from the second input file, and the remaining two columns are the chromosome name and chromosome coordinate, respectively. Example input/output is provided below.

## Solution

To solve this problem I have implemented a script in Python 3, that takes two inputs, the first one specifying the Transcript to Chromosome alignments and the second one are a list of queries in Transcript coordiantes, where the program will then transform them from these transcript coordinates to chromosome choordinates based on the alignments provided by the first input file.

It should be noted that the script does not assume that a transcript is only mapped to one genome, it also takes into account multi-mapped transcripts, i.e. if transcript TR1 is mapped to both chromosomes CHR1 and CHR2, then it will list choordinates on the two different chromosomes respectively. 

The script also handles the following exceptions and warns the user about them:
 * If an invalid CIGAR alignment string is provided by the user, the script warns the user about this and skips it.
 * If an invalid query is asked by the user over the provided alignments, i.e. if the coordinate for the transcript is out of bounds, then the script will warn the user and will skip this query.
 
 By the end of the exectution the script will list the different transcript to chromosome alignments, as well construct their alignments based on the CIGAR string and chromosome start information and output a tab separated file listing these transformations for each of the queries, from transcript cooridinates to chromosome coordinates. The output will be saved in the same directory as the input files, with the suffix: 
