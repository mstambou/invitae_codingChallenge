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
 
 By the end of the exectution the script will list the different transcript to chromosome alignments, as well construct their alignments based on the CIGAR string and chromosome start information and output a tab separated file listing these transformations for each of the queries, from transcript cooridinates to chromosome coordinates. The output will be saved in the same directory as the input files, with the suffix: _T2Gcoords.tsv
 
 Example to run the script:
 
 ```python3
 python3 Tr2Chr.py -i test_files/T2GAlignment_1.txt -q test_files/query_1.txt
 
 python3 Tr2Chr.py -i test_files/T2GAlignment_2.txt -q test_files/query_2.txt
```


```python3
 python3 Tr2Chr.py -i test_files/T2GAlignment_1.txt -q test_files/query_1.txt
    0) Transcript : TR1	Chromosome : CHR1	AlnStart : 3	CIGAR : 8M7D6M2I2M11D7M
    Chromosome:	GGGGGGGGGGGGGGGGGGGGGGGG--GGGGGGGGGGGGGGGGGGGG...
    Alignment :	...||||||||.......||||||..||...........|||||||
    Transcript:	   TTTTTTTT-------TTTTTTTTTT-----------TTTTTTT
    
    1) Transcript : TR2	Chromosome : CHR2	AlnStart : 10	CIGAR : 20M
    Chromosome:	GGGGGGGGGGGGGGGGGGGGGGGGGGGGGG...
    Alignment :	..........||||||||||||||||||||
    Transcript:	          TTTTTTTTTTTTTTTTTTTT
    
    python3 Tr2Chr.py -i test_files/T2GAlignment_2.txt -q test_files/query_2.txt
    0) Transcript : TR1	Chromosome : CHR1	AlnStart : 0	CIGAR : 3I2M2D3I3M
    Chromosome:	---GGGG---GGG...
    Alignment :	...||.....|||
    Transcript:	TTTTT--TTTTTT
    
    1) Transcript : TR1	Chromosome : CHR2	AlnStart : 1	CIGAR : 3M5I1M1D2M
    Chromosome:	GGGG-----GGGG...
    Alignment :	.|||.....|.||
    Transcript:	 TTTTTTTTT-TT
```



```python
global consumes_both, consumes_query, consumes_reference
consumes_both = ['M', '=', 'X']
consumes_query = ['I', 'S']
consumes_reference = ['D', 'N', 'H', 'P']


def get_T2Gidx(CHR_start, CIGAR):
    """
    Function that will take 0-indexed chromosome start position
    and a CIGAR string, and will parse the CIGAR string and will
    create a dictionary mapping between the 0-indexed Transcript 
    index and 0-indexed chromosome coordiantes, and will also return
    an alignment scheme, where G denotes chromosome nucleotides, T deenotes
    transcript nucleotides, - denoting gaps, | denoting match states
    and . denoting indels.
    """
    tr2g_dic = {}
    n, o = '', ''
    genome, transcript, alignment = 'G'*CHR_start, ' '*CHR_start, '.'*CHR_start    
    TR_PTR = 0
    CHR_PTR = CHR_start
    if CIGAR.replace('=', '').isalnum(): #make sure first it only contains CIGAR defined operations
        for ch in CIGAR:
            if ch.isnumeric():
                n+= ch
            if ch.isalpha() or ch == '=':
                o = ch
                n = int(n)
                if o in consumes_both:
                    for i in range(n):
                        tr2g_dic[TR_PTR] = CHR_PTR
                        TR_PTR += 1
                        CHR_PTR += 1
                        genome += 'G'
                        transcript += 'T'
                        alignment += '|'
                elif o in consumes_query:
                    for i in range(n):
                        tr2g_dic[TR_PTR] = CHR_PTR
                        TR_PTR += 1
                        genome += '-'
                        transcript += 'T'
                        alignment += '.'
                elif o in consumes_reference:
                    for i in range(n):
                        tr2g_dic[TR_PTR] = CHR_PTR                    
                        CHR_PTR += 1
                        genome += 'G'
                        transcript += '-'
                        alignment += '.'
                n, o = '', ''
        genome += '...'
        return tr2g_dic, (genome, alignment, transcript)
    else:
        return 'err', '.'


def get_alignments(alignment_f):
    """
    Function that will take the first input file as argument, 
    will parse it and return a mapping scheme between the transcript
    coordiantes to chromosome coordinates of the different transcript 
    to chromosome alignment schemes listed in the first input file
    """
    with open(alignments_f, 'r') as in_f:
        T2Galignment_dic = {}
        for c, line in enumerate(in_f):
            line = line.strip().split('\t')
            Tr, Chr, start, CIGAR = line[0], line[1], int(line[2]), line[3]
            if Tr not in T2Galignment_dic:
                T2Galignment_dic[Tr] = {}
            T2Galignment_dic[Tr][Chr], aln = get_T2Gidx(start, CIGAR)
            if T2Galignment_dic[Tr][Chr] == 'err':
                print(f'Warning {CIGAR} is not a proper CIGAR string!!!')
            print(f'{c}) Transcript : {line[0]}\tChromosome : {line[1]}\tAlnStart : {line[2]}\tCIGAR : {line[3]}')            
            print(f'Chromosome:\t{aln[0]}\nAlignment :\t{aln[1]}\nTranscript:\t{aln[2]}\n')
        return T2Galignment_dic
        

def get_T2Gcoords(query_f, T2Galignment_dic):
    """
    Function that will take in the second input file (list of queries)
    and output a tab separated file denoting the translated coordinates
    for the queries from transcript coordinates to chromosoem coordinates
    """
    out_fname = query_f.rsplit('.')[0]+'_T2Gcoords.tsv'    
    with open(query_f, 'r') as in_f, open(out_fname, 'w') as out_f:
        for line in in_f:
            line = line.strip().split('\t')
            Tr, p = line[0], int(line[1])
            if Tr not in T2Galignment_dic:
                print(f'Warning, skipping transcript {Tr}, not defined in the alignment file')
            else:
                for Chr in T2Galignment_dic[Tr]:
                    if T2Galignment_dic[Tr][Chr] == 'err':
                        print(f'Warning, skipping {Tr}:{p} to {Chr} mapping, transcript {Tr} to Chromosome {Chr} alignment is not well defined!!!')
                    elif p not in T2Galignment_dic[Tr][Chr]:
                        print(f'Warning, skipping {Tr}:{p} to {Chr} mapping, coordinate {p} is not defined in transcript {Tr}')
                    else:
                        out_f.write(f'{Tr}\t{p}\t{Chr}\t{T2Galignment_dic[Tr][Chr][p]}\n')


alignments_f = '/data/mstambou/invitae/T2GAlignment_1.txt'
query_f = '/data/mstambou/invitae/query_1.txt'

print('python3 Tr2Chr.py -i test_files/T2GAlignment_1.txt -q test_files/query_1.txt')
T2Galignment_dic = get_alignments(alignments_f)
get_T2Gcoords(query_f, T2Galignment_dic)


alignments_f = '/data/mstambou/invitae/T2GAlignment_2.txt'
query_f = '/data/mstambou/invitae/query_2.txt'

print('python3 Tr2Chr.py -i test_files/T2GAlignment_2.txt -q test_files/query_2.txt')
T2Galignment_dic = get_alignments(alignments_f)
get_T2Gcoords(query_f, T2Galignment_dic)


```

   
    



```python

```

