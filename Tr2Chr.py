#!/usr/bin/python3                                                                                                                                                               
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


if __name__ == "__main__":

    import sys, getopt

    alignments_f = ''
    query_f = ''

    argv = sys.argv[1:]
    try:
        opts, args = getopt.getopt(argv,"hi:q:",["ifile=","query="])
    except getopt.GetoptError:
        print ('Tr2Chr.py -i <inputFile> -q <queryFile>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print ('please specify two command line arguments to run this scipt, i.e.\nTr2Chr.py -i <inputFile> -q <queryFile>')
            sys.exit()
        elif opt in ("-i", "--ifile"):
            alignments_f = arg
        elif opt in ("-q", "--query"):
            query_f = arg
    print ('Input file selected:', alignments_f)
    print ('Query file selected:', query_f)

    T2Galignment_dic = get_alignments(alignments_f)
    get_T2Gcoords(query_f, T2Galignment_dic)
