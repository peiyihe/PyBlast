__author__ = 'Jiashun'
import re
import numpy as np
from collections import Counter
from math import ceil
from math import floor

# Global variables
chr_names = None
chrom_seek_index = None
fasta_file = None

def init_blast(data_dir='dataset'):
    """Initialize global variables required for BLAST"""
    global chr_names, chrom_seek_index, fasta_file
    chr_names = np.load(f'{data_dir}/sarscov2_chr_names.npy')
    chrom_seek_index = np.load(f'{data_dir}/sarscov2_chrom_seek_index.npy')
    fasta_file = open(f"{data_dir}/sarscov2.fasta")
    print(f"BLAST initialization complete! Loaded {len(chr_names)} sequence(s).")
    return chr_names, chrom_seek_index, fasta_file

# compare single base
def SingleBaseCompare(seq1,seq2,i,j):
    if seq1[i] == seq2[j]:
        return 2
    else:
        return -1
    
# Smith–Waterman Alignment 
def SMalignment(seq1, seq2):
    m = len(seq1)
    n = len(seq2)
    g = -3
    matrix = []
    for i in range(0, m):
        tmp = []
        for j in range(0, n):
            tmp.append(0)
        matrix.append(tmp)
    for sii in range(0, m):
        matrix[sii][0] = sii*g
    for sjj in range(0, n):
        matrix[0][sjj] = sjj*g
    for siii in range(1, m):
        for sjjj in range(1, n):
            matrix[siii][sjjj] = max(matrix[siii-1][sjjj] + g, matrix[siii - 1][sjjj - 1] + SingleBaseCompare(seq1,seq2,siii, sjjj), matrix[siii][sjjj-1] + g)
    sequ1 = [seq1[m-1]]
    sequ2 = [seq2[n-1]]
    while m > 1 and n > 1:
        if max(matrix[m-1][n-2], matrix[m-2][n-2], matrix[m-2][n-1]) == matrix[m-2][n-2]:
            m -= 1
            n -= 1
            sequ1.append(seq1[m-1])
            sequ2.append(seq2[n-1])
        elif max(matrix[m-1][n-2], matrix[m-2][n-2], matrix[m-2][n-1]) == matrix[m-1][n-2]:
            n -= 1
            sequ1.append('-')
            sequ2.append(seq2[n-1])
        else:
            m -= 1
            sequ1.append(seq1[m-1])
            sequ2.append('-')
    sequ1.reverse()
    sequ2.reverse()
    align_seq1 = ''.join(sequ1)
    align_seq2 = ''.join(sequ2)
    align_score = 0.
    for k in range(0, len(align_seq1)):
        if align_seq1[k] == align_seq2[k]:
            align_score += 1
    align_score = float(align_score)/len(align_seq1)
    return align_seq1, align_seq2, align_score

# Display BlAST result
def Display(seque1, seque2):
    le = 60
    while len(seque1)-le >= 0:
        print('sequence1: ',end='')
        for a in list(seque1)[le-40:le]:
            print(a,end='')
        print("\n")
        print('           ',end='')
        for k in range(le-40, le):
            if seque1[k] == seque2[k]:
                print('|',end='')
            else:
                print(' ',end='')
        print("\n")
        print('sequence2: ',end='')
        for b in list(seque2)[le-40:le]:
            print(b,end='')
        print("\n")
        le += 40
    if len(seque1) > le-40:
        print('sequence1: ',end='')
        for a in list(seque1)[le-40:len(seque1)]:
            print(a,end='')
        print("\n")
        print('           ',end='')
        for k in range(le-40, len(seque1)):
            if seque1[k] == seque2[k]:
                print('|',end='')
            else:
                print(' ',end='')
        print("\n")
        print('sequence2: ',end='')
        for b in list(seque2)[le-40:len(seque2)]:
            print(b,end='')
        print("\n")

# transform base to numeric value
def WordToNum(word):
    tmp = []
    trans = {'A':1,'C':2,'G':3,'T':4}
    for w in word:
        tmp.append(trans[w])
    return tmp

# transform word with 11 bases to its index
def WordToIndex(word,word_len):
    tmp = 0
    word_num = WordToNum(word)
    for i,v in enumerate(word_num):
        tmp += (v-1)*4**(word_len-i)
    return tmp   

# Get word's postion in genome from library
def GetWordPos(word):
    assert len(word)== 11
    seek_index = WordToIndex(word,11-1)
    positions = []
    # Only need to open file once and load index once (for single sequence)
    with open('dataset/sarscov2.txt','r') as chr_seq:
        seeks = np.load("dataset/sarscov2_library_seeks.npy")
        chr_seq.seek(seeks[seek_index,0])
        position = chr_seq.read(seeks[seek_index,1])
        try:
            positions.append(list(map(int, position[:-1].split(","))))
        except:
            positions.append([])
    return positions

# Extract subsequence from fasta file
def ExtractSeq(chr_index,pos,length):
    # Fix: fasta file has 70 characters per line (not 60)
    pos = pos+floor(pos/70)
    fasta_file.seek(chrom_seek_index[chr_index,1]+pos-1)
    return re.sub(r'\n', '', fasta_file.read(length))

# main blast function
def Blast(query_seq):
    i = 0
    query_words = []
    query_seq_length = len(query_seq)
    words_length = query_seq_length-11+1
    while i < words_length:
        query_words.append(query_seq[i:i+11])
        i += 1
    words_positions = []
    for word in query_words:
        words_positions.append(GetWordPos(word))
    # Dynamically iterate through all sequences
    for chr_index in range(len(chr_names)):
        for word_index in range(words_length):
            for pos in range(len(words_positions[word_index][chr_index])):
                words_positions[word_index][chr_index][pos] += words_length - word_index - 1
        
        words_positions_corrects = []
        for word_index in range(words_length):
            words_positions_corrects += words_positions[word_index][chr_index]
        
        words_positions_corrects_count = Counter(words_positions_corrects)
        
        # Debug info: display matching statistics
        if words_positions_corrects_count:
            top_matches = words_positions_corrects_count.most_common(5)
            print(f"Top 5 matching positions found in sequence {chr_names[chr_index]}:")
            for pos, count in top_matches:
                print(f"  Position {pos}: {count} 11-mer matches")
        
        finded_postions = []
        for count_ in words_positions_corrects_count:
            # we can select the bigger threshold of words_positions_corrects_count[count_] just 
            # like we select the highly similar sequence in NCBI BLAST
            if words_positions_corrects_count[count_] > 5:
                finded_postions.append(count_)
        
        if not finded_postions and words_positions_corrects_count:
            print(f"⚠️ Warning: Found matches but all below threshold (>5), max match count: {max(words_positions_corrects_count.values())}")
        if finded_postions:
            for finded_postion in finded_postions:
                candidate_seq_pos = finded_postion - query_seq_length + 11 - 5
                candidate_seq_length = query_seq_length + 11
                candidate_sequence = ExtractSeq(chr_index,candidate_seq_pos,candidate_seq_length)
                i_start_indexs = []
                for i_start in range(15):
                    _,_,score = SMalignment(candidate_sequence[i_start:],query_seq)
                    i_start_indexs.append(score)
                i_start = np.array(i_start_indexs).argmax()
                i_end_indexs = []
                for i_end in range(1,16):
                    _,_,score = SMalignment(candidate_sequence[:-i_end],query_seq)
                    i_end_indexs.append(score)
                i_end = np.array(i_end_indexs).argmax()+1
                candidate_sequence = candidate_sequence[i_start:-i_end]
                align_seq1,align_seq2,align_score = SMalignment(candidate_sequence,query_seq)
                if align_score>0.8:
                    print("find in chromosome "+chr_names[chr_index]+": "+str(candidate_seq_pos+i_start)+' ---> '+str(candidate_seq_pos+i_start+len(candidate_sequence)-1)+", align score: "+str(align_score))
                    Display(align_seq1, align_seq2)
    return None

if __name__ == "__main__":
    # Initialize BLAST
    init_blast()
    
    # Execute query
    query_sequence = 'TAACCAGAATGGAGAACGCAGTGGGGCGCGATCAAAACAACGTCGGCCCCAAGGTTTACCCAATAATACT'
    Blast(query_sequence)