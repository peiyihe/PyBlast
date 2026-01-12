import numpy as np
import re
from math import ceil

def BaseToNum(chr_seq):
    chr_seq = re.sub(r'A', '1', chr_seq)
    chr_seq = re.sub(r'C', '2', chr_seq)
    chr_seq = re.sub(r'G', '3', chr_seq)
    chr_seq = re.sub(r'T', '4', chr_seq)
    return chr_seq

def BaseToIndex(word,word_len):
    tmp = 0
    for i,v in enumerate(word):
        tmp += (int(v)-1)*4**(word_len-i)
    return tmp    

def GenSeek(library,word_len):
    seeks = np.zeros((4**word_len,2),dtype=int)
    tmp = 0
    for i,l in enumerate(library):
        seeks[i,0] = tmp
        seeks[i,1] = len(l)
        tmp += len(l)
    return seeks

def BuildLibrary(chr_name):
    word_len = 11
    chr_seq = chrom_dict[chr_name]
    chr_seq = BaseToNum(chr_seq)
    chr_len = len(chr_seq)
    library = np.zeros(4**word_len,dtype=str).tolist()
    ii = 0
    while ii<chr_len-word_len:
        w = chr_seq[ii:ii+word_len]
        ii += 1
        if 'N' in w:
            continue
        try:
            library[BaseToIndex(w,word_len-1)] += str(ii)+","
        except:
            pass
    
    seeks = GenSeek(library,word_len)
    lib_seq = ''.join(library)
    with open('dataset/sarscov2.txt', 'w') as f:
        f.write(lib_seq)
        f.close()
    np.save('dataset/sarscov2_library_seeks.npy',seeks)
    

if __name__ == '__main__':
    hg19 = open("dataset/sarscov2.fasta")
    head = True
    chrom_dict = {}
    head_line = []
    chr_names = []
    for line in hg19:
        if line.startswith(">"):
            head_line.append(line)
            if head:
                head = False
            else:
                chr_seq = re.sub(r'\n', '', chr_seq)
                chr_seq = chr_seq.upper()
                chrom_dict[chr_name] = chr_seq
            chr_name = line.split()[0][1:]
            chr_names.append(chr_name)
            chr_seq = ''
            print(chr_name,end=",")
        else:
            chr_seq += line
    chr_seq = re.sub(r'\n', '', chr_seq)
    chr_seq = chr_seq.upper()
    chrom_dict[chr_name] = chr_seq
    
    # Build index based on actual sequence length (not parsed from header)
    chrom_seek_index = []
    for i, (chr_name_item, line) in enumerate(zip(chr_names, head_line)):
        seq_length = len(chrom_dict[chr_name_item])
        header_length = len(line)
        chrom_seek_index.append([seq_length, header_length])
    chrom_seek_index = np.array(chrom_seek_index)
    
    # Dynamically process index for multiple sequences
    for i in range(1, len(head_line)):
        chrom_seek_index[i,1]=chrom_seek_index[i,1]+chrom_seek_index[i-1,1]+chrom_seek_index[i-1,0]+ceil(chrom_seek_index[i-1,0]/70)
    np.save('dataset/sarscov2_chrom_seek_index.npy',chrom_seek_index)
    np.save('dataset/sarscov2_chr_names.npy',np.array(chr_names))
    print(chr_names)
    
    # For single sequence, process directly in loop (no multiprocessing needed)
    print("\nStarting to build index library...")
    for chr_name in chr_names:
        print(f"Processing: {chr_name}")
        BuildLibrary(chr_name)
    print("Index library construction completed!")
