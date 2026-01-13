# Python for BLAST input

## 1. Edit distance main function

```python
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
    g = -3 # this is gap penalty, meet gap, -3
    
    # initialize all zero matrix
    matrix = []
    for i in range(0, m):
        tmp = []
        for j in range(0, n):
            tmp.append(0)
        matrix.append(tmp)
    
    # 初始化第一列：每个位置都是上一个位置加一个空位罚分
    for sii in range(0, m):
        matrix[sii][0] = sii*g
    # 初始化第一行：同理
    for sjj in range(0, n):
        matrix[0][sjj] = sjj*g
        
    
    for siii in range(1, m):
        for sjjj in range(1, n):
    # 状态转移方程：当前格子的值取以下三个来源的最大值：
            # 1. matrix[siii-1][sjjj] + g : 来自上方（序列 2 产生空位）
            # 2. matrix[siii-1][sjjj-1] + SingleBaseCompare(...) : 来自左上角（碱基匹配或错配）
            # 3. matrix[siii][sjjj-1] + g : 来自左方（序列 1 产生空位
            
            matrix[siii][sjjj] = max(matrix[siii-1][sjjj] + g, matrix[siii - 1][sjjj - 1] + SingleBaseCompare(seq1,seq2,siii, sjjj), matrix[siii][sjjj-1] + g)
    sequ1 = [seq1[m-1]]
    sequ2 = [seq2[n-1]]
    
    
    # Reverse step
		sequ1 = [seq1[m-1]] # 存放比对后的序列 1（从末尾字符开始）
    sequ2 = [seq2[n-1]] # 存放比对后的序列 2
    
    while m > 1 and n > 1:
        # 寻找当前格子得分的来源（上、左、左上哪个分最高）
        # 1. 如果来自左上角 (matrix[m-2][n-2])：说明两个碱基是对齐的
        if max(matrix[m-1][n-2], matrix[m-2][n-2], matrix[m-2][n-1]) == matrix[m-2][n-2]:
            m -= 1
            n -= 1
            sequ1.append(seq1[m-1])
            sequ2.append(seq2[n-1])
        # 2. 如果来自左边 (matrix[m-1][n-2])：序列 1 需要插入空位 '-'
        elif max(matrix[m-1][n-2], matrix[m-2][n-2], matrix[m-2][n-1]) == matrix[m-1][n-2]:
            n -= 1
            sequ1.append('-')
            sequ2.append(seq2[n-1])
        # 3. 如果来自上边 (matrix[m-2][n-1])：序列 2 需要插入空位 '-'
        else:
            m -= 1
            sequ1.append(seq1[m-1])
            sequ2.append('-')
            
    sequ1.reverse() # 因为是从后往前找的，所以需要反转列表
    sequ2.reverse()
    align_seq1 = ''.join(sequ1) # 转成字符串
    align_seq2 = ''.join(sequ2)
    
    # Calculate score
    align_score = 0.
    # 遍历比对后的序列，统计碱基完全一致的位置
    for k in range(0, len(align_seq1)):
        if align_seq1[k] == align_seq2[k]:
            align_score += 1
            
    # 计算一致性比例 (Identity)
    align_score = float(align_score)/len(align_seq1)
    
    # 返回比对后的两个字符串以及相似度得分
    return align_seq1, align_seq2, align_score
```

## 2. Display

```python
# Display BlAST result
def Display(seque1, seque2):
    le = 40
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
```

## 3. Reference and seeding

Similar to reference

```python
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
```

从索引库获取位置

```python
# Get word's postion in genome from library
def GetWordPos(word):
    assert len(word)== 11
    seek_index = WordToIndex(word,11-1)
    positions = []
    # Only need to open file once and load index once (for single sequence)
    with open('dataset/sarscov2.txt','r') as chr_seq:
        seeks = np.load("dataset/sarscov2_library_seeks.npy")
        # read location, seeks[seek_index, 0] is offset
        chr_seq.seek(seeks[seek_index,0])
        # read, seeks[seek_index, 1] is length
        position = chr_seq.read(seeks[seek_index,1])
        try:
						# position[:-1] 是去掉最后的逗号，split(",") 将其切分为列表
            # map(int, ...) 将字符串列表转为整数列表
            positions.append(list(map(int, position[:-1].split(","))))
        except:
            positions.append([])
    return positions
```

当你在索引库里找到了一个位置（比如 5000），这个函数负责去原始 `.fasta` 文件中把那一块的序列“抠”出来。chrom_seek_index是fasta文件的，而sarscov2_library_seeks是k-mer的。

```python
 def ExtractSeq(chr_index, pos, length):
    # 自动计算由于换行符导致的偏移补偿
    # pos / fasta_line_width 计算出当前位置前有多少行，即有多少个换行符
    pos = pos + floor(pos / fasta_line_width)
    
    # 计算最终在文件中的跳转位置：
    # chrom_seek_index[chr_index, 1] 是该染色体标题行结束的位置
    # build 里面初始化:chrom_seek_index.append([seq_length, header_length])
    # pos-1 是从标题结束处开始计算的相对偏移
    fasta_file.seek(chrom_seek_index[chr_index, 1] + pos - 1)
    
    # 读取指定长度的内容，并用正则去掉中间可能混入的换行符
    return re.sub(r'\n', '', fasta_file.read(length))
```

## 4. Main function

GetWordPos返回的是 [ [位置1, 位置2...] ]，是一个列表嵌套，因为只有一个reference所以只有一个列表[位置1, 位置2...]。

**代码里的解决方法：**`pos += words_length - word_index - 1`

这个公式的作用是**把所有种子的坐标都推算到同一个“终点线”上**。

- 对于种子 0：`1000 + (3 - 0 - 1) = 1002`
- 对于种子 1：`1001 + (3 - 1 - 1) = 1002`
- 对于种子 2：`1002 + (3 - 2 - 1) = 1002`

**奇迹发生了：** 经过计算，这三个种子的坐标都变成了 **1002**！

```python
def Blast(query_seq):
    i = 0
    query_words = [] # k-mer
    query_seq_length = len(query_seq)
    words_length = query_seq_length-11+1
    while i < words_length:
        query_words.append(query_seq[i:i+11])
        i += 1
    words_positions = []
    # seeding
    for word in query_words:
        words_positions.append(GetWordPos(word)) #每个k-mer在基因组里的位置
    
    # voting
    # Dynamically iterate through all sequences
    for chr_index in range(len(chr_names)):
        for word_index in range(words_length): # 一共有几个k-mer
            for pos in range(len(words_positions[word_index][chr_index])):
                words_positions[word_index][chr_index][pos] += words_length - word_index - 1
        
        words_positions_corrects = []
        for word_index in range(words_length):
		        # 把校准后的所有位置，不分种子，全部倒进一个大桶里。现在我们不再关心是哪个种子匹配上的，我们只关心校准后的坐标。
            # word_index = 0的时候 words_positions[word_index][chr_index]就代表第一个seed匹配的位置[100, 500]
            words_positions_corrects += words_positions[word_index][chr_index]
        
        words_positions_corrects_count = Counter(words_positions_corrects)
        
        # Debug info: display matching statistics
        # if words_positions_corrects_count:
        #     top_matches = words_positions_corrects_count.most_common(5)
        #     print(f"Top 5 matching positions found in sequence {chr_names[chr_index]}:")
        #     for pos, count in top_matches:
        #         print(f"  Position {pos}: {count} 11-mer matches")
        
        # Filtering
        finded_postions = []
        for count_ in words_positions_corrects_count:
            # we can select the bigger threshold of words_positions_corrects_count[count_] just 
            # like we select the highly similar sequence in NCBI BLAST
            if words_positions_corrects_count[count_] > 5:
                finded_postions.append(count_)
        
        # Alignment
        if not finded_postions and words_positions_corrects_count:
            print(f"⚠️ Warning: Found matches but all below threshold (>5), max match count: {max(words_positions_corrects_count.values())}")
        if finded_postions:
            for finded_postion in finded_postions:
		            # 1. 计算要从基因组提取的范围（前后多预留一点空间，防止截断）
                candidate_seq_pos = finded_postion - query_seq_length + 11 - 5
                candidate_seq_length = query_seq_length + 11
                # 提取候选序列
                candidate_sequence = ExtractSeq(chr_index,candidate_seq_pos,candidate_seq_length)
                
                # 2. 寻找最佳开头：尝试微调起始位置（0-14位），看哪种对齐得分最高
                i_start_indexs = []
                for i_start in range(15):
                    _,_,score = SMalignment(candidate_sequence[i_start:],query_seq)
                    i_start_indexs.append(score)
                i_start = np.array(i_start_indexs).argmax()
                
                # 3. 寻找最佳结尾：同理，尝试微调结束位置
                i_end_indexs = []
                for i_end in range(1,16):
                    _,_,score = SMalignment(candidate_sequence[:-i_end],query_seq)
                    i_end_indexs.append(score)
                i_end = np.array(i_end_indexs).argmax()+1
                
                # 4. 最终裁剪并执行最后一次比对
                candidate_sequence = candidate_sequence[i_start:-i_end]
                align_seq1,align_seq2,align_score = SMalignment(candidate_sequence,query_seq)
                # print(align_seq1, align_seq2, align_score)
                if align_score>0.8:
                    print("find in chromosome "+chr_names[chr_index]+": "+str(candidate_seq_pos+i_start)+' ---> '+str(candidate_seq_pos+i_start+len(candidate_sequence)-1)+", align score: "+str(align_score))
                    Display(align_seq1, align_seq2)
    return None
```