# Python for BLAST reference

# Reference

1. 将 DNA 序列中的碱基字符（A, C, G, T）替换为数字字符串（1, 2, 3, 4）。

```python
def BaseToNum(chr_seq):
    chr_seq = re.sub(r'A', '1', chr_seq)
    chr_seq = re.sub(r'C', '2', chr_seq)
    chr_seq = re.sub(r'G', '3', chr_seq)
    chr_seq = re.sub(r'T', '4', chr_seq)
    return chr_seq
```

1. 计算 $k$-mer 的哈希索引：$\sum (digit - 1) \times 4^{position}$

```python
def BaseToIndex(word,word_len):
    tmp = 0
    for i,v in enumerate(word):
        tmp += (int(v)-1)*4**(word_len-i)
    return tmp
```

1. 生成寻址信息
****假设 $k$-mer 只有四种：AA(0), AC(1), AG(2), AT(3)。
扫描完基因组后，library 可能是这样的：
1. `library[0]` (AA): `"1,5,20,"` （长度 7）
2. `library[1]` (AC): `"2,"` （长度 2）
3. `library[2]` (AG): `""` （长度 0）
4. `library[3]` (AT): `"3,15,"` （长度 5）
**`GenSeek` 执行过程：**
• `i=0`: `seeks[0,0]=0`, `seeks[0,1]=7`, `tmp` 变成 7。
• `i=1`: `seeks[1,0]=7`, `seeks[1,1]=2`, `tmp` 变成 9。
• `i=2`: `seeks[2,0]=9`, `seeks[2,1]=0`, `tmp` 保持 9。
• `i=3`: `seeks[3,0]=9`, `seeks[3,1]=5`, `tmp` 变成 14。

`GenSeek` 的作用就是记录“书的目录”：
• **`seeks[i, 0]` (Offset)**：第 $i$ 个 $k$-mer 的位置信息从文件的第几个字符开始。
• **`seeks[i, 1]` (Length)**：往后读多少个字符能读完这个 $k$-mer 的所有位置。

```python
def GenSeek(library,word_len):
    seeks = np.zeros((4**word_len,2),dtype=int)
    tmp = 0
    for i,l in enumerate(library):
        seeks[i,0] = tmp         # 存储在文件中的起始偏移量（Offset）
        seeks[i,1] = len(l)       # 存储该段信息的长度
        tmp += len(l)
    return seeks
```

1. 构建索引库

```python
def BuildLibrary(chr_name):
    word_len = 11  # 定义 k-mer 长度为 11
    chr_seq = chrom_dict[chr_name] # 从全局词典中取出对应染色体
    chr_seq = BaseToNum(chr_seq) # 转为数字
    chr_len = len(chr_seq)
    # 创建一个长度为 4^11 的空列表，用来存每个 k-mer 出现的位置
    library = np.zeros(4**word_len,dtype=str).tolist()
    
    ii = 0
    while ii < chr_len - word_len:
        w = chr_seq[ii : ii + word_len] # 从当前位置 ii 开始，切出长度为 11 的片段
        ii += 1
        if 'N' in w: # 跳过含有未知碱基 'N' 的片段
            continue
        try:
            # 将当前位置 ii 记录到对应的索引位置，用逗号隔开
            # 在library对应的位置加上这个ii的索引
            library[BaseToIndex(w, word_len-1)] += str(ii) + ","
        except:
            pass
    
    seeks = GenSeek(library, word_len) # 生成指针表
    lib_seq = ''.join(library)         # 将所有位置信息合并成一个巨大的字符串
    
    # 存入文件
    with open('dataset/sarscov2.txt', 'w') as f:
        f.write(lib_seq)
    np.save('dataset/sarscov2_library_seeks.npy', seeks)
```

1. main function

```python
if __name__ == '__main__':
    hg19 = open("dataset/sarscov2.fasta") # 打开原始基因组文件
    head = True                          # 标志位：判断是否是文件的第一个标题行
    chrom_dict = {}                      # 字典：存放 {序列名: 完整序列}
    head_line = []                       # 列表：存放所有的标题行文本
    chr_names = []                       # 列表：存放所有的序列名
    
    for line in hg19:                    # 逐行读取文件
        if line.startswith(">"):         # 如果这一行以 '>' 开头（说明是新序列的开始）
            head_line.append(line)       # 记录原始标题行
            if head:
                head = False             # 如果是第一个标题，不保存之前的序列（因为还没读到序列）
            else:
		            # 这个地方并没有用到 可能是多个reference的时候才需要
                # 当遇到下一个 '>' 时，说明上一个序列读完了，进行清理和保存
                chr_seq = re.sub(r'\n', '', chr_seq) # 去掉序列里的换行符
                chr_seq = chr_seq.upper()            # 全部转为大写
                chrom_dict[chr_name] = chr_seq       # 存入字典
            
            # 解析当前标题行，提取名字（去掉 '>' 和空格后的第一部分）
            chr_name = line.split()[0][1:] 
            chr_names.append(chr_name)
            chr_seq = ''                 # 重置序列字符串，准备接收新序列
            print(chr_name, end=",")      # 打印进度
        else:
		        # 非标题行都是做这个，主要执行的这一行
            chr_seq += line              # 如果不是标题行，就把这一行拼接到当前序列中
            
# 循环结束后，最后一个序列后面没有新的 '>' 了，所以需要手动去掉换行符号
    chr_seq = re.sub(r'\n', '', chr_seq)
    chr_seq = chr_seq.upper()
    chrom_dict[chr_name] = chr_seq
    
# 初始化一个列表，记录每个序列的长度和标题行长度
    chrom_seek_index = []
    for i, (chr_name_item, line) in enumerate(zip(chr_names, head_line)):
        seq_length = len(chrom_dict[chr_name_item]) # 碱基数量
        header_length = len(line)                   # 标题行字符数
        chrom_seek_index.append([seq_length, header_length])
    
    chrom_seek_index = np.array(chrom_seek_index)

for i in range(1, len(head_line)):
        # This part is also for multipe reference sequence not checked yet
        # 公式 = 上一个序列的起始位置 
        #        + 上一个序列标题长度 
        #        + 上一个序列碱基总数 
        #        + 换行符的数量 (ceil(长度/70))
        chrom_seek_index[i,1] = chrom_seek_index[i,1] + chrom_seek_index[i-1,1] + \
                               chrom_seek_index[i-1,0] + ceil(chrom_seek_index[i-1,0]/70)

# 将计算好的偏移量索引和名字保存为 numpy 文件
    np.save('dataset/sarscov2_chrom_seek_index.npy', chrom_seek_index)
    np.save('dataset/sarscov2_chr_names.npy', np.array(chr_names))
    print(chr_names)
    
    print("\nStarting to build index library...")
    # 遍历所有读入的序列，正式开始构建位置索引库
    for chr_name in chr_names:
        print(f"Processing: {chr_name}")
        BuildLibrary(chr_name) # 调用核心函数，生成 .txt 和 .npy 指针表
    print("Index library construction completed!")
```