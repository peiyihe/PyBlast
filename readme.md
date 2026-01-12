# PyBLAST

conda create -n python_blast python=3.10 -y

## Install

conda activate python_blast

conda install numpy numba -y

python -c "import numpy; import numba; print('Installation successful!')"

conda install jupyter -y

## Usage

python build_library.py

python blast.py

Or use main.ipynb

## Warning

- **ExtractSeq** function take care 70 need to match with input fasta length each row

- Similar in build_library.py


```
chrom_seek_index[i,1]=chrom_seek_index[i,1]+chrom_seek_index[i-1,1]+chrom_seek_index[i-1,0]+ceil(chrom_seek_index[i-1,0]/70)
```

Still have some bugs, will return to this project soon.

Modified from: https://github.com/JiaShun-Xiao/BLAST-bioinfor-tool