# PyBLAST

Python implementation for BLAST sequence alignment. Mainly for learning.

## Install

```
conda create -n python_blast python=3.10 -y

conda activate python_blast

conda install numpy numba -y

python -c "import numpy; import numba; print('Installation successful!')"

conda install jupyter -y
```

## Usage

```
python build_library.py

python blast.py
```

## Warning

Still have some bugs mainly for multiple reference, will return to this project soon.

Modified from: https://github.com/JiaShun-Xiao/BLAST-bioinfor-tool