# PyBLAST

An educational Python implementation of seed-based DNA sequence alignment. PyBLAST combines exact 11-mer lookup with Smith–Waterman local alignment to demonstrate reference indexing, candidate selection, and alignment traceback.

The repository includes a runnable example and a Jupyter notebook. The bundled SARS-CoV-2 reference sequence (`NC_045512.2`) serves as example input for demonstrating the algorithm. The project focuses on programming and algorithm exploration; it does not implement the full NCBI BLAST algorithm or its statistical scoring.

## Features

- Exact 11-mer indexing with independent position lists for each FASTA record.
- Local alignment with substitutions, insertions, and deletions.
- Reference extraction across variable line widths and LF/CRLF line endings.
- Case-insensitive DNA queries with support for IUPAC ambiguity symbols.
- Console output showing reference coordinates, aligned sequences, and sequence identity.

## Installation

Create and activate the Conda environment:

```bash
conda create -n python_blast python=3.10 -y
conda activate python_blast
conda install numpy numba -y
```

Verify the installation:

```bash
python -c "import numpy; import numba; print('Installation successful!')"
```

NumPy is required by the current implementation. Numba is retained in the environment setup, but the alignment code does not currently use JIT acceleration.

For the interactive notebook, also install Jupyter and register the environment as a kernel:

```bash
conda install jupyter -y
python -m ipykernel install --user --name python_blast --display-name "Python (python_blast)"
```

## Quick Start

Run all commands from the repository root with `python_blast` activated.

### 1. Build the reference indexes

```bash
python build_library.py
```

This reads `dataset/sarscov2.fasta` and regenerates the position library and NumPy indexes in `dataset/`. Existing generated files are overwritten. Rebuild after changing the reference FASTA or updating the indexing implementation.

### 2. Run the example query

```bash
python blast.py
```

The script initializes the database and aligns its embedded 70-base query. The bundled reference produces the following result summary, followed by the aligned sequences:

```text
BLAST initialization complete! Loaded 1 sequence(s).
find in chromosome NC_045512.2: 28351 ---> 28420, align score: 1.0
```

Coordinates are **one-based and inclusive**. Here, `align score: 1.0` means 100% sequence identity across the alignment columns.

### 3. Run a custom query

Use the Python interface from a script or notebook:

```python
from blast import Blast, init_blast

init_blast()
query = (
    'TAACCAGAATGGAGAACGCAGTGGGGCGCGATCAAAACAA'
    'CGTCGGCCCCAAGGTTTACCCAATAATACT'
)
Blast(query)
```

Call `init_blast()` before searching. Subsequent calls to `Blast()` reuse the loaded reference and index. Results are printed to the console; `Blast()` returns `None`.

To explore the supplied notebook:

```bash
jupyter notebook main.ipynb
```

Select the **Python (python_blast)** kernel and run the cells in order.

## Custom Reference Data

Create a directory containing a FASTA file named `sarscov2.fasta`. The filename is currently fixed, even when the file contains sequences from another organism. Each record must have a unique identifier—the first whitespace-separated token after `>`—and a nonempty sequence.

For example:

```text
my_dataset/
└── sarscov2.fasta
```

Build and load the database using the same directory:

```python
from build_library import build_libraries
from blast import Blast, init_blast

build_libraries('my_dataset')
init_blast('my_dataset')
Blast('ACGTTGCAAGT')  # Replace with a query appropriate for your reference.
```

Multiple records can share one FASTA file. Their coordinates are relative to their individual reference sequences. Keep the FASTA and all generated files together; rebuild and reinitialize after editing the reference.

| Generated file | Contents |
| --- | --- |
| `sarscov2.txt` | Comma-separated, one-based seed positions |
| `sarscov2_chr_names.npy` | Reference identifiers in FASTA order |
| `sarscov2_chrom_seek_index.npy` | Reference lengths and FASTA byte offsets |
| `sarscov2_library_seeks.npy` | Byte offsets and lengths for seed-position lookups |

## Alignment Method

1. **Normalize the query.** Remove whitespace, convert to uppercase, and validate DNA IUPAC symbols. Queries shorter than 11 bases raise `ValueError`.
2. **Look up seeds.** Search overlapping 11-mers containing only `A`, `C`, `G`, and `T`. Seeds containing ambiguous bases are skipped in both queries and references.
3. **Select candidate regions.** Group seed hits by their implied reference start. A candidate requires at least `min(6, valid_query_seed_count)` supporting hits. Extend its reference window by up to five bases on either side.
4. **Align locally.** Apply Smith–Waterman alignment within each candidate window, using the scoring scheme below.
5. **Filter and report.** Report alignments with identity strictly greater than 80% and query coverage of at least 80%.

| Alignment operation | Score |
| --- | ---: |
| Matching canonical bases | +2 |
| Mismatch, including ambiguous bases | −1 |
| Gap, per base | −3 |

The printed `align score` is an **identity fraction**, not the dynamic-programming score:

```text
identity = matching canonical-base columns / total alignment columns
query coverage = aligned query bases / total query bases
```

Gap columns contribute to the identity denominator. Ambiguous symbols count as mismatches, including identical symbols such as `N` aligned with `N`. Local alignment may omit query ends, so identity and query coverage describe different properties of a hit.

For direct pairwise alignment without a reference database, `SMalignment(reference, query)` returns `(aligned_reference, aligned_query, identity)`.

## Repository Structure

```text
PyBlast/
├── blast.py             # Database loading, seed lookup, alignment, and output
├── build_library.py     # FASTA parsing and reference index generation
├── dataset/             # Example FASTA input; indexes are generated locally
├── main.ipynb           # Interactive alignment example
└── readme.md
```

## Limitations

- Only the supplied strand is searched; reverse-complement searching is not implemented.
- A hit requires an exact 11-base seed. Seed thresholds and the small candidate window can miss divergent sequences or larger indels.
- Scoring parameters and reporting thresholds are fixed in the source code.
- Reference sequences are loaded into memory. Dense seek tables require approximately 64 MiB per reference on disk, regardless of sequence length; this layout is intended for small databases.
- Results are console output only. E-values, bit scores, protein alignment, and structured result export are not implemented.

## Acknowledgments

Adapted from [JiaShun-Xiao/BLAST-bioinfor-tool](https://github.com/JiaShun-Xiao/BLAST-bioinfor-tool).
