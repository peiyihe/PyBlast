from collections import Counter
from pathlib import Path

import numpy as np

from build_library import WORD_LENGTH, read_fasta


# Global variables initialized by init_blast().
chr_names = None
chrom_seek_index = None
fasta_file = None
fasta_line_width = None
reference_sequences = None
library_seeks = None
library_path = None


def init_blast(data_dir='dataset'):
    """Load reference sequences and memory-map the seed index once."""
    global chr_names, chrom_seek_index, fasta_file, fasta_line_width
    global reference_sequences, library_seeks, library_path
    data_dir = Path(data_dir)
    records = read_fasta(data_dir / 'sarscov2.fasta')
    names = np.load(data_dir / 'sarscov2_chr_names.npy')
    chrom_index = np.load(data_dir / 'sarscov2_chrom_seek_index.npy')
    seeks = np.load(data_dir / 'sarscov2_library_seeks.npy', mmap_mode='r')
    if seeks.ndim == 2:
        seeks = seeks[np.newaxis, :, :]
    expected_chrom_index = np.array(
        [(len(sequence), offset) for _, sequence, offset in records])
    if (names.tolist() != [name for name, _, _ in records]
            or not np.array_equal(chrom_index, expected_chrom_index)
            or seeks.shape != (len(records), 4 ** WORD_LENGTH, 2)):
        raise ValueError('Reference indexes do not match the FASTA. Rebuild with build_libraries(data_dir).')
    new_library_path = data_dir / 'sarscov2.txt'
    if not new_library_path.is_file():
        raise FileNotFoundError(new_library_path)
    new_fasta_file = open(data_dir / 'sarscov2.fasta')
    if fasta_file is not None:
        fasta_file.close()
    chr_names, chrom_seek_index = names, chrom_index
    reference_sequences = [sequence for _, sequence, _ in records]
    library_seeks, library_path = seeks, new_library_path
    fasta_file = new_fasta_file
    fasta_file.readline()
    fasta_line_width = len(fasta_file.readline().strip())
    fasta_file.seek(0)
    print(f"BLAST initialization complete! Loaded {len(chr_names)} sequence(s).")
    return chr_names, chrom_seek_index, fasta_file


def _require_initialized():
    if reference_sequences is None:
        raise RuntimeError('Call init_blast() before querying the reference.')


def SingleBaseCompare(seq1, seq2, i, j):
    # Ambiguous bases are mismatches, even when the symbols are identical.
    return 2 if seq1[i] == seq2[j] and seq1[i] in 'ACGT' else -1


def _local_alignment(seq1, seq2):
    """Smith-Waterman with +2 match, -1 mismatch, and -3 linear gap cost.

    Return aligned strings, identity, and zero-based half-open intervals.
    """
    seq1, seq2 = seq1.upper(), seq2.upper()
    m, n = len(seq1), len(seq2)
    matrix = [[0] * (n + 1) for _ in range(m + 1)]
    best_score, end_i, end_j = 0, 0, 0
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            matrix[i][j] = max(
                0,
                matrix[i - 1][j - 1] + SingleBaseCompare(seq1, seq2, i - 1, j - 1),
                matrix[i - 1][j] - 3,
                matrix[i][j - 1] - 3,
            )
            if matrix[i][j] > best_score:
                best_score, end_i, end_j = matrix[i][j], i, j

    i, j = end_i, end_j
    aligned1, aligned2 = [], []
    while i > 0 and j > 0 and matrix[i][j] > 0:
        # Follow a transition that produced this cell, including its cost.
        if matrix[i][j] == matrix[i - 1][j - 1] + SingleBaseCompare(seq1, seq2, i - 1, j - 1):
            aligned1.append(seq1[i - 1])
            aligned2.append(seq2[j - 1])
            i, j = i - 1, j - 1
        elif matrix[i][j] == matrix[i - 1][j] - 3:
            aligned1.append(seq1[i - 1])
            aligned2.append('-')
            i -= 1
        else:
            aligned1.append('-')
            aligned2.append(seq2[j - 1])
            j -= 1
    align_seq1 = ''.join(reversed(aligned1))
    align_seq2 = ''.join(reversed(aligned2))
    matches = sum(a == b and a in 'ACGT' for a, b in zip(align_seq1, align_seq2))
    identity = matches / len(align_seq1) if align_seq1 else 0.0
    return align_seq1, align_seq2, identity, i, end_i, j, end_j


def SMalignment(seq1, seq2):
    """Return a local alignment and its identity fraction (0 to 1)."""
    return _local_alignment(seq1, seq2)[:3]


# Display BlAST result
def Display(seque1, seque2):
    # fix bug the first le - 40 is not 0, previous code is le = 60
    le = 40
    while len(seque1)-le >= 0:
        print('sequence1: ',end='')
        for a in list(seque1)[le-40:le]:
            print(a,end='')
        print("\n")
        print('           ',end='')
        for k in range(le-40, le):
            if seque1[k] == seque2[k] and seque1[k] in 'ACGT':
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
            if seque1[k] == seque2[k] and seque1[k] in 'ACGT':
                print('|',end='')
            else:
                print(' ',end='')
        print("\n")
        print('sequence2: ',end='')
        for b in list(seque2)[le-40:len(seque2)]:
            print(b,end='')
        print("\n")

# Transform a canonical DNA word to its base-4 index.
def WordToNum(word):
    trans = {'A': 1, 'C': 2, 'G': 3, 'T': 4}
    return [trans[base] for base in word.upper()]


def WordToIndex(word, word_len):
    return sum((value - 1) * 4 ** (word_len - i)
               for i, value in enumerate(WordToNum(word)))


def GetWordPos(word):
    """Return one list of one-based seed positions per reference."""
    _require_initialized()
    word = word.upper()
    if len(word) != WORD_LENGTH:
        raise ValueError(f'Seed words must contain {WORD_LENGTH} bases.')
    if not set(word) <= set('ACGT'):
        return [[] for _ in chr_names]
    seek_index = WordToIndex(word, WORD_LENGTH - 1)
    positions = []
    with open(library_path, 'rb') as handle:
        for seeks in library_seeks:
            offset, length = map(int, seeks[seek_index])
            handle.seek(offset)
            entry = handle.read(length).decode('ascii').rstrip(',')
            positions.append([int(value) for value in entry.split(',')] if entry else [])
    return positions


def ExtractSeq(chr_index, pos, length):
    """Extract bases using a one-based start, clipping at the reference end."""
    _require_initialized()
    if not 0 <= chr_index < len(reference_sequences):
        raise IndexError('Reference index out of range.')
    if pos < 1 or length < 0:
        raise ValueError('Sequence position must be >= 1 and length must be >= 0.')
    # Slice normalized bases, so newlines never consume requested length.
    return reference_sequences[chr_index][pos - 1:pos - 1 + length]


def Blast(query_seq):
    _require_initialized()
    query_seq = ''.join(query_seq.split()).upper()
    if len(query_seq) < WORD_LENGTH:
        raise ValueError(f'Query must contain at least {WORD_LENGTH} bases.')
    if not set(query_seq) <= set('ACGTRYSWKMBDHVN'):
        raise ValueError('Query must contain DNA IUPAC bases only.')

    words_positions = []
    for word_index in range(len(query_seq) - WORD_LENGTH + 1):
        word = query_seq[word_index:word_index + WORD_LENGTH]
        if set(word) <= set('ACGT'):
            words_positions.append((word_index, GetWordPos(word)))
    # Retain the six-seed threshold for longer queries; allow 11-15 bases.
    threshold = min(6, len(words_positions))
    found = False
    for chr_index, reference in enumerate(reference_sequences):
        starts = Counter(
            pos - word_index
            for word_index, positions in words_positions
            for pos in positions[chr_index]
        )
        reported = set()
        for start, count in starts.items():
            if count < threshold:
                continue
            # A small extension allows local alignment to recover nearby indels.
            candidate_start = max(1, start - 5)
            candidate_end = min(len(reference), start + len(query_seq) - 1 + 5)
            candidate = ExtractSeq(chr_index, candidate_start,
                                   max(0, candidate_end - candidate_start + 1))
            aligned1, aligned2, identity, ref_start, ref_end, query_start, query_end = _local_alignment(candidate, query_seq)
            coverage = (query_end - query_start) / len(query_seq)
            # Identity alone would accept tiny local matches to a long query.
            if identity <= 0.8 or coverage < 0.8:
                continue
            first = candidate_start + ref_start
            last = candidate_start + ref_end - 1
            key = (first, last, aligned1, aligned2)
            if key in reported:
                continue
            reported.add(key)
            found = True
            print(f'find in chromosome {chr_names[chr_index]}: {first} ---> {last}, align score: {identity}')
            Display(aligned1, aligned2)
    if not found:
        print('No alignments found (requires an exact 11-base seed).')
    return None


if __name__ == "__main__":
    # Initialize BLAST
    init_blast()
    
    # Execute query
    query_sequence = 'TAACCAGAATGGAGAACGCAGTGGGGCGCGATCAAAACAACGTCGGCCCCAAGGTTTACCCAATAATACT'
    Blast(query_sequence)
