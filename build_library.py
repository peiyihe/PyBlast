from collections import defaultdict
from pathlib import Path

import numpy as np


WORD_LENGTH = 11
chrom_dict = {}


def read_fasta(path):
    """Return (name, sequence, byte offset) records, independent of line width."""
    records = []
    seen = set()
    name = None
    chunks = []
    offset = 0
    # Binary offsets account for CRLF and different wrapping in each record.
    with open(path, 'rb') as handle:
        for line in iter(handle.readline, b''):
            if line.startswith(b'>'):
                if name is not None:
                    records.append((name, ''.join(chunks), offset))
                header = line[1:].decode('utf-8').split()
                if not header or header[0] in seen:
                    raise ValueError('FASTA headers must have unique, nonempty names.')
                name = header[0]
                seen.add(name)
                chunks = []
                offset = handle.tell()
            elif line.strip():
                if name is None:
                    raise ValueError('FASTA sequence encountered before its header.')
                chunks.append(b''.join(line.split()).decode('ascii').upper())
    if name is not None:
        records.append((name, ''.join(chunks), offset))
    if not records or any(not sequence for _, sequence, _ in records):
        raise ValueError('FASTA must contain nonempty reference sequences.')
    return records


def BaseToNum(chr_seq):
    return chr_seq.upper().translate(str.maketrans('ACGT', '1234'))


def BaseToIndex(word, word_len):
    return sum((int(value) - 1) * 4 ** (word_len - i)
               for i, value in enumerate(word))


def GenSeek(library, word_len):
    seeks = np.zeros((4 ** word_len, 2), dtype=np.int64)
    offset = 0
    for i, entry in enumerate(library):
        seeks[i] = offset, len(entry)
        offset += len(entry)
    return seeks


def BuildLibrary(chr_name, data_dir='dataset', append=False):
    """Write one reference's positions and return its absolute seek offsets."""
    sequence = BaseToNum(chrom_dict[chr_name])
    positions = defaultdict(list)
    # Include the final word; stored sequence coordinates are one-based.
    for start in range(len(sequence) - WORD_LENGTH + 1):
        word = sequence[start:start + WORD_LENGTH]
        if set(word) <= set('1234'):
            positions[BaseToIndex(word, WORD_LENGTH - 1)].append(str(start + 1))

    seeks = np.zeros((4 ** WORD_LENGTH, 2), dtype=np.int64)
    with open(Path(data_dir) / 'sarscov2.txt', 'ab' if append else 'wb') as handle:
        for index in sorted(positions):
            entry = (','.join(positions[index]) + ',').encode('ascii')
            seeks[index] = handle.tell(), len(entry)
            handle.write(entry)
    if not append:
        np.save(Path(data_dir) / 'sarscov2_library_seeks.npy', seeks)
    return seeks


def build_libraries(data_dir='dataset'):
    """Regenerate all indexes for sarscov2.fasta in the selected directory."""
    global chrom_dict
    data_dir = Path(data_dir)
    records = read_fasta(data_dir / 'sarscov2.fasta')
    chrom_dict = {name: sequence for name, sequence, _ in records}
    names = [name for name, _, _ in records]
    np.save(data_dir / 'sarscov2_chr_names.npy', np.array(names))
    np.save(data_dir / 'sarscov2_chrom_seek_index.npy',
            np.array([(len(sequence), offset) for _, sequence, offset in records],
                     dtype=np.int64))

    print('Starting to build index library...')
    for i, name in enumerate(names):
        print(f'Processing: {name}')
        seeks = BuildLibrary(name, data_dir, append=i > 0)
        if i == 0 and len(names) > 1:
            # Keep the original 2-D format for single-reference databases.
            all_seeks = np.lib.format.open_memmap(
                data_dir / 'sarscov2_library_seeks.npy', mode='w+',
                dtype=np.int64, shape=(len(names), 4 ** WORD_LENGTH, 2))
        if len(names) > 1:
            all_seeks[i] = seeks
    if len(names) > 1:
        all_seeks.flush()
        del all_seeks
    print('Index library construction completed!')


if __name__ == '__main__':
    build_libraries()
