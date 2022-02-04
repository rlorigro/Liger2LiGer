from subprocess import run
import sys
import os


class FastqIndexElement:
    def __init__(self, sequence_offset, quality_offset, length):
        self.sequence_offset = sequence_offset
        self.quality_offset = quality_offset
        self.length = length

    def __str__(self):
        return str(self.sequence_offset) + "\t" + str(self.quality_offset) + "\t" + str(self.length)


class FastaIndexElement:
    def __init__(self, sequence_offset, length):
        self.sequence_offset = sequence_offset
        self.length = length

    def __str__(self):
        return str(self.sequence_offset) + "\t" + str(self.length)


'''
faidx format:
-------------
NAME	Name of this reference sequence
LENGTH	Total length of this reference sequence, in bases
OFFSET	Offset in the FASTA/FASTQ file of this sequence's first base
LINEBASES	The number of bases on each line
LINEWIDTH	The number of bytes in each line, including the newline
QUALOFFSET	Offset of sequence's first quality within the FASTQ file
'''
def load_fastq_index(faidx_path):
    name_to_offset = dict()
    index_elements = list()

    with open(faidx_path, 'r') as file:
        for line in file:
            data = line.strip().split('\t')

            name = data[0]
            length = int(data[1])
            sequence_offset = int(data[2])
            quality_offset = int(data[5])

            index_element = FastqIndexElement(sequence_offset, quality_offset, length)

            index_elements.append(index_element)

            name_to_offset[name] = len(index_elements) - 1

    return name_to_offset, index_elements


def load_fasta_index(faidx_path):
    name_to_offset = dict()
    index_elements = list()

    with open(faidx_path, 'r') as file:
        for line in file:
            data = line.strip().split('\t')

            name = data[0]
            length = int(data[1])
            sequence_offset = int(data[2])

            index_element = FastaIndexElement(sequence_offset, length)

            index_elements.append(index_element)

            name_to_offset[name] = len(index_elements) - 1

    return name_to_offset, index_elements


def extract_bytes_from_file(mmap_file_object, offset, n_bytes):
    s = None

    data = mmap_file_object[offset:offset+n_bytes]

    s = bytes(data)

    return s


# Use a system call to samtools faidx to build the index
def build_index(path):
    index_path = path + ".fai"
    if os.path.exists(index_path):
        sys.stderr.write("Index exists\n")
        sys.stderr.flush()

        if os.path.getmtime(index_path) < os.path.getmtime(path):
            sys.stderr.write("WARNING: fasta/q index path is older than file itself, may be out of date: " + index_path + "\n")
    else:
        sys.stderr.write("No index found, indexing... ")
        sys.stderr.flush()

        arguments = ["samtools", "faidx", path]
        run(arguments, check=True)
        sys.stderr.write("Done\n")

    return index_path

