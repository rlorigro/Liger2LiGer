from modules.Fastx import *

from subprocess import run
import argparse
import struct
import mmap
import sys
import os


def load_query_ids(query_ids_path):
    ids = list()
    with open(query_ids_path, 'r') as file:
        for line in file:
            if line == '\n':
                continue

            ids.append(line.strip())

    return ids


def main(fastq_path, query_ids_path):
    faidx_path = build_index(fastq_path)

    name_to_offset, index_elements = load_fastq_index(faidx_path=faidx_path)

    queries = load_query_ids(query_ids_path=query_ids_path)

    output_path = os.path.splitext(fastq_path)[0] + "_" + os.path.splitext(os.path.basename(query_ids_path))[0] + ".fastq"
    sys.stderr.write("Writing to: " + output_path + '\n')

    with open(fastq_path, 'rb') as input_file, open(output_path, 'wb') as output_file:
        mm = mmap.mmap(input_file.fileno(), 0, prot=mmap.PROT_READ)

        for name in queries:
            print("fetching %s" % name)

            if name in name_to_offset:
                offset_index = name_to_offset[name]
            else:
                exit("ERROR: read name not found in fastq index")

            if offset_index < len(index_elements):
                index_element = index_elements[offset_index]
            else:
                exit("ERROR: attempted to access fastq index element " + offset_index + " which is greater than the "
                     "size of the list of indexes")


            s = extract_bytes_from_file(mmap_file_object=mm,
                                        offset=index_element.sequence_offset,
                                        n_bytes=index_element.length)

            q = extract_bytes_from_file(mmap_file_object=mm,
                                        offset=index_element.quality_offset,
                                        n_bytes=index_element.length)

            output_file.write(b'@')
            output_file.write(name.encode('utf-8'))
            output_file.write(b'\n')
            output_file.write(s)
            output_file.write(b'\n')
            output_file.write(b'+')
            output_file.write(b'\n')
            output_file.write(q)
            output_file.write(b'\n')

    return


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--fastq",
        type=str,
        required=True,
        help="path of file containing FASTA/FASTQ sequence"
    )
    parser.add_argument(
        "--ids",
        type=str,
        required=True,
        help="path of file containing 1 id per line to be queried"
    )

    args = parser.parse_args()

    main(
        fastq_path=args.fastq,
        query_ids_path=args.ids,
    )
