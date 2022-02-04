import argparse
import sys


def filter_by_id(paf_path, names):
    """
    PAF has columns (0-based):
        0 = read_name
    """

    with open(paf_path, 'r') as file:
        for l,line in enumerate(file):
            tokens = line.strip().split()

            read_name = tokens[0]

            if read_name in names:
                print(line.strip())

    return


def load_names_as_set(path):
    names = set()
    with open(path, 'r') as file:
        for line in file:
            name = line.strip()

            names.add(name)

    return names


def main(paf_path, names_path):
    names = load_names_as_set(path=names_path)
    filter_by_id(paf_path=paf_path, names=names)

    return


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--paf",
        type=str,
        required=True,
        help="path of PAF file to be filtered"
    )
    parser.add_argument(
        "--names",
        type=str,
        required=True,
        help="path of file containing a list of read ids, one on each line"
    )

    args = parser.parse_args()

    main(
        paf_path=args.paf,
        names_path=args.names
    )
