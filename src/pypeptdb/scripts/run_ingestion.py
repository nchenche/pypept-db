import argparse
from pathlib import Path

from pypeptdb.ingestion.ingest import ingest_data_to_pypeptdb


def get_parser():
    parser = argparse.ArgumentParser(
        prog='pypept-db data ingestion',
        description='Command line to ingest data from an SDF file to the pypept database.'
    )

    parser.add_argument(
        '--sdf-file', '-F',
        type=str,
        required=True,
        help='Monomers SDF file. See expected format [here](https://github.com/Boehringer-Ingelheim/pyPept/blob/master/src/pyPept/data/monomers.sdf)'
    )

    args = parser.parse_args()
    return args


def main():
    args = get_parser()
    sdf_file = args.sdf_file
    ingest_data_to_pypeptdb(source=sdf_file)


if __name__ == '__main__':
    main()
