import argparse
from datetime import datetime
from pathlib import Path

from utils.data_pypept import get_combined_sdf




def get_parser():
    timestamp = datetime.now().strftime("%Y%m%d%H%M%S")

    parser = argparse.ArgumentParser(
        prog='pypept-db data conversion to an SDF file',
        description='Command line to convert data from the pypept database to an SDF file.'
    )

    parser.add_argument(
        '--outdir', '-O',
        type=str,
        default='.',
        required=False,
        help='Path where to output the monomers SDF file. Default to ".".)'
    )

    parser.add_argument(
        '--outname', '-N',
        type=str,
        default=f"monomers_{timestamp}.sdf",
        required=False,
        help='Name of the monomers SDF file.)'
    )

    args = parser.parse_args()
    return args


def main():
    args = get_parser()
    outdir = args.outdir
    outname = args.outname
    output = Path(outdir) / Path(outname)

    combined_sdf = get_combined_sdf()
    
    output.parent.mkdir(exist_ok=True, parents=True)
    with open(output, 'w') as outfile:
        outfile.write(combined_sdf)

if __name__ == '__main__':
    main()
