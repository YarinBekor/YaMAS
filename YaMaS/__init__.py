import argparse

import pkg_resources

from .dataset_downloding import download
from .export_data import export

def main():
    parser = argparse.ArgumentParser(description='YaMaS package')
    parser.add_argument('-v', '--version', action='version',
                        version='%(prog)s {version}'.format(version=pkg_resources.require("YaMaS")[0].version))

    parser.add_argument('--download', nargs='+', help='Add datasets to be downloaded')
    parser.add_argument('--verbose', action='store_true', help='Enable verbose mode')
    parser.add_argument('--export', nargs='+', help="output_dir, trim, trunc, classifier_file, otu_output_file, taxonomy_output_file, threads")
    args = parser.parse_args()

    if args.export:
        export(args.export)

    if args.download:
        for dataset_name in args.download:
            download(dataset_name, args.verbose)



