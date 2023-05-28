import argparse
import json
import pkg_resources
from .dataset_downloading import download
from .export_data import export


def main():
    parser = argparse.ArgumentParser(description='YMS package')
    parser.add_argument('-v', '--version', action='version',
                        version='%(prog)s {version}'.format(version=pkg_resources.require("YMS")[0].version))

    parser.add_argument('--download', nargs='+', help='Add datasets to be downloaded')
    parser.add_argument('--type', nargs=1, choices=['16S', 'Shotgun'], help='Type of data to be downloaded')
    parser.add_argument('--export', nargs=7,
                        help="output_dir, trim, trunc, classifier_file, otu_output_file, taxonomy_output_file, threads")
    parser.add_argument('--config', help='Path to config file')
    parser.add_argument('--verbose', action='store_true', help='Enable verbose mode')
    args = parser.parse_args()

    if args.config:
        with open(args.config) as f:
            config = json.load(f)
        specific_location = config.get('specific_location')
    else:
        config_path = pkg_resources.resource_filename(__name__, "config.json")
        with open(config_path) as f:
            config = json.load(f)
        specific_location = config.get('specific_location')

    if args.export:
        try:
            output_dir = args.export[0]
            trim = int(args.export[1])
            trunc = int(args.export[2])
            classifier_file = args.export[3]
            export(output_dir, trim, trunc, classifier_file)
        except IndexError:
            print(f"missing {len(args.export)-1} arguments")

    if args.download:
        if not(args.type):
            raise ValueError("Missing dataset type. Use --type 16S/Shotgun")
        else:
            data_type = args.type[0]
            for dataset_name in args.download:
                download(dataset_name,data_type, args.verbose, specific_location)
