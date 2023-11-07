import argparse
import json
import pkg_resources
from .dataset_downloading import download
from .export_data import export
from .prerun_configs import set_environment

def main():
    # Initialize the argument parser with a description.
    parser = argparse.ArgumentParser(description='YMS package')

    # Add an argument for displaying the version of the YMS package.
    parser.add_argument('-v', '--version', action='version',
                        version='%(prog)s {version}'.format(version=pkg_resources.require("YMS")[0].version))

    parser.add_argument('--ready',nargs=1,choices=['Ubuntu', 'CentOS'], help='Type of Operating system')
    
    # Add an argument for specifying datasets to be downloaded.
    parser.add_argument('--download', nargs='+', help='Add datasets to be downloaded')

    # Add an argument for specifying the type of data to be downloaded (16S or Shotgun).
    parser.add_argument('--type', nargs=1, choices=['16S', 'Shotgun'], help='Type of data to be downloaded')

    # Add an argument for specifying export parameters.
    parser.add_argument('--export', nargs=6,
                        help="origin_dir_path, data_type, start, end, classifier_file, threads")

    # Add an argument for specifying the path to a configuration file.
    parser.add_argument('--config', help='Path to config file')

    # Add a flag for enabling verbose mode.
    parser.add_argument('--verbose', action='store_true', help='Enable verbose mode')

    # Parse the command line arguments.
    args = parser.parse_args()


    if args.ready:
        set_environment(args.ready[0])

    if args.config:
        # If a config file path is provided, load the configuration from the file.
        with open(args.config) as f:
            config = json.load(f)
        specific_location = config.get('specific_location')
    else:
        # If no config file path is provided, use the default configuration bundled with the package.
        config_path = pkg_resources.resource_filename(__name__, "config.json")
        with open(config_path) as f:
            config = json.load(f)
        specific_location = config.get('specific_location')

    if args.export:
        try:
            # Extract export parameters from the command line arguments.
            origin_dir = args.export[0]
            data_type= args.export[1]
            trim = int(args.export[2])
            trunc = int(args.export[3])
            classifier_file = args.export[4]
            threads = args.export[5]

            # Call the export function with the specified parameters.
            export(origin_dir,data_type,trim, trunc, classifier_file, threads)
        except IndexError:
            # Handle the case where the number of export arguments is insufficient.
            print(f"missing {len(args.export)-1} arguments")

    if args.download:
        if not(args.type):
            # Ensure that a dataset type is specified when downloading datasets.
            raise ValueError("Missing dataset type. Use --type 16S/Shotgun")
        else:
            # Extract the dataset type and iterate over the specified dataset names for downloading.
            data_type = args.type[0]
            for dataset_name in args.download:
                download(dataset_name, data_type, args.verbose, specific_location)
