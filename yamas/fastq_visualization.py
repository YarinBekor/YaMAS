import os

import pandas as pd
import csv
import gzip
import os.path
import pickle
import datetime
from tqdm import tqdm
from metaphlan.utils.merge_metaphlan_tables import merge
from .utilities import run_cmd, ReadsData, check_conda_qiime2
import json
import shutil
import tarfile
import os
import yaml


def check_metadata(metadata_path):
    metadata = pd.read_csv(metadata_path, sep='\t')
    if 'barcode' in metadata.columns:
        return "yes"
    return "no"


def qiime_import(dir_path):
    qza_path = os.path.join(dir_path, "qza")
    multiplexed_qza_file_path = os.path.join(qza_path, "multiplexed-seqs.qza")

    command = [
        "qiime", "tools", "import",
        "--type", "EMPSingleEndSequences",
        "--input-path", dir_path,
        "--output-path", multiplexed_qza_file_path
    ]

    run_cmd(command)
    return multiplexed_qza_file_path


def qiime_demux(dir_path, multiplexed_qza_file_path, metadata_path):
    demux_qza_file_path = os.path.join(dir_path, "qza", "demux-single-end.qza")
    demux_details_single = os.path.join(dir_path, "qza", "demux_details_single.qza")

    command = [
        "qiime", "demux", "emp-single",
        "--m-barcodes-file", metadata_path,
        "--m-barcodes-column", "barcode",
        "--i-seqs", multiplexed_qza_file_path,
        "--o-per-sample-sequences", demux_qza_file_path,
        "--o-error-correction-details", demux_details_single,
        "--p-rev-comp-barcodes",
        "--p-rev-comp-mapping-barcodes"
    ]
    run_cmd(command)
    return demux_qza_file_path


def qiime_summarize(dir_path, demux_qza_file_path):
    vis_path = os.path.join(dir_path, "vis")
    vis_file_path = os.path.join(vis_path, "dataset_view.qzv")
    command = [
        "qiime", "demux", "summarize",
        "--i-data", demux_qza_file_path,
        "--o-visualization", vis_file_path
    ]
    run_cmd(command)
    return vis_file_path


def get_reads_data(dir_path, demux_qza_file_path):
    output_dir_path = os.path.join(dir_path, "extracted-reads")

    # extract data from demultiplexed-seqs.qza file
    command = [
        "qiime", "tools", "extract",
        "--input-path", demux_qza_file_path,
        "--output-path", output_dir_path
    ]
    run_cmd(command)

    os.chdir(output_dir_path)
    subdirectories = [d for d in os.listdir() if os.path.isdir(d)]
    if not subdirectories:
        print("Something wrong with the data. Its doesn't follow the rules of qiita.")
        return
    else:
        subdirectory_path = subdirectories[0]

    os.chdir(subdirectory_path)

    metadata_file_path = 'metadata.yaml'
    if os.path.isfile(metadata_file_path):
        with open(metadata_file_path, 'r') as metadata_file:
            metadata_content = yaml.safe_load(metadata_file)

            # Check if 'single' is in the metadata content
            if 'Single' in str(metadata_content):
                return ReadsData(dir_path, fwd=True, rev=False)

            else:
                return ReadsData(dir_path, fwd=True, rev=True)
    else:
        print("Something wrong with the data. Its doesn't follow the rules of qiita")


def fastq_visualization(fastq_path, barcode_path, metadata_path, data_type, verbose_print):
    verbose_print("\n")
    verbose_print(datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S'))
    dir_path = os.path.commonpath([os.path.dirname(fastq_path), os.path.dirname(metadata_path), os.path.dirname(barcode_path)])

    verbose_print("\n")
    verbose_print('Checking environment...', end=" ")
    check_conda_qiime2()
    verbose_print('Done.')
    verbose_print("\n")
    verbose_print('Checking metadata...', end=" ")
    if check_metadata(metadata_path) == "no":
        verbose_print("The 'barcode' column does not exist in metadata.tsv. check and try again.")
        return
    verbose_print('Done.')
    verbose_print("\n")
    run_cmd(["mkdir", os.path.join(dir_path, "qza")])
    run_cmd(["mkdir", os.path.join(dir_path, "vis")])

    verbose_print("Find ALL NEW data in the directory you created:", dir_path)

    if data_type == '16S' or data_type == '18S':

        verbose_print("\n")
        verbose_print(
            f"{datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S')} -- Start 'Import the multiplexed sequences' (1/4)")
        multiplexed_qza_file_path = qiime_import(dir_path)
        verbose_print(
            f"{datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S')} -- Finish 'Import the multiplexed sequences' (1/4)")
        verbose_print("\n")

        verbose_print(f"{datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S')} -- Start 'Demultiplex the reads' (2/4)")
        demux_qza_file_path = qiime_demux(dir_path, multiplexed_qza_file_path, metadata_path)
        verbose_print(
            f"{datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S')} -- Finish 'Demultiplex the reads' (2/4)")
        verbose_print("\n")

        verbose_print(
            f"{datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S')} -- Start 'Summarize demultiplexed and trimmed reads' (3/3)")
        vis_file_path = qiime_summarize(dir_path, demux_qza_file_path)
        verbose_print(
            f"{datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S')} -- Finish 'Summarize demultiplexed and trimmed reads' (3/3)")

        # getting values about fwd and rev
        reads_data = get_reads_data(dir_path, demux_qza_file_path)
        verbose_print(f"{datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S')} -- Finish creating visualization\n")

        pickle.dump(reads_data, open(os.path.join(reads_data.dir_path, "reads_data.pkl"), "wb"))

        print(f"Visualization file is located in {vis_file_path}\n"
              f"Please drag this file to https://view.qiime2.org/ and continue.\n")
        if reads_data.fwd and reads_data.rev:
            print(f"Note: The data has both forward and reverse reads.\n"
                  f"Therefore, you must give the parameters 'trim' and 'trunc' of export() "
                  f"as a tuple of two integers."
                  f"The first place related to the forward read and the second to the reverse.")
        else:
            print(f"Note: The data has only a forward read.\n"
                  f"Therefore, you must give the parameters 'trim' and 'trunc' of export() "
                  f"exactly one integers value which is related to the forward read.")

        return reads_data.dir_path

    else:
        print("YaMAS doesnt support downloading shotgun, yet.")
