from __future__ import annotations

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

def metaphlan_extraction(reads_data):
    paired = reads_data.rev and reads_data.fwd
    fastq_path = os.path.join(reads_data.dir_path, "fastq")
    export_path = os.path.join(reads_data.dir_path, "export")
    run_cmd([f"mkdir {export_path}"])
    final_output_path = os.path.join(export_path, 'final.txt')
    run_cmd([f"touch {final_output_path}"])
    fastq_files = [a for a in os.listdir(fastq_path) if a.split(".")[-1] == "fastq"]
    args = []
    if paired:
        print("paired")
        for i in tqdm(range(0,len(fastq_files),2)):
            fastq_name = fastq_files[i].split('_')[0]
            fastq_1 = os.path.join(fastq_path, fastq_files[i])
            fastq_2 = os.path.join(fastq_path, fastq_files[i+1])
            output = os.path.join(fastq_path,f"{fastq_name}.bowtie2.bz2")
            command = [f"metaphlan {fastq_1},{fastq_2} --input_type fastq --bowtie2out {output} --nproc 24"]
            run_cmd(command)
            final_output_file = os.path.join(os.path.join(reads_data.dir_path, 'qza'), f'{fastq_name}_profile.txt')
            command = [f"metaphlan {output} --input_type bowtie2out --nproc 24 > {final_output_file}"]
            run_cmd(command)
            args.append(final_output_file)
    else:
        print("not paired")
        for fastq in tqdm(fastq_files):
            output = os.path.join(os.path.join(reads_data.dir_path, 'qza'), f'{fastq}_profile.txt')
            fastq = os.path.join(fastq_path,fastq)
            command = [f"metaphlan {fastq} --input_type fastq --nproc 24 > {output}"]
            run_cmd(command)
            args.append(output)
    merge(args, open(final_output_path, 'w'), gtdb=False)

def metaphlan_txt_csv(reads_data, dataset_id):
    export_path = os.path.join(reads_data.dir_path, "export")
    input_file = os.path.join(export_path,f"{dataset_id}_final.txt")
    output_file = os.path.join(export_path,f"{dataset_id}_final_table.csv")
    with open(input_file, 'r') as txt_file:
        lines = txt_file.readlines()

    # Extract data from the text file
    headers = lines[0].strip().split('\t')
    data = [line.strip().split('\t') for line in lines[1:]]

    # Replace "|" with ","
    headers = [header.replace('|', ',') for header in headers]
    data = [[entry.replace('|', ',') for entry in row] for row in data]

    # Transpose the data
    transposed_data = list(map(list, zip(*data)))

    # Write the transposed data to a CSV file
    with open(output_file, 'w', newline='') as csv_file:
        writer = csv.writer(csv_file)
        writer.writerow(headers)
        writer.writerows(transposed_data)

    print(f"CSV file '{output_file}' has been created.")

def qiime_import(dir_path,fastq_path):
    qza_path = os.path.join(dir_path, "qza")

    multiplexed_qza_file_path = os.path.join(qza_path, "multiplexed-seqs.qza")

    command = [
        "qiime", "tools", "import",
        "--type", "MultiplexedSingleEndBarcodeInSequence",
        "--input-path", fastq_path,
        "--output-path", multiplexed_qza_file_path,
    ]
    run_cmd(command)
    return multiplexed_qza_file_path


def qiime_demux(dir_path,multiplexed_qza_file_path, metadata_path):

    demux_qza_file_path = os.path.join(dir_path,"qza", "demux-single-end.qza")
    untrim_qza_file_path = os.path.join(dir_path, "qza","untrimmed.qza")


    command = [
        "qiime", "cutadapt", "demux-single",
        "--i-seqs", multiplexed_qza_file_path,
        "--m-barcodes-file", metadata_path,
        "--m-barcodes-column", "barcode",
        "--p-error-rate", "0",
        "--o-per-sample-sequences", demux_qza_file_path,
        "--o-untrimmed-sequences", untrim_qza_file_path,
        "--verbose"
    ]
    run_cmd(command)
    return demux_qza_file_path

def check_metadata(metadata_path):

    metadata = pd.read_csv(metadata_path, sep='\t')
    if 'barcode' in metadata.columns:
        return "yes"
    return "no"


def trim_single(dir_path,demux_qza_file_path):
    #Trim adapters from demultiplexed reads
    #If there are sequencing adapters or PCR primers in the reads which you'd like to remove, you can do that next as follows.

    trimmed_seqs_file_path = os.path.join(dir_path,"qza", "trimmed-seqs.qza")
    command = [
        "qiime", "cutadapt", "trim-single",
        "--i-demultiplexed-sequences", demux_qza_file_path,
        "--p-front", "GCTACGGGGGG",
        "--p-error-rate", "0",
        "--o-trimmed-sequences", trimmed_seqs_file_path,
        "--verbose"
    ]
    run_cmd(command)
    return trimmed_seqs_file_path
def qiime_summarize(dir_path,trimmed_seqs_file_path):
    #Summarize demultiplexed and trimmed reads
    vis_path = os.path.join(dir_path, "vis")
    vis_file_path = os.path.join(vis_path, "trimmed-seqs.qzv")

    command=[
        "qiime", "demux", "summarize",
        "--i-data", trimmed_seqs_file_path,
        "--o-visualization", vis_file_path
    ]
    run_cmd(command)
    return vis_file_path


def get_reads_data(dir_path,demux_qza_file_path):

    output_dir_path = os.path.join(dir_path,"extracted-reads")

    #extract data from demultiplexed-seqs.qza file
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




def qiita_visualization(fastq_path,metadata_path,data_type, verbose_print):

    verbose_print("\n")
    verbose_print(datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S'))
    dir_path = os.path.commonpath([os.path.dirname(fastq_path), os.path.dirname(metadata_path)])

    verbose_print("\n")
    verbose_print('Checking environment...', end=" ")
    check_conda_qiime2()
    verbose_print('Done.')
    verbose_print("\n")
    verbose_print('Checking metadata...', end= " ")
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
        verbose_print(f"{datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S')} -- Start 'Import the multiplexed sequences' (1/4)")
        multiplexed_qza_file_path= qiime_import(dir_path,fastq_path)
        verbose_print(f"{datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S')} -- Finish 'Import the multiplexed sequences' (1/4)")
        verbose_print("\n")

        verbose_print(f"{datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S')} -- Start 'Demultiplex the reads' (2/4)")
        demux_qza_file_path=qiime_demux(dir_path, multiplexed_qza_file_path,metadata_path)
        verbose_print(f"{datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S')} -- Finish 'Demultiplex the reads' (2/4)")
        verbose_print("\n")

        verbose_print(f"{datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S')} -- Start 'Trim adapters from demultiplexed reads' (3/4)")
        trimmed_seqs_file_path=trim_single(dir_path,demux_qza_file_path)
        verbose_print(f"{datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S')} -- Finish 'rim adapters from demultiplexed reads' (3/4)")
        verbose_print("\n")

        verbose_print(f"{datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S')} -- Start 'Summarize demultiplexed and trimmed reads' (4/4)")
        vis_file_path= qiime_summarize(dir_path,trimmed_seqs_file_path)
        verbose_print(f"{datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S')} -- Finish 'Summarize demultiplexed and trimmed reads' (4/4)")

        #getting values about fwd and rev
        reads_data= get_reads_data(dir_path,demux_qza_file_path)
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