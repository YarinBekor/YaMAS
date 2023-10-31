from __future__ import annotations

import csv
import os.path
import pickle
import datetime
from tqdm import tqdm
from metaphlan.utils.merge_metaphlan_tables import merge
from .utilities import run_cmd, ReadsData, check_conda_qiime2

CONDA_PREFIX = os.environ.get("CONDA_PREFIX", None)


def check_input(acc_list: str):
    # flag = False

    print("Input path:", acc_list, "... ", end=" ")
    if not (os.path.exists(acc_list)):
        FileNotFoundError("Invalid. File does not exist.")
    elif not (os.path.isfile(acc_list)):
        FileNotFoundError("Invalid. Not a file.")
    else:
        print("Valid.")


import os

def create_dir(dir_name, specific_location):
    dir_path = os.path.join(os.path.abspath(specific_location), dir_name)

    os.makedirs(dir_path, exist_ok=True)
    print(f"{dir_path} created.")

    sra_path = os.path.join(dir_path, 'sra')
    os.makedirs(sra_path, exist_ok=True)
    print(f"{sra_path} created.")

    fastq_path = os.path.join(dir_path, 'fastq')
    os.makedirs(fastq_path, exist_ok=True)
    print(f"{fastq_path} created.")

    qza_path = os.path.join(dir_path, 'qza')
    os.makedirs(qza_path, exist_ok=True)
    print(f"{qza_path} created.")

    vis_path = os.path.join(dir_path, 'vis')
    os.makedirs(vis_path, exist_ok=True)
    print(f"{vis_path} created.")

    return dir_path


def download_data_from_sra(dir_path: str, acc_list: str = ""):
    run_cmd(['prefetch',
             "--option-file", acc_list,
             "--output-directory", os.path.join(dir_path, "sra"),
             "--max-size", "u"])


def sra_to_fastq(dir_path: str,data_type):
    print(f"converting files from .sra to .fastq.")
    for sra_dir in tqdm(os.listdir(os.path.join(dir_path, "sra")), desc="converted files"):
        sra_file = os.listdir(os.path.join(dir_path, "sra", sra_dir))[0]
        sra_path = os.path.join(dir_path, "sra", sra_dir, sra_file)
        fastq_path = os.path.join(dir_path, "fastq")
        run_cmd(["fasterq-dump", "--split-files", sra_path, "-O", fastq_path])
    if data_type == '16S':
        # check if reads include fwd and rev
        fastqs = sorted(os.listdir(os.path.join(dir_path, "fastq")))[:3]
        if len(set([fastq.split("_")[0] for fastq in fastqs])) == 1:
            return ReadsData(dir_path, fwd=True, rev=True)
        return ReadsData(dir_path, fwd=True, rev=False)
    else:
        # check if reads include fwd and rev
        fastqs = os.listdir(os.path.join(dir_path, "fastq"))
        for fastq in fastqs:
            if '_' in fastq:
                return ReadsData(dir_path, fwd=True, rev=True)
        return ReadsData(dir_path, fwd=True, rev=False)


def create_manifest(reads_data: ReadsData):
    fastq_path = os.path.join(reads_data.dir_path, "fastq")

    if not reads_data.rev:
        files = [os.path.join(fastq_path, f) for f in os.listdir(fastq_path)
                 if os.path.isfile(os.path.join(fastq_path, f))]
        names = [f.split('/')[-1].split('.')[0] for f in files]

        with open(os.path.join(reads_data.dir_path, 'manifest.tsv'), 'w') as manifest:
            tsv_writer = csv.writer(manifest, delimiter='\t')
            tsv_writer.writerow(["SampleID", "absolute-filepath"])
            for n, f in zip(*(names, files)):
                tsv_writer.writerow([n, f])
        return

    files_fwd = sorted([os.path.join(fastq_path, f) for f in os.listdir(fastq_path)
                        if os.path.isfile(os.path.join(fastq_path, f)) and "_1" in f])
    files_rev = sorted([os.path.join(fastq_path, f) for f in os.listdir(fastq_path)
                        if os.path.isfile(os.path.join(fastq_path, f)) and "_2" in f])
    names = sorted([f.split('.')[0] for f in os.path.join(reads_data.dir_path, "sra")])

    with open(os.path.join(reads_data.dir_path, 'manifest.tsv'), 'w') as manifest:
        tsv_writer = csv.writer(manifest, delimiter='\t')
        tsv_writer.writerow(["SampleID", "forward-absolute-filepath", "reverse-absolute-filepath"])
        for n, ff, fr in zip(*(names, files_fwd, files_rev)):
            tsv_writer.writerow([n, ff, fr])


def qiime_import(reads_data: ReadsData):
    qza_path = os.path.join(reads_data.dir_path, "qza")
    paired = reads_data.rev and reads_data.fwd

    qza_file_path = os.path.join(qza_path, f"demux-{'paired' if paired else 'single'}-end.qza")
    command = [
        "qiime", "tools", "import",
        "--type", f"SampleData[{'PairedEndSequencesWithQuality' if paired else 'SequencesWithQuality'}]",
        "--input-path", f"{os.path.join(reads_data.dir_path, 'manifest.tsv')}",
        "--input-format", "PairedEndFastqManifestPhred33V2" if paired else "SingleEndFastqManifestPhred33V2",
        "--output-path", qza_file_path
    ]
    run_cmd(command)

    return qza_file_path


def qiime_demux(reads_data: ReadsData, qza_file_path: str, dataset_id):
    vis_file_path = os.path.join(reads_data.dir_path, "vis", dataset_id + ".qzv")

    command = [
        "qiime", "demux", "summarize",
        "--i-data", qza_file_path,
        "--o-visualization", vis_file_path
    ]
    run_cmd(command)
    return vis_file_path

def metaphlan_extraction(reads_data, dataset_id):
    paired = reads_data.rev and reads_data.fwd
    fastq_path = os.path.join(reads_data.dir_path, "fastq")
    export_path = os.path.join(reads_data.dir_path, "export")
    run_cmd([f"mkdir {export_path}"])
    final_output_path = os.path.join(export_path, f'{dataset_id}_final.txt')
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

# This function is the main function to download the project. It Hnadles all the download flow for 16S and Shotgun.
def visualization(acc_list, dataset_id, data_type, verbose_print, specific_location):
    verbose_print("\n")
    verbose_print(datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S'))

    # Decide directory name
    dir_name = f"{dataset_id}-{datetime.datetime.now().strftime('%d-%m-%Y_%H-%M-%S')}"

    # Check if the enviorment is qiime2 type.
    verbose_print("\n")
    verbose_print('Checking environment...', end=" ")
    check_conda_qiime2()
    verbose_print('Done.')
    verbose_print("Checking inputs:")
    check_input(acc_list)

    verbose_print("\n")
    verbose_print("Creating a new directory for this dataset import:'", dir_name, "'")
    # Path to working dir:
    dir_path = create_dir(dir_name,specific_location)
    verbose_print("Find ALL NEW data in:", dir_path)

    # This two stages are for all data type. Prefetching the data from the relevant database, and converting it to fatsqs.
    verbose_print("\n")
    verbose_print(f"{datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S')} -- Start prefetch (1/5)")
    download_data_from_sra(dir_path, acc_list)
    verbose_print(f"{datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S')} -- Finish prefetch (1/5)")

    verbose_print("\n")
    verbose_print(f"{datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S')} -- Start conversion (2/5)")
    reads_data = sra_to_fastq(dir_path,data_type)
    verbose_print(f"{datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S')} -- Finish conversion (2/5)")



    if data_type == '16S':

        verbose_print("\n")
        verbose_print(f"{datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S')} -- Start creating manifest (3/5)")
        create_manifest(reads_data)
        verbose_print(f"{datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S')} -- Finish creating manifest (3/5)")

        verbose_print("\n")
        verbose_print(f"{datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S')} -- Start 'qiime import' (4/5)")
        qza_file_path = qiime_import(reads_data)
        verbose_print(f"{datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S')} -- Finish 'qiime import' (4/5)")

        verbose_print("\n")
        verbose_print(f"{datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S')} -- Start 'qiime demux' (5/5)")
        vis_file_path = qiime_demux(reads_data, qza_file_path, dataset_id)
        verbose_print(f"{datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S')} -- Finish 'qiime demux' (5/5)")

        pickle.dump(reads_data, open(os.path.join(reads_data.dir_path, "reads_data.pkl"), "wb"))
        verbose_print(f"{datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S')} -- Finish creating visualization\n")

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
        verbose_print("\n")
        verbose_print(f"{datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S')} -- Start metaphlan extraction (4/5)")
        metaphlan_extraction(reads_data, dataset_id)
        verbose_print(f"{datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S')} -- Finish metaphlan extraction (4/5)")
        
        verbose_print("\n")
        verbose_print(f"{datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S')} -- Start converting resualts to CSV (5/5)")
        metaphlan_txt_csv(reads_data, dataset_id)
        verbose_print(f"{datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S')} -- finish converting resualts to CSV (5/5)")
        
        verbose_print("\n")
        verbose_print(f"{datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S')} -- Finished downloading.\n")


