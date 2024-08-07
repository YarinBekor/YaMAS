from .create_visualization import visualization
from .create_visualization import visualization_continue
from .create_visualization import visualization_continue_fastq
from .qiita_visualization import qiita_visualization
from .fastq_visualization import fastq_visualization

import os


# This function downloads the accession list for the specified project (Using the SRA databse, which holds all the necessery metadata.)
def get_acc_list(bio_project_name, verbose_print):
    get_run_info_command = f'esearch -db sra -query {bio_project_name} | efetch -format runinfo > {bio_project_name}_run_info.csv'
    os.system(get_run_info_command)
    verbose_print(f"downloaded the run info at {bio_project_name}_run_info.csv")

    get_acc_info_command = f"cat {bio_project_name}_run_info.csv | cut -f 1 -d ',' | grep -e ERR -e SRR > {bio_project_name}_acc_info.txt"
    os.system(get_acc_info_command)
    verbose_print(f"downloaded the accession list at {bio_project_name}_acc_info.txt")

    return f"{bio_project_name}_acc_info.txt"


def download(dataset_name, data_type, acc_list, verbose, specific_location):
    verbose_print = print if verbose else lambda *a, **k: None

    verbose_print("\n")
    verbose_print("download starts.")

    # checking if the acc_list is provided
    acc_list_path= acc_list if acc_list else None
    acc_list_path = f"{acc_list_path}"
    if acc_list_path is None:
        acc_list_path = get_acc_list(dataset_name, verbose_print)
    visualization(acc_list_path, dataset_name, data_type, verbose_print, specific_location)


def continue_from(dataset_id,continue_path, data_type, verbose, specific_location):
    verbose_print = print if verbose else lambda *a, **k: None

    verbose_print("\n")
    verbose_print(f"Continue downloading from {continue_path}.")
    visualization_continue(dataset_id,continue_path, data_type, verbose_print, specific_location)


def continue_from_fastq(dataset_id,continue_path, data_type, verbose, specific_location):
    verbose_print = print if verbose else lambda *a, **k: None

    verbose_print("\n")
    verbose_print(f"Continue downloading from {continue_path}.")
    visualization_continue_fastq(dataset_id,continue_path, data_type, verbose_print, specific_location)


def download_qiita(fastq_path,metadata_path,data_type, verbose):
    verbose_print = print if verbose else lambda *a, **k: None

    verbose_print("\n")
    verbose_print("download starts.")

    qiita_visualization(fastq_path,metadata_path,data_type, verbose_print)


def download_fastq(fastq_path,barcode_path,metadata_path,data_type, verbose):
    verbose_print = print if verbose else lambda *a, **k: None

    verbose_print("\n")
    verbose_print("download starts.")

    fastq_visualization(fastq_path,barcode_path, metadata_path, data_type, verbose_print)