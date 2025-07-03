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


def download(dataset_name, data_type, acc_list, verbose, specific_location,as_single, 
              threads: int = 8, pathways: str = "no"):
    
    verbose_print = print if verbose else lambda *a, **k: None
    as_single= True if as_single else False

    verbose_print("\n")
    verbose_print("download starts.")

    # checking if the acc_list is provided
    acc_list_path = acc_list if acc_list else None

    if acc_list_path is None:
        acc_list_path = get_acc_list(dataset_name, verbose_print)
    else:
        acc_list_path=get_project_list(dataset_name,acc_list_path, verbose_print)
    visualization(acc_list_path, dataset_name, data_type, verbose_print, specific_location,as_single, 
                   threads=threads, pathways=pathways)


def get_project_list(bio_project_name, acc_list_path, verbose_print):
    # if we are given a list of accession numbers, we will use them and produce the run info file of each sample and its metadata
    with open(acc_list_path, 'r') as f:
        acc_list = f.read().splitlines()
        get_run_info_command = f'esearch -db sra -query {acc_list[0]} | efetch -format runinfo > {bio_project_name}.csv'
        os.system(get_run_info_command)

        # Loop through the rest of the accession numbers and append the results
        for acc in acc_list[1:]:
            get_run_info_command = f'esearch -db sra -query {acc} | efetch -format runinfo >> {bio_project_name}.csv'
            os.system(get_run_info_command)
    verbose_print(f"downloaded the run info at {bio_project_name}.csv")
    return f"{acc_list_path}"


def continue_from(dataset_id,continue_path, data_type, verbose, specific_location):
    verbose_print = print if verbose else lambda *a, **k: None

    verbose_print("\n")
    verbose_print(f"Continue downloading from {continue_path}.")
    
    visualization_continue(dataset_id,continue_path, data_type, verbose_print, specific_location)

    
# This function is used to continue the download from a specific path
def continue_from_fastq(dataset_id, continue_path, data_type, verbose, specific_location, 
                        run_humann: bool = False, threads: int = 8, pathways: str = "no"):
    verbose_print = print if verbose else lambda *a, **k: None
    verbose_print("\n")
    verbose_print(f"Continue downloading from {continue_path}.")
    visualization_continue_fastq(dataset_id,continue_path, data_type, verbose_print, specific_location, 
                                  threads=threads, pathways=pathways)


# This function is used to download the qiita data
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