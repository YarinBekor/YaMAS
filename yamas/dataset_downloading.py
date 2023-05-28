from .create_visualization import visualization
import os


def get_acc_list(bio_project_name, verbose_print):
    get_run_info_command = f'esearch -db sra -query {bio_project_name} | efetch -format runinfo > {bio_project_name}_run_info.csv'
    os.system(get_run_info_command)
    verbose_print(f"downloaded the run info at {bio_project_name}_run_info.csv")

    get_acc_info_command = f"cat {bio_project_name}_run_info.csv | cut -f 1 -d ',' | grep -e ERR -e SRR > {bio_project_name}_acc_info.txt"
    os.system(get_acc_info_command)
    verbose_print(f"downloaded the accession list at {bio_project_name}_acc_info.txt")

    return f"{bio_project_name}_acc_info.txt"


def download(dataset_name,data_type, verbose,specific_location):
    verbose_print = print if verbose else lambda *a, **k: None

    verbose_print("\n")
    verbose_print("download starts.")

    acc_list_path = get_acc_list(dataset_name, verbose_print)
    visualization(acc_list_path, dataset_name,data_type, verbose_print,specific_location)
