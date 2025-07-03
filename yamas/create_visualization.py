from __future__ import annotations
import os
import csv
import os.path
import pickle
import datetime
from tqdm import tqdm
from metaphlan.utils.merge_metaphlan_tables import merge
from .utilities import run_cmd, ReadsData, check_conda_qiime2
import json
import shutil

from .generate_pathways import run_humann_pipeline
from pathlib import Path

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
             "--max-size", "100G"])
    repo_root = Path(os.environ.get("NCBI_VDB_REPOSITORY_ROOT",
                                    Path.home() / "ncbi"))
    src_dir   = repo_root / "public" / "sra"
    dst_dir   = Path(dir_path) / "sra"
    dst_dir.mkdir(parents=True, exist_ok=True)

    for sra_file in src_dir.glob("*.sra"):
        shutil.move(sra_file, dst_dir / sra_file.name)
        print(f"moved {sra_file} → {dst_dir/sra_file.name}")


def sra_to_fastq(dir_path: str, as_single):
    print(f"converting files from .sra to .fastq.")
    for sra_dir in tqdm(os.listdir(os.path.join(dir_path, "sra")), desc="converted files"):
        sra_file = os.listdir(os.path.join(dir_path, "sra", sra_dir))[0]
        sra_path = os.path.join(dir_path, "sra", sra_dir, sra_file)
        fastq_path = os.path.join(dir_path, "fastq")
        run_cmd(["fasterq-dump", "--split-files", sra_path, "-O", fastq_path])

    # check if reads include fwd and rev
    fastqs = sorted(os.listdir(os.path.join(dir_path, "fastq")))[:3]
    if len(set([fastq.split("_")[0] for fastq in fastqs])) == 2:
        if as_single:
            delete__2_files = f"rm {os.path.join(dir_path, '*_2.fastq')}"
            run_cmd([delete__2_files])
            print("Single reads- only forward reads are kept, reverse reads are deleted.")
            return ReadsData(dir_path, fwd=True, rev=False)
        else:
            return ReadsData(dir_path, fwd=True, rev=True)


    return ReadsData(dir_path, fwd=True, rev=False)


def create_manifest(reads_data: ReadsData):
    base_dir = os.path.abspath(reads_data.dir_path)
    fastq_path = os.path.join(base_dir, "fastq")

    #not paired reads
    if not reads_data.rev:
        files = [
            os.path.join(fastq_path, f) 
            for f in os.listdir(fastq_path)
                 if os.path.isfile(os.path.join(fastq_path, f))
                 ]
        names = [os.path.splitext(os.path.basename(f))[0] for f in files]
        
        manifest_path = os.path.join(base_dir, 'manifest.tsv')
        with open(manifest_path, 'w', newline='') as manifest:
            tsv_writer = csv.writer(manifest, delimiter='\t')
            tsv_writer.writerow(["SampleID", "absolute-filepath"])
            for n, f in zip(*(names, files)):
                abs_fp = os.path.abspath(f)
                tsv_writer.writerow([n, abs_fp])
        return
    #paired-end reads
    files_fwd = sorted([
        os.path.join(fastq_path, f)
        for f in os.listdir(fastq_path)
        if os.path.isfile(os.path.join(fastq_path, f)) and "_1" in f])
    files_rev = sorted([
        os.path.join(fastq_path, f)
        for f in os.listdir(fastq_path)
        if os.path.isfile(os.path.join(fastq_path, f)) and "_2" in f])
    
    names = sorted([
        os.path.splitext(os.path.basename(f))[0]
        for f in os.listdir(os.path.join(base_dir, "sra"))])
    
    manifest_path = os.path.join(base_dir, 'manifest.tsv')

    with open(manifest_path, 'w', newline='') as manifest:
        tsv_writer = csv.writer(manifest, delimiter='\t')
        tsv_writer.writerow([
            "SampleID",
            "forward-absolute-filepath",
            "reverse-absolute-filepath"
        ])
        for n, ff, fr in zip(names, files_fwd, files_rev):
            abs_ff = os.path.abspath(ff)
            abs_fr = os.path.abspath(fr)
            tsv_writer.writerow([n, abs_ff, abs_fr])
    
    # remove sra folder
    shutil.rmtree(os.path.join(base_dir, "sra"))



def qiime_import(reads_data: ReadsData):
    qza_path = os.path.join(reads_data.dir_path, "qza")
    paired = reads_data.rev and reads_data.fwd
    if paired: print("Paired reads")

    qza_file_path = os.path.join(qza_path, f"demux-{'paired' if paired else 'single'}-end.qza")
    command = [
        "qiime", "tools", "import",
        "--type", f"SampleData[{'PairedEndSequencesWithQuality' if paired else 'SequencesWithQuality'}]",
        "--input-path", f"{os.path.join(reads_data.dir_path, 'manifest.tsv')}",
        "--input-format", "PairedEndFastqManifestPhred33V2" if paired else "SingleEndFastqManifestPhred33V2",
        "--output-path", qza_file_path,

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

def get_files_in_directory(directory, extension=""):
    return [os.path.join(directory, f) for f in os.listdir(directory) if f.endswith(extension)]


def metaphlan_extraction(reads_data, dataset_id):
    paired = reads_data.rev and reads_data.fwd
    fastq_path = os.path.join(reads_data.dir_path, "fastq")
    export_path = os.path.join(reads_data.dir_path, "export")
    run_cmd([f"mkdir {export_path}"])
    final_output_path = os.path.join(export_path, f'{dataset_id}_final.txt')
    run_cmd([f"touch {final_output_path}"])
    fastq_files = [a for a in os.listdir(fastq_path) if a.split(".")[-1] == "fastq"]

    if paired:
        print("paired")
        for i in tqdm(range(0, len(fastq_files), 2)):
            fastq_name = fastq_files[i].split('_')[0]
            fastq_1 = os.path.join(fastq_path, fastq_files[i])
            fastq_2 = os.path.join(fastq_path, fastq_files[i + 1])
            output = os.path.join(fastq_path, f"{fastq_name}.bowtie2.bz2")
            command = [f"metaphlan {fastq_1},{fastq_2} --input_type fastq --bowtie2out {output} --nproc 24"]
            run_cmd(command)


            final_output_file = os.path.join(os.path.join(reads_data.dir_path, 'qza'), f'{fastq_name}_profile.txt')
            command = [f"metaphlan {output} --input_type bowtie2out --nproc 24 > {final_output_file}"]
            run_cmd(command)
            # after converting the fastq to bowtie and bowtie to profile we can delete these files
            run_cmd([f"rm {fastq_1} {fastq_2} {output}"])
    else:
        print("not paired")
        for fastq in tqdm(fastq_files):
            output = os.path.join(os.path.join(reads_data.dir_path, 'qza'), f'{fastq}_profile.txt')
            fastq = os.path.join(fastq_path, fastq)
            command = [f"metaphlan {fastq} --input_type fastq --nproc 24 > {output}"]
            run_cmd(command)
            # after converting the fastq to profile we can delete the fastq files
            run_cmd([f"rm {fastq}"])

    qza_dir= os.path.join(reads_data.dir_path, 'qza')
    # Gather all profile files from the directory
    profile_files = get_files_in_directory(qza_dir, extension="_profile.txt")

    # Merge the profile files
    merge(profile_files, open(final_output_path, 'w'), gtdb=False)
    # delete the fastq dir, we convert all the fastq to profile
    shutil.rmtree(fastq_path)

def metaphlan_txt_csv(reads_data, dataset_id):
    export_path = os.path.join(reads_data.dir_path, "export")
    input_file = os.path.join(export_path, f"{dataset_id}_final.txt")
    output_file = os.path.join(export_path, f"{dataset_id}_final_table.csv")
    with open(input_file, 'r') as txt_file:
        lines = txt_file.readlines()

    # Extract data  the text file
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
def visualization(acc_list, dataset_id, data_type, verbose_print, specific_location,as_single, 
                  run_humann: bool = False, threads: int = 8, pathways: str = "no"):
    
    verbose_print("\n")
    verbose_print(datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S'))
    data_json = {}
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
    dir_path = create_dir(dir_name, specific_location)
    verbose_print("Find ALL NEW data in:", dir_path)
    json_file_path = f"{dir_path}/metadata.json"

    # This two stages are for all data type. Prefetching the data from the relevant database, and converting it to fatsqs.
    verbose_print("\n")
    verbose_print(f"{datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S')} -- Start prefetch (1/6)")
    download_data_from_sra(dir_path, acc_list)

    data_json["dir_path"] = dir_path
    data_json["dataset_id"] = dataset_id
    with open(json_file_path, "w") as json_file:
        json.dump(data_json, json_file)
    verbose_print(f"{datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S')} -- Finish prefetch (1/6)")

    verbose_print("\n")
    verbose_print(f"{datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S')} -- Start conversion (2/6)")
    reads_data = sra_to_fastq(dir_path, as_single)
    verbose_print(f"{datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S')} -- Finish conversion (2/6)")

    verbose_print(f"{datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S')} -- Start creating metadata.json (3/6)")

    # Store dir_path and reads_data in the data_json dictionary
    data_json["type"] = data_type
    data_json["read_data_fwd"] = reads_data.fwd
    data_json["read_data_rev"] = reads_data.rev

    with open(json_file_path, "w") as json_file:
        json.dump(data_json, json_file)

    verbose_print(f"{datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S')} -- Finish creating metadata.json (3/6)")

    if data_type == '16S' or data_type == '18S':

        verbose_print("\n")
        verbose_print(f"{datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S')} -- Start creating manifest (4/6)")
        create_manifest(reads_data)
        verbose_print(f"{datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S')} -- Finish creating manifest (4/6)")

        verbose_print("\n")
        verbose_print(f"{datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S')} -- Start 'qiime import' (5/6)")
        qza_file_path = qiime_import(reads_data)
        verbose_print(f"{datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S')} -- Finish 'qiime import' (5/6)")

        verbose_print("\n")
        verbose_print(f"{datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S')} -- Start 'qiime demux' (6/6)")
        vis_file_path = qiime_demux(reads_data, qza_file_path, dataset_id)
        verbose_print(f"{datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S')} -- Finish 'qiime demux' (6/6)")

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
    
    else:  # shotgun
        verbose_print("\n")
        verbose_print(f"{datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S')} -- Start metaphlan extraction (4/6)")
        metaphlan_extraction(reads_data, dataset_id)
        verbose_print(f"{datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S')} -- Finish metaphlan extraction (4/6)")

        verbose_print("\n")
        verbose_print(
            f"{datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S')} -- Start converting resualts to CSV (5/6)")
        metaphlan_txt_csv(reads_data, dataset_id)
        verbose_print(
            f"{datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S')} -- finish converting resualts to CSV (5/6)")

        verbose_print("\n")
        verbose_print(f"{datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S')} -- Finished downloading. 6/6 \n")
        
        if pathways == "yes":  # If pathways is set to "yes", run the HUMAnN pipeline
            fastq_dir   = Path(dir_path) / "fastq"
            fastqs      = sorted(fastq_dir.glob("*.fastq"))
            fastq_input = ",".join(map(str, fastqs))          # "f1.fastq,f2.fastq,…"
            output_dir  = Path(dir_path) / "humann_results"

            meta_profile = Path(dir_path) / "export" / f"{dataset_id}_final.txt"
            
            run_humann_pipeline(fastq_input, output_dir, meta_profile=meta_profile)


def visualization_continue_fastq(dataset_id, continue_path, data_type, verbose_print, specific_location, 
                                  threads: int = 8, pathways: str = "no"):
    continue_path = Path(continue_path)
    verbose_print("\n")
    verbose_print('Checking environment...', end=" ")
    check_conda_qiime2()
    verbose_print('Done.')

    verbose_print("\n")
    data_json = {}
    verbose_print(f"{datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S')} -- Start conversion (1/5)")
    reads_data = sra_to_fastq(continue_path, as_single=False)
    verbose_print(f"{datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S')} -- Finish conversion (1/5)")

    verbose_print(f"{datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S')} -- Start creating metadata.json (2/5)")

    verbose_print(f"{datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S')} -- Finish creating metadata.json (2/5)")

    if data_type == '16S' or data_type == '18S':

        verbose_print("\n")
        verbose_print(f"{datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S')} -- Start creating manifest (3/5)")
        create_manifest(reads_data)
        verbose_print(f"{datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S')} -- Finish creating manifest (3/5)")
        sra_dir = continue_path / "sra"
        if sra_dir.exists():          # only delete if it still exists
            shutil.rmtree(sra_dir)
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
    
    else:  # shotgun
        verbose_print("\n")
        verbose_print(f"{datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S')} -- Start metaphlan extraction (3/5)")
        metaphlan_extraction(reads_data, dataset_id)
        verbose_print(f"{datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S')} -- Finish metaphlan extraction (3/5)")

        verbose_print("\n")
        verbose_print(
            f"{datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S')} -- Start converting resualts to CSV (4/5)")
        metaphlan_txt_csv(reads_data, dataset_id)
        verbose_print(f"{datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S')} -- finish converting resualts to CSV (4/5)")

        if pathways == "yes":  # If pathways is set to "yes", run the HUMAnN pipeline
            fastq_dir   = continue_path / "fastq"
            fastqs      = sorted(fastq_dir.glob("*.fastq"))
            fastq_input = ",".join(map(str, fastqs))
            output_dir  = continue_path / "humann_results"

            meta_profile = Path(continue_path) / "export" / f"{dataset_id}_final.txt"
            
            run_humann_pipeline(fastq_input, output_dir, meta_profile=meta_profile)
                
        verbose_print("\n")
        verbose_print(f"{datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S')} -- Finished downloading. 5/5 \n")


def visualization_continue(dataset_id, continue_path, data_type, verbose_print, specific_location):
    verbose_print("\n")
    verbose_print('Checking environment...', end=" ")
    check_conda_qiime2()
    verbose_print('Done.')

    json_file_path = f"{continue_path}/metadata.json"
    try:
        with open(json_file_path, 'r') as json_file:
            data_json = json.load(json_file)
            reads_data_fwd = data_json.get("read_data_fwd")
            reads_data_rev = data_json.get("read_data_rev")
            if reads_data_fwd == "true":
                print("yes")
            print(reads_data_fwd, reads_data_rev)
            reads_data = ReadsData(continue_path, fwd=reads_data_fwd, rev=reads_data_rev)

    except FileNotFoundError:
        print(
            f"Error: Metadata file not found at {json_file_path}. Please check the path and try again, or try --download command to start a new download of this data.")
        return

    if data_type == '16S' or data_type == '18S':

        verbose_print("\n")
        verbose_print(f"{datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S')} -- Start creating manifest (1/3)")
        create_manifest(reads_data)
        verbose_print(f"{datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S')} -- Finish creating manifest (1/3)")
        sra_path = os.path.join(continue_path, "sra")
        if os.path.isdir(sra_path):
            shutil.rmtree(sra_path)        
        verbose_print("\n")
        verbose_print(f"{datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S')} -- Start 'qiime import' (2/3)")
        qza_file_path = qiime_import(reads_data)
        verbose_print(f"{datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S')} -- Finish 'qiime import' (2/3)")

        verbose_print("\n")
        verbose_print(f"{datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S')} -- Start 'qiime demux' (3/3)")
        vis_file_path = qiime_demux(reads_data, qza_file_path, dataset_id)
        verbose_print(f"{datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S')} -- Finish 'qiime demux' (3/3)")

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
        verbose_print(f"{datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S')} -- Start metaphlan extraction (1/2)")
        metaphlan_extraction(reads_data, dataset_id)
        verbose_print(f"{datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S')} -- Finish metaphlan extraction (1/2)")

        verbose_print("\n")
        verbose_print(
            f"{datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S')} -- Start converting resualts to CSV (2/2)")
        metaphlan_txt_csv(reads_data, dataset_id)
        verbose_print(
            f"{datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S')} -- finish converting resualts to CSV (2/2)")

        verbose_print("\n")
        verbose_print(f"{datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S')} -- Finished downloading.\n")
