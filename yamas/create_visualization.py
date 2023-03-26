from __future__ import annotations

import csv
import os.path
import pickle
import datetime
from tqdm import tqdm

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


def create_dir(dir_name, specific_location):
    dir_path = os.path.join(os.path.abspath(specific_location), dir_name)

    os.system(f"mkdir {dir_name}")
    print(f"{dir_name} created.")

    os.system(f"mkdir {os.path.join(dir_path, 'sra')}")
    print(f"{dir_name}/sra created.")

    os.system(f"mkdir {os.path.join(dir_path, 'fastq')}")
    print(f"{dir_name}/fastq created.")

    os.system(f"mkdir {os.path.join(dir_path, 'qza')}")
    print(f"{dir_name}/qza created.")

    os.system(f"mkdir {os.path.join(dir_path, 'vis')}")
    print(f"{dir_name}/vis created.")

    return dir_path

def download_data_from_sra(dir_path: str, acc_list: str = ""):
    run_cmd(['prefetch',
             "--option-file", acc_list,
             "--output-directory", os.path.join(dir_path, "sra"),
             "--max-size", "u"])


def sra_to_fastq(dir_path: str):
    print(f"converting files from .sra to .fastq.")
    for sra_dir in tqdm(os.listdir(os.path.join(dir_path, "sra")), desc="converted files"):
        sra_file = os.listdir(os.path.join(dir_path, "sra", sra_dir))[0]
        sra_path = os.path.join(dir_path, "sra", sra_dir, sra_file)
        fastq_path = os.path.join(dir_path, "fastq")
        run_cmd(["fasterq-dump", "--split-files", sra_path, "-O", fastq_path])

    # check if reads include fwd and rev
    fastqs = sorted(os.listdir(os.path.join(dir_path, "fastq")))[:3]
    if len(set([fastq.split("_")[0] for fastq in fastqs])) == 1:
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

def visualization(acc_list, dataset_id, verbose_print, specific_location):
    verbose_print("\n")
    verbose_print(datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S'))
    dir_name = f"{dataset_id}-{datetime.datetime.now().strftime('%d-%m-%Y_%H-%M-%S')}"
    verbose_print("Starting conversion to VIS file")

    verbose_print("\n")
    verbose_print('Checking environment...', end=" ")
    check_conda_qiime2()
    verbose_print('Done.')
    verbose_print("Checking inputs:")
    check_input(acc_list)

    verbose_print("\n")
    verbose_print("Creating a new directory for this dataset import:'", dir_name, "'")
    dir_path = create_dir(dir_name,specific_location)
    verbose_print("Find ALL NEW data in:", dir_path)

    verbose_print("\n")
    verbose_print(f"{datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S')} -- Start prefetch (1/5)")
    download_data_from_sra(dir_path, acc_list)
    verbose_print(f"{datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S')} -- Finish prefetch (1/5)")

    verbose_print("\n")
    verbose_print(f"{datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S')} -- Start conversion (2/5)")
    reads_data = sra_to_fastq(dir_path)
    verbose_print(f"{datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S')} -- Finish conversion (2/5)")

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
