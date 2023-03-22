from __future__ import annotations

import datetime
import os
import pickle

from .utilities import ReadsData, run_cmd, download_classifier_url, check_conda_qiime2


def trim_trunc_check(reads_data: ReadsData, trim: int | tuple[int, int], trunc: int | tuple[int, int]):
    if reads_data.fwd and reads_data.rev:
        if not isinstance(trim, tuple) or not isinstance(trunc, tuple):
            raise TypeError("The read consist of both forward and reverse, "
                            "so 'trim' and 'trunc' must be tuples of integers.\n"
                            f"Got {type(trim)}, {type(trunc)}.")

        if len(trim) != 2 or len(trunc) != 2:
            raise ValueError("The read consist of both forward and reverse, "
                             "so 'trim' and 'trunc' must be tuples of 2 integers.\n"
                             f"Got tuples of length {len(trim)}, {len(trunc)}.")
    if not isinstance(trim, int) or not isinstance(trunc, int):
        raise TypeError("The read consist of only forward, "
                        "so 'trim' and 'trunc' must be integers.\n"
                        f"Got {type(trim)}, {type(trunc)}.")


def classifier_exists(classifier_path: str):
    if not (os.path.exists(classifier_path) and os.path.isfile(classifier_path)):
        raise FileNotFoundError("Classifier not found! Please give the right path to the classifier.\n"
                                f"Download it from: {download_classifier_url()}")


def qiime_dada2(reads_data: ReadsData, input_path: str,
                left: int | tuple[int, int], right: int | tuple[int, int], threads: int = 12):
    paired = reads_data.fwd and reads_data.rev

    trim_range = ["--p-trim-left-f", str(left[0]), "--p-trim-left-r", str(left[1])] if paired \
        else ["--p-trim-left", str(left)]
    trunc_range = ["--p-trunc-len-f", str(right[0]), "--p-trunc-len-r", str(right[1])] if paired \
        else ["--p-trunc-len", str(right)]

    command = [
                  "qiime", "dada2", "denoise-paired" if paired else "denoise-single",
                  "--i-demultiplexed-seqs", input_path,
              ] + trim_range + trunc_range + [
                  "--o-table", os.path.join(reads_data.dir_path, "qza", "dada2_table.qza"),
                  "--p-n-threads", str(threads),
                  "--p-chimera-method", "consensus",
                  "--o-representative-sequences", os.path.join(reads_data.dir_path, "qza", "dada2_rep-seqs.qza"),
                  "--o-denoising-stats", os.path.join(reads_data.dir_path, "qza", "dada2_denoising-stats.qza"),
              ]
    run_cmd(command)


def cluster_features(reads_data: ReadsData):
    qza_path = lambda filename: os.path.join(reads_data.dir_path, "qza", filename)
    command = [
        "qiime", "vsearch", "cluster-features-de-novo",
        "--i-table", qza_path("dada2_table.qza"),
        "--i-sequences", qza_path("dada2_rep-seqs.qza"),
        "--p-perc-identity", "0.99",
        "--o-clustered-table", qza_path("table-dn-99.qza"),
        "--o-clustered-sequences", qza_path("rep-seqs-dn-99.qza")
    ]
    run_cmd(command)


def assign_taxonomy(reads_data: ReadsData, classifier_path: str):
    qza_path = lambda filename: os.path.join(reads_data.dir_path, "qza", filename)
    command = [
        "qiime", "feature-classifier", "classify-sklearn",
        "--i-reads", qza_path("rep-seqs-dn-99.qza"),
        "--i-classifier", classifier_path,
        "--o-classification", qza_path("gg-13-8-99-nb-classified.qza")
    ]
    run_cmd(command)


def clean_taxonomy1(reads_data: ReadsData):
    qza_path = lambda filename: os.path.join(reads_data.dir_path, "qza", filename)
    command = [
        "qiime", "taxa", "filter-table",
        "--i-table", qza_path("table-dn-99.qza"),
        "--i-taxonomy", qza_path("gg-13-8-99-nb-classified.qza"),
        "--p-exclude", "mitochondria,chloroplast",
        "--o-filtered-table", qza_path("clean_table.qza")
    ]
    run_cmd(command)


def clean_taxonomy2(reads_data: ReadsData):
    qza_path = lambda filename: os.path.join(reads_data.dir_path, "qza", filename)
    command = [
        "qiime", "feature-table", "filter-features",
        "--i-table", qza_path("clean_table.qza"),
        "--p-min-samples", "3",
        "--p-min-frequency", "10",
        "--o-filtered-table", qza_path("feature-frequency-filtered-table.qza")
    ]
    run_cmd(command)


def export_otu(reads_data: ReadsData):
    output_file = os.path.join(reads_data.dir_path, "exports", "otu.tsv")
    # export
    command = [
        "qiime", "tools", "export",
        "--input-path", os.path.join(reads_data.dir_path, "qza", "feature-frequency-filtered-table.qza"),
        "--output-path", os.path.join(reads_data.dir_path, "exports")
    ]
    run_cmd(command)

    # convert
    command = [
        "biom", "convert",
        "-i", os.path.join(reads_data.dir_path, "exports", "feature-table.biom"),
        "-o", output_file,
        "--to-tsv"
    ]
    run_cmd(command)


def export_taxonomy(reads_data: ReadsData, classifier_file_path):
    output_file = os.path.join(reads_data.dir_path, "exports", "tax.tsv")
    # export
    command = [
        "qiime", "tools", "export",
        "--input-path", os.path.join(reads_data.dir_path, "qza", "gg-13-8-99-nb-classified.qza"),
        "--output-path", output_file
    ]
    run_cmd(command)


def export_phylogeny(reads_data: ReadsData):
    # sequence alignment using mafft
    input_file_path = os.path.join(reads_data.dir_path,'qza', 'rep-seqs-dn-99.qza')
    output_file_path = os.path.join(reads_data.dir_path,'qza', 'aligned-rep-seqs.qza')
    command = ["qiime", "alignment", "mafft", "--i-sequences", input_file_path, "--o-alignment", output_file_path]
    run_cmd(command)

    # Construct a phylogeny using fastree:
    input_file_path = output_file_path
    output_file_path = os.path.join(reads_data.dir_path,'exports', 'fasttree-tree.qza')

    command = ["qiime", "phylogeny", "fasttree", "--i-alignment", input_file_path, "--o-tree", output_file_path,
               "--verbose"]
    run_cmd(command)

    # Root the phylogeny:
    input_file_path = output_file_path
    output_file_path = os.path.join(reads_data.dir_path,'exports', "fasttree-tree-rooted.qza")

    command = ["qiime", "phylogeny", "midpoint-root", "--i-tree", input_file_path, "--o-rooted-tree", output_file_path]
    run_cmd(command)


def export(output_dir: str, trim, trunc, classifier_file_path: str, threads: int = 12):
    print("\n")
    print(datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S'))
    print("Starting OTU & TAXONOMY files extraction")

    check_conda_qiime2()
    reads_data: ReadsData = pickle.load(open(os.path.join(output_dir, "reads_data.pkl"), "rb"))
    trim_trunc_check(reads_data, trim, trunc)
    classifier_exists(classifier_file_path)

    paired = reads_data.rev and reads_data.fwd
    output_path = os.path.join(reads_data.dir_path, "qza", f"demux-{'paired' if paired else 'single'}-end.qza")

    print("\n")
    print(f"{datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S')} -- Start dada2 (1/7)")
    qiime_dada2(reads_data, output_path, left=trim, right=trunc, threads=threads)
    print(f"{datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S')} -- Finish dada2 (1/7)")

    print("\n")
    print(f"{datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S')} -- Start clustering features (2/7)")
    cluster_features(reads_data)
    print(f"{datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S')} -- Finish clustering features (2/7)")

    print("\n")
    print(f"{datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S')} -- Start assigning taxonomy (3/7)")
    assign_taxonomy(reads_data, classifier_file_path)
    print(f"{datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S')} -- Finish assigning taxonomy (3/7)")

    run_cmd(["mkdir", os.path.join(reads_data.dir_path, "exports")])

    print("\n")
    print(f"{datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S')} -- Start cleaning taxonomy (4/7)")
    clean_taxonomy1(reads_data)
    clean_taxonomy2(reads_data)
    print(f"{datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S')} -- Finish cleaning taxonomy (4/7)")

    print(f"{datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S')} -- Start exporting OTU (5/7)")
    export_otu(reads_data)
    print(f"{datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S')} -- Finish exporting OTU (5/7)")

    print(f"{datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S')} -- Start exporting taxonomy (6/7)")
    export_taxonomy(reads_data, classifier_file_path)
    print(f"{datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S')} -- Finish exporting taxonomy (6/7)")

    print(f"{datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S')} -- Start exporting phylogeny (7/7)")
    export_phylogeny(reads_data)
    print(f"{datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S')} -- Finish exporting phylogeny (7/7)")
