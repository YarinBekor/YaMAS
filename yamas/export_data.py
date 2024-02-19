from __future__ import annotations

import csv
import datetime
import os
import pickle
import re

import pandas as pd
from biom import Table, load_table
import qiime2 as q2
from skbio import TreeNode
from Bio import Phylo

from .utilities import ReadsData, run_cmd, download_classifier_url, check_conda_qiime2

nodes_names = []


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


# TODO add to download 18S
def classifier_exists(classifier_path: str):
    if not (os.path.exists(classifier_path) and os.path.isfile(classifier_path)):
        raise FileNotFoundError("Classifier not found! Please give the right path to the classifier.\n"
                                f"Download it from: 16S: {download_classifier_url()}\n")


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
                  "--verbose"
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


def assign_taxonomy(reads_data: ReadsData, data_type, classifier_path: str):
    if data_type == '16S':
        qza_path = lambda filename: os.path.join(reads_data.dir_path, "qza", filename)
        command = [
            "qiime", "feature-classifier", "classify-sklearn",
            "--i-reads", qza_path("rep-seqs-dn-99.qza"),
            "--i-classifier", classifier_path,
            "--o-classification", qza_path("gg-13-8-99-nb-classified.qza")
        ]
        run_cmd(command)

    if data_type == '18S':
        qza_path = lambda filename: os.path.join(reads_data.dir_path, "qza", filename)
        command = [
            "qiime", "feature-classifier", "classify-sklearn",
            "--i-reads", qza_path("rep-seqs-dn-99.qza"),
            "--i-classifier", classifier_path,
            "--o-classification", qza_path("silva-132-99-nb-classifier.qza")
        ]
        run_cmd(command)


def clean_taxonomy1(reads_data: ReadsData, data_type):
    if data_type == '16S':
        qza_path = lambda filename: os.path.join(reads_data.dir_path, "qza", filename)
        command = [
            "qiime", "taxa", "filter-table",
            "--i-table", qza_path("table-dn-99.qza"),
            "--i-taxonomy", qza_path("gg-13-8-99-nb-classified.qza"),
            "--p-exclude", "mitochondria,chloroplast",
            "--o-filtered-table", qza_path("clean_table.qza")
        ]
        run_cmd(command)

    if data_type == '18':
        qza_path = lambda filename: os.path.join(reads_data.dir_path, "qza", filename)
        command = [
            "qiime", "taxa", "filter-table",
            "--i-table", qza_path("table-dn-99.qza"),
            "--i-taxonomy", qza_path("silva-132-99-nb-classifier.qza"),
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


def export_taxonomy(reads_data: ReadsData, data_type, classifier_file_path):
    output_file = os.path.join(reads_data.dir_path, "exports", "tax.tsv")
    # export
    if data_type == '16S':
        command = [
            "qiime", "tools", "export",
            "--input-path", os.path.join(reads_data.dir_path, "qza", "gg-13-8-99-nb-classified.qza"),
            "--output-path", output_file
        ]
        run_cmd(command)
    if data_type == '18S':
        command = [
            "qiime", "tools", "export",
            "--input-path", os.path.join(reads_data.dir_path, "qza", "silva-132-99-nb-classifier.qza"),
            "--output-path", output_file
        ]
        run_cmd(command)


def export_phylogeny(reads_data: ReadsData):
    # sequence alignment using mafft
    input_file_path = os.path.join(reads_data.dir_path, 'qza', 'rep-seqs-dn-99.qza')
    output_file_path = os.path.join(reads_data.dir_path, 'qza', 'aligned-rep-seqs.qza')
    command = ["qiime", "alignment", "mafft", "--i-sequences", input_file_path, "--o-alignment", output_file_path]
    run_cmd(command)

    # Construct a phylogeny using fastree:
    input_file_path = output_file_path
    output_file_path = os.path.join(reads_data.dir_path, 'exports', 'fasttree-tree.qza')

    command = ["qiime", "phylogeny", "fasttree", "--i-alignment", input_file_path, "--o-tree", output_file_path,
               "--verbose"]
    run_cmd(command)

    # Root the phylogeny:
    input_file_path = output_file_path
    output_file_path = os.path.join(reads_data.dir_path, 'exports', "fasttree-tree-rooted.qza")

    command = ["qiime", "phylogeny", "midpoint-root", "--i-tree", input_file_path, "--o-rooted-tree", output_file_path]
    run_cmd(command)


def export_tree(reads_data: ReadsData):
    # convert the rooted tree to newick format
    input_file_path = os.path.join(reads_data.dir_path, 'exports', 'fasttree-tree-rooted.qza')
    output_file_path = os.path.join(reads_data.dir_path, 'exports', 'tree.nwk')
    tree = q2.Artifact.load(input_file_path)
    tntree = tree.view(TreeNode)
    tntree.prune()
    tntree.write(output_file_path)


def convert_to_csv(reads_data: ReadsData):
    # convert otu.tsv to otu.csv and taxonomy.tsv to taxonomy.csv
    otu_tsv = os.path.join(reads_data.dir_path, 'exports', 'otu.tsv')
    otu_csv = os.path.splitext(otu_tsv)[0] + '.csv'
    with open(otu_tsv, 'r', newline='') as tsv_in_file, open(otu_csv, 'w', newline='') as csv_out_file:
        csv_writer = csv.writer(csv_out_file)
        tsv_reader = csv.reader(tsv_in_file, delimiter='\t')
        # skipping the first line-not relevant
        next(tsv_reader)
        for row in tsv_reader:
            csv_writer.writerow(row)

    tax_tsv = os.path.join(reads_data.dir_path, 'exports', 'tax.tsv', 'taxonomy.tsv')
    tax_csv = os.path.splitext(tax_tsv)[0] + '.csv'
    with open(tax_tsv, 'r', newline='') as tsv_in_file, open(tax_csv, 'w', newline='') as csv_out_file:
        csv_writer = csv.writer(csv_out_file)
        tsv_reader = csv.reader(tsv_in_file, delimiter='\t')
        for row in tsv_reader:
            csv_writer.writerow(row)


def export_otu_padding_for_tree(reads_data: ReadsData):
    tree_file = os.path.join(reads_data.dir_path, 'exports', 'tree.nwk')
    otu_path = os.path.join(reads_data.dir_path, 'exports', 'otu.csv')
    otu_padding_path = os.path.join(reads_data.dir_path, 'exports', 'otu_padding.csv')
    otu = pd.read_csv(otu_path)
    #getting the first column of the otu.csv- the ASV list
    asv_list = otu['#OTU ID'].tolist()
    tree = Phylo.read(tree_file, "newick")

    append_nodes_names(tree.root)

    # getting all the asv that are in the otu but not in the tree
    not_in_tree = [i for i in asv_list if i not in nodes_names]
    if not not_in_tree:
        print("All the ASVs in the OTU are in the tree.nwk")
    else:
        print("Some ASVs in the OTU are not in the tree.nwk\n The following ASVs are not in the tree:\n ", not_in_tree)

    # getting all the nodes that are in the tree but not in the ASV
    in_tree_not_asv = [i for i in nodes_names if i not in asv_list]
    print(f"There is {len(in_tree_not_asv)} ASVs that in tree.nwk but not in otu.csv.\n Starting creating otu_padding.csv.")

    # making a copy of the otu.csv and adding the ASVs that are in the tree but not in the otu with 0 values
    with open(otu_path, 'r', newline='') as otu_original, open(otu_padding_path, 'w', newline='') as otu_pad:
        reader = csv.reader(otu_original)
        writer = csv.writer(otu_pad)
        for row in reader:
            writer.writerow(row)

    with open(otu_padding_path, 'a', newline='') as otu_pad:
        writer = csv.writer(otu_pad)
        for i in in_tree_not_asv:
            new_row = [i] + [0] * len(otu.columns[1:])
            writer.writerow(new_row)

    print("otu_padding.csv created successfully.")


def append_nodes_names(clade):
    # Function to recursively append names of terminal nodes to the list from the tree.nwk
    if clade.is_terminal():
        nodes_names.append(clade.name)
    else:
        for sub_clade in clade.clades:
            append_nodes_names(sub_clade)


def export(output_dir: str, data_type, trim, trunc, classifier_file_path: str, threads: int = 12):
    print("\n")
    print(datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S'))
    print(f"### Exporting {data_type} ###")
    print("Starting OTU & TAXONOMY files extraction")

    check_conda_qiime2()
    reads_data: ReadsData = pickle.load(open(os.path.join(output_dir, "reads_data.pkl"), "rb"))
    # trim_trunc_check(reads_data, trim, trunc)
    classifier_exists(classifier_file_path)

    paired = reads_data.rev and reads_data.fwd
    output_path = os.path.join(reads_data.dir_path, "qza", f"demux-{'paired' if paired else 'single'}-end.qza")

    print("\n")
    print(f"{datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S')} -- Start dada2 (1/8)")
    qiime_dada2(reads_data, output_path, left=trim, right=trunc, threads=threads)
    print(f"{datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S')} -- Finish dada2 (1/8)")

    print("\n")
    print(f"{datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S')} -- Start clustering features (2/8)")
    cluster_features(reads_data)
    print(f"{datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S')} -- Finish clustering features (2/8)")

    print("\n")
    print(f"{datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S')} -- Start assigning taxonomy (3/8)")
    assign_taxonomy(reads_data, data_type, classifier_file_path)
    print(f"{datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S')} -- Finish assigning taxonomy (3/8)")

    run_cmd(["mkdir", os.path.join(reads_data.dir_path, "exports")])

    print("\n")
    print(f"{datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S')} -- Start cleaning taxonomy (4/8)")
    clean_taxonomy1(reads_data, data_type)
    clean_taxonomy2(reads_data)
    print(f"{datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S')} -- Finish cleaning taxonomy (4/8)")

    print(f"{datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S')} -- Start exporting OTU (5/8)")
    export_otu(reads_data)
    print(f"{datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S')} -- Finish exporting OTU (5/8)")

    print(f"{datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S')} -- Start exporting taxonomy (6/8)")
    export_taxonomy(reads_data, data_type, classifier_file_path)
    print(f"{datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S')} -- Finish exporting taxonomy (6/8)")

    print(f"{datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S')} -- Start exporting phylogeny (7/8)")
    export_phylogeny(reads_data)
    print(f"{datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S')} -- Finish exporting phylogeny (7/8)")

    print(f"{datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S')} -- Start exporting tree (7/8)")
    export_tree(reads_data)
    convert_to_csv(reads_data)
    export_otu_padding_for_tree(reads_data)
    print(f"{datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S')} -- Finish exporting tree (7/8)")
