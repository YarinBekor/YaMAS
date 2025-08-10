import os
from metaphlan import install_metaphlan, install_metaphlan_database, run_metaphlan
import multiprocessing
from findspecies import search_species
from extract_file_KEGG import download_file
from combine_fna import combine_fna_files
import rpy2.robjects as robjects
from reference_data2 import reference_data
from translate_to_K0 import translate_genes_to_K0s
from find_pathways_no_enrich import filter_genes
from more_pathways import fetch_more_gene_pathway_data
from remove_zero import remove_zeros_from_count_matrix
from bacteria_spreadsheet import create_spreadsheet
import argparse

def install_metaphlan(metaphlan_database_folder):
    install_metaphlan()
    install_metaphlan_database(metaphlan_database_folder)


def main():
    # Get the directory where the script is located
    script_directory = os.path.dirname(os.path.abspath(__file__))

    parser = argparse.ArgumentParser(description="Your script description here")

    parser.add_argument("--main-folder", required=True, help="Path to main folder")
    parser.add_argument("--reference-folder", required=True, help="Path to reference folder")
    parser.add_argument("--metaphlan-database", required=True, help="Path to metaphlan database folder")
    parser.add_argument("--fastq-file-1", required=True, help="Path to first fastq file")
    parser.add_argument("--fastq-file-2", required=True, help="Path to second fastq file")
    parser.add_argument("--install-metaphlan", action="store_true", help="Install metaphlan")
    parser.add_argument("--run-metaphlan-from-code", action="store_true", help="Run metaphlan from code")
    parser.add_argument("--run-metaphlan-from-terminal", action="store_true", help="Run metaphlan from terminal")
    parser.add_argument("--metaphlan-output", help="Path to metaphlan output file")

    args = parser.parse_args()

    main_folder_path = args.main_folder
    reference_folder_path = args.reference_folder
    metaphlan_database_folder = args.metaphlan_database
    fastq_file_1 = args.fastq_file_1
    fastq_file_2 = args.fastq_file_2

    dictionary_path = os.path.join(script_directory, "species_codes.json")
    main_csv = os.path.join(script_directory, "species2.csv")
    kegg_species = os.path.join(script_directory, "all_KEGG_species.txt")
    pathways_genes = os.path.join(script_directory, "dictionary.txt")

    if args.install_metaphlan:
        install_metaphlan(metaphlan_database_folder)

    if args.run_metaphlan_from_code:
        command = f"metaphlan {fastq_file_1} --bowtie2out metagenome.bowtie2.bz2 --nproc 2 --input_type fastq -o {main_folder_path}/profiled_metagenome.txt --bowtie2db {metaphlan_database_folder}"
        os.system(command)
        metaphlan_output = f"{main_folder_path}/profiled_metagenome.txt"

    if args.run_metaphlan_from_terminal:
        print("Ran metaphlan from terminal before starting program")
        if args.metaphlan_output:
            metaphlan_output = args.metaphlan_output
        else:
            print("Error: Please provide the --metaphlan-output argument.")

    # Create a list containing the names of the unique species identified
    species_list = search_species(metaphlan_output, main_folder_path, kegg_species)
    species_list.sort()

    # Build a reference file for reading and alignment
    # The FNA (fasta) and GTF file of each species will be downloaded from the NCBI database
    # If they have been downloaded from a previous session, they will not be downloaded again
    for species in species_list:
        species_directory = os.path.join(reference_folder_path, species)
        if not os.path.exists(species_directory):
            # Download the fna/gtf file only if the species folder doesn't exist
            download_file(species, "fna", reference_folder_path, main_csv)
            download_file(species, "gtf", reference_folder_path, main_csv)
        else:
            print(f"Skipping {species}. Folder already exists in the output directory.")

    # The FNA files of each species will be combined into one large annotation file for alignment
    annotation_file = combine_fna_files(species_list, reference_folder_path, main_folder_path)

    # The final steps of the processing are done in R
    r_code = f'''
    # List of packages to check and install
    packages_to_install <- c("Rsubread", "GenomicRanges", "GenomicAlignments", "rtracklayer", "clusterProfiler")

    # Check if BiocManager is installed, if not, install it
    if (!requireNamespace("BiocManager", quietly = TRUE)) {{
        install.packages("BiocManager")
    }}

    # Load the BiocManager package
    library(BiocManager)

    # Loop through the list of packages
    for (pkg in packages_to_install) {{
        # Check if the package is already installed
        if (!requireNamespace(pkg, quietly = TRUE)) {{
        # Install the package from Bioconductor
            BiocManager::install(pkg)
            cat(paste(pkg, "has been installed.\n"))
        }} else {{
            cat(paste(pkg, "is already installed.\n"))
        }}
    }}

    library(Rsubread)
    library(GenomicRanges)
    library(GenomicAlignments)
    library(rtracklayer)
    library(clusterProfiler)

    # Set a folder path on the desktop
    folder_path <- "{main_folder_path}"

    # Set the paths to your input files
    buildindex(basename = "my_index", reference = "{annotation_file}")
    readfile1 <- "{fastq_file_1}"  # Path to the first FASTQ file
    readfile2 <- "{fastq_file_2}"  # Path to the second FASTQ file
    annotation_directory <- "{reference_folder_path}"  # Path to the annotation directory
    dictionary_path <- "{dictionary_path}"  # Path to the json dictionary file

    # Align reads (execute only once)
    align(index = "my_index", readfile1 = readfile1, readfile2 = readfile2, type = "rna", output_file = file.path(folder_path, "aligned_reads.bam"))

    # Load species_list from a .txt file in the folder path
    species_list <- scan(file.path(folder_path, "species_list.txt"), what = "character")

    # Loop through species_list
    for (species_name in species_list) {{
        # Create a subfolder for the current species
        species_folder <- file.path(folder_path, species_name)
        dir.create(species_folder, showWarnings = FALSE)

        # Build the gtf_file path for each species
        gtf_file_subdirectory <- file.path(annotation_directory, species_name)
        gtf_file <- file.path(gtf_file_subdirectory, paste0(species_name, ".gtf"))

        # Load the aligned reads and create gene id list
        gr <- import(gtf_file, format = "gtf")
        bam_file <- file.path(folder_path, "aligned_reads.bam")  # Use the aligned BAM file saved previously
        bam_data <- readGAlignmentPairs(bam_file)
        overlaps <- findOverlaps(bam_data, gr)
        gene_ids <- mcols(gr)$gene_id[subjectHits(overlaps)]
        gene_ids <- gene_ids[!is.na(gene_ids)]

        # Get unique gene_ids
        unique_genes <- unique(gene_ids)

        # Count the occurrences of each gene_id
        gene_counts <- table(gene_ids)

        # Save the unique genes to a file
        unique_genes_output_file <- file.path(species_folder, paste0(species_name, "_unique_genes.txt"))
        write.table(unique_genes, file = unique_genes_output_file, col.names = FALSE, row.names = FALSE, quote = FALSE)

        # Save the gene counts to a file
        gene_counts_output_file <- file.path(species_folder, paste0(species_name, "_gene_counts.txt"))
        write.table(as.data.frame(gene_counts), file = gene_counts_output_file, col.names = TRUE, row.names = FALSE, quote = FALSE)

        # Perform pathway analysis
        # Read the json dictionary file to get the organism code
        json_data <- jsonlite::fromJSON(dictionary_path)
        organism <- json_data[[species_name]]

        pwRes <- enrichKEGG(gene = unique_genes,  # Use the unique gene list for enrichKEGG
                            organism = organism,
                            keyType = "kegg",
                            pvalueCutoff = 0.05,
                            qvalueCutoff = 0.2)

        # Get the enriched pathways and associated genes
        pathways <- pwRes$ID
        genes <- pwRes$geneID

        output_file <- file.path(species_folder, paste0(species_name, "_pathway_data.txt"))
        cat("Enriched Pathways:\\n", file = output_file)
        for (i in seq_along(pathways)) {{
            cat("\\nPathway:", pathways[i], "\\n", "Associated Genes:", genes[i], "\\n", file = output_file, append = TRUE)
        }}

        # Save the enriched pathways to a file
        output_file_pwRes <- file.path(species_folder, paste0(species_name, "_pwRes_output.txt"))
        pwRes_text <- capture.output(print(pwRes))
        writeLines(pwRes_text, con = output_file_pwRes)
    }}
    '''
    robjects.r(r_code)

    # The following section will use parallel processing to quickly download the reference species data faster
    # This data gives a list of pathways and the genes they contain
    def process_species_subset(species_list, annotation_folder_path, dictionary_path, start_idx, end_idx):
        subset_species_list = species_list[start_idx:end_idx]
        reference_data(subset_species_list, annotation_folder_path, dictionary_path)

    if __name__ == "__main__":
        num_processes = 3  # Number of parallel processes you want to run

        # Divide the species list into equal parts for each process
        chunk_size = len(species_list) // num_processes
        processes = []

        for i in range(num_processes):
            start_idx = i * chunk_size
            end_idx = (i + 1) * chunk_size if i < num_processes - 1 else len(species_list)

            # Create a process for each subset of the species list
            p = multiprocessing.Process(target=process_species_subset,
                                        args=(species_list, reference_folder_path, main_csv, start_idx, end_idx))
            processes.append(p)
            p.start()
        # Wait for all processes to finish
        for p in processes:
            p.join()

    # Translate the unique gene names of each count matrix to KO codes, which are recognized for all species in the KEGG database
    print("Translated gene ID's")

    # Use all the data extracted to create a spreadsheet with bacteria names and their genes counts
    print("Made spreadsheet")

    # Translate the unique gene names of each count matrix to KO codes, which are recognized for all species in the KEGG database
    for species in species_list:
        count_matrix = f"{main_folder_path}/{species}/{species}_gene_counts.txt"
        translate_genes_to_K0s(species, count_matrix, main_csv)

    # Use all the data extracted to create a spreadsheet with bacteria names and their genes counts
    create_spreadsheet(species_list, main_folder_path)

if __name__ == "__main__":
    main()