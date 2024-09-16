import subprocess
import os

def install_metaphlan():
    try:
        # Use subprocess to run the pip install command for MetaPhlAn
        subprocess.check_call(["pip", "install", "metaphlan"])
        print("MetaPhlAn installed successfully.")
    except subprocess.CalledProcessError as e:
        print(f"Error installing MetaPhlAn: {e}")

def install_metaphlan_database(database_folder):
    try:
        # Use subprocess to run the command to install the MetaPhlAn database
        subprocess.check_call(["metaphlan", "--install", "--bowtie2db", database_folder])
        print("MetaPhlAn database installed successfully.")
    except subprocess.CalledProcessError as e:
        print(f"Error installing MetaPhlAn database: {e}")

def run_metaphlan(metagenome_fastq, database_folder, output_folder):
    try:
        # Construct the command as a list of individual arguments
        command = [
            "metaphlan",
            metagenome_fastq,
            "--bowtie2out", os.path.join(output_folder, "metagenome.bowtie2.bz2"),
            "--nproc", "5",
            "--input_type", "fastq",
            "-o", os.path.join(output_folder, "profiled_metagenome.txt"),
            "--bowtie2db", database_folder
        ]

        # Use subprocess to run the command
        subprocess.check_call(command)
        print("MetaPhlAn analysis completed successfully.")

        # Return the path to the output file
        output_file_path = os.path.join(output_folder, "profiled_metagenome.txt")
        return output_file_path

    except subprocess.CalledProcessError as e:
        print(f"Error running MetaPhlAn: {e}")
        return None