import os
import wget
import re
from bs4 import BeautifulSoup
import json
import requests
import find

# This program extracts FNA, GTF, and general table files from the NCBI database
# It extracts only the reference files used by KEGG and that contain the local KEGG gene IDs
def find_gtf_link(url):
    try:
        response = requests.get(url)
        response.raise_for_status()

        soup = BeautifulSoup(response.text, "html.parser")
        gtf_link = soup.find(href=re.compile(r'^GCA_.*?\.gtf\.gz$'))

        if gtf_link:
            return gtf_link['href']

    except requests.exceptions.RequestException as e:
        print(f"Error fetching the page: {e}")

    return None

def find_gff_link(url):
    try:
        response = requests.get(url)
        response.raise_for_status()

        soup = BeautifulSoup(response.text, "html.parser")
        gff_link = soup.find(href=re.compile(r'^GCA_.*?\.gff\.gz$'))

        if gff_link:
            return gff_link['href']

    except requests.exceptions.RequestException as e:
        print(f"Error fetching the page: {e}")

    return None

def find_table_link(url):
    try:
        response = requests.get(url)
        response.raise_for_status()

        soup = BeautifulSoup(response.text, "html.parser")
        table_link = soup.find(href=re.compile(r'^GCA_.*?_table\.txt\.gz$'))
        if table_link:
            return table_link['href']

    except requests.exceptions.RequestException as e:
        print(f"Error fetching the page: {e}")

    return None


def find_second_fna_link(url):
    try:
        response = requests.get(url)
        response.raise_for_status()

        soup = BeautifulSoup(response.text, "html.parser")
        fna_links = soup.find_all(href=re.compile(r'genomic\.fna\.gz$'))

        if len(fna_links) >= 2:
            return fna_links[1]['href']

    except requests.exceptions.RequestException as e:
        print(f"Error fetching the page: {e}")

    return None

def download_file(species_name, file_type, output_directory, main_csv):
    base_link = find.extract_csv(main_csv, species_name, 'genome_link')

    # Step 4: Extract the file links and download the file
    if file_type == "fna":
        file_link = find_second_fna_link(base_link)
    elif file_type == "table":
        file_link = find_table_link(base_link)
    elif file_type == "gtf":
        file_link = find_gtf_link(base_link)
    elif file_type == "gff":
        file_link = find_gff_link(base_link)
    else:
        print("Invalid file type. Use 'fna' or 'gtf' or 'table'.")
        return None

    # Step 5: Create the species folder if it doesn't exist
    species_folder = os.path.join(output_directory, species_name)
    if not os.path.exists(species_folder):
        os.makedirs(species_folder)

    # Step 6: Download the file and extract it
    if file_type == "fna":
        output_file = f"{species_name}.fna.gz"  # Set the filename with .fna.gz suffix
    elif file_type == "table":
        output_file = f"{species_name}.txt.gz"  # Set the filename with .txt.gz suffix
    elif file_type == "gtf":
        output_file = f"{species_name}.gtf.gz"  # Set the filename with .gtf.gz suffix
    elif file_type == "gff":
        output_file = f"{species_name}.gff.gz"  # Set the filename with .gtf.gz suffix
    else:
        print("Invalid file type. Use 'fna' or 'gtf' or 'table.")
        return None

    file_path = os.path.join(species_folder, output_file)
    file_url = f"{base_link}/{file_link}"

    wget.download(file_url, out=file_path)

    # Use the full filename for the gunzip command
    gz_file_path = file_path

    # Change working directory to the species folder before extraction
    os.chdir(species_folder)

    # Extract the gz file using gunzip command
    os.system(f"gunzip {gz_file_path}")

    if file_type == "fna":
        print(f"{species_name} FNA file downloaded and extracted")
    elif file_type == "table":
        print(f"{species_name} table file downloaded and extracted")
    elif file_type == "gtf":
        print(f"{species_name} GTF file downloaded and extracted")
    elif file_type == "gff":
        print(f"{species_name} GFF file downloaded and extracted")

    if file_link:
        return f"{base_link}/{file_link}"
