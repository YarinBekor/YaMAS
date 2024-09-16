import requests
from bs4 import BeautifulSoup
import re
import wget
import gzip
import os

# This program extracts species information directly from the NCBI database
# The reference genomes are up to date and accurate, but often don't contain the gene ID's used by KEGG
def add_underscore_if_needed(species, dashlist_path):
    try:
        with open(dashlist_path) as f:
            dashlist = f.read().splitlines()
            if species in dashlist:
                return f"_{species}"
            else:
                return species
    except FileNotFoundError:
        return species


def find_gtf_link(url):
    try:
        response = requests.get(url)
        response.raise_for_status()

        soup = BeautifulSoup(response.text, "html.parser")
        gtf_link = soup.find(href=re.compile(r'^GCF_.*?\.gtf\.gz$'))

        if gtf_link:
            return gtf_link['href']

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


def get_url(species_name, file_type):
    base_url = "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/"
    species_name = add_underscore_if_needed(species_name, '/home/yisrael/Desktop/Species/dashlist.txt')

    # Construct the URL based on the species name
    species_url = f"{base_url}{species_name}/all_assembly_versions"

    try:
        response = requests.get(species_url)
        response.raise_for_status()

        soup = BeautifulSoup(response.text, "html.parser")
        gcf_link = soup.find(href=re.compile(r'^GCF_.*?/'))

        if gcf_link:
            gcf_url = f"{species_url}/{gcf_link['href']}"
            if file_type == "fna":
                file_link = find_second_fna_link(gcf_url)
            elif file_type == "gtf":
                file_link = find_gtf_link(gcf_url)
            else:
                print("Invalid file type. Use 'fna' or 'gtf'.")
                return None

            if file_link:
                return f"{gcf_url}/{file_link}"

    except requests.exceptions.RequestException as e:
        print(f"Error fetching the page: {e}")

    return None


def download_file_NCBI(species_name, destination_folder, file_type):
    try:
        url = get_url(species_name, file_type)
        filename = url.split("/")[-1]

        # Create the directory for the species file if it doesn't exist
        species_file_path = os.path.join(destination_folder, species_name)
        os.makedirs(species_file_path, exist_ok=True)

        gz_file_path = os.path.join(species_file_path, filename)
        wget.download(url, gz_file_path)

        # Rename the downloaded file with appropriate extension
        if file_type == "fna":
            extracted_file_path = os.path.join(species_file_path, f"GCF_{species_name}.fna")
        elif file_type == "gtf":
            extracted_file_path = os.path.join(species_file_path, f"GCF_{species_name}.gtf")
        else:
            print("Invalid file type. Use 'fna' or 'gtf'.")
            return None

        with gzip.open(gz_file_path, 'rb') as f_in:
            with open(extracted_file_path, 'wb') as f_out:
                f_out.write(f_in.read())

        print(f"{file_type.upper()} File for {species_name} downloaded and extracted successfully!")

        return extracted_file_path

    except Exception as e:
        print(f"Error downloading or extracting the {file_type.upper()} file for {species_name}: {e}")
        return None
