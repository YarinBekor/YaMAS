import os
import requests
from bs4 import BeautifulSoup
import time
import csv
import random

def fetch_more_gene_pathway_data(output_directory, species_code, gene_list, output_file):
    # Create the output directory if it doesn't exist
    os.makedirs(output_directory, exist_ok=True)

    # Create a new CSV file to store the pathway data
    output_file_path = os.path.join(output_directory, output_file)

    existing_paths = set()

    if os.path.exists(output_file_path):
        # Read the existing CSV file and extract the pathways
        with open(output_file_path, 'r') as csv_file:
            csv_reader = csv.reader(csv_file, delimiter=',')
            try:
                next(csv_reader)  # Skip the header row
            except StopIteration:
                # CSV file is empty, and there's no header row
                pass
            else:
                existing_paths = set(row[0] for row in csv_reader)

    with open(output_file_path, 'a', newline='') as output_file:
        csv_writer = csv.writer(output_file, delimiter=',')

        # If the file is empty or there's no header row, write the header row
        if output_file.tell() == 0:
            csv_writer.writerow(["Gene", "Pathways"])  # Updated header row

        error_links = []

        for gene in gene_list:
            if gene in existing_paths:
                continue

            retries = 3
            for attempt in range(retries):
                url = f"https://www.genome.jp/dbget-bin/get_linkdb?-t+pathway+{species_code}:{gene}"
                headers = {'User-Agent': 'Mozilla/5.0'}  # Use a popular browser's User-Agent
                response = requests.get(url, headers=headers)

                if response.status_code == 200:
                    soup = BeautifulSoup(response.content, 'html.parser')
                    path_info = soup.get_text()

                    separator = "--------------------------------------------------"
                    index = path_info.find(separator)
                    new_text = path_info[index + len(separator):]

                    lines = [line for line in new_text.splitlines()]

                    extracted_genes = []

                    for line in lines:
                        start_index_path = line.find(":") + 1
                        end_index_path = line.find('  ', start_index_path)
                        if start_index_path >= 0 and end_index_path > start_index_path:
                            extracted_gene = line[start_index_path:end_index_path]
                            extracted_gene = extracted_gene.strip()
                            extracted_genes.append(extracted_gene)

                    # Write pathway, functions, and genes to the CSV file
                    csv_writer.writerow([gene, ", ".join(extracted_genes)])
                    break
                elif response.status_code == 403:
                    if attempt < retries - 1:
                        # Retry with exponential backoff: wait longer after each attempt
                        print(f"retrying: {url}")
                        wait_time = 2 ** attempt + random.random()  # Add random jitter
                        time.sleep(wait_time)
                        print(f"Wait time: {wait_time}")
                        continue
                    else:
                        print(f"Request for {gene} failed after {retries} attempts.")
                        error_links.append(gene)
                else:
                    error_links.append(gene)

        if error_links:
            print("Some links encountered errors while fetching data:", error_links)

    return output_file_path
