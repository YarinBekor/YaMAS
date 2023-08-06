import requests
from bs4 import BeautifulSoup
import math
import re
from find import extract_csv


def extract_genes_and_KO(t_code, species_code, info_type):
    code = f"{species_code}"
    code_length = len(code)

    # Step 1: Create the base link using the species code
    genes_link = f"https://www.genome.jp/dbget-bin/get_linkdb?-t+genes+gn:{t_code}"

    response = requests.get(genes_link)
    soup = BeautifulSoup(response.text, 'html.parser')
    hits_tag = soup.find(string=lambda text: 'Hits:' in text)
    num_hits = int(re.search(r'\d+', hits_tag).group())
    num_pages = math.ceil(num_hits / 1000)

    extracted_genes = []
    extracted_K0s = []

    # Step 4 & 6: Extract the text from all pages and save it
    for page_number in range(1, num_pages + 1):
        url = f"https://www.genome.jp/dbget-bin/get_linkdb?-t+genes+-p+{page_number}+gn:{t_code}"
        response = requests.get(url)
        response.raise_for_status()  # Check for any errors in the request

        # Assuming the website uses UTF-8 encoding for text
        text_content = response.text

        # Remove lines containing "no KO assigned"
        lines = [line for line in text_content.splitlines() if "<a href=\"/entry/" in line
                 and "no KO assigned" not in line]

        for line in lines:
            if code_length == 3:
                start_index_gene = line.find(f"{species_code}:") + 4
            else:
                start_index_gene = line.find(f"{species_code}:") + 5
            end_index_gene = line.find('"', start_index_gene)
            if start_index_gene >= 0 and end_index_gene > start_index_gene:
                extracted_gene = line[start_index_gene:end_index_gene]
                extracted_gene = extracted_gene.strip()
                extracted_genes.append(extracted_gene)

            start_index_K0 = line.find("  K") + 2
            end_index_K0 = line.find(" ", start_index_K0)
            if start_index_K0 >= 0 and end_index_K0 > start_index_K0:
                extracted_K0 = line[start_index_K0:end_index_K0]
                extracted_K0 = extracted_K0.strip()
                extracted_K0s.append(extracted_K0)

    if info_type == "genes":
        return extracted_genes
    if info_type == "K0s":
        return extracted_K0s
    if info_type == "dictionary":
        gene_K0_dict = {}  # Initialize an empty dictionary

        if len(extracted_genes) == len(extracted_K0s):
            for i in range(len(extracted_genes)):
                gene = extracted_genes[i]
                K0 = extracted_K0s[i]
                gene_K0_dict[gene] = K0  # Map gene to K0 in the dictionary
            return gene_K0_dict
        else:
            print("Lists of genes and K0s have different lengths.")
            return None

def translate_genes_to_K0s(species_name, count_matrix_filename, csv_filename):
    t_code = extract_csv(csv_filename, species_name, "t_number")
    species_code = extract_csv(csv_filename, species_name, "code")

    gene_K0_dict = extract_genes_and_KO(t_code, species_code, "dictionary")

    if gene_K0_dict is None:
        print("Failed to create the gene-K0 dictionary.")
        return None

    translated_lines = []
    with open(count_matrix_filename, 'r') as matrix_file:
        gene_K0_counts = {}

        for line in matrix_file:
            gene_name, count = line.strip().split()
            if gene_name in gene_K0_dict:
                K0_name = gene_K0_dict[gene_name]
                if K0_name in gene_K0_counts:
                    gene_K0_counts[K0_name] += int(count)
                else:
                    gene_K0_counts[K0_name] = int(count)

        for K0_name, total_count in gene_K0_counts.items():
            translated_line = f"{K0_name} {total_count}"
            translated_lines.append(translated_line)

    translated_filename = count_matrix_filename.replace(".txt", "_translated.txt")

    with open(translated_filename, 'w') as translated_file:
        for line in translated_lines:
            translated_file.write(line + '\n')

    print(f"Translated count matrix saved to: {translated_filename}")

    return translated_lines

