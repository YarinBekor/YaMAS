import os
import pandas as pd
from findspecies import search_species

def create_spreadsheet(bacteria_list, main_folder_path):
    species_data = {}  # Dictionary to store species data

    # Loop through each species
    for species_name in bacteria_list:
        species_folder = os.path.join(main_folder_path, species_name)
        count_matrix_filename = os.path.join(species_folder, f"{species_name}_gene_counts_translated.txt")

        if os.path.exists(count_matrix_filename):
            # Read the count matrix file for the current species
            with open(count_matrix_filename, 'r') as matrix_file:
                genes_counts = [line.strip().split() for line in matrix_file]

            # Extract genes and counts for the current species
            genes = [gene_count[0] for gene_count in genes_counts]
            counts = [int(gene_count[1]) for gene_count in genes_counts]

            # Add species data to the dictionary
            species_data[species_name] = dict(zip(genes, counts))

    # Create a DataFrame from the collected species data
    main_spreadsheet = pd.DataFrame(species_data)

    # Transpose the DataFrame to have genes as rows and bacteria as columns
    main_spreadsheet = main_spreadsheet.transpose()

    # Save the main spreadsheet to the main folder path
    spreadsheet_filename = os.path.join(main_folder_path, "large_spreadsheet.csv")
    main_spreadsheet.to_csv(spreadsheet_filename)
    print(f"Spreadsheet saved to: {spreadsheet_filename}")
