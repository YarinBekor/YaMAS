import csv

# This function extracts specified information for a given species from the main csv file

def extract_csv(csv_file, species_name, info_type):
    with open(csv_file, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            if row['species_name'].lower() == species_name.lower():
                if info_type == 't_number':
                    return row['t_number']
                elif info_type == 'code':
                    return row['code']
                elif info_type == 'genome_link':
                    return row['genome_link']
                else:
                    return "Invalid info_type provided. Please choose from 't_number', 'species_kegg_code', or 'species_archive_link'."

