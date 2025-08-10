import csv

# This code removes all zero entries from a count matrix
def remove_zeros_from_count_matrix(file_path):
    remaining_genes = []

    # Read the count matrix from the CSV file
    with open(file_path, 'r') as csv_file:
        csv_reader = csv.reader(csv_file)
        header = next(csv_reader)  # Assuming the first row contains the header
        for row in csv_reader:
            gene = row[0]  # Assuming the first column contains the gene names
            counts = list(map(int, row[1:]))  # Convert the count values to integers
            if any(counts):  # Check if any count value is not zero
                remaining_genes.append(gene)

    return remaining_genes
