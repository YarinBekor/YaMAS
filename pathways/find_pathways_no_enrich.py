import pandas as pd

def filter_genes(input_csv, genes_list_file, output_csv):
    import pandas as pd

    # Read the CSV file without header
    df = pd.read_csv(input_csv, sep='\t', header=None, names=["Pathway Name", "Pathway Function", "Genes"], skiprows=1)

    # Read the gene list from the text file and convert it to a set for faster membership checking
    with open(genes_list_file, "r") as file:
        genes_list = {gene.strip() for gene in file}

    # Function to filter genes for each row
    def filter_genes_row(row):
        genes = row["Genes"]
        if isinstance(genes, str):
            genes = genes.split(", ")
            filtered_genes = [gene for gene in genes if gene in genes_list]
            return ", ".join(filtered_genes)
        return genes

    # Apply the filtering function to the dataframe
    df["Genes"] = df.apply(filter_genes_row, axis=1)

    # Save the filtered dataframe to a new CSV file
    df.to_csv(output_csv, sep=',', index=False)
