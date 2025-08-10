import os
import pandas as pd
import numpy as np

# This function take spreadsheets containing gene counts and uses statistical analysis to give them relative pathway scores
def calculate_log_adjusted(actual, average):
    if actual == 0 or average == 0:
        return 0  # Avoid division by zero
    return np.log(actual / average)

def calculate_pathway_sums(output_folder, pathway_dict_file, selected_pathways):
    # Step 1: Combine all spreadsheets into a single combined spreadsheet with averages
    combined_averages_df = None
    for number in range(10):
        spreadsheet_path = os.path.join(output_folder, f"fake_spreadsheet_{number}.csv")
        if os.path.exists(spreadsheet_path):
            df = pd.read_csv(spreadsheet_path, index_col=0)
            if combined_averages_df is None:
                combined_averages_df = df
            else:
                combined_averages_df += df

    # Calculate averages
    combined_averages_df /= 10  # Assuming 10 spreadsheets were used

    # Save the combined averages spreadsheet
    combined_averages_output_path = os.path.join(output_folder, "combined_spreadsheet_with_averages.csv")
    combined_averages_df.to_csv(combined_averages_output_path)
    print("Combined averages spreadsheet saved:", combined_averages_output_path)

    # Iterate over the range of spreadsheet numbers (0 to 9)
    for number in range(10):
        spreadsheet_path = os.path.join(output_folder, f"fake_spreadsheet_{number}.csv")
        log_adjusted_output_path = os.path.join(output_folder, f"log_adjusted_{number}.csv")

        if os.path.exists(spreadsheet_path):
            df = pd.read_csv(spreadsheet_path, index_col=0)
            log_adjusted_df = df.copy()

            for species in df.index:
                for gene in df.columns:
                    if species in combined_averages_df.index and gene in combined_averages_df.columns:
                        average_value = combined_averages_df.loc[species, gene]
                        actual_value = df.loc[species, gene]
                        log_adjusted_value = calculate_log_adjusted(actual_value, average_value)
                        log_adjusted_df.loc[species, gene] = log_adjusted_value

            log_adjusted_df.to_csv(log_adjusted_output_path)
            print("Log-adjusted spreadsheet saved:", log_adjusted_output_path)

    # Load pathway dictionary
    with open(pathway_dict_file, 'r') as f:
        pathway_dict = eval(f.read())

    # Iterate over the range of spreadsheet numbers (0 to 9)
    for number in range(10):
        spreadsheet_path = os.path.join(output_folder, f"log_adjusted_{number}.csv")
        pathway_sums_output_path = os.path.join(output_folder, f"pathway_sums_{number}.csv")

        if os.path.exists(spreadsheet_path):
            log_adjusted_df = pd.read_csv(spreadsheet_path, index_col=0)
            pathway_sums_df = pd.DataFrame(index=log_adjusted_df.index)

            for pathway in selected_pathways:
                pathway_genes = pathway_dict.get(pathway, [])
                pathway_sum = log_adjusted_df[pathway_genes].sum(axis=1)
                pathway_sums_df[pathway] = pathway_sum

            pathway_sums_df.to_csv(pathway_sums_output_path)
            print("Pathway sums spreadsheet saved:", pathway_sums_output_path)

output_folder = "/home/yisrael/Desktop/sample_csvs"
pathways_dict = "/home/yisrael/Desktop/pathways_files/pathway_dictionary.txt"
pathways = [ "00010" , "00020" , "00030"]
calculate_pathway_sums(output_folder, pathways_dict, pathways)