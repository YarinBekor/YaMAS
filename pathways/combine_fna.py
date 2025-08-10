import os

# This function takes all the FNA files and combines them into one file for alignment
def combine_fna_files(species_list, annotation_directory, output_file_path):
    try:
        combined_fna_path = os.path.join(output_file_path, "combined_fna.fna")

        with open(combined_fna_path, "w") as output_file:
            for species in species_list:
                species_directory = os.path.join(annotation_directory, species)
                species_fna_file_path = os.path.join(species_directory, f"{species}.fna")

                # Check if the fna file for the species exists
                if os.path.exists(species_fna_file_path):
                    with open(species_fna_file_path, "r") as fna_file:
                        # Read the contents of the fna file and write it to the output file
                        output_file.write(fna_file.read())
                    print(f"FNA file for {species} added to the output file.")
                else:
                    print(f"Warning: FNA file for {species} not found.")

        print(f"All fna files have been combined into {combined_fna_path}.")
        return combined_fna_path

    except Exception as e:
        print(f"Error combining fna files: {e}")
        return None
