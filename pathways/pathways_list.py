import requests

# This function extracts a list of pathways for a given species using its t_number
def extract_paths_from_pages(t_number, info_type):
    url = f"https://www.genome.jp/dbget-bin/get_linkdb?-t+pathway+gn:{t_number}"

    try:
        response = requests.get(url)
        response.raise_for_status()  # Check for any errors in the request

        # Assuming the website uses UTF-8 encoding for text
        text_content = response.text

        # Remove lines containing "no KO assigned"
        lines = [line for line in text_content.splitlines() if "<a href=\"/pathway/" in line]

        extracted_paths = []
        extracted_functions = []

        for line in lines:
            start_index_path = line.find("y/") + 2
            end_index_path = line.find('"', start_index_path)
            if start_index_path >= 0 and end_index_path > start_index_path:
                extracted_path = line[start_index_path:end_index_path]
                extracted_path = extracted_path.strip()
                extracted_paths.append(extracted_path)

            start_index_function = line.find("            ") + 12
            end_index_function = line.find(" -", start_index_function)
            if start_index_function >= 0 and end_index_function > start_index_function:
                extracted_function = line[start_index_function:end_index_function]
                extracted_function = extracted_function.strip()
                extracted_functions.append(extracted_function)


        if not extracted_paths:
            print("No data found for the given t_number.")
            return

        if info_type == "pathways":
            return extracted_paths
        elif info_type == "functions":
            return extracted_functions
        else:
            print("Invalid info_type provided. Please use 'pathways' or 'functions'.")

    except requests.exceptions.RequestException as e:
        print("Error occurred while fetching the URL:", e)
        return None

    except requests.exceptions.RequestException as e:
        print("Error occurred while fetching the URL:", e)

    if info_type == "name":
        return extracted_paths

    if info_type == "function":
        return extracted_functions
